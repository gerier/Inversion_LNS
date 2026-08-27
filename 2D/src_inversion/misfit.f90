subroutine f(m,fm)
!-----------------------------------------------------------------------
! Compute the objective function.
!
! The objective function consists of:
!   (i) a data misfit
!   (ii) an optional regularization term.
!
! The forward problem is solved for the current model and the synthetic
! pressure field is compared with the observed pressure at each receiver.
!-----------------------------------------------------------------------
use MPI
use parameters
implicit none
 double precision :: fm
 double precision, dimension(1:Nflat) :: m
 integer :: irec
 double precision, dimension(NREC) :: fm_local_per_rec
 double precision :: regul_term_p0, regul_term_rho0,regul_term_windx
 double precision, dimension(1:Nflat) :: normsq_m, normsq_mm0
 
 ! Initialize variables
 fm_local_per_rec(:) = 0.0d0
 factor_regul_SRdist(:) = 1.0

 call MPI_BARRIER(MPI_COMM_WORLD, code)
 
 call flatmodel2priormodel(m)

! solve the forward problem
 call reset_forward() 
 call forwardproblem(p0_prior, rho0_prior, windx_prior, windy_prior,  1, NSTEP, 1) 
 sispressure_prior(:,:) = sispressure(:,:)
 
! misfit of pressure waveform
 do irec=1,NREC
   if (i_rank == (ix_rec(irec)-1)/NX_LOCAL .and. j_rank == (iy_rec(irec)-1)/NY_LOCAL) then
   if ( normsq_pressure_true_per_rec(irec) > TINYVAL) then
   fm_local_per_rec(irec) = sum((sispressure_prior(:,irec) - sispressure_true(:,irec))**2) &
                    / normsq_pressure_true_per_rec(irec) 
   endif
   endif
 enddo

 call MPI_BARRIER(MPI_COMM_WORLD, code)
 call MPI_ALLREDUCE(sum(fm_local_per_rec), fm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)

 fm = DELTAT * DELTAX * DELTAY * fm / 2
 fx_data = fm 

 count_f = count_f + 1
 
 ! add a regularization term
 if (type_regul_term == 1) then
  normsq_m(:) = 0.0d0
  normsq_mm0= factor_regul_SRdist*(m-m0)**2
  normsq_mm0(1:2*NY_LOCAL) = normsq_mm0(1:2*NY_LOCAL) / m0(1:2*NY_LOCAL)**2
  
  call MPI_ALLREDUCE( sum(normsq_mm0(1:NY_LOCAL)), regul_term_rho0, &
              1,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, code)
  call MPI_ALLREDUCE( sum(normsq_mm0(NY_LOCAL+1:2*NY_LOCAL)), regul_term_p0, &
              1,MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,code)
  call MPI_ALLREDUCE( sum(normsq_mm0(2*NY_LOCAL+1:Nflat)), regul_term_windx, &
              1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD, code)
                                                                                    
  fx_regul =  regul_weight * 0.5d0 * DELTAY * &
         (regul_term_p0 + regul_term_rho0 + regul_term_windx)
  fm = fm + fx_regul
 
 else
  if (type_regul_term /= 0) then 
   print *, "ERROR: Unknown regularization type"
   stop
   endif
 endif

endsubroutine f


!!!! Note: the regularization term requires further validation (TODO).
! Special care is required when m = m0
! The implementation depends on the selected parameterization.

subroutine df(m,flat_grad)
!-----------------------------------------------------------------------
! Compute the gradient of the objective function.
!
! The gradient is obtained from the sensitivity kernels and optionally
! augmented with the gradient of the regularization term.
!-----------------------------------------------------------------------
use parameters
implicit none 
double precision, dimension(1:Nflat) :: m,flat_grad,reg_grad
double precision :: ONE_OVER_DX2, ONE_OVER_DY2, ONE_OVER_DXDY

  ! initialization
  ONE_OVER_DX2 = 1.0d0 / DELTAX**4
  ONE_OVER_DY2 = 1.0d0 / DELTAY**4
  ONE_OVER_DXDY = 1.0d0 / DELTAX**2 / DELTAY**2

  factor_regul_SRdist(:) = 1.0
 
  ! sensitivity kernel computations
  call MPI_BARRIER(MPI_COMM_WORLD, code)
  call flatmodel2priormodel(m)
  call compute_kernel()
  call kernelparam2inversionparam(flat_grad)

  ! regularization term
  if (type_regul_term == 1) then
 
    reg_grad(:) = 0.0d0
    reg_grad(1:NY_LOCAL) = (m(1:NY_LOCAL) - m0(1:NY_LOCAL)) /m0(1:NY_LOCAL)**2
    reg_grad(NY_LOCAL+1:2*NY_LOCAL) = (m(NY_LOCAL+1:2*NY_LOCAL) - m0(NY_LOCAL+1:2*NY_LOCAL)) / m0(NY_LOCAL+1:2*NY_LOCAL)**2
    reg_grad(2*NY_LOCAL+1:3*NY_LOCAL) = (m(2*NY_LOCAL+1:3*NY_LOCAL) - m0(2*NY_LOCAL+1:3*NY_LOCAL)) 

    flat_grad = flat_grad + DELTAY * regul_weight * reg_grad * factor_regul_SRdist
   
  else 
    if (type_regul_term /= 0) then 
      print *, "ERROR: Type of the regularisation term: Unknown"
      stop
    endif
  endif

  count_grad = count_grad + 1
endsubroutine df



subroutine kernelparam2inversionparam(flat_grad)
!-----------------------------------------------------------------------
! Convert sensitivity kernels into the inversion parameterization.
!
! Depending on the selected parameterization, the physical kernels
! (density, pressure, horizontal wind) are transformed into the
! variables optimized by the inversion.
!-----------------------------------------------------------------------
use mpi
use parameters, only : p0_prior, rho0_prior, c0_prior, gamma_chemestry, &
                       K_rho0, K_p0, K_windx, &
                       NX_LOCAL, NY_LOCAL, Nflat, ZERO, JSOURCE, &
                       parametrisation, scale_model, &
                       code, ierr, j_rank, row_Comm         
implicit none
double precision, dimension (1:NX_LOCAL,1:NY_LOCAL) :: grad_c0, grad_rho0, grad_lnc0, grad_lnrho0, grad_lnp0
double precision, dimension(1:Nflat) :: flat_grad
double precision :: scalar_grad_lnrho0, scalar_grad_lnc0, scalar_grad_lnp0, &
                    scalar_grad_rho0, scalar_grad_c0, scalar_grad_p0, &
                    scalar_grad_windx
integer :: j

  if (parametrisation == 1) then
    ! Density, wave speed, horizontal wind
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL) &
                  * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_c0(:,:) = 2 * c0_prior(1:NX_LOCAL,1:NY_LOCAL) * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                    / gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL) &
                    * K_p0(1:NX_LOCAL,1:NY_LOCAL)
    grad_rho0(:,:) = p0_prior(1:NX_LOCAL,1:NY_LOCAL) / rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL) &
                     + K_rho0(1:NX_LOCAL,1:NY_LOCAL) 

    do j=1,NY_LOCAL
      ! Sum contributions from all MPI processes 
      call MPI_ALLREDUCE( sum(grad_c0(:,j)), scalar_grad_c0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(grad_rho0(:,j)), scalar_grad_rho0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_windx(:,j)), scalar_grad_windx,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM, row_Comm,ierr)
      ! Initialize the flattened gradient.
      flat_grad(j) = ZERO
      flat_grad(NY_LOCAL+j) = ZERO
      ! Store the gradient components in the flattened gradient vector
      flat_grad(j) = scale_model(1) * scalar_grad_rho0
      flat_grad(NY_LOCAL + j) = scale_model(2) * scalar_grad_c0
      flat_grad(2*NY_LOCAL + j)= scale_model(3) * scalar_grad_windx
    enddo

    ! Interpolate the gradient at the source location
    if (j_rank == JSOURCE / NY_LOCAL) then
      flat_grad(jsource) = (flat_grad(jsource-1) + flat_grad(jsource+1)) *0.5
      flat_grad(NY_LOCAL + jsource) = (flat_grad(NY_LOCAL+ jsource-1) + flat_grad(NY_LOCAL+jsource+1)) * 0.5
      flat_grad(2*NY_LOCAL+jsource) = (flat_grad(2*NY_LOCAL + jsource-1) + flat_grad(2*NY_LOCAL + jsource+1)) *0.5
    endif
    
  elseif (parametrisation == 2) then
    ! Density, pressure, horizontal wind
    do j=1,NY_LOCAL
      ! Sum contributions from all MPI processes
      call MPI_ALLREDUCE( sum(K_rho0(:,j)), scalar_grad_rho0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_p0(:,j)), scalar_grad_p0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_windx(:,j)), scalar_grad_windx,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM, row_Comm,ierr)
      ! Initialize
      flat_grad(j) = ZERO
      flat_grad(NY_LOCAL+j) = ZERO
      ! Store the gradient components in the flattened gradient vector
      flat_grad(j) = scale_model(1) * scalar_grad_rho0
      flat_grad(NY_LOCAL + j) = scale_model(2) * scalar_grad_p0
      flat_grad(2*NY_LOCAL + j) = scale_model(3) * scalar_grad_windx
    enddo
    
  elseif (parametrisation == 3) then
   ! Log-density, log-wave speed, horizontal wind
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL) &
                     * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_lnc0(:,:) = 2 * c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2 * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                      / gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL) &
                       * K_p0(1:NX_LOCAL,1:NY_LOCAL)
    grad_lnrho0(:,:) =  p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL) &
                       + rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL) 
  
    call MPI_BARRIER(MPI_COMM_WORLD, code)
  
    do j=1,NY_LOCAL
      ! Sum contributions from all MPI processes
      call MPI_ALLREDUCE( sum(grad_lnrho0(:,j)), scalar_grad_lnrho0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(grad_lnc0(:,j)), scalar_grad_lnc0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_windx(:,j)), scalar_grad_windx,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM, row_Comm,ierr)

      flat_grad(j) = scale_model(1) * scalar_grad_lnrho0
      flat_grad(NY_LOCAL + j) = scale_model(2) * scalar_grad_lnc0
      flat_grad(2*NY_LOCAL + j) = scale_model(3) * scalar_grad_windx
    enddo
  

   ! Interpolate the gradient at the source location
   if (j_rank == JSOURCE / NY_LOCAL) then
     flat_grad(jsource) = (flat_grad(jsource-1) + flat_grad(jsource+1)) *0.5 
     flat_grad(NY_LOCAL + jsource) = (flat_grad(NY_LOCAL+ jsource-1) + flat_grad(NY_LOCAL+jsource+1)) * 0.5
     flat_grad(2*NY_LOCAL+jsource) = (flat_grad(2*NY_LOCAL + jsource-1) + flat_grad(2*NY_LOCAL + jsource+1)) *0.5
   endif

   call MPI_BARRIER(MPI_COMM_WORLD, code)    

    
  elseif (parametrisation == 35) then
   print *, "Selected parameterization: ", parametrisation
   ! Log-density, wave speed, horizontal wind
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL) &
                     * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_c0(:,:) = 2 * c0_prior(1:NX_LOCAL,1:NY_LOCAL) * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                      / gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL) &
                       * K_p0(1:NX_LOCAL,1:NY_LOCAL)
    grad_lnrho0(:,:) =  p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL) &
                       + rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL) 
  
    call MPI_BARRIER(MPI_COMM_WORLD, code)
  
    do j=1,NY_LOCAL
      call MPI_ALLREDUCE( sum(grad_lnrho0(:,j)), scalar_grad_lnrho0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(grad_c0(:,j)), scalar_grad_c0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_windx(:,j)), scalar_grad_windx,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM, row_Comm,ierr)

      flat_grad(j) = scale_model(1) * scalar_grad_lnrho0
      flat_grad(NY_LOCAL + j) = scale_model(2) * scalar_grad_c0 
      flat_grad(2*NY_LOCAL + j) = scale_model(3) * scalar_grad_windx
    enddo
  

   ! Interpolate the gradient at the source location
   if (j_rank == JSOURCE / NY_LOCAL) then
     flat_grad(jsource) = (flat_grad(jsource-1) + flat_grad(jsource+1)) *0.5 
     flat_grad(NY_LOCAL + jsource) = (flat_grad(NY_LOCAL+ jsource-1) + flat_grad(NY_LOCAL+jsource+1)) * 0.5
     flat_grad(2*NY_LOCAL+jsource) = (flat_grad(2*NY_LOCAL + jsource-1) + flat_grad(2*NY_LOCAL + jsource+1)) *0.5
   endif

   call MPI_BARRIER(MPI_COMM_WORLD, code)    


  elseif (parametrisation == 36) then
   print *, "Selected parameterization: ", parametrisation
   ! Log-normalised-density, log-wave speed, horizontal wind
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL) &
                     * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_c0(:,:) = 2 * c0_prior(1:NX_LOCAL,1:NY_LOCAL) * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                      / gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL) &
                       * K_p0(1:NX_LOCAL,1:NY_LOCAL)
    grad_lnrho0(:,:) =  p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL) &
                       + rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL) 
  
    call MPI_BARRIER(MPI_COMM_WORLD, code)
  
    do j=1,NY_LOCAL
      call MPI_ALLREDUCE( sum(grad_lnrho0(:,j)), scalar_grad_lnrho0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(grad_c0(:,j)), scalar_grad_c0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_windx(:,j)), scalar_grad_windx,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM, row_Comm,ierr)


      flat_grad(j) = scale_model(1) * scalar_grad_lnrho0
      flat_grad(NY_LOCAL + j) = scale_model(2) * scalar_grad_c0 
      flat_grad(2*NY_LOCAL + j) = scale_model(3) * scalar_grad_windx
    enddo
  

    ! Interpolate the gradient at the source location
    if (j_rank == JSOURCE / NY_LOCAL) then
      flat_grad(jsource) = (flat_grad(jsource-1) + flat_grad(jsource+1)) *0.5 
      flat_grad(NY_LOCAL + jsource) = (flat_grad(NY_LOCAL+ jsource-1) + flat_grad(NY_LOCAL+jsource+1)) * 0.5
      flat_grad(2*NY_LOCAL+jsource) = (flat_grad(2*NY_LOCAL + jsource-1) + flat_grad(2*NY_LOCAL + jsource+1)) *0.5
   endif

   call MPI_BARRIER(MPI_COMM_WORLD, code)


  elseif (parametrisation == 4) then
    ! Log-density, log-pressure, horizontal wind
    
    grad_lnp0(:,:) = p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL)
    grad_lnrho0(:,:) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL)
    
    do j=1,NY_LOCAL
      ! Sum contributions from all MPI processes
      call MPI_ALLREDUCE( sum(grad_lnp0(:,j)), scalar_grad_lnp0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(grad_lnrho0(:,j)), scalar_grad_lnrho0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_windx(:,j)), scalar_grad_windx,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM, row_Comm,ierr)
      ! Initialize
      flat_grad(j) = ZERO
      flat_grad(NY_LOCAL+j) = ZERO
      ! Store the gradient components in the flattened gradient vector
      flat_grad(j) = scale_model(1) * scalar_grad_lnrho0
      flat_grad(NY_LOCAL + j) = scale_model(2) * scalar_grad_lnp0
      flat_grad(2*NY_LOCAL + j) = scale_model(3) * scalar_grad_windx
    enddo
    
  elseif (parametrisation == 5) then
    ! Log-wave speed, log-pressure, horizontal wind
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL) &
                          * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_lnc0(:,:) = - 2 * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL)
    grad_lnp0(:,:) = p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL) &
                       + rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL) 
    
    do j=1,NY_LOCAL
      ! Log-wave speed, log-pressure, horizontal wind
      call MPI_ALLREDUCE( sum(grad_lnp0(:,j)), scalar_grad_lnp0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(grad_lnc0(:,j)), scalar_grad_lnc0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_windx(:,j)), scalar_grad_windx,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM, row_Comm,ierr)
      ! Initialize
      flat_grad(j) = ZERO
      flat_grad(NY_LOCAL+j) = ZERO
      ! Store the gradient components in the flattened gradient vector
      flat_grad(j) = scale_model(1) * scalar_grad_lnc0
      flat_grad(NY_LOCAL + j) = scale_model(2) * scalar_grad_lnp0
      flat_grad(2*NY_LOCAL + j) = scale_model(3) * scalar_grad_windx
    enddo 
  else
    print *, "ERROR: Unknown parameterization"
    stop  
  endif
  

endsubroutine kernelparam2inversionparam


subroutine flatmodel2priormodel(flat_model)
!-----------------------------------------------------------------------
! Convert the flattened inversion vector into the physical model.
!
! The flattened vector contains one vertical profile per parameter:
!
!   [ parameter_1(y),
!     parameter_2(y),
!     horizontal_wind(y) ]
!
! The meaning of each parameter depends on the selected
! inversion parameterization.
!
! The profiles are expanded over the computational domain and ghost
! cells are updated through MPI communications.
!-----------------------------------------------------------------------
use parameters
implicit none
double precision, dimension(1:Nflat) :: flat_model
double precision, dimension(1:NY_LOCAL) :: flat_rho0, flat_c0, flat_p0,flat_windx
integer :: j

  if (parametrisation == 1) then
   ! Density, wave speed, horizontal wind
   flat_rho0 = scale_model(1)  * flat_model(:NY_LOCAL)
   flat_c0 = scale_model(2)  * flat_model(NY_LOCAL+1:2*NY_LOCAL)
   flat_windx = scale_model(3) * flat_model(1+2*NY_LOCAL:Nflat)
   
   do j=1,NY_LOCAL
      rho0_prior(:,j) = flat_rho0(j)
      c0_prior(:,j) = flat_c0(j)
      windx_prior(:,j) = flat_windx(j)
   enddo
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                     * c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2&
                      / gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL)
   
  else if (parametrisation == 2) then
   ! Density, pressure, horizontal wind
   flat_rho0 = scale_model(1)  * flat_model(:NY_LOCAL)
   flat_p0 = scale_model(2)  * flat_model(NY_LOCAL+1:2*NY_LOCAL)
   flat_windx = scale_model(3) * flat_model(1+2*NY_LOCAL:)
   
   do j=1,NY_LOCAL
      rho0_prior(:,j) = flat_rho0(j)
      p0_prior(:,j) = flat_p0(j)
      windx_prior(:,j) = flat_windx(j)
   enddo
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt( p0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                          * gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL) &
                             /rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
   
  else if (parametrisation == 35) then
   ! Log-density, wave speed, horizontal wind
   flat_rho0(1:NY_LOCAL) = exp(scale_model(1)  * flat_model(1:NY_LOCAL))
   flat_c0(1:NY_LOCAL) = scale_model(2)  * flat_model(NY_LOCAL+1:2*NY_LOCAL)
   flat_windx(1:NY_LOCAL) = scale_model(3) * flat_model(1+2*NY_LOCAL:Nflat)
   
   do j=1,NY_LOCAL
      rho0_prior(:,j) = flat_rho0(j)
      c0_prior(:,j) = flat_c0(j)
      windx_prior(:,j) = flat_windx(j)
   enddo
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * &
                                   c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2 &
                                   / gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL)
  

  else if (parametrisation == 36) then
   ! Log-normalised-density, log-wave speed, horizontal wind
   flat_rho0(1:NY_LOCAL) = density_prior * exp(scale_model(1)  * flat_model(1:NY_LOCAL))
   flat_c0(1:NY_LOCAL) = scale_model(2)  * flat_model(NY_LOCAL+1:2*NY_LOCAL)
   flat_windx(1:NY_LOCAL) = scale_model(3) * flat_model(1+2*NY_LOCAL:Nflat)
   
   do j=1,NY_LOCAL
      rho0_prior(:,j) = flat_rho0(j)
      c0_prior(:,j) = flat_c0(j)
      windx_prior(:,j) = flat_windx(j)
   enddo
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * &
                                   c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2 &
                                   / gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL)


  else if (parametrisation == 3) then
   ! Log-density, log-wave speed, horizontal wind
   flat_rho0(1:NY_LOCAL) = exp(scale_model(1)  * flat_model(1:NY_LOCAL))
   flat_c0(1:NY_LOCAL) = exp(scale_model(2)  * flat_model(NY_LOCAL+1:2*NY_LOCAL))
   flat_windx(1:NY_LOCAL) = scale_model(3) * flat_model(1+2*NY_LOCAL:Nflat)
   
   do j=1,NY_LOCAL
      rho0_prior(:,j) = flat_rho0(j)
      c0_prior(:,j) = flat_c0(j)
      windx_prior(:,j) = flat_windx(j)
   enddo
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * &
                                   c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2 &
                                   / gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL)

  else if (parametrisation == 4) then
   ! Log-density, log-pressure, horizontal wind
  
   flat_rho0 = exp(scale_model(1)  * flat_model(1:NY_LOCAL))
   flat_p0 = exp(scale_model(2)  * anint(1e16*flat_model(NY_LOCAL+1:2*NY_LOCAL))/1e16)
   flat_windx = scale_model(3) * flat_model(1+2*NY_LOCAL:)

   do j=1,NY_LOCAL
      rho0_prior(:,j) = flat_rho0(j)  
      p0_prior(:,j) = flat_p0(j)
      windx_prior(:,j) = flat_windx(j)
   enddo
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt( p0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                             * gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL) &
                             /rho0_prior(1:NX_LOCAL,1:NY_LOCAL))

 
  elseif (parametrisation == 5) then
  ! Log-wave speed, log-pressure, horizontal wind
  
   flat_c0 = exp(scale_model(1)  * flat_model(1:NY_LOCAL))
   flat_p0 = exp(scale_model(2)  * flat_model(NY_LOCAL+1:2*NY_LOCAL))
   flat_windx = scale_model(3) * flat_model(1+2*NY_LOCAL:)
   
   do j=1,NY_LOCAL
      c0_prior(:,j) = flat_c0(j)  
      p0_prior(:,j) = flat_p0(j)
      windx_prior(:,j) = flat_windx(j)
   enddo
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) =  p0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                             * gamma_chemestry(1:NX_LOCAL,1:NY_LOCAL) &
                             /c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2
                             
  else
    print *, "ERROR: Unknown parameterization"
    stop 
    
  endif
  
  call MPI_BARRIER(MPI_COMM_WORLD, code)    
                                    
  call send_receive_rightleft(windx_prior)
  call send_receive_rightleft(windy_prior)
  call send_receive_rightleft(rho0_prior)
  call send_receive_rightleft(p0_prior)
  
  call send_receive_topbottom(windx_prior)
  call send_receive_topbottom(windy_prior)
  call send_receive_topbottom(rho0_prior)
  call send_receive_topbottom(p0_prior)
     
  call send_receive_corners(windx_prior) 
  call send_receive_corners(windy_prior)
  call send_receive_corners(rho0_prior)
  call send_receive_corners(p0_prior)


  if (USE_PML_XMIN .and. i_rank == 0) then
    windx_prior(-1,:) = windx_prior(1,:) 
    windx_prior(0,:) = windx_prior(1,:) 
    windy_prior(-1,:) = windy_prior(1,:) 
    windy_prior(0,:) = windy_prior(1,:) 
    rho0_prior(-1,:) = rho0_prior(1,:) 
    rho0_prior(0,:) = rho0_prior(1,:)
    p0_prior(-1,:) = p0_prior(1,:) 
    p0_prior(0,:) = p0_prior(1,:) 
  endif
  
  if (USE_PML_XMAX .and. i_rank == NPROC_X-1) then

    windx_prior(NX_LOCAL+1,:) = windx_prior(NX_LOCAL,:) 
    windx_prior(NX_LOCAL+2,:) = windx_prior(NX_LOCAL,:) 
    windy_prior(NX_LOCAL+1,:) = windy_prior(NX_LOCAL,:) 
    windy_prior(NX_LOCAL+2,:) = windy_prior(NX_LOCAL,:) 
    rho0_prior(NX_LOCAL+1,:) = rho0_prior(NX_LOCAL,:) 
    rho0_prior(NX_LOCAL+2,:) = rho0_prior(NX_LOCAL,:)
    p0_prior(NX_LOCAL+1,:) = p0_prior(NX_LOCAL,:) 
    p0_prior(NX_LOCAL+2,:) = p0_prior(NX_LOCAL,:) 
    
  endif
  
  if ( .not. USE_PML_YMIN .and. USE_PML_YMAX .and. j_rank == 0) then
    do j=-1,1
      p0_prior(:,j) = p0_prior(:,3-j)
      rho0_prior(:,j) = rho0_prior(:,3-j)
      windx_prior(:,j) = - windx_prior(:,3-j)
    enddo    
     
  endif

  if (USE_PML_YMAX .and. j_rank == NPROC_Y-1) then
    windx_prior(:,NY_LOCAL+1) = windx_prior(:,NY_LOCAL) 
    windx_prior(:,NY_LOCAL+2) = windx_prior(:,NY_LOCAL) 
    rho0_prior(:,NY_LOCAL+1) = rho0_prior(:,NY_LOCAL) 
    rho0_prior(:,NY_LOCAL+2) = rho0_prior(:,NY_LOCAL)
    p0_prior(:,NY_LOCAL+1) = p0_prior(:,NY_LOCAL) 
    p0_prior(:,NY_LOCAL+2) = p0_prior(:,NY_LOCAL) 
  endif
  
  if (USE_PML_YMAX .and. USE_PML_XMAX .and. i_rank == NPROC_X-1 .and. j_rank == NPROC_Y-1) then
    windx_prior(NX_LOCAL+1:NX_LOCAL+2,NY_LOCAL+1:NY_LOCAL+2) = windx_prior(NX_LOCAL,NY_LOCAL) 
    rho0_prior(NX_LOCAL+1:NX_LOCAL+2,NY_LOCAL+1:NY_LOCAL+2) = rho0_prior(NX_LOCAL,NY_LOCAL) 
    p0_prior(NX_LOCAL+1:NX_LOCAL+2,NY_LOCAL+1:NY_LOCAL+2) = p0_prior(NX_LOCAL,NY_LOCAL) 
  endif
  
    if (USE_PML_YMAX .and. USE_PML_XMIN .and. i_rank == 0 .and. j_rank == NPROC_Y-1) then
    windx_prior(-1:0,NY_LOCAL+1:NY_LOCAL+2) = windx_prior(1,NY_LOCAL) 
    rho0_prior(-1:0,NY_LOCAL+1:NY_LOCAL+2) = rho0_prior(1,NY_LOCAL) 
    p0_prior(-1:0,NY_LOCAL+1:NY_LOCAL+2) = p0_prior(1,NY_LOCAL) 
  endif

endsubroutine flatmodel2priormodel



subroutine priormodel2flatmodel(flat_model)
!-----------------------------------------------------------------------
! Convert the current physical model into the flattened inversion vector.
! The inversion vector contains one-dimensional vertical profiles:
!
!   flat_model =
!     [ rho(y),
!       c(y),
!       wind_x(y) ]
!
! depending on the selected parameterization.
!-----------------------------------------------------------------------
use parameters, only : NY_LOCAL, Nflat,density_prior, &
                       scale_model, rho0_prior, c0_prior, p0_prior, windx_prior, parametrisation
implicit none
double precision, dimension(1:Nflat):: flat_model
 
  if (parametrisation == 1) then
  ! Density, wave speed, horizontal wind
  flat_model(1:NY_LOCAL) = rho0_prior(25,1:NY_LOCAL) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = c0_prior(25,1:NY_LOCAL)/ scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)
    
  elseif (parametrisation == 2) then
  ! Density, pressure, horizontal wind
  flat_model(1:NY_LOCAL) = rho0_prior(25,1:NY_LOCAL) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = p0_prior(25,1:NY_LOCAL)/ scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)
  
  elseif (parametrisation == 3) then
  ! Log-density, log-wave speed, horizontal wind
  flat_model(1:NY_LOCAL) = log(rho0_prior(25,1:NY_LOCAL)) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = log(c0_prior(25,1:NY_LOCAL))/ scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)
  
  elseif (parametrisation == 35) then
  ! Log-density, wave speed, horizontal wind
  flat_model(1:NY_LOCAL) = log(rho0_prior(25,1:NY_LOCAL)) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = c0_prior(25,1:NY_LOCAL)/ scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)
  
 elseif (parametrisation == 36) then
  ! Log-normalised-density, log-wave speed, horizontal wind
  flat_model(1:NY_LOCAL) = log(rho0_prior(25,1:NY_LOCAL)/density_prior) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = c0_prior(25,1:NY_LOCAL)/ scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)

  elseif (parametrisation == 4) then 
  ! Log-density, log-pressure, horizontal wind
  flat_model(1:NY_LOCAL) = log(rho0_prior(25,1:NY_LOCAL)) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = log(p0_prior(25,1:NY_LOCAL)) / scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)
 
  
  elseif (parametrisation == 5) then
  ! Log-wave speed, log pressure, horizontal wind
  flat_model(1:NY_LOCAL) = log(c0_prior(25,1:NY_LOCAL)) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = log(p0_prior(25,1:NY_LOCAL))/ scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)
  

  else
    print *, "ERROR: Unknown parameterization"
    stop 
    
  endif

endsubroutine priormodel2flatmodel

