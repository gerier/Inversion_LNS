subroutine f(m,fm)
use MPI
use parameters
implicit none
 double precision :: fm
 double precision, dimension(1:Nflat) :: m
 integer :: irec
 double precision, dimension(NREC) :: fm_local_per_rec
 double precision :: regul_term_p0, regul_term_rho0,regul_term_windx
 double precision, dimension(1:Nflat) :: normsq_m, normsq_mm0
 fm_local_per_rec(:) = 0.0d0

 factor_regul_SRdist(:) = 1.0

 call MPI_BARRIER(MPI_COMM_WORLD, code)
 
 call flatmodel2priormodel(m)

 call reset_forward() 
 call forwardproblem(p0_prior, rho0_prior, windx_prior, windy_prior,  1, NSTEP, 1) 
 sispressure_prior(:,:) = sispressure(:,:)
 
 !
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
 
 if (type_regul_term == 1) then
normsq_m(:) = 0.0d0
  normsq_mm0= factor_regul_SRdist*(m-m0)**2
  normsq_mm0(1:2*NY_LOCAL) = normsq_mm0(1:2*NY_LOCAL) / m0(1:2*NY_LOCAL)**2
  
  call MPI_ALLREDUCE( sum(normsq_mm0(1:NY_LOCAL)), regul_term_rho0, 1, MPI_DOUBLE_PRECISION,&
                                                                                    MPI_SUM,  MPI_COMM_WORLD, code)
  call MPI_ALLREDUCE( sum(normsq_mm0(NY_LOCAL+1:2*NY_LOCAL)), regul_term_p0,1,MPI_DOUBLE_PRECISION,&
                                                                                    MPI_SUM, MPI_COMM_WORLD,code)
  call MPI_ALLREDUCE( sum(normsq_mm0(2*NY_LOCAL+1:Nflat)), regul_term_windx, 1, MPI_DOUBLE_PRECISION, &
                                                                                    MPI_SUM,  MPI_COMM_WORLD, code)
                                                                                    
  fx_regul =  regul_weight * 0.5d0 * &
         (regul_term_p0 + regul_term_rho0 + regul_term_windx)
  fm = fm + fx_regul
 
 else
  if (type_regul_term /= 0) then 
   print *, "ERROR: Type of the regularisation term: Unknown"
   stop
   endif
 endif

endsubroutine f


!!!! ! regularisation
! attention au cas m = m0
! attention à la parametrisation choisie

subroutine df(m,flat_grad)

use parameters
implicit none 
double precision, dimension(1:Nflat) :: m,flat_grad,reg_grad
double precision :: ONE_OVER_DX2, ONE_OVER_DY2, ONE_OVER_DXDY


ONE_OVER_DX2 = 1 / DELTAX**4
ONE_OVER_DY2 = 1 / DELTAY**4
ONE_OVER_DXDY = 1 / DELTAX**2 / DELTAY**2

  call MPI_BARRIER(MPI_COMM_WORLD, code)
  call flatmodel2priormodel(m)
  call compute_kernel()
  call kernelparam2inversionparam(flat_grad)

 if (type_regul_term == 1) then
 
    reg_grad(:) = 0.0d0
    !reg_grad(1:NY_LOCAL) = (m(1:NY_LOCAL) - m0(1:NY_LOCAL)) /m0(1:NY_LOCAL)**2
    !reg_grad(NY_LOCAL+1:2*NY_LOCAL) = (m(NY_LOCAL+1:2*NY_LOCAL) - m0(NY_LOCAL+1:2*NY_LOCAL)) / m0(NY_LOCAL+1:2*NY_LOCAL)**2
    reg_grad(2*NY_LOCAL+1:3*NY_LOCAL) = (m(2*NY_LOCAL+1:3*NY_LOCAL) - m0(2*NY_LOCAL+1:3*NY_LOCAL)) 

   flat_grad = flat_grad +  regul_weight * reg_grad * factor_regul_SRdist
   
   
 
 else 
   if (type_regul_term /= 0) then 
   print *, "ERROR: Type of the regularisation term: Unknown"
   stop
   endif
 endif

   
  count_grad = count_grad + 1
endsubroutine df



subroutine kernelparam2inversionparam(flat_grad)

use parameters
implicit none
double precision, dimension (1:NX_LOCAL,1:NY_LOCAL) :: grad_c0, grad_rho0, grad_lnc0, grad_lnrho0, grad_lnp0
double precision, dimension(1:Nflat) :: flat_grad
double precision :: scalar_grad_lnrho0, scalar_grad_lnc0, scalar_grad_lnp0, &
                scalar_grad_rho0, scalar_grad_c0, scalar_grad_p0, &
                scalar_grad_windx
integer :: j

  K_rho0(:,:) = K_rho0(:,:) * DELTAX * DELTAY
  K_p0(:,:) = K_p0(:,:) * DELTAX * DELTAY
  K_windx(:,:) = K_windx(:,:) * DELTAX * DELTAY

  if (parametrisation == 1) then
    ! density, wind, velocity
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                  * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_c0(:,:) = 2 * c0_prior(1:NX_LOCAL,1:NY_LOCAL) * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                    / gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                    * K_p0(1:NX_LOCAL,1:NY_LOCAL)
    grad_rho0(:,:) = p0_prior(1:NX_LOCAL,1:NY_LOCAL) / rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL) &
                     + K_rho0(1:NX_LOCAL,1:NY_LOCAL) 

    do j=1,NY_LOCAL
      ! gather information of different mpi processus 
      call MPI_ALLREDUCE( sum(grad_c0(:,j)), scalar_grad_c0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(grad_rho0(:,j)), scalar_grad_rho0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_windx(:,j)), scalar_grad_windx,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM, row_Comm,ierr)
      ! init
      flat_grad(j) = ZERO
      flat_grad(NY_LOCAL+j) = ZERO
      ! load information in gradient vector
      flat_grad(j) = scale_model(1) * scalar_grad_rho0
      flat_grad(NY_LOCAL + j) = scale_model(2) * scalar_grad_c0
      flat_grad(2*NY_LOCAL + j)= scale_model(3) * scalar_grad_windx
    enddo
    
  elseif (parametrisation == 2) then
    ! density, wind, pressure
    do j=1,NY_LOCAL
      ! gather information of different mpi processus 
      call MPI_ALLREDUCE( sum(K_rho0(:,j)), scalar_grad_rho0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_p0(:,j)), scalar_grad_p0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_windx(:,j)), scalar_grad_windx,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM, row_Comm,ierr)
      ! init
      flat_grad(j) = ZERO
      flat_grad(NY_LOCAL+j) = ZERO
      ! load information in gradient vector
      flat_grad(j) = scale_model(1) * scalar_grad_rho0
      flat_grad(NY_LOCAL + j) = scale_model(2) * scalar_grad_p0
      flat_grad(2*NY_LOCAL + j) = scale_model(3) * scalar_grad_windx
    enddo
    
  elseif (parametrisation == 3) then
   ! log density, wind, log velocity
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                     * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_lnc0(:,:) = 2 * c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2 * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                      / gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                       * K_p0(1:NX_LOCAL,1:NY_LOCAL)
    grad_lnrho0(:,:) =  p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL) &
                       + rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL) 
  
    call MPI_BARRIER(MPI_COMM_WORLD, code)
  
    do j=1,NY_LOCAL
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
  

   ! Interpolation
   if (j_rank == JSOURCE / NY_LOCAL) then
     flat_grad(jsource) = (flat_grad(jsource-1) + flat_grad(jsource+1)) *0.5 
     flat_grad(NY_LOCAL + jsource) = (flat_grad(NY_LOCAL+ jsource-1) + flat_grad(NY_LOCAL+jsource+1)) * 0.5
     flat_grad(2*NY_LOCAL+jsource) = (flat_grad(2*NY_LOCAL + jsource-1) + flat_grad(2*NY_LOCAL + jsource+1)) *0.5
   endif

   call MPI_BARRIER(MPI_COMM_WORLD, code)    

    
  elseif (parametrisation == 35) then
   print *, "chosen param is ", parametrisation
   ! log density, wind, log velocity
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                     * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_c0(:,:) = 2 * c0_prior(1:NX_LOCAL,1:NY_LOCAL) * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                      / gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
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
  

   ! Interpolation
   if (j_rank == JSOURCE / NY_LOCAL) then
     flat_grad(jsource) = (flat_grad(jsource-1) + flat_grad(jsource+1)) *0.5 
     flat_grad(NY_LOCAL + jsource) = (flat_grad(NY_LOCAL+ jsource-1) + flat_grad(NY_LOCAL+jsource+1)) * 0.5
     flat_grad(2*NY_LOCAL+jsource) = (flat_grad(2*NY_LOCAL + jsource-1) + flat_grad(2*NY_LOCAL + jsource+1)) *0.5
   endif

   call MPI_BARRIER(MPI_COMM_WORLD, code)    

  elseif (parametrisation == 4) then
    ! log density, wind, log pressure
    
    grad_lnp0(:,:) = p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL)
    grad_lnrho0(:,:) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL)
    
    do j=1,NY_LOCAL
      ! gather information of different mpi processus 
      call MPI_ALLREDUCE( sum(grad_lnp0(:,j)), scalar_grad_lnp0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(grad_lnrho0(:,j)), scalar_grad_lnrho0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_windx(:,j)), scalar_grad_windx,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM, row_Comm,ierr)
      ! init
      flat_grad(j) = ZERO
      flat_grad(NY_LOCAL+j) = ZERO
      ! load information in gradient vector
      flat_grad(j) = scale_model(1) * scalar_grad_lnrho0
      flat_grad(NY_LOCAL + j) = scale_model(2) * scalar_grad_lnp0
      flat_grad(2*NY_LOCAL + j) = scale_model(3) * scalar_grad_windx
    enddo
    
  elseif (parametrisation == 5) then
    ! log celerity, wind, log pressure
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                          * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_lnc0(:,:) = - 2 * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL)
    grad_lnp0(:,:) = p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL) &
                       + rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL) 
    
    do j=1,NY_LOCAL
      ! log celerity, wind, log pressure
      call MPI_ALLREDUCE( sum(grad_lnp0(:,j)), scalar_grad_lnp0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(grad_lnc0(:,j)), scalar_grad_lnc0,1,MPI_DOUBLE_PRECISION,& 
                        MPI_SUM, row_Comm,ierr)
      call MPI_ALLREDUCE( sum(K_windx(:,j)), scalar_grad_windx,1,MPI_DOUBLE_PRECISION,&
                        MPI_SUM, row_Comm,ierr)
      ! init
      flat_grad(j) = ZERO
      flat_grad(NY_LOCAL+j) = ZERO
      ! load information in gradient vector
      flat_grad(j) = scale_model(1) * scalar_grad_lnc0
      flat_grad(NY_LOCAL + j) = scale_model(2) * scalar_grad_lnp0
      flat_grad(2*NY_LOCAL + j) = scale_model(3) * scalar_grad_windx
    enddo 
  else
    print *, "ERROR: parametrisation unknown"
    stop  
  endif
  

endsubroutine kernelparam2inversionparam


subroutine flatmodel2priormodel(flat_model)

use parameters
implicit none
double precision, dimension(1:Nflat) :: flat_model
double precision, dimension(1:NY_LOCAL) :: flat_rho0, flat_c0, flat_p0,flat_windx
integer :: j

  if (parametrisation == 1) then
   ! density, wind, velocity
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
                      / gamma_chimie(1:NX_LOCAL,1:NY_LOCAL)
   
  else if (parametrisation == 2) then
   ! density, wind, pressure
   flat_rho0 = scale_model(1)  * flat_model(:NY_LOCAL)
   flat_p0 = scale_model(2)  * flat_model(NY_LOCAL+1:2*NY_LOCAL)
   flat_windx = scale_model(3) * flat_model(1+2*NY_LOCAL:)
   
   do j=1,NY_LOCAL
      rho0_prior(:,j) = flat_rho0(j)
      p0_prior(:,j) = flat_p0(j)
      windx_prior(:,j) = flat_windx(j)
   enddo
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt( p0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                          * gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                             /rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
   
  else if (parametrisation == 35) then
   ! log density, wind, log velocity
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
                                   / gamma_chimie(1:NX_LOCAL,1:NY_LOCAL)
  

  else if (parametrisation == 3) then
   ! log density, wind, log velocity
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
                                   / gamma_chimie(1:NX_LOCAL,1:NY_LOCAL)

  else if (parametrisation == 4) then
   ! log density, wind, log pressure
  
   flat_rho0 = exp(scale_model(1)  * flat_model(1:NY_LOCAL))
   flat_p0 = exp(scale_model(2)  * anint(1e16*flat_model(NY_LOCAL+1:2*NY_LOCAL))/1e16)
   flat_windx = scale_model(3) * flat_model(1+2*NY_LOCAL:)

   do j=1,NY_LOCAL
      rho0_prior(:,j) = flat_rho0(j)  
      p0_prior(:,j) = flat_p0(j)
      windx_prior(:,j) = flat_windx(j)
   enddo
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt( p0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                             * gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                             /rho0_prior(1:NX_LOCAL,1:NY_LOCAL))

 
  elseif (parametrisation == 5) then
  ! log celerity, log pressure, wind
  
   flat_c0 = exp(scale_model(1)  * flat_model(1:NY_LOCAL))
   flat_p0 = exp(scale_model(2)  * flat_model(NY_LOCAL+1:2*NY_LOCAL))
   flat_windx = scale_model(3) * flat_model(1+2*NY_LOCAL:)
   
   do j=1,NY_LOCAL
      c0_prior(:,j) = flat_c0(j)  
      p0_prior(:,j) = flat_p0(j)
      windx_prior(:,j) = flat_windx(j)
   enddo
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) =  p0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                             * gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                             /c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2
                             
  else
    print *, "ERROR: parametrisation unknown"
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

use parameters, only : NX_LOCAL, NY_LOCAL, Nflat, &
                       scale_model, rho0_prior, c0_prior, p0_prior, windx_prior, parametrisation
implicit none
double precision, dimension(1:Nflat):: flat_model
 
  if (parametrisation == 1) then
  ! density, wind, velocity
  flat_model(1:NY_LOCAL) = rho0_prior(25,1:NY_LOCAL) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = c0_prior(25,1:NY_LOCAL)/ scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)
    
  elseif (parametrisation == 2) then
  ! density, wind, pressure
  flat_model(1:NY_LOCAL) = rho0_prior(25,1:NY_LOCAL) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = p0_prior(25,1:NY_LOCAL)/ scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)
  
  elseif (parametrisation == 3) then
  ! log density, wind, log velocity
  flat_model(1:NY_LOCAL) = log(rho0_prior(25,1:NY_LOCAL)) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = log(c0_prior(25,1:NY_LOCAL))/ scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)
  
  elseif (parametrisation == 35) then
  ! log density, wind, log velocity
  flat_model(1:NY_LOCAL) = log(rho0_prior(25,1:NY_LOCAL)) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = c0_prior(25,1:NY_LOCAL)/ scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)
  

  elseif (parametrisation == 4) then 
  ! log density, wind, log pressure
  flat_model(1:NY_LOCAL) = log(rho0_prior(25,1:NY_LOCAL)) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = log(p0_prior(25,1:NY_LOCAL)) / scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)
  !flat_model(NY_LOCAL+1:2*NY_LOCAL) = log(anint(1e16*p0_prior(10,1:NY_LOCAL))/1e16) / scale_model(2)

 
  
  elseif (parametrisation == 5) then
  ! log celerity, wind, log pressure
  flat_model(1:NY_LOCAL) = log(c0_prior(25,1:NY_LOCAL)) / scale_model(1)
  flat_model(NY_LOCAL+1:2*NY_LOCAL) = log(p0_prior(25,1:NY_LOCAL))/ scale_model(2)
  flat_model(2*NY_LOCAL+1:Nflat)= windx_prior(25,1:NY_LOCAL)/ scale_model(3)
  

  else
    print *, "ERROR: parametrisation unknown"
    stop 
    
  endif

endsubroutine priormodel2flatmodel

