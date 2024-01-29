subroutine f(m,fm)
use MPI
use parameters
implicit none
 double precision :: fm
 double precision, dimension(Nflat) :: m
 integer :: it,irec
 double precision, dimension(NREC) :: fm_local_per_rec
 character(len=100) :: file_name
 integer :: ii,jj
 double precision :: regul_term
 
 fm_local_per_rec(:) = 0.0d0
 
 
 call flatmodel2priormodel(m)
 
  !   write(file_name, "('./TEST/windx_',i6.6,'_',i6.6,'.txt')") count_f,rank
  !  OPEN(UNIT=12, FILE=file_name, ACTION="write")
  !  DO ii=1,NX_LOCAL
  !    WRITE(12,*) (windx_prior(ii,jj), jj=1,NY_LOCAL)
  !  END DO
  !  CLOSE(12)
 
  !  write(file_name, "('./TEST/windy_',i6.6,'_',i6.6,'.txt')") count_f,rank
  !  OPEN(UNIT=12, FILE=file_name, ACTION="write")
  !  DO ii=1,NX_LOCAL
  !    WRITE(12,*) (windy_prior(ii,jj), jj=1,NY_LOCAL)
  !  END DO
  !  CLOSE(12)
    
  !   write(file_name, "('./TEST/p0_',i6.6,'_',i6.6,'.txt')") count_f,rank
  !  OPEN(UNIT=12, FILE=file_name, ACTION="write")
  !  DO ii=1,NX_LOCAL
  !    WRITE(12,*) (p0_prior(ii,jj), jj=1,NY_LOCAL)
  !  END DO
  !  CLOSE(12)
 
  !  write(file_name, "('./TEST/rho0_',i6.6,'_',i6.6,'.txt')") count_f,rank
  !  OPEN(UNIT=12, FILE=file_name, ACTION="write")
  !  DO ii=1,NX_LOCAL
  !    WRITE(12,*) (rho0_prior(ii,jj), jj=1,NY_LOCAL)
  !  END DO
  !  CLOSE(12)
 
 call forwardproblem(p0_prior, rho0_prior, windx_prior, windy_prior, kappa_unrelaxed_prior,  1, NSTEP, 1) 
 sispressure_prior(:,:) = sispressure(:,:)
 
 !
 do irec=1,NREC
   if (i_rank == ix_rec(irec)/NX_LOCAL .and. j_rank == iy_rec(irec)/NY_LOCAL) then
   if ( norm_pressure_true_per_rec(irec) > TINYVAL) then
   fm_local_per_rec(irec) = sum((sispressure_prior(:,irec) - sispressure_true(:,irec))**2) &
                    / norm_pressure_true_per_rec(irec) 
   endif
   endif
enddo

print*, ""
 call MPI_BARRIER(MPI_COMM_WORLD, code)
 call MPI_ALLREDUCE(sum(fm_local_per_rec), fm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)

 fm = DELTAT * fm / 2

 count_f = count_f + 1
 
 if (type_regul_term == 1) then
    
  call MPI_ALLREDUCE( sum(factor_regul_SRdist*(m-m0)**2), regul_term, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)
  fm = fm + regul_weight * regul_term * 0.5d0
 
 endif

endsubroutine f




subroutine df(m,flat_grad)

use parameters
implicit none 
double precision, dimension(Nflat) :: m,flat_grad

  call flatmodel2priormodel(m)
  call compute_kernel()
  call kernelparam2inversionparam(flat_grad)

 if (type_regul_term == 1) then 
   flat_grad = flat_grad + regul_weight * factor_regul_SRdist * (m - m0)
 endif

   
  count_grad = count_grad + 1
endsubroutine df



subroutine kernelparam2inversionparam(flat_grad)

use parameters
implicit none
double precision, dimension (NX_LOCAL,NY_LOCAL) :: grad_c0, grad_rho0, grad_lnc0, grad_lnrho0, grad_lnp0
double precision, dimension(Nflat) :: flat_grad

  if (parametrisation == 1) then
    ! density, wind, velocity
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chimie * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_c0(:,:) = 2 * c0_prior(1:NX_LOCAL,1:NY_LOCAL) * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) / gamma_chimie &
                    * K_p0(1:NX_LOCAL,1:NY_LOCAL)
    grad_rho0(:,:) = p0_prior(1:NX_LOCAL,1:NY_LOCAL) / rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL) &
                     + K_rho0(1:NX_LOCAL,1:NY_LOCAL) 
    
    flat_grad(1:prod_NXNY_LOCAL) = scale_model(1) * reshape(grad_rho0,[prod_NXNY_LOCAL])
    flat_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = scale_model(3) * reshape(grad_c0,[prod_NXNY_LOCAL])
    flat_grad(2*prod_NXNY_LOCAL+1:Nflat)= scale_model(4) * K_windx(1,1:NY_LOCAL)
    
  elseif (parametrisation == 2) then
    ! density, wind, pressure
    flat_grad(:prod_NXNY_LOCAL) = scale_model(1) * reshape(K_rho0(1:NX_LOCAL,1:NY_LOCAL),[prod_NXNY_LOCAL])
    flat_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = scale_model(2) * reshape(K_p0(1:NX_LOCAL,1:NY_LOCAL),[prod_NXNY_LOCAL])
    flat_grad(2*prod_NXNY_LOCAL+1:Nflat)= scale_model(4) * K_windx(1,1:NY_LOCAL)
    
  elseif (parametrisation == 3) then
   ! log density, wind, log velocity
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chimie * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_lnc0(:,:) = 2 * c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2 * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) / gamma_chimie &
                       * K_p0(1:NX_LOCAL,1:NY_LOCAL)
    grad_lnrho0(:,:) =  p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL) &
                       + rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL) 
    
    flat_grad(:prod_NXNY_LOCAL) = scale_model(1) * reshape(grad_lnrho0,[prod_NXNY_LOCAL])
    flat_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = scale_model(3) * reshape(grad_lnc0,[prod_NXNY_LOCAL])
    flat_grad(2*prod_NXNY_LOCAL+1:Nflat)= scale_model(4) * K_windx(1,1:NY_LOCAL)
  
  elseif (parametrisation == 4) then
    ! log density, wind, log pressure
    !if (i_rank == ISOURCE / NX_LOCAL .and. j_rank == JSOURCE / NY_LOCAL) then
    !  K_p0(ISOURCE-offset_i,JSOURCE-offset_j) = 0.0d0 
    !   K_rho0(ISOURCE-offset_i,JSOURCE-offset_j) = 0.0d0
    !endif
    grad_lnp0(:,:) = p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL)
    grad_lnrho0(:,:) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL)
    
    flat_grad(:prod_NXNY_LOCAL) = scale_model(1) * reshape(grad_lnrho0(1:NX_LOCAL,1:NY_LOCAL),[prod_NXNY_LOCAL])
    flat_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = scale_model(2) * reshape(grad_lnp0(1:NX_LOCAL,1:NY_LOCAL),[prod_NXNY_LOCAL])
  
    flat_grad(2*prod_NXNY_LOCAL+1:Nflat)= scale_model(4) * K_windx(1,1:NY_LOCAL)
    
  elseif (parametrisation == 5) then
    ! log celerity, wind, log pressure
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chimie * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_lnc0(:,:) = - 2 * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL)
    grad_lnp0(:,:) = p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL) &
                       + rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL) 
    
    flat_grad(:prod_NXNY_LOCAL) = scale_model(3) * reshape(grad_lnc0,[prod_NXNY_LOCAL])
    flat_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = scale_model(2) * reshape(grad_lnp0,[prod_NXNY_LOCAL])
    flat_grad(2*prod_NXNY_LOCAL+1:Nflat)= scale_model(4) * K_windx(1,1:NY_LOCAL)
     
  else
    print *, "ERROR: parametrisation unknown"
    stop  
  endif

endsubroutine kernelparam2inversionparam


subroutine flatmodel2priormodel(flat_model)

use parameters
implicit none
double precision, dimension(Nflat) :: flat_model
double precision, dimension(prod_NXNY_LOCAL) :: flat_rho0, flat_c0, flat_p0
double precision, dimension(NY_LOCAL) :: flat_windx
double precision, dimension(1:NX_LOCAL, 1:NY_LOCAL) :: test
integer :: i

  if (parametrisation == 1) then
   ! density, wind, velocity
   flat_rho0 = scale_model(1)  * flat_model(:prod_NXNY_LOCAL)
   flat_c0 = scale_model(3)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL)
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_rho0, (/NX_LOCAL, NY_LOCAL/))
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_c0, (/NX_LOCAL, NY_LOCAL/))

   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                     * c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2 / gamma_chimie
   
  else if (parametrisation == 2) then
   ! density, wind, pressure
   flat_rho0 = scale_model(1)  * flat_model(:prod_NXNY_LOCAL)
   flat_p0 = scale_model(2)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL)
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_rho0, (/NX_LOCAL, NY_LOCAL/))
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_p0, (/NX_LOCAL, NY_LOCAL/))
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt( p0_prior(1:NX_LOCAL,1:NY_LOCAL) * gamma_chimie &
                             /rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
   
  else if (parametrisation == 3) then
   ! log density, wind, log velocity
   flat_rho0 = exp(scale_model(1)  * flat_model(:prod_NXNY_LOCAL))
   flat_c0 = exp(scale_model(3)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL))
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_rho0, (/NX_LOCAL, NY_LOCAL/))
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_c0, (/NX_LOCAL, NY_LOCAL/))

   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * &
                                   c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2 / gamma_chimie
   
  else if (parametrisation == 4) then
   ! log density, wind, log pressure
  
   flat_rho0 = exp(scale_model(1)  * flat_model(1:prod_NXNY_LOCAL))
   flat_p0 = exp(scale_model(2)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL))
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_rho0, (/NX_LOCAL, NY_LOCAL/))  
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_p0, (/NX_LOCAL, NY_LOCAL/))

    
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt( p0_prior(1:NX_LOCAL,1:NY_LOCAL) * gamma_chimie &
                             /rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
   
  elseif (parametrisation == 5) then
  ! log celerity, wind, log pressure
  
   flat_c0 = exp(scale_model(3)  * flat_model(1:prod_NXNY_LOCAL))
   flat_p0 = exp(scale_model(2)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL))
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_c0, (/NX_LOCAL, NY_LOCAL/))  
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_p0, (/NX_LOCAL, NY_LOCAL/))
    
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) =  p0_prior(1:NX_LOCAL,1:NY_LOCAL) * gamma_chimie &
                             /c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2
                             
  else
    print *, "ERROR: parametrisation unknown"
    stop 
    
  endif
  
  kappa_unrelaxed_prior(1:NX_LOCAL,1:NY_LOCAL) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * &
                                     c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2
      
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
  
  if (USE_PML_YMIN .and. j_rank == 0) then
    windx_prior(:,-1) = windx_prior(:,1) 
    windx_prior(:,0) = windx_prior(:,1) 
    windy_prior(:,-1) = windy_prior(:,1) 
    windy_prior(:,0) = windy_prior(:,1) 
    rho0_prior(:,-1) = rho0_prior(:,1) 
    rho0_prior(:,0) = rho0_prior(:,1)
    p0_prior(:,-1) = p0_prior(:,1) 
    p0_prior(:,0) = p0_prior(:,1) 
  endif
  
  if (USE_PML_YMAX .and. j_rank == NPROC_Y-1) then
    windx_prior(:,NY_LOCAL+1) = windx_prior(:,NY_LOCAL) 
    windx_prior(:,NY_LOCAL+2) = windx_prior(:,NY_LOCAL) 
    windy_prior(:,NY_LOCAL+1) = windy_prior(:,NY_LOCAL) 
    windy_prior(:,NY_LOCAL+2) = windy_prior(:,NY_LOCAL) 
    rho0_prior(:,NY_LOCAL+1) = rho0_prior(:,NY_LOCAL) 
    rho0_prior(:,NY_LOCAL+2) = rho0_prior(:,NY_LOCAL)
    p0_prior(:,NY_LOCAL+1) = p0_prior(:,NY_LOCAL) 
    p0_prior(:,NY_LOCAL+2) = p0_prior(:,NY_LOCAL) 
  endif
  
  if (USE_PML_YMAX .and. USE_PML_XMAX .and. i_rank == NPROC_X-1 .and. j_rank == NPROC_Y-1) then
    windx_prior(NX_LOCAL+1:NX_LOCAL+2,NY_LOCAL+1:NY_LOCAL+2) = windx_prior(NX_LOCAL,NY_LOCAL) 
    windy_prior(NX_LOCAL+1:NX_LOCAL+2,NY_LOCAL+1:NY_LOCAL+2) = windy_prior(NX_LOCAL,NY_LOCAL) 
    rho0_prior(NX_LOCAL+1:NX_LOCAL+2,NY_LOCAL+1:NY_LOCAL+2) = rho0_prior(NX_LOCAL,NY_LOCAL) 
    p0_prior(NX_LOCAL+1:NX_LOCAL+2,NY_LOCAL+1:NY_LOCAL+2) = p0_prior(NX_LOCAL,NY_LOCAL) 
  endif
  
    if (USE_PML_YMAX .and. USE_PML_XMIN .and. i_rank == 0 .and. j_rank == NPROC_Y-1) then
    windx_prior(-1:0,NY_LOCAL+1:NY_LOCAL+2) = windx_prior(1,NY_LOCAL) 
    windy_prior(-1:0,NY_LOCAL+1:NY_LOCAL+2) = windy_prior(1,NY_LOCAL) 
    rho0_prior(-1:0,NY_LOCAL+1:NY_LOCAL+2) = rho0_prior(1,NY_LOCAL) 
    p0_prior(-1:0,NY_LOCAL+1:NY_LOCAL+2) = p0_prior(1,NY_LOCAL) 
  endif
  
    if (USE_PML_YMIN .and. USE_PML_XMAX .and. i_rank == NPROC_X-1 .and. j_rank == 0) then
    windx_prior(NX_LOCAL+1:NX_LOCAL+2,-1:0) = windx_prior(NX_LOCAL,1) 
    windy_prior(NX_LOCAL+1:NX_LOCAL+2,-1:0) = windy_prior(NX_LOCAL,1) 
    rho0_prior(NX_LOCAL+1:NX_LOCAL+2,-1:0) = rho0_prior(NX_LOCAL,1) 
    p0_prior(NX_LOCAL+1:NX_LOCAL+2,-1:0) = p0_prior(NX_LOCAL,1) 
  endif

    if (USE_PML_YMIN .and. USE_PML_XMIN .and. i_rank == 0 .and. j_rank == 0) then
    windx_prior(-1:0,-1:0) = windx_prior(1,1) 
    windy_prior(-1:0,-1:0) = windy_prior(1,1) 
    rho0_prior(-1:0,-1:0) = rho0_prior(1,1) 
    p0_prior(-1:0,-1:0) = p0_prior(1,1) 
  endif
    
endsubroutine flatmodel2priormodel



subroutine priormodel2flatmodel(flat_model)

use parameters, only : NX_LOCAL, NY_LOCAL, Nflat,prod_NXNY_LOCAL, &
                       scale_model, rho0_prior, c0_prior, p0_prior, windx_prior, parametrisation
implicit none
double precision, dimension(Nflat), intent(out) :: flat_model

  if (parametrisation == 1) then
  ! density, wind, velocity
  flat_model(:prod_NXNY_LOCAL) = reshape(rho0_prior(1:NX_LOCAL,1:NY_LOCAL) / scale_model(1),[prod_NXNY_LOCAL])
  flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = reshape(c0_prior(1:NX_LOCAL,1:NY_LOCAL)/ scale_model(3),[prod_NXNY_LOCAL])
  flat_model(2*prod_NXNY_LOCAL+1:Nflat)= windx_prior(1,1:NY_LOCAL)/ scale_model(4)
    
  elseif (parametrisation == 2) then
  ! density, wind, pressure
  flat_model(:prod_NXNY_LOCAL) = reshape(rho0_prior(1:NX_LOCAL,1:NY_LOCAL) / scale_model(1),[prod_NXNY_LOCAL])
  flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = reshape(p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ scale_model(2),[prod_NXNY_LOCAL])
  flat_model(2*prod_NXNY_LOCAL+1:Nflat)= windx_prior(1,1:NY_LOCAL)/ scale_model(4)
  
  elseif (parametrisation == 3) then
  ! log density, wind, log velocity
  flat_model(:prod_NXNY_LOCAL) = reshape( log(rho0_prior(1:NX_LOCAL,1:NY_LOCAL) / scale_model(1)),[prod_NXNY_LOCAL])
  flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = reshape(log(c0_prior(1:NX_LOCAL,1:NY_LOCAL)/ scale_model(3)),[prod_NXNY_LOCAL])
  flat_model(2*prod_NXNY_LOCAL+1:Nflat)= windx_prior(1,1:NY_LOCAL)/ scale_model(4)
  
  elseif (parametrisation == 4) then
  ! log density, wind, log pressure
  flat_model(1:prod_NXNY_LOCAL) = reshape( log(rho0_prior(1:NX_LOCAL,1:NY_LOCAL) / scale_model(1)),[prod_NXNY_LOCAL])
  flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = reshape(log(p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ scale_model(2)),[prod_NXNY_LOCAL])
  flat_model(2*prod_NXNY_LOCAL+1:Nflat)= windx_prior(1,1:NY_LOCAL)/ scale_model(4)
  
  elseif (parametrisation == 5) then
  ! log celerity, wind, log pressure
  flat_model(1:prod_NXNY_LOCAL) = reshape( log(c0_prior(1:NX_LOCAL,1:NY_LOCAL) / scale_model(3)),[prod_NXNY_LOCAL])
  flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = reshape(log(p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ scale_model(2)),[prod_NXNY_LOCAL])
  flat_model(2*prod_NXNY_LOCAL+1:Nflat)= windx_prior(1,1:NY_LOCAL)/ scale_model(4)
  
  
  else
    print *, "ERROR: parametrisation unknown"
    stop 
    
  endif

endsubroutine priormodel2flatmodel

