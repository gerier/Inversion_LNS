subroutine f(m,fm)
use MPI
use parameters
implicit none
 double precision :: fm, fm_local
 double precision, dimension(Nflat) :: m
 integer :: it
 
 call flatmodel2priormodel(m)

 call forwardproblem(p0_prior, rho0_prior, windx_prior, windy_prior, kappa_unrelaxed_prior,  1, NSTEP, 1) 
 sispressure_prior(:,:) = sispressure(:,:)
 
 ! TODO : revoir a cause du mpi
 fm_local = sum((sispressure_prior(:,:) - sispressure_true(:,:))**2)
 
 call MPI_BARRIER(MPI_COMM_WORLD, code)
 call MPI_ALLREDUCE(fm_local, fm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)

 fm = DELTAT * fm / (2 * norm_pressure_true)

endsubroutine f




subroutine df(m,flat_grad)

use parameters
implicit none 
double precision, dimension(Nflat) :: m,flat_grad

  call flatmodel2priormodel(m)
  call compute_kernel()
  call kernelparam2inversionparam(flat_grad)
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
    
    flat_grad(1:prod_NXNY_LOCAL) = scale_model(1) !* reshape(grad_rho0,[prod_NXNY_LOCAL])
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
    grad_lnrho0(:,:) = p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL) &
                       + K_rho0(1:NX_LOCAL,1:NY_LOCAL) 
    
    flat_grad(:prod_NXNY_LOCAL) = scale_model(1) * reshape(grad_rho0,[prod_NXNY_LOCAL])
    flat_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = scale_model(3) * reshape(grad_c0,[prod_NXNY_LOCAL])
    flat_grad(2*prod_NXNY_LOCAL+1:Nflat)= scale_model(4) * K_windx(1,1:NY_LOCAL)
  
  elseif (parametrisation == 4) then
    ! log density, wind, log pressure
    grad_lnp0(:,:) = p0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_p0(1:NX_LOCAL,1:NY_LOCAL)
    grad_lnrho0(:,:) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * K_rho0(1:NX_LOCAL,1:NY_LOCAL)
    
    flat_grad(:prod_NXNY_LOCAL) = scale_model(1) * reshape(K_rho0(1:NX_LOCAL,1:NY_LOCAL),[prod_NXNY_LOCAL])
    flat_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = scale_model(2) * reshape(K_p0(1:NX_LOCAL,1:NY_LOCAL),[prod_NXNY_LOCAL])
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
double precision, dimension(prod_NXNY_LOCAL) :: flat_rho0, flat_c0, flat_windx, flat_p0
integer :: i

  if (parametrisation == 1) then
   ! density, wind, velocity
   flat_rho0 = scale_model(1)  * flat_model(:prod_NXNY_LOCAL)
   flat_c0 = scale_model(3)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL)
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_rho0, (/NX_LOCAL, NY_LOCAL/))
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_c0, (/NX_LOCAL, NY_LOCAL/))
   do i=1,NX
     windx_prior(i,1:NY_LOCAL) = flat_windx(:)
   enddo
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                     * c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2 / gamma_chimie
   
  else if (parametrisation == 2) then
   ! density, wind, pressure
   flat_rho0 = scale_model(1)  * flat_model(:prod_NXNY_LOCAL)
   flat_p0 = scale_model(2)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL)
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_rho0, (/NX_LOCAL, NY_LOCAL/))
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_p0, (/NX_LOCAL, NY_LOCAL/))
   do i=1,NX
     windx_prior(i,1:NY_LOCAL) = flat_windx(:)
   enddo
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt( p0_prior(1:NX_LOCAL,1:NY_LOCAL) * gamma_chimie &
                             /rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
   
  else if (parametrisation == 3) then
   ! log density, wind, log velocity
   flat_rho0 = exp(scale_model(1)  * flat_model(:prod_NXNY_LOCAL))
   flat_c0 = exp(scale_model(3)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL))
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_rho0, (/NX_LOCAL, NY_LOCAL/))
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_c0, (/NX_LOCAL, NY_LOCAL/))
   do i=1,NX
     windx_prior(i,1:NY_LOCAL) = flat_windx(:)
   enddo
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * &
                                   c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2 / gamma_chimie
   
  else if (parametrisation == 4) then
   ! log density, wind, log pressure
   flat_rho0 = exp(scale_model(1)  * flat_model(:prod_NXNY_LOCAL))
   flat_p0 = exp(scale_model(2)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL))
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_rho0, (/NX_LOCAL, NY_LOCAL/))
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_p0, (/NX_LOCAL, NY_LOCAL/))
   do i=1,NX
     windx_prior(i,1:NY_LOCAL) = flat_windx(:)
   enddo
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt( p0_prior(1:NX_LOCAL,1:NY_LOCAL) * gamma_chimie &
                             /rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
   
  else
    print *, "ERROR: parametrisation unknown"
    stop 
    
  endif
  
  kappa_unrelaxed_prior(1:NX_LOCAL,1:NY_LOCAL) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * &
                                     c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2
                                     
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
  
endsubroutine flatmodel2priormodel



subroutine priormodel2flatmodel(flat_model)

use parameters
implicit none
double precision, dimension(Nflat) :: flat_model
print *, "1#################################### Rank ", rank
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
  flat_model(:prod_NXNY_LOCAL) = reshape( log(rho0_prior(1:NX_LOCAL,1:NY_LOCAL) / scale_model(1)),[prod_NXNY_LOCAL])
  flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = reshape(log(p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ scale_model(2)),[prod_NXNY_LOCAL])
  flat_model(2*prod_NXNY_LOCAL+1:Nflat)= windx_prior(1,1:NY_LOCAL)/ scale_model(4)
  
  else
    print *, "ERROR: parametrisation unknown"
    stop 
    
  endif
  print *, "2#################################### Rank ", rank
endsubroutine priormodel2flatmodel

