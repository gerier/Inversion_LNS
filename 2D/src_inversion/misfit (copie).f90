subroutine f(m,fm)
use MPI
use parameters
implicit none
 double precision :: fm
 double precision, dimension(Nflat) :: m
 integer :: irec
 double precision, dimension(NREC) :: fm_local_per_rec
 integer :: i,j
 double precision :: regul_term
 double precision :: regul_term_p0
 double precision :: regul_term_rho0
 double precision :: regul_term_windx
 double precision, dimension(Nflat) :: norm2_m
 fm_local_per_rec(:) = 0.0d0

 call flatmodel2priormodel(m)
  
 call forwardproblem(p0_prior, rho0_prior, windx_prior, windy_prior,  1, NSTEP, 1) 
 sispressure_prior(:,:) = sispressure(:,:)
 
 !
 do irec=1,NREC
   if (i_rank == (ix_rec(irec)-1)/NX_LOCAL .and. j_rank == (iy_rec(irec)-1)/NY_LOCAL) then
   if ( norm_pressure_true_per_rec(irec) > TINYVAL) then
   fm_local_per_rec(irec) = sum((sispressure_prior(:,irec) - sispressure_true(:,irec))**2) &
                    / norm_pressure_true_per_rec(irec) 
   endif
   endif
enddo

 call MPI_BARRIER(MPI_COMM_WORLD, code)
 call MPI_ALLREDUCE(sum(fm_local_per_rec), fm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)

 fm = DELTAT * fm / 2

 count_f = count_f + 1
 
 if (type_regul_term == 1) then
norm2_m(:) = 0.0d0
  !norm2_m(1:2*prod_NXNY_LOCAL) = factor_regul_SRdist*(m - m0)**2 ! (exp(m)-exp(m0))**2
  norm2_m= factor_regul_SRdist&
  *(m-m0)**2 !*(exp(m(1:prod_NXNY_LOCAL))/exp(m0(1:prod_NXNY_LOCAL)) - 1)**2
  
  !norm2_m(1+2*prod_NXNY_LOCAL:Nflat) = factor_regul_SRdist*(m - m0)**2
  
  call MPI_ALLREDUCE( sum(norm2_m(1:prod_NXNY_LOCAL)), regul_term_rho0, 1, MPI_DOUBLE_PRECISION,&
                                                                                    MPI_SUM,  MPI_COMM_WORLD, code)
  call MPI_ALLREDUCE( sum(norm2_m(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL)), regul_term_p0,1,MPI_DOUBLE_PRECISION,&
                                                                                    MPI_SUM, MPI_COMM_WORLD,code)
  call MPI_ALLREDUCE( sum(norm2_m(2*prod_NXNY_LOCAL+1:Nflat)), regul_term_windx, 1, MPI_DOUBLE_PRECISION, &
                                                                                    MPI_SUM,  MPI_COMM_WORLD, code)
                                                                                    
  fm = fm + regul_weight * 0.5d0 * &
         (regul_term_p0 + regul_term_rho0 + regul_term_windx)
 
 elseif (type_regul_term == 2) then
   regul_term_p0 = 0.0d0 
   regul_term_rho0 = 0.0d0 
   regul_term_windx = 0.0d0 
   
   do j=1,NY_LOCAL
     do i=1,NX_LOCAL
       regul_term_rho0 = regul_term_rho0 +  &
        ((rho0_prior(i,j) - rho0_prior(i-1,j))/ DELTAX + (rho0_prior(i,j) - rho0_prior(i,j-1))/DELTAY )**2
                        
       regul_term_p0 = regul_term_p0 + & 
                ((p0_prior(i,j) - p0_prior(i-1,j))/DELTAX + (p0_prior(i,j)  - p0_prior(i,j-1))/DELTAY )**2 
                   
     enddo
     regul_term_windx = regul_term_windx + ((windx_prior(1,j) - windx_prior(1,j-1))/DELTAY )**2
   enddo
   print *, "Regul = ", regul_term
   fm = fm + regul_weight * 0.5d0 * &
          (1e5 * regul_term_p0/norm_p0_prior + regul_term_rho0/norm_rho0_prior + regul_term_windx/norm_windx_prior)

   
 elseif (type_regul_term == 3) then
   regul_term_p0 = 0.0d0 
   regul_term_rho0 = 0.0d0 
   regul_term_windx = 0.0d0 
   
   do j=1,NY_LOCAL
     do i=1,NX_LOCAL
       regul_term_rho0 = regul_term_rho0 + &
        ((-rho0_prior(i+1,j) + 2*rho0_prior(i,j) - rho0_prior(i-1,j))/ DELTAX**2 &
          + (- rho0_prior(i,j+1) + 2*rho0_prior(i,j) - rho0_prior(i,j-1))/DELTAY**2 )**2 
                        
      regul_term_p0 = regul_term_p0 + &
        ((-p0_prior(i+1,j) + 2*rho0_prior(i,j) - p0_prior(i-1,j))/ DELTAX**2 &
          + (- p0_prior(i,j+1) + 2*p0_prior(i,j) - p0_prior(i,j-1))/DELTAY**2 )**2
        
     enddo      
	regul_term_windx = regul_term_windx + &
          (( windx_prior(i,j+1) - 2*windx_prior(i,j) + windx_prior(i,j-1))/DELTAY**2 )**2
   enddo
   
   
   
   fm = fm + regul_weight * 0.5d0 *&
           (regul_term_p0/norm_p0_prior +regul_term_rho0/norm_rho0_prior + regul_term_windx/norm_windx_prior)


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
double precision, dimension(Nflat) :: m,flat_grad,reg_grad
double precision :: ONE_OVER_DX2, ONE_OVER_DY2, ONE_OVER_DXDY
integer :: i,j


ONE_OVER_DX2 = 1 / DELTAX**4
ONE_OVER_DY2 = 1 / DELTAY**4
ONE_OVER_DXDY = 1 / DELTAX**2 / DELTAY**2


  call flatmodel2priormodel(m)
  call compute_kernel()
  call kernelparam2inversionparam(flat_grad)

 if (type_regul_term == 1) then
 
    !call MPI_ALLREDUCE(sum(m0(1:prod_NXNY_LOCAL)**2), norm_rho0_prior,&
    !1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)
   !call MPI_ALLREDUCE(sum(m0(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL)**2), norm_p0_prior,&
   ! 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)
    
    reg_grad(:) = 0.0d0
    reg_grad(1:2*prod_NXNY_LOCAL) = regul_weight * factor_regul_SRdist(1:2*prod_NXNY_LOCAL)  &
                                    * (m(1:2*prod_NXNY_LOCAL) - m0(1:2*prod_NXNY_LOCAL))
                                    !(exp(m(1:2*prod_NXNY_LOCAL)) / exp(m0(1:2*prod_NXNY_LOCAL))-1)
     	
     
   reg_grad(1+2*prod_NXNY_LOCAL:Nflat) = regul_weight * factor_regul_SRdist(1+2*prod_NXNY_LOCAL:Nflat) *&
                  (m(1+2*prod_NXNY_LOCAL:Nflat) - m0(1+2*prod_NXNY_LOCAL:Nflat))
                 
   reg_grad(:prod_NXNY_LOCAL) = reg_grad(:prod_NXNY_LOCAL) !* exp(m(:prod_NXNY_LOCAL)) 
   reg_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = reg_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) !&
                            ! * exp(m(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL))
   reg_grad(2*prod_NXNY_LOCAL+1:Nflat) = reg_grad(2*prod_NXNY_LOCAL+1:Nflat)
   
   
   print *, flat_grad(100),reg_grad(100),flat_grad(prod_NXNY_LOCAL+100),reg_grad(prod_NXNY_LOCAL+100), &
    norm_p0_prior, norm_rho0_prior
    
    
   flat_grad = flat_grad + reg_grad 
   
   
   
   
 elseif (type_regul_term == 2) then
   do j=1,NY_LOCAL
     do i=1,NX_LOCAL
     reg_grad(:) = 0.0d0
     
       reg_grad(i+(j-1)*NX_LOCAL) =  factor_regul_SRdist(i+(j-1)*NX_LOCAL) &
        * rho0_prior(i,j) *  &
        ((-rho0_prior(i+1,j) + 2*rho0_prior(i,j) - rho0_prior(i-1,j))/ DELTAX**2 &
          + (- rho0_prior(i,j+1) + 2*rho0_prior(i,j) - rho0_prior(i,j-1))/DELTAY**2 )
                        
      reg_grad(prod_NXNY_LOCAL+i+(j-1)*NX_LOCAL) = factor_regul_SRdist(prod_NXNY_LOCAL+i+(j-1)*NX_LOCAL) &
       * p0_prior(i,j) *  &
        ((-p0_prior(i+1,j) + 2*p0_prior(i,j) - p0_prior(i-1,j))/ DELTAX**2 &
          + (- p0_prior(i,j+1) + 2*p0_prior(i,j) - p0_prior(i,j-1))/DELTAY**2 )
      
     enddo      
      reg_grad(2*prod_NXNY_LOCAL+j) = factor_regul_SRdist(2*prod_NXNY_LOCAL+j) &
        * (- windx_prior(i,j+1) + 2*windx_prior(i,j) - windx_prior(i,j-1))/DELTAY**2
   enddo
   ! penser au cas ou m = m0
   
   reg_grad(:prod_NXNY_LOCAL) = reg_grad(:prod_NXNY_LOCAL) /  norm_rho0_prior
   reg_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = reg_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) /  norm_p0_prior
   reg_grad(2*prod_NXNY_LOCAL+1:Nflat) = reg_grad(2*prod_NXNY_LOCAL+1:Nflat) /  norm_windx_prior
   
   flat_grad = flat_grad + reg_grad
   
 elseif (type_regul_term == 3) then
    do j=1,NY_LOCAL
     do i=1,NX_LOCAL
     
       reg_grad(:) = 0.0d0
     
       reg_grad(i+(j-1)*NX_LOCAL) = factor_regul_SRdist(i+(j-1)*NX_LOCAL) * ( &
             rho0_prior(i,j) * (ONE_OVER_DX2 + ONE_OVER_DY2 + ONE_OVER_DXDY ) &
             - 8 * (rho0_prior(i+1,j) + rho0_prior(i-1,j)) * (ONE_OVER_DX2 + ONE_OVER_DXDY) &
             - 8 * (rho0_prior(i,j+1) + rho0_prior(i,j-1)) * (ONE_OVER_DY2 + ONE_OVER_DXDY) &
             + 4 * (rho0_prior(i+1,j+1) + rho0_prior(i+1,j-1) + rho0_prior(i-1,j+1) + rho0_prior(i-1,j-1)) * ONE_OVER_DXDY &
             + 2 * (rho0_prior(i+2,j) + rho0_prior(i-2,j)) * ONE_OVER_DX2 &
             + 2 * (rho0_prior(i,j+2) + rho0_prior(i,j-2)) * ONE_OVER_DY2 )
                              
      reg_grad(prod_NXNY_LOCAL+i+(j-1)*NX_LOCAL) = factor_regul_SRdist(prod_NXNY_LOCAL+i+(j-1)*NX_LOCAL) * ( &
              p0_prior(i,j) * (ONE_OVER_DX2 + ONE_OVER_DY2 + ONE_OVER_DXDY ) &
             - 8 * (p0_prior(i+1,j) + p0_prior(i-1,j)) * (ONE_OVER_DX2 + ONE_OVER_DXDY) &
             - 8 * (p0_prior(i,j+1) + p0_prior(i,j-1)) * (ONE_OVER_DY2 + ONE_OVER_DXDY) &
             + 4 * (p0_prior(i+1,j+1) + p0_prior(i+1,j-1) + p0_prior(i-1,j+1) + p0_prior(i-1,j-1)) * ONE_OVER_DXDY &
             + 2 * (p0_prior(i+2,j) + p0_prior(i-2,j)) * ONE_OVER_DX2 &
             + 2 * (p0_prior(i,j+2) + p0_prior(i,j-2)) * ONE_OVER_DY2 )
        
     enddo      
       reg_grad(2*prod_NXNY_LOCAL+j) = reg_grad(2*prod_NXNY_LOCAL+j) + factor_regul_SRdist(2*prod_NXNY_LOCAL+j) * &
         (windx_prior(1,j+2) - 4*windx_prior(1,j+1) + 6*windx_prior(1,j) - 4*windx_prior(1,j-1) + windx_prior(1,j-2))/DELTAY**4
   enddo
   
      
   reg_grad(:prod_NXNY_LOCAL) = reg_grad(:prod_NXNY_LOCAL) /  norm_rho0_prior
   reg_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) = reg_grad(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL) /  norm_p0_prior
   reg_grad(2*prod_NXNY_LOCAL+1:Nflat) = reg_grad(2*prod_NXNY_LOCAL+1:Nflat) /  norm_windx_prior
   
   flat_grad = flat_grad + reg_grad
   
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
double precision, dimension (NX_LOCAL,NY_LOCAL) :: grad_c0, grad_rho0, grad_lnc0, grad_lnrho0, grad_lnp0
double precision, dimension(Nflat) :: flat_grad

  if (parametrisation == 1) then
    ! density, wind, velocity
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                  * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_c0(:,:) = 2 * c0_prior(1:NX_LOCAL,1:NY_LOCAL) * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                    / gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
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
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                     * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
    grad_lnc0(:,:) = 2 * c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2 * rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                      / gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
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
    c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt(gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                          * p0_prior(1:NX_LOCAL,1:NY_LOCAL)/ rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
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

  if (parametrisation == 1) then
   ! density, wind, velocity
   flat_rho0 = scale_model(1)  * flat_model(:prod_NXNY_LOCAL)
   flat_c0 = scale_model(3)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL)
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_rho0, (/NX_LOCAL, NY_LOCAL/))
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_c0, (/NX_LOCAL, NY_LOCAL/))

   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                     * c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2&
                      / gamma_chimie(1:NX_LOCAL,1:NY_LOCAL)
   
  else if (parametrisation == 2) then
   ! density, wind, pressure
   flat_rho0 = scale_model(1)  * flat_model(:prod_NXNY_LOCAL)
   flat_p0 = scale_model(2)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL)
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_rho0, (/NX_LOCAL, NY_LOCAL/))
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_p0, (/NX_LOCAL, NY_LOCAL/))
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt( p0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                          * gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                             /rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
   
  else if (parametrisation == 3) then
   ! log density, wind, log velocity
   flat_rho0 = exp(scale_model(1)  * flat_model(:prod_NXNY_LOCAL))
   flat_c0 = exp(scale_model(3)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL))
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_rho0, (/NX_LOCAL, NY_LOCAL/))
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_c0, (/NX_LOCAL, NY_LOCAL/))

   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = rho0_prior(1:NX_LOCAL,1:NY_LOCAL) * &
                                   c0_prior(1:NX_LOCAL,1:NY_LOCAL)**2 &
                                   / gamma_chimie(1:NX_LOCAL,1:NY_LOCAL)
   
  else if (parametrisation == 4) then
   ! log density, wind, log pressure
  
   flat_rho0 = exp(scale_model(1)  * flat_model(1:prod_NXNY_LOCAL))
   flat_p0 = exp(scale_model(2)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL))
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_rho0, (/NX_LOCAL, NY_LOCAL/))  
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_p0, (/NX_LOCAL, NY_LOCAL/))

    
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = sqrt( p0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                             * gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
                             /rho0_prior(1:NX_LOCAL,1:NY_LOCAL))
   
  elseif (parametrisation == 5) then
  ! log celerity, wind, log pressure
  
   flat_c0 = exp(scale_model(3)  * flat_model(1:prod_NXNY_LOCAL))
   flat_p0 = exp(scale_model(2)  * flat_model(prod_NXNY_LOCAL+1:2*prod_NXNY_LOCAL))
   flat_windx = scale_model(4) * flat_model(1+2*prod_NXNY_LOCAL:)
   
   c0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_c0, (/NX_LOCAL, NY_LOCAL/))  
   p0_prior(1:NX_LOCAL,1:NY_LOCAL) = reshape(flat_p0, (/NX_LOCAL, NY_LOCAL/))
    
   rho0_prior(1:NX_LOCAL,1:NY_LOCAL) =  p0_prior(1:NX_LOCAL,1:NY_LOCAL) &
                             * gamma_chimie(1:NX_LOCAL,1:NY_LOCAL) &
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

