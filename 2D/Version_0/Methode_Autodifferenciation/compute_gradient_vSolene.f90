

subroutine compute_gradient()

  use parameters
  implicit none

  integer :: i,j, ii,jj
  double precision :: grad, delta_p0, delta_rho, delta_v0x, obj_function_apriori
  double precision, dimension(0:NX+1,0:NY+1) :: rho_test, p0_test,kappa_unrelaxed_test, v0x_test, v0y_test

   delta_p0 =  0.01 * p0(1,1)! v0x(1,1)
   delta_rho =  0.01 * rho(1,1)! v0x(1,1)
   delta_v0x = 0.01 

   print *, "DELTA = ", delta_v0x

  ! INITIALISATION
  sisvx(:,:) = ZERO
  sisvy(:,:) = ZERO
  sispressure(:,:) = ZERO
  sisrhop(:,:) = ZERO
  
  pressure(:,:)    = ZERO
  rhop(:,:)        = ZERO
  vx(:,:)          = ZERO
  vy(:,:)          = ZERO
  
  ! Compute the reality
  call forwardproblem(pressure, rhop, vx, vy, p0_obs, rho_obs, v0x_obs, v0y_obs, kappa_unrelaxed_obs,  1, NSTEP) 
  obspressure(:,:) = sispressure(:,:)
  norm_obs = maxval(abs(obspressure))**2 

  ! INITIALISATION
  sisvx(:,:) = ZERO
  sisvy(:,:) = ZERO
  sispressure(:,:) = ZERO
  sisrhop(:,:) = ZERO
  
  pressure(:,:)    = ZERO
  rhop(:,:)        = ZERO
  vx(:,:)          = ZERO
  vy(:,:)          = ZERO
  ! compute the a priori
  call forwardproblem(pressure, rhop, vx, vy, p0, rho, v0x, v0y, kappa_unrelaxed,  1, NSTEP) 

  ! compute the objective function for the a priori model
  call compute_objective_function(obj_function_apriori, sispressure, obspressure)

  ! make perturbations on the a priori model in (i,j), compute the wavefront for this new model and compute the associated objective function.
  ! then compute the derivative for the parameter p(i,j) 
  do j= 70,130
   do i=80,150
    
    ! INITIALISATION
    sisvx(:,:) = ZERO
    sisvy(:,:) = ZERO
    sispressure(:,:) = ZERO
    sisrhop(:,:) = ZERO
  
    pressure(:,:)    = ZERO
    rhop(:,:)        = ZERO
    vx(:,:)          = ZERO
    vy(:,:)          = ZERO
  
  
    ! pressure test
    !p0_test(:,:) = p0(:,:) 
    !p0_test(i,j) = p0(i,j) + delta_p0
    !kappa_unrelaxed_test = kappa_unrelaxed(:,:)
    !kappa_unrelaxed_test(i,j) =p0_test(i,j) * gamma_chimie
    ! density test
    !rho_test(:,:) = rho(:,:) 
    !rho_test(i,j) = rho(i,j) + delta_rho
    ! wind test
    !v0x_test(:,:) = v0x(:,:)
    !v0x_test(i,j) = v0x(i,j) + delta_v0x
    v0y_test(:,:) = v0y(:,:)
    v0y_test(i,j) = v0y(i,j) + delta_v0x
 
 
    call forwardproblem(pressure, rhop, vx, vy, p0, rho, v0x, v0y_test, kappa_unrelaxed,  1, NSTEP)
    
    call compute_objective_function(grad, sispressure, obspressure)
    Kvx(i,j) = (grad  -  obj_function_apriori) / delta_v0x
    print *, grad
     
    if ((mod(j,10) == 0) .or. i == 120) then
     call create_color_image(Kvx,NX,NY,0,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              0,.FALSE.,.FALSE.,.FALSE.,.FALSE.,6)
    OPEN(UNIT=12, FILE="OUTPUT/Kernel.txt", ACTION="write")
  DO ii=1,NX
    WRITE(12,*) (Kvx(ii,jj), jj=1,NY)
  END DO
	close(12)  
    endif
  enddo
 enddo
    
 call create_color_image(Kvx,NX,NY,0,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              0,.FALSE.,.FALSE.,.FALSE.,.FALSE.,6)

       OPEN(UNIT=12, FILE="aoutput.txt", ACTION="write")
	  DO ii=1,NX
	    WRITE(12,*) (Kvx(ii,jj), jj=1,NY)
	  END DO
       close(12)

endsubroutine compute_gradient



subroutine compute_objective_function(grad, signal1, signal2)

 use parameters
  implicit none

  integer :: it,irec
  double precision :: diff,grad
  double precision, dimension(NSTEP,NREC) :: signal1, signal2

  grad = 0.d0
  do it=1,NSTEP
   do irec=1,nrec
   	diff = signal1(it,irec) - signal2(it,irec)
   	if (abs(diff) < 1e-16) then
   	  diff = 0.0d0
   	endif 
     grad = grad + (signal1(it,irec) - signal2(it,irec))**2 
   enddo
  enddo


     grad = grad * 0.5d0 * DELTAT / norm_obs
endsubroutine compute_objective_function
