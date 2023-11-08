

subroutine compute_gradient()

  use parameters
  implicit none

  integer :: i,j, ii,jj
  double precision :: grad, delta, obj_function_apriori
  double precision, dimension(0:NX+1,0:NY+1) :: rho0_perturb, p0_perturb,kappa_unrelaxed_perturb, windx_perturb, windy_perturb

  if (which_kernel == 1) then
    delta =  0.0001 * p0_prior(1,1)
  elseif (which_kernel == 2) then
   delta =  0.0001 * rho0_prior(1,1)
  elseif (which_kernel == 3 .or. which_kernel == 4) then
   delta = 0.0001 
  endif
   
  ! Initialisation
  call reset_forward()
  
  ! compute the solution associated to the a priori model
  call forwardproblem(p0_prior,rho0_prior,windx_prior,windy_prior,kappa_unrelaxed_prior, 1,NSTEP, 3) 
  
  ! compute the objective function evaluated at the a priori model
  call compute_objective_function(obj_function_apriori, sispressure, sispressure_true)

  ! initialise variable recording the perturbated model
  rho0_perturb(:,:) = rho0_prior(:,:)
  p0_perturb(:,:) = p0_prior(:,:)
  kappa_unrelaxed_perturb(:,:) = kappa_unrelaxed_prior(:,:)
  windx_perturb(:,:) = windx_prior(:,:) 
  windy_perturb(:,:) = windy_prior(:,:)  


  ! make perturbations on the a priori model in (i,j) 
  ! compute the wavefront for this new model and compute the associated objective function.
  ! then compute the derivative for the parameter p(i,j) 
  do j= 70,130
   do i=80,150
    
    ! Initialisation
    call reset_forward()
    
    ! add the perturbation
    if (which_kernel == 1) then
      p0_perturb(:,:) = p0_prior(:,:) 
      p0_perturb(i,j) = p0_prior(i,j) + delta
      kappa_unrelaxed_perturb      = kappa_unrelaxed_prior(:,:)
      kappa_unrelaxed_perturb(i,j) = p0_perturb(i,j) * gamma_chimie
    elseif (which_kernel == 2) then
      rho0_perturb(:,:) = rho0_prior(:,:) 
      rho0_perturb(i,j) = rho0_prior(i,j) + delta
    elseif (which_kernel == 3) then
      windx_perturb(:,:) = windx_prior(:,:)
      windx_perturb(i,j) = windx_prior(i,j) + delta
     elseif (which_kernel == 4) then
      windy_perturb(:,:) = windy_prior(:,:)
      windy_perturb(i,j) = windy_prior(i,j) + delta
    endif
 
    ! compute the waveform with the perturbated model 
    call forwardproblem(p0_perturb,rho0_perturb,windx_perturb,windy_perturb,kappa_unrelaxed_perturb, 1,NSTEP,3)

    ! compute the objective function evaluated at the perturbation model
    call compute_objective_function(grad, sispressure, sispressure_true)
    K(i,j) = (grad  -  obj_function_apriori) / delta
     
    ! intermediate saving
    if ((mod(j,10) == 0) .or. i == 120) then
     call create_color_image(K,NX,NY,0,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
             NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,6)
     open(unit=12, file="OUTPUT/Kernel.txt", action="write")
     do ii=1,NX
       write(12,*) (K(ii,jj), jj=1,NY)
     enddo
     close(12)  
    endif
    
  enddo
 enddo
    
 ! save the kernel information and associated pictures
 call create_color_image(K,NX,NY,0,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,6)

 open(unit=12, file="OUTPUT/Kernel.txt", action="write")
 do ii=1,NX
  write(12,*) (K(ii,jj), jj=1,NY)
 enddo
 close(12)

endsubroutine compute_gradient



subroutine compute_objective_function(grad, signal1, signal2)

 use parameters
  implicit none

  integer :: it,irec
  double precision :: diff,grad
  double precision, dimension(NSTEP,NREC) :: signal1, signal2

  !!!
  ! Compute the objective function : the norm L2 of the difference between signal1 and signal 2
  ! Signal 2 is the list of signals observed on each receiver. 
  ! Signal 1 is the list of signals, solution of the direct problem computing with a model m, at each receiver. 
  !!!
  
  grad = 0.d0
  do it=1,NSTEP
   do irec=1,nrec
   	diff = signal1(it,irec) - signal2(it,irec)
   	if (abs(diff) < TINYVAL) then
   	  diff = 0.0d0
   	endif 
     grad = grad + (signal1(it,irec) - signal2(it,irec))**2 
   enddo
  enddo

  grad = grad * 0.5d0 * DELTAT / norm_pressure_true

endsubroutine compute_objective_function
