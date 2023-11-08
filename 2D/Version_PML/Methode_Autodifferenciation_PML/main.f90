

program main

 use parameters
 implicit none
 
 integer i,j, irec, ii, jj
 double precision :: Courant_number
 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                    INITIALISATION                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  print *
  print *,'2D acoustic finite-difference code in density, velocity and pressure'
  print *

! display size of the model
  print *
  print *,'NX = ',NX
  print *,'NY = ',NY
  print *
  print *,'size of the model along X = ',(NX - 1) * DELTAX
  print *,'size of the model along Y = ',(NY - 1) * DELTAY
  print *
  print *,'Total number of grid points = ',NX * NY
  print *




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                       PML INIT                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call computePML()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                Background domain                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! compute the Lame parameter and density


!!! Compute the true model
  do j = 1,NY
    do i = 1,NX
      if ((i < IObs_start) .or. (i >= IObs_end) .or. (j < JObs_start) .or. (j >= JObs_end)) then 
        ! to have adiscontinuity in the background model
        ! attention, pour la comparaison
        !avec le code python, les bornes sont placés en Iobs-1 et non en Iobs
        rho0_true(i,j)            = density
        p0_true(i,j)              = density * cp_unrelaxed * cp_unrelaxed / gamma_chimie
        kappa_unrelaxed_true(i,j) = gamma_chimie * p0_true(i,j)
        windx_true(i,j)           = 0
        windy_true(i,j)           = 0
      else
        rho0_true(i,j)            = density  * obstacle_factor_rho
        p0_true(i,j)              = rho0_true(i,j) * (cp_unrelaxed)**2 * obstacle_factor_c2 / gamma_chimie
        kappa_unrelaxed_true(i,j) = gamma_chimie * p0_true(i,j)
        windx_true(i,j)           = 0
        windy_true(i,j)           = 0
       endif
       
      ! to have a wind profil  
      if (add_wind_profile .and. (j > 25 .and. j <= 175)) then 
        windy_true(i,j) = exp(- ((j-1)*0.1d0 - 10)**2 / 5 ) * 150
      endif
      
    enddo
  enddo


!!! Compute the a priori model
  do j = 1,NY
    do i = 1,NX
        ! to have adiscontinuity in the background model
        ! attention, pour la comparaison
        !avec le code python, les bornes sont placés en Iobs-1 et non en Iobs
        rho0_prior(i,j)            = density
        kappa_unrelaxed_prior(i,j) = density * cp_unrelaxed * cp_unrelaxed
        p0_prior(i,j)              = density * cp_unrelaxed * cp_unrelaxed / gamma_chimie
        windx_prior(i,j)           = 0
        windy_prior(i,j)           = 0
    enddo
  enddo
        
  ! write the background model
  call write_background(rho0_true, kappa_unrelaxed_true, p0_true, windx_true, windy_true, gamma_chimie)
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          source / receivers INFORMATION          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! print position of the source
  print *,'Position of the source:'
  print *
  print *,'x = ',xsource
  print *,'y = ',ysource
  print *

! define location of receivers
  print *,'There are ',nrec,' receivers'
  print *
  if (NREC > 1) then
! this is to avoid a warning with GNU gfortran at compile time about division by zero when NREC = 1
    myNREC = NREC
    xspacerec = (xfin-xdeb) / dble(myNREC-1)
    yspacerec = (yfin-ydeb) / dble(myNREC-1)
  else
    xspacerec = 0.d0
    yspacerec = 0.d0
  endif
  do irec=1,nrec
    xrec(irec) = xdeb + dble(irec-1)*xspacerec
    yrec(irec) = ydeb + dble(irec-1)*yspacerec
  enddo

! find closest grid point for each receiver
  do irec=1,nrec
    dist = HUGEVAL
    do j = 1,NY
    do i = 1,NX
      distval = sqrt((DELTAX*dble(i-1) - xrec(irec))**2 + (DELTAY*dble(j-1) - yrec(irec))**2)
      if (distval < dist) then
        dist = distval
        ix_rec(irec) = i
        iy_rec(irec) = j
      endif
    enddo
    enddo
    print *,'receiver ',irec,' x_target,y_target = ',xrec(irec),yrec(irec)
    print *,'closest grid point found at distance ',dist,' in i,j = ',ix_rec(irec),iy_rec(irec)
    print *
  enddo
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                   CFL Condition                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! check the Courant stability condition for the explicit time scheme
! R. Courant, K. O. Friedrichs and H. Lewy (1928)
! For this O(2,4) scheme, when DELTAX == DELTAY the Courant number is given by Levander, Geophysics, vol. 53(11), p. 1427,
! equation (7) and is equal to 0.606 (it is thus smaller than that of the O(2,2) scheme, which is 1/sqrt(2) = 0.707,
! i.e. when switching to a fourth-order spatial scheme one needs a time step that is about 0.707 / 0.606 = 1.167 times smaller.
  if (DELTAX == DELTAY) then
    Courant_number = cp_unrelaxed * DELTAT / DELTAX
    print *,'Courant number is ',Courant_number
    print *,' (the maximum possible value is 0.606; in practice for accuracy reasons a value not larger than 0.30 is recommended)'
    print *
    if (Courant_number > 0.606) stop 'time step is too large, simulation will be unstable'
  endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!             Retrieve the observations                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  ! INITIALISATION
  call reset_forward()
  save_sismos = .True.
  
  call forwardproblem(p0_true, rho0_true, windx_true, windy_true, kappa_unrelaxed_true,  1, NSTEP, 3) 
  sispressure_true(:,:) = sispressure(:,:)
  norm_pressure_true = maxval(abs(sispressure_true))**2
  
  OPEN(UNIT=12, FILE="OUTPUT/p_true.txt", ACTION="write")
  DO ii=1,NX
    WRITE(12,*) (pressure(ii,jj), jj=1,NY)
  END DO
  CLOSE(12)

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                  Compute gradient                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  call compute_gradient()
    
 end program  main
