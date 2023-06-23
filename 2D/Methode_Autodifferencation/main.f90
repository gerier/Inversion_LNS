

program main

 use parameters
 implicit none
 
 integer i,j, irec
 double precision :: Courant_number
 

integer :: it, it_last_frame


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
        rho_obs(i,j) = density
        p0_obs(i,j) = density * cp_unrelaxed * cp_unrelaxed / gamma_chimie
        kappa_unrelaxed_obs(i,j) = gamma_chimie * p0_obs(i,j)
        v0x_obs(i,j) = 0
        v0y_obs(i,j) = 0
      else
        rho_obs(i,j) = density  * obsfactor_rho
        p0_obs(i,j) = rho_obs(i,j) * (cp_unrelaxed)**2 * obsfactor_c2 / gamma_chimie
        kappa_unrelaxed_obs(i,j) = gamma_chimie * p0_obs(i,j)
        v0x_obs(i,j) = 0
        v0y_obs(i,j) = 0
       endif
       
      ! to have a wind profil  
      if (.False. .and. j > 25 .and. j <= 175) then 
        v0y(i,j) = exp(- ((j-1)*0.1d0 - 10)**2 / 5 ) * 15
      endif
    enddo
  enddo


!!! Compute the a priori model
  do j = 1,NY
    do i = 1,NX
        ! to have adiscontinuity in the background model
        ! attention, pour la comparaison
        !avec le code python, les bornes sont placés en Iobs-1 et non en Iobs
        rho(i,j) = density
        kappa_unrelaxed(i,j) = density*cp_unrelaxed*cp_unrelaxed
        p0(i,j) = density * cp_unrelaxed * cp_unrelaxed / gamma_chimie
        v0x(i,j) = 0
        v0y(i,j) = 0
     enddo
   enddo
        
        
  ! write the background model
  open(unit=201,file='./OUTPUT/true_background_model.dat',status='unknown')
  do j = 1,NY
    write(201,*) (j-1)*0.1d0, rho_obs(1,j), kappa_unrelaxed_obs(1,j), p0_obs(1,j), v0x_obs(1,j), v0y_obs(1,j), gamma_chimie 
  enddo
  close(201)
  
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
  
  call compute_gradient()

 end program  main
