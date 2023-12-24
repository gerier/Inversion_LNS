
subroutine init_backgrounds()

use parameters

integer :: i,j


!!! Compute the true model
  do j = 1,NY_LOCAL
    do i = 1,NX_LOCAL
    
      i_global = i + offset_i
      j_global = j + offset_j 
      
      if ((i_global < IObs_start) .or. (i_global >= IObs_end) .or. (j_global < JObs_start) .or. (j_global >= JObs_end)) then 
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
      if (add_wind_profile .and. (j_global > ymin_wind .and. j_global <= ymax_wind)) then 
        windy_true(i,j) = exp(- ((j_global-1)*0.1d0 - mean_gauss_wind)**2 / sigma2_gauss_wind ) * max_wind_factor
      endif
      
    enddo
  enddo

  call send_receive_rightleft(windx_true)
  call send_receive_rightleft(windy_true)
  call send_receive_rightleft(rho0_true)
  call send_receive_rightleft(p0_true)
     
  call send_receive_topbottom(windx_true)
  call send_receive_topbottom(windy_true)
  call send_receive_topbottom(rho0_true)
  call send_receive_topbottom(p0_true)
     
  call send_receive_corners(windx_true) 
  call send_receive_corners(windy_true)
  call send_receive_corners(rho0_true)

  ! write the background model
  call write_background(rho0_true, kappa_unrelaxed_true, p0_true, windx_true, windy_true, gamma_chimie)


!!! Compute the a priori model
  do j = 1,NY_LOCAL
    do i = 1,NX_LOCAL
        ! to have a discontinuity in the background model
        ! attention, pour la comparaison
        !avec le code python, les bornes sont placés en Iobs-1 et non en Iobs
        rho0_prior(i,j)            = density
        kappa_unrelaxed_prior(i,j) = density * cp_unrelaxed * cp_unrelaxed
        p0_prior(i,j)              = density * cp_unrelaxed * cp_unrelaxed / gamma_chimie
        windx_prior(i,j)           = 0
        windy_prior(i,j)           = 0
    enddo
  enddo

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
       
    ! write the background model
  call write_background(rho0_prior, kappa_unrelaxed_prior, p0_prior, windx_prior, windy_prior, gamma_chimie)


endsubroutine init_backgrounds
 


subroutine init_source_recvrs() 

use parameters

  if (rank == 0) then
    ! print position of the source
    print *,'Position of the source:'
    print *
    print *,'x = ',xsource
    print *,'y = ',ysource
    print *

    ! define location of receivers
    print *,'There are ',nrec,' receivers'
    print *
  endif
  
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
    
    if (rank ==0) then
    print *,'receiver ',irec,' x_target,y_target = ',xrec(irec),yrec(irec)
    print *,'closest grid point found at distance ',dist,' in i,j = ',ix_rec(irec),iy_rec(irec)
    print *
    endif
  enddo
  
  
endsubroutine init_source_recvrs
