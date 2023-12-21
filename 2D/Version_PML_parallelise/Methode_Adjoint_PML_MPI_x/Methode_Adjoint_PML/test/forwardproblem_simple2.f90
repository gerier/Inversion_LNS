subroutine forwardproblem(CI_pressure,CI_rhop,CI_vx,CI_vy,  p0,rho0,windx,windy,kappa_unrelaxed, it_start, it_end, field_number)
 !

use parameters , only : pressure, rhop, vx, vy, &
                       sispressure, sisrhop, sisvx, sisvy, NREC,IT_DISPLAY, ix_rec, iy_rec, &
                       NX, NY, NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX, &
                       NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,ONE_OVER_SIX_DELTAX,ONE_OVER_SIX_DELTAY, &
                       DELTAX, DELTAY, DELTAT, NSTEP, t,   &
                       a, f0, t0, pi, factor, ISOURCE, JSOURCE, source_term, &
                       distance2, factor_ssf, SSF_Sigma, &
                       ZERO, gamma_chimie, stability_threshold, save_sismos,&
                       a_x, a_y, a_x_half, a_y_half, &
                       b_x, b_y, b_x_half, b_y_half, &
                       c_x, c_y, c_x_half, c_y_half, &
                       one_over_K_x, one_over_K_x_half, one_over_K_y, one_over_K_y_half
implicit none



! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here

  integer :: it_start, it_end

  integer :: field_number

  double precision, dimension(0:NX+1,0:NY+1) :: &
    CI_pressure, CI_rhop, CI_vx, CI_vy, &
    p0, rho0, windx, windy, kappa_unrelaxed

  double precision :: &
      value_dvx_dx, &
      value_dvy_dy, &
      value_dpressure_dx, &
      value_dpressure_dy, &
      value_dp0_dx, &
      value_dp0_dy, &
      value_drho0_dx, &
      value_drho0_dy, &
      value_drhop_dx, &
      value_drhop_dy, &
      value_dwindx_dx, &
      value_dwindy_dy, &
      value_vdrho0,    &
      value_v0drhop, &
      value_vdp0 , &
      value_v0dp

  double precision :: &
    value_rho0dv, &
    value_dvy_dx, &
    value_dvx_dy, &
    value_dwindx_dx_prec, &
    value_dwindx_dx_next, &
    value_dwindy_dy_prec, &
    value_dwindy_dy_next, &
    value_drho0_dy_prec, &
    value_drho0_dy_next, &
    value_drho0_dx_prec, &
    value_drho0_dx_next

  double precision :: &
      vx_half_x, &
      vy_half_y,&
      windx_half_x, &
      windy_half_y, &
      windx_half_x_half_y, &
      windy_half_x_half_y, &
      vy_half_x_half_y, &
      vx_half_x_half_y

   double precision :: &
      value_vdwindy, &
      value_vdwindx, &
      value_v0dvy, &
      value_v0dvx, &
      value_v0dwindy, &
      value_v0dwindx, &
      value_dwindy_dx, &
      value_dwindx_dy

  double precision :: &
      value_dwindx_dy_prec,&
      value_dwindy_dx_prec,&
      value_dwindx_dy_next,&
      value_dwindy_dx_next,&
      value_windx_dwindy_dx,&
      value_windy_dwindx_dy,&
      value_windx_dwindx_dx,&
      value_windy_dwindy_dy

  ! memory variables 
  double precision, dimension(NX,NY) ::        &
  ! for equation on rhop and pressure
      eq1_memory_dp0_dx, eq1_memory_dp0_dy,            &
      eq1_memory_drho0_dx, eq1_memory_drho0_dy,        &
      eq1_memory_dpressure_dx, eq1_memory_dpressure_dy,&
      eq1_memory_drhop_dx, eq1_memory_drhop_dy,        &
      eq1_memory_dvx_dx, eq1_memory_dvy_dy,            &
      eq1_memory_dwindx_dx, eq1_memory_dwindy_dy,      &
  ! for equation on vx
      eq2_memory_dpressure_dx,                         &
      eq2_memory_drho0_dx, eq2_memory_drho0_dy,        &
      eq2_memory_dvx_dx, eq2_memory_dvx_dy,            &
      eq2_memory_dwindx_dx, eq2_memory_dwindx_dy,      &
      eq2_memory_dwindy_dy,                            &
  ! for equation on vy
      eq3_memory_dpressure_dy,                         &
      eq3_memory_drho0_dy, eq3_memory_drho0_dx,        &
      eq3_memory_dvy_dy, eq3_memory_dvy_dx,            &
      eq3_memory_dwindy_dy, eq3_memory_dwindy_dx,      &
      eq3_memory_dwindx_dx
     
  
        
        

  double precision  :: u_mm, u_m, u_p, u_pp


      
        
  integer :: Ip1,Im1, Ip2, Im2, Jp1,Jm1, Jp2, Jm2


  integer :: i,j,it,irec

  double precision :: velocnorm,pressurenorm

  double precision, dimension(0:NX+1,0:NY+1) :: rhop_old, p_old, vx_old, vy_old

 ! to interpolate material parameters or velocity at the right location in the staggered grid cell
  double precision :: rho0_half_x,rho0_half_y,rhop_half_x,rhop_half_y

!---
!--- program starts here
!---


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!           Variable initialisation                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! suppress old files (can be commented out if "call system" is missing in your compiler)
  call system('rm -f density*.dat pressure_*.dat Vx_*.dat Vy_*.dat image*.pnm image*.gif source*.dat')

! initialize arrays
  vx(:,:) = CI_vx              ! ZERO
  vy(:,:) = CI_vy              ! ZERO
  pressure(:,:) = CI_pressure  ! ZERO
  rhop(:,:) = CI_rhop          ! ZERO

! PML
  eq1_memory_dp0_dx(:,:)       = ZERO
  eq1_memory_dp0_dy(:,:)       = ZERO
  eq1_memory_drho0_dx(:,:)     = ZERO
  eq1_memory_drho0_dy(:,:)     = ZERO
  eq1_memory_dpressure_dx(:,:) = ZERO
  eq1_memory_dpressure_dy(:,:) = ZERO
  eq1_memory_drhop_dx(:,:)     = ZERO
  eq1_memory_drhop_dy(:,:)     = ZERO
  eq1_memory_dvx_dx(:,:)       = ZERO
  eq1_memory_dvy_dy(:,:)       = ZERO
  eq1_memory_dwindx_dx(:,:)    = ZERO
  eq1_memory_dwindy_dy(:,:)    = ZERO
  ! for equation on vx
  eq2_memory_dpressure_dx(:,:) = ZERO
  eq2_memory_drho0_dx(:,:)     = ZERO
  eq2_memory_drho0_dy(:,:)     = ZERO
  eq2_memory_dvx_dx(:,:)       = ZERO
  eq2_memory_dvx_dy(:,:)       = ZERO
  eq2_memory_dwindx_dx(:,:)    = ZERO
  eq2_memory_dwindx_dy(:,:)    = ZERO
  eq2_memory_dwindy_dy(:,:)    = ZERO
  ! for equation on vy
  eq3_memory_dpressure_dy(:,:) = ZERO
  eq3_memory_drho0_dy(:,:)     = ZERO
  eq3_memory_drho0_dx(:,:)     = ZERO
  eq3_memory_dvy_dy(:,:)       = ZERO
  eq3_memory_dvy_dx(:,:)       = ZERO
  eq3_memory_dwindy_dy(:,:)    = ZERO
  eq3_memory_dwindy_dx(:,:)    = ZERO
  eq3_memory_dwindx_dx(:,:)    = ZERO
      

!---
!---  beginning of time loop
!---

  do it = it_start, it_end

!-----------------------------------------------------------------------
! compute pressure and update memory variables for C-PML
! also update memory variables for viscoacoustic attenuation if needed
!-----------------------------------------------------------------------

  p_old(:,:) = pressure(:,:)
  rhop_old(:,:) = rhop(:,:)

  do j = 2,NY
    do i = 1,NX-1

       Im1 = i-1
       Im2 = i-2
       Ip1 = i+1
       Ip2 = i+2
       Jm1 = j-1
       Jm2 = j-2
       Jp1 = j+1
       Jp2 = j+2
      

    if (j == 1) then
      Jm1 = 1
      Jm2 = 1
    elseif (j == 2) then
      Jm2 = 1
    else if (j == NY) then
      Jp1 = NY
      Jp2 = NY
    else if (j == NY-1) then
      Jp2 = NY
    endif
    if (i == 1) then
      Im1 = 1
      Im2 = 1
    elseif (i == 2) then
      Im2 = 1
    else if (i == NX) then
       Ip1 = NX
       Ip2 = NX
    else if (i == NX-1) then
      Ip2 = NX 
    endif


        ! derivative computations
	u_mm = vx(Im1,j); u_m = vx(i,j); u_p = vx(Ip1,j); u_pp = vx(Ip2,j)
        call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
        
        u_mm = vy(i,Jm2); u_m = vy(i,Jm1); u_p = vy(i,j); u_pp = vy(i,Jp1)
        call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      
    
        eq1_memory_dvx_dx(i,j) = b_x_half(i) * eq1_memory_dvx_dx(i,j) + a_x_half(i) * value_dvx_dx
        eq1_memory_dvy_dy(i,j) = b_y(j) * eq1_memory_dvy_dy(i,j) + a_y(j) * value_dvy_dy 
                
        ! 
        value_dvx_dx = value_dvx_dx * one_over_K_x_half(i) + eq1_memory_dvx_dx(i,j)
        value_dvy_dy = value_dvy_dy * one_over_K_y(j) + eq1_memory_dvy_dy(i,j) 
      
        
    
        
        ! updateT
        pressure(i,j) = pressure(i,j) - (kappa_unrelaxed(i,j) * (value_dvx_dx + value_dvy_dy)) * DELTAT
     

      enddo
    enddo


  

! add source
! add the source (pressure located at a given grid point)
  a = pi*pi*f0*f0
  t = dble(it-1)*DELTAT

! Gaussian
!  source_term = factor * exp(- 4 * a*(t-t0)*(t-t0))  ! / (2.d0 * a)
! derivative of guassian
  source_term = -8 * a* (t-t0) *  factor * exp(- 4 * a*(t-t0)*(t-t0))    ! / (2.d0 * a)
  !source_term = (8 * a* (t-t0))**2  *  factor * exp(- 4 * a*(t-t0)*(t-t0))  ! / (2.d0 * a)
  !source_term = factor * exp(- 4 * a*(t-t0)*(t-t0))  ! / (2.d0 * a)

! define location of the source
  !i = ISOURCE
  j = JSOURCE
  i = ISOURCE
! the pressure source is added to d(pressure)/dt in this split pressure / velocity scheme
! and that is why we need to select the first derivative of a Gaussian as a source time wavelet
! above instead of a Ricker (i.e. a second derivative) added to d2(pressure)/dt2
! as in the unsplit equation written in pressure only.
! Since the formula is d(pressure)/dt = (pressure_new - pressure_old) / DELTAT = pressure_source_term
! we also need to multiply by DELTAT here to avoid having an amplitude of the seismogram
! that varies when one changes the time step, i.e. we write:
! pressure_new = pressure_old + pressure_source_term * DELTAT at the source grid point
    !do j = 1,NY
      !rho0_half_y = 0.5d0 * (rho0(i,j+1) + rho0(i,j))
      !vy(i,j) = vy(i,j) + source_term * DELTAT / rho0_half_y * 0.5d0
      !distance2 = ((i - Isource) * DELTAX)**2 + ((j - Jsource) * DELTAY)**2
      factor_ssf = 1 !exp( - distance2 / SSF_Sigma**2 )

      pressure(i,j) = pressure(i,j) + source_term * factor_ssf * DELTAT
      !rhop(i,j) = rhop(i,j) + source_term * factor_ssf * DELTAT * (rho0(i,j)/gamma_chimie / p0(i,j))
    !enddo

  ! write the source
  open(unit=211,file='./OUTPUT/source_time_function_model.dat',status='unknown',position="append")
    write(211,*) (t-t0), source_term
  close(211)




!--------------------------------------------------------
! compute velocity and update memory variables for C-PML
!--------------------------------------------------------

  vx_old(:,:) = vx(:,:)
  vy_old(:,:) = vy(:,:)


  do j = 2,NY
    do i = 2,NX

       Im1 = i-1
       Im2 = i-2
       Ip1 = i+1
       Ip2 = i+2
       Jm1 = j-1
       Jm2 = j-2
       Jp1 = j+1
       Jp2 = j+2
      
    if (j == 1) then
      Jm1 = 1
      Jm2 = 1
    elseif (j == 2) then
      Jm2 = 1
    else if (j == NY) then
      Jp1 = NY
      Jp2 = NY
    else if (j == NY-1) then
      Jp2 = NY
    endif
    if (i == 1) then
      Im1 = 1
      Im2 = 1
    elseif (i == 2) then
      Im2 = 1
    else if (i == NX) then
       Ip1 = NX
       Ip2 = NX
    else if (i == NX-1) then
      Ip2 = NX 
    endif


      ! some interpolations
      rho0_half_x  = 0.5d0 * (rho0(i,j)  + rho0(Im1,j))


     
      u_mm = pressure(Im2,j); u_m = pressure(Im1,j); u_p = pressure(i,j); u_pp = pressure(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dpressure_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      
      eq2_memory_dpressure_dx(i,j) = b_x(i) * eq2_memory_dpressure_dx(i,j) + a_x(i) * value_dpressure_dx
      
     
     !
     value_dpressure_dx = value_dpressure_dx * one_over_K_x(i) + eq2_memory_dpressure_dx(i,j)
    
     ! update
     vx(i,j) = vx(i,j) - value_dpressure_dx * DELTAT / rho0(i,j)
     
   enddo
  enddo


  do j = 1,NY-1
    do i = 1,NX-1

       Im1 = i-1
       Im2 = i-2
       Ip1 = i+1
       Ip2 = i+2
       Jm1 = j-1
       Jm2 = j-2
       Jp1 = j+1
       Jp2 = j+2
      
      if (j == 1) then
      Jm1 = 1
      Jm2 = 1
    elseif (j == 2) then
      Jm2 = 1
    else if (j == NY) then
      Jp1 = NY
      Jp2 = NY
    else if (j == NY-1) then
      Jp2 = NY
    endif
    if (i == 1) then
      Im1 = 1
      Im2 = 1
    elseif (i == 2) then
      Im2 = 1
    else if (i == NX) then
       Ip1 = NX
       Ip2 = NX
    else if (i == NX-1) then
      Ip2 = NX 
    endif


      ! compute derivative of the pressure and density according to y
      u_mm = pressure(i,Jm1); u_m = pressure(i,j); u_p = pressure(i,Jp1); u_pp = pressure(i,Jp2)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dpressure_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
     
      eq3_memory_dpressure_dy(i,j) = b_y_half(j) * eq3_memory_dpressure_dy(i,j) + a_y_half(j) * value_dpressure_dy


      !
      value_dpressure_dy = value_dpressure_dy * one_over_K_y_half(j) + eq3_memory_dpressure_dy(i,j)


      ! update
      vy(i,j) = vy(i,j) - value_dpressure_dy * DELTAT / rho0(i,j)


    enddo
  enddo

  ! Dircihlet conditions
  vx(1,:) = ZERO
  vx(NX,:) = ZERO
  
  vx(:,1) = ZERO
  vx(:,NY) = ZERO
  
  vy(1,:) = ZERO
  vy(NX,:) = ZERO
  
  vy(:,1) = ZERO
  vy(:,NY) = ZERO
  
  
! store seismograms
  do irec = 1,NREC

! beware here that the two components of the velocity vector are not defined at the same point
! in a staggered grid, and thus the two components of the velocity vector are recorded at slightly different locations,
! vy is staggered by half a grid cell along X and along Y with respect to vx
    sisvx(it,irec) = vx(ix_rec(irec),iy_rec(irec))
    sisvy(it,irec) = vy(ix_rec(irec),iy_rec(irec))
    sispressure(it,irec) = pressure(ix_rec(irec),iy_rec(irec))
    sisrhop(it,irec) = rhop(ix_rec(irec),iy_rec(irec))

  enddo


! output information
  if ((mod(it,IT_DISPLAY) == 0 .or. it == 5) .and. (save_sismos)) then


! print maximum of pressure and of norm of velocity
    pressurenorm = maxval(abs(pressure))
    velocnorm = maxval(sqrt(vx**2 + vy**2))
    print *,'Time step # ',it,' out of ',NSTEP
    print *,'Time: ',sngl((it-1)*DELTAT),' seconds'
    print *,'Max absolute value of pressure = ',pressurenorm
    print *,'Max norm velocity vector V (m/s) = ',velocnorm

    print *
! check stability of the code, exit if unstable
    if (pressurenorm > STABILITY_THRESHOLD .or. velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

    call create_color_image(pressure,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
                         0,.FALSE.,.FALSE.,.FALSE.,.FALSE.,field_number)

! save the part of the seismograms that has been computed so far, so that users can monitor the progress of the simulation

    call write_seismograms(sisvx,sisvy,sispressure,sisrhop,NSTEP,NREC,DELTAT,t0)

  endif


  enddo   ! end of time loop

! save seismograms
  if (save_sismos) then
    call write_seismograms(sisvx,sisvy,sispressure,sisrhop,NSTEP,NREC,DELTAT,t0)
  endif


  if (it == NSTEP) then
! create script for Gnuplot
  open(unit=20,file='./OUTPUT/plotgnu',status='unknown')
    write(20,*) 'set term x11'
    write(20,*) '# set term postscript landscape monochrome dashed "Helvetica" 22'
    write(20,*)
    write(20,*) 'set xlabel "Time (s)"'
    write(20,*) 'set ylabel "Amplitude (m / s)"'
    write(20,*)

    write(20,*) 'set output "v_sigma_Vx_receiver_001.eps"'
    write(20,*) 'plot "Vx_file_001.dat" t ''Vx C-PML'' w l lc 1'
    write(20,*) 'pause -1 "Hit any key..."'
    write(20,*)

    write(20,*) 'set output "v_sigma_Vy_receiver_001.eps"'
    write(20,*) 'plot "Vy_file_001.dat" t ''Vy C-PML'' w l lc 1'
    write(20,*) 'pause -1 "Hit any key..."'
    write(20,*)

    write(20,*) 'set output "v_sigma_Vx_receiver_002.eps"'
    write(20,*) 'plot "Vx_file_002.dat" t ''Vx C-PML'' w l lc 1'
    write(20,*) 'pause -1 "Hit any key..."'
    write(20,*)

    write(20,*) 'set output "v_sigma_Vy_receiver_002.eps"'
    write(20,*) 'plot "Vy_file_002.dat" t ''Vy C-PML'' w l lc 1'
    write(20,*) 'pause -1 "Hit any key..."'
    write(20,*)

    close(20)

    print *
    print *,'End of the simulation'
    print *
  endif



  endsubroutine forwardproblem
