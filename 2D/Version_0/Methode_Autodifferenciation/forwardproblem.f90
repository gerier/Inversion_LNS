subroutine forwardproblem(CI_pressure, CI_rhop, CI_vx, CI_vy,  p0, rho, v0x, v0y, kappa_unrelaxed, time_start, time_end) 
 !

use parameters , only : pressure, rhop, vx, vy, &
                       sispressure, sisrhop, sisvx, sisvy, NREC,IT_DISPLAY, ix_rec, iy_rec, &
                       NX, NY, NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX, &
                       NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,ONE_OVER_SIX_DELTAX,ONE_OVER_SIX_DELTAY, &
                       DELTAX, DELTAY, DELTAT, NSTEP, t,   &
                       a, f0, t0, pi, factor, ISOURCE, JSOURCE, source_term, &
                       distance2, factor_ssf, SSF_Sigma, & 
                       ZERO, gamma_chimie, stability_threshold, use_checkpoint
implicit none



! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here

  integer :: time_start, time_end



  double precision, dimension(0:NX+1,0:NY+1) :: &
    CI_pressure, CI_rhop, CI_vx, CI_vy, &
    p0, rho, v0x, v0y, kappa_unrelaxed                                           

  double precision :: &
      value_dvx_dx, &
      value_dvy_dy, &
      value_dpressure_dx, &
      value_dpressure_dy, &
      value_dp0_dx, &
      value_dp0_dy, &
      value_drho_dx, &
      value_drho_dy, &
      value_drhop_dx, &
      value_drhop_dy, &
      value_dv0x_dx, &
      value_dv0y_dy, &
      value_vdrho,&
      value_v0drhop, &
      value_vdp0 , &
      value_v0dp

  double precision :: &
    value_rhodv, &
    value_dvy_dx, &
    value_dvx_dy, &
    value_dv0x_dx_prec, &
    value_dv0x_dx_next, &
    value_dv0y_dy_prec, &
    value_dv0y_dy_next, &
    value_drho_dy_prec, &
    value_drho_dy_next, &
    value_drho_dx_prec, &
    value_drho_dx_next

  double precision :: &
      vx_half_x, &
      vy_half_y,&
      v0x_half_x, &
      v0y_half_y, &
      v0x_half_x_half_y, &
      v0y_half_x_half_y, &
      vy_half_x_half_y, &
      vx_half_x_half_y

   double precision :: &
      value_vdv0y, &
      value_vdv0x, &
      value_v0dvy, &
      value_v0dvx, &
      value_v0dv0y, &
      value_v0dv0x, &
      value_dv0y_dx, &
      value_dv0x_dy,&
      dkappax,&
      dkappay,&
      dkappa_v

  double precision :: &
      value_dv0x_dy_prec,&
      value_dv0y_dx_prec,&
      value_dv0x_dy_next,&
      value_dv0y_dx_next,&
      value_v0x_dv0y_dx,&
      value_v0y_dv0x_dy,&
      value_v0x_dv0x_dx,&
      value_v0y_dv0y_dy

  double precision :: &
      value_drhovxv0x_dx,&
      value_drhovxv0y_dy,&
      value_drhovyv0x_dx,&
      value_drhovyv0y_dy
 
  double precision, dimension (0:NX+1,0:NY+1) :: & 
      rhovxv0x,&
      rhovxv0y,&
      rhovyv0x,&
      rhovyv0y

  integer :: Ip1,Im1, Jp1,Jm1


  integer :: i,j,it,irec

  double precision :: velocnorm,pressurenorm

  double precision, dimension(0:NX+1,0:NY+1) :: aux_rhop, aux_p, aux_vx, aux_vy

 ! to interpolate material parameters or velocity at the right location in the staggered grid cell
  double precision kappa_half_x,rho_half_x_half_y,vy_interpolated,&
                 rho_half_x,rho_half_y,rhop_half_x,rhop_half_y
 
double precision :: derivative

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



!---
!---  beginning of time loop
!---

  do it = time_start, time_end

!-----------------------------------------------------------------------
! compute pressure and update memory variables for C-PML
! also update memory variables for viscoacoustic attenuation if needed
!-----------------------------------------------------------------------
  do j=1,NY
    do i=1,NX
      aux_p(i,j) = pressure(i,j)
      aux_rhop(i,j) = rhop(i,j)
    enddo
  enddo


    do j = 1,NY   
      do i = 1,NX
      
      Im1 = i-1
      Ip1 = i+1
      Jm1 = j-1
      Jp1 = j+1

    if (j == 1) then
       Jm1 = NY-1  
    else if (j == NY) then
       Jp1 = 2
    endif 
    if (i == 1) then
       Im1 = NX-1
    else if (i == NX) then
       Ip1 = 2
    endif
    

! interpolate material parameters at the right location in the staggered grid cell
        vx_half_x = (vx(Ip1,j) + vx(i,j)) * 0.5d0
        vy_half_y = (vy(i,j) + vy(i,Jm1)) * 0.5d0
        v0x_half_x = (v0x(Ip1,j) + v0x(i,j)) * 0.5d0
        v0y_half_y = (v0y(i,j) + v0y(i,Jm1)) * 0.5d0
        
        call compute_centered_dU_dx_in_i(vx, i, j, value_dvx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
        call compute_centered_dU_dx_in_i(v0x, i, j, value_dv0x_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)

        call compute_centered_dU_dz_in_i(vy, i, Jm1, value_dvy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)  
        call compute_centered_dU_dz_in_i(v0y, i, Jm1, value_dv0y_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
        
        ! decentered derivative in x                     
        call compute_decentered_dU_dx_in_i(p0, i, j, value_dp0_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
        call compute_decentered_dU_dx_in_i(aux_p, i, j, value_dpressure_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
        call compute_decentered_dU_dx_in_i(rho, i, j, value_drho_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
        call compute_decentered_dU_dx_in_i(aux_rhop, i, j, value_drhop_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)

        ! decentered derivative in y                                
        call compute_decentered_dU_dz_in_i(p0, i, j, value_dp0_dy, ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y)   
        call compute_decentered_dU_dz_in_i(aux_p, i, j, value_dpressure_dy, ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y) 
        call compute_decentered_dU_dz_in_i(rho, i, j, value_drho_dy, ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y)   
        call compute_decentered_dU_dz_in_i(aux_rhop, i, j, value_drhop_dy, ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y) 

        ! TEST TO UNDERSTAND WHY WE DO NOT HAVE THE SAME SOLUTION WITH
        ! SPECFEM2D-DG-LNS
        call compute_decentered_dU_dx_in_i(v0y, i, Jm1, value_dv0y_dx_prec,ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x_half_y)
        call compute_decentered_dU_dz_in_i(v0x, i, j, value_dv0x_dy_prec,ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_x_half_y)
        call compute_decentered_dU_dx_in_i(v0y, i, j, value_dv0y_dx_next,ONE_OVER_SIX_DELTAX,NX,NY,v0y_half_x_half_y)
        call compute_decentered_dU_dz_in_i(v0x, Ip1, j, value_dv0x_dy_next,ONE_OVER_SIX_DELTAY,NX,NY,v0x_half_x_half_y)
        value_dv0y_dx = 0.5d0 * (value_dv0y_dx_prec + value_dv0y_dx_next)  
        value_dv0x_dy = 0.5d0 * (value_dv0x_dy_prec + value_dv0x_dy_next)
 


        ! intermediate computations
        value_vdp0   = (vx_half_x  * value_dp0_dx)       + (vy_half_y  * value_dp0_dy)
        value_v0dp   = (v0x_half_x * value_dpressure_dx) + (v0y_half_y * value_dpressure_dy)

        value_vdrho   = (vx_half_x  * value_drho_dx)   + (vy_half_y  * value_drho_dy)
        value_v0drhop = (v0x_half_x * value_drhop_dx)  + (v0y_half_y * value_drhop_dy)
 
        value_v0x_dv0x_dx = v0x_half_x * value_dv0x_dx 
        value_v0y_dv0x_dy = v0y_half_y * value_dv0x_dy 
        value_v0x_dv0y_dx = v0x_half_x * value_dv0y_dx
        value_v0y_dv0y_dy = v0y_half_y * value_dv0y_dy

        ! updateT
        pressure(i,j) = pressure(i,j) - (kappa_unrelaxed(i,j) * (value_dvx_dx + value_dvy_dy)) * DELTAT
        pressure(i,j) = pressure(i,j) - value_vdp0 * DELTAT ! * gamma_chimie
        pressure(i,j) = pressure(i,j) - gamma_chimie * pressure(i,j) * (value_dv0x_dx + value_dv0y_dy) * DELTAT
        pressure(i,j) = pressure(i,j) - value_v0dp * DELTAT 

        rhop(i,j) = rhop(i,j) - rho(i,j) * (value_dvx_dx + value_dvy_dy) * DELTAT
        rhop(i,j) = rhop(i,j) - value_vdrho * DELTAT
        rhop(i,j) = rhop(i,j) - rhop(i,j) * (value_dv0x_dx + value_dv0y_dy) * DELTAT
        rhop(i,j) = rhop(i,j) - value_v0drhop * DELTAT

      enddo
    enddo

  !do i=1,NX
  !  pressure(i,1) = pressure(i,NY)
  !  rhop(i,1) = rhop(i,NY)
  !enddo

  !do i=2,NY
  !  pressure(NX,i) = pressure(1,i)
  !  rhop(NX,i) = rhop(1,i)
  !enddo

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
      !rho_half_y = 0.5d0 * (rho(i,j+1) + rho(i,j))
      !vy(i,j) = vy(i,j) + source_term * DELTAT / rho_half_y * 0.5d0
      distance2 = ((i - Isource) * DELTAX)**2 + ((j - Jsource) * DELTAY)**2
      factor_ssf = exp( - distance2 / SSF_Sigma**2 )
      
      pressure(i,j) = pressure(i,j) + source_term * factor_ssf * DELTAT
      rhop(i,j) = rhop(i,j) + source_term * factor_ssf * DELTAT * (rho(i,j)/gamma_chimie / p0(i,j)) 
    !enddo

  ! write the source
  open(unit=211,file='./OUTPUT/source_time_function_model.dat')
    write(211,*) (t-t0), source_term
  close(211)


!--------------------------------------------------------
! compute velocity and update memory variables for C-PML
!--------------------------------------------------------

  do j=1,NY
   do i=1,NX
    aux_vx(i,j) = vx(i,j)
    aux_vy(i,j) = vy(i,j)

    Im1 = i-1
    Ip1 = i+1
    Jm1 = j-1
    Jp1 = j+1

    if (j == 1) then
       Jm1 = NY-1  
    else if (j == NY) then
       Jp1 = 2
    endif 
    if (i == 1) then
       Im1 = NX-1
    else if (i == NX) then
       Ip1 = 2
    endif

    v0y_half_x_half_y = (v0y(i,j)   + v0y(i,Jm1) + v0y(Im1,j) + v0y(Im1,Jm1)) * 0.25d0
    rho_half_x        = 0.5d0 * (rho(i,j)  + rho(Im1,j))
    v0x_half_x_half_y = (v0x(i,j)   + v0x(i,Jp1) + v0x(Ip1,j) + v0x(Ip1,Jp1)) * 0.25d0
    rho_half_y        = 0.5d0 * (rho(i,j) + rho(i,Jp1))
    
    rhovxv0x(i,j) = rho_half_x * aux_vx(i,j) * v0x(i,j) 
    rhovxv0y(i,j) = rho_half_x * aux_vx(i,j) * v0y_half_x_half_y
    rhovyv0x(i,j) = rho_half_y * aux_vy(i,j) * v0x_half_x_half_y
    rhovyv0y(i,j) = rho_half_y * aux_vy(i,j) * v0y(i,j)
   enddo
  enddo

  do j = 1,NY
    do i = 1,NX

    Im1 = i-1
    Ip1 = i+1
    Jm1 = j-1
    Jp1 = j+1

    if (j == 1) then
       Jm1 = NY-1  
    else if (j == NY) then
       Jp1 = 2
    endif 
    if (i == 1) then
       Im1 = NX-1
    else if (i == NX) then
       Ip1 = 2
    endif
    
      ! some interpolations
      rho_half_x  = 0.5d0 * (rho(i,j)  + rho(Im1,j))
      rhop_half_x = 0.5d0 * (rhop(i,j) + rhop(Im1,j))

      vy_half_x_half_y  = (vy(i,j)    + vy(i,Jm1)  + vy(Im1,j) + vy(Im1,Jm1))   * 0.25d0
      v0y_half_x_half_y = (v0y(i,j)   + v0y(i,Jm1) + v0y(Im1,j) + v0y(Im1,Jm1)) * 0.25d0

      ! compute derivative of the pressure, density according to x
      call compute_centered_dU_dx_in_i(rho, Im1, j, value_drho_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
      call compute_centered_dU_dx_in_i(pressure, Im1, j, value_dpressure_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)

      ! compute derivative of v0x, vx according to x
      call compute_decentered_dU_dx_in_i(v0x, i, j, value_dv0x_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x(i,j))
      call compute_decentered_dU_dx_in_i(aux_vx, i, j, value_dvx_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x(i,j))
     
      ! compute derivative of v0x, v0y, vx according to y
      call compute_decentered_dU_dz_in_i(v0x, i, j, value_dv0x_dy, ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_x_half_y)
      call compute_decentered_dU_dz_in_i(aux_vx, i, j, value_dvx_dy, ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_x_half_y)
      
      call compute_decentered_dU_dz_in_i(rho, Im1, j, value_drho_dy_prec, ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y)
      call compute_decentered_dU_dz_in_i(rho, i, j, value_drho_dy_next, ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y)
      value_drho_dy = 0.5d0 * (value_drho_dy_next + value_drho_dy_prec)
      
      call compute_centered_dU_dz_in_i(v0y, Im1, Jm1, value_dv0y_dy_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
      call compute_centered_dU_dz_in_i(v0y, i, Jm1, value_dv0y_dy_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
      value_dv0y_dy = 0.5d0 * (value_dv0y_dy_next + value_dv0y_dy_prec)
      
      !
      call compute_decentered_dU_dx_in_i(rhovxv0x, i, j, value_drhovxv0x_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x(i,j))
      call compute_decentered_dU_dz_in_i(rhovxv0y, i, j, value_drhovxv0y_dy,ONE_OVER_SIX_DELTAY,NX,NY,v0y(i,j))

      ! intermediate computations: (v0 . nabla) v0_x ; (v' . nabla) v0_x ; (v0 . nabla) v0_x
      value_v0dv0x = v0x(i,j)     * value_dv0x_dx + v0y_half_x_half_y * value_dv0x_dy
      value_vdv0x  = aux_vx(i,j)  * value_dv0x_dx + vy_half_x_half_y  * value_dv0x_dy
      value_v0dvx  = v0x(i,j)     * value_dvx_dx  + v0y_half_x_half_y * value_dvx_dy

      ! compute rho div v0 ; v0 div rho
      value_rhodv = rho_half_x * value_dv0x_dx + rho_half_x        * value_dv0y_dy
      value_vdrho = v0x(i,j)   * value_drho_dx + v0y_half_x_half_y * value_drho_dy

      ! update
      vx(i,j) = vx(i,j) - value_dpressure_dx * DELTAT / rho_half_x
      !  VERSION NON CONSERVATIVE 
      vx(i,j) = vx(i,j) - aux_vx(i,j)     * (value_vdrho + value_rhodv) * DELTAT / rho_half_x
      vx(i,j) = vx(i,j) - value_v0dvx  * DELTAT 
      ! VERSION CONSERVATIVE
      !vx(i,j) = vx(i,j) - value_drhovxv0x_dx * DELTAT / rho_half_x
      !vx(i,j) = vx(i,j) - value_drhovxv0y_dy * DELTAT / rho_half_x

 
      vx(i,j) = vx(i,j) - rhop_half_x * value_v0dv0x * DELTAT / rho_half_x
      vx(i,j) = vx(i,j) - value_vdv0x  * DELTAT
    enddo
  enddo

  do j = 1,NY
    do i = 1,NX


    Im1 = i-1
    Ip1 = i+1
    Jm1 = j-1
    Jp1 = j+1

    if (j == 1) then
       Jm1 = NY-1  
    else if (j == NY) then
       Jp1 = 2
    endif 
    if (i == 1) then
       Im1 = NX-1
    else if (i == NX) then
       Ip1 = 2
    endif
    
!     interpolate density at the right location in the staggered grid cell
      rho_half_y  = 0.5d0 * (rho(i,j)  + rho(i,Jp1))
      rhop_half_y = 0.5d0 * (rhop(i,j) + rhop(i,Jp1))

      vx_half_x_half_y  = (aux_vx(Ip1,Jp1)  + aux_vx(i,Jp1)  +aux_vx(Ip1,j)  + aux_vx(i,j))  * 0.25d0
      vy_half_y         = (aux_vy(i,j)      + aux_vy(i,Jp1))                         * 0.5d0
      v0x_half_x_half_y = (v0x(Ip1,Jp1) + v0x(i,Jp1) +v0x(Ip1,j) + v0x(i,j)) * 0.25d0
      v0y_half_y        = (v0y(i,j)     + v0y(i,Jm1))                        * 0.5d0
      
      ! compute derivative of the pressure and density according to y
      call compute_centered_dU_dz_in_i(pressure, i, j, value_dpressure_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
      call compute_centered_dU_dz_in_i(rho, i, j, value_drho_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)

      ! compute derivative of v0x, vx according to x
      call compute_decentered_dU_dx_in_i(v0y, i, j, value_dv0y_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x_half_y)
      call compute_decentered_dU_dx_in_i(aux_vy, i, j, value_dvy_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x_half_y)
      
      call compute_decentered_dU_dx_in_i(rho, i, j, value_drho_dx_prec, ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
      call compute_decentered_dU_dx_in_i(rho, i, Jp1, value_drho_dx_next, ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
      value_drho_dx = 0.5d0 * (value_drho_dx_next + value_drho_dx_prec)
      
      call compute_centered_dU_dx_in_i(v0x, i, j, value_dv0x_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
      call compute_centered_dU_dx_in_i(v0x, i, Jp1, value_dv0x_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
      value_dv0x_dx = 0.5d0 * (value_dv0x_dx_next + value_dv0x_dx_prec)
      
      ! compute derivative of v0y, vy according to y
      call compute_decentered_dU_dz_in_i(v0y, i, j, value_dv0y_dy, ONE_OVER_SIX_DELTAX,NX,NY,v0y(i,j))
      call compute_decentered_dU_dz_in_i(aux_vy, i, j, value_dvy_dy, ONE_OVER_SIX_DELTAX,NX,NY,v0y(i,j))

      ! 
      call compute_decentered_dU_dx_in_i(rhovyv0x, i, j, value_drhovyv0x_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x_half_y)
      call compute_decentered_dU_dz_in_i(rhovyv0y, i, j, value_drhovyv0y_dy, ONE_OVER_SIX_DELTAY,NX,NY, v0y(i,j))

      ! intermediate computations: (v0 . nabla) v0_y ; (v' . nabla) v0_y ; (v0 . nabla) v0_y
      value_v0dv0y = v0x_half_x_half_y * value_dv0y_dx + v0y(i,j)     * value_dv0y_dy 
      value_vdv0y  = vx_half_x_half_y  * value_dv0y_dx + aux_vy(i,j)  * value_dv0y_dy 
      value_v0dvy  = v0x_half_x_half_y * value_dvy_dx  + v0y(i,j)     * value_dvy_dy 

      ! compute rho div v0 ; v0 div rho
      value_rhodv = rho_half_y        * value_dv0x_dx  + rho_half_y * value_dv0y_dy
      value_vdrho = v0x_half_x_half_y * value_drho_dx  + v0y(i,j)   * value_drho_dy

      ! update
      vy(i,j) = vy(i,j) - value_dpressure_dy * DELTAT / rho_half_y
      ! VERSION NON CONSERVATIVE
      vy(i,j) = vy(i,j) - aux_vy(i,j)     * (value_vdrho + value_rhodv) * DELTAT / rho_half_y
      vy(i,j) = vy(i,j) - value_v0dvy  * DELTAT
      ! VERSION CONSERVATIVE
      !vy(i,j) = vy(i,j) - value_drhovyv0x_dx * DELTAT / rho_half_y      
      !vy(i,j) = vy(i,j) - value_drhovyv0y_dy * DELTAT / rho_half_y
      
      vy(i,j) = vy(i,j) - rhop_half_y * value_v0dv0y * DELTAT / rho_half_y
      vy(i,j) = vy(i,j) - value_vdv0y  * DELTAT

    enddo
  enddo


! Dirichlet conditions (rigid boundaries) on the edges or at the bottom of the PML layers
 ! do i=1,NX
 !   vx(i,1) = vx(i,NY)
 !   vy(i,NY) = vy(i,1)
 ! enddo

 ! do i=1,NY
 !   vy(NX,i) = vy(1,i)
 !   vx(1,i) =  vx(NX,i)
 ! enddo


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
  if ((mod(it,IT_DISPLAY) == 0 .or. it == 0) .and. (.not. use_checkpoint)) then


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
                         0,.FALSE.,.FALSE.,.FALSE.,.FALSE.,3)

! save the part of the seismograms that has been computed so far, so that users can monitor the progress of the simulation

    call write_seismograms(sisvx,sisvy,sispressure,sisrhop,NSTEP,NREC,DELTAT,t0)

  endif
  

  enddo   ! end of time loop

! save seismograms
  if (.not. use_checkpoint) then
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





