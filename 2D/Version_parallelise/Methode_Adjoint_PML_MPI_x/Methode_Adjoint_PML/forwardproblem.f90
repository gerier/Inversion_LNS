subroutine forwardproblem(p0,rho0,windx,windy,kappa_unrelaxed, it_start, it_end, field_number)
 !

use parameters !, only : pressure, rhop, vx, vy, &
               !        sispressure, sisrhop, sisvx, sisvy, NREC,IT_DISPLAY, ix_rec, iy_rec, &
               !        NX, NY, NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX, &
               !        NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,ONE_OVER_SIX_DELTAX,ONE_OVER_SIX_DELTAY, &
               !        DELTAX, DELTAY, DELTAT, NSTEP, t,   &
               !        a, f0, t0, pi, factor, ISOURCE, JSOURCE, source_term, &
               !        distance2, factor_ssf, SSF_Sigma, &
               !        ZERO, gamma_chimie, stability_threshold, save_sismos,&
               !        a_x, a_y, a_x_half, a_y_half, &
               !        b_x, b_y, b_x_half, b_y_half, &
               !        c_x, c_y, c_x_half, c_y_half, &
               !        one_over_K_x, one_over_K_x_half, one_over_K_y, one_over_K_y_half, &
               !        NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX
implicit none



! arrays for the memory variables
! could declare these arrays in PML only to save a lot of memory, but proof of concept only here

  integer :: it_start, it_end

  integer :: field_number

  double precision, dimension(-1:NX_LOCAL+2,0:NY+1) :: &
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

        
  double precision  :: u_mm, u_m, u_p, u_pp
  
  integer :: Ip1,Im1, Ip2, Im2, Jp1,Jm1, Jp2, Jm2


  integer :: i,j,it,irec

  double precision :: velocnorm,pressurenorm

  double precision, dimension(-1:NX_LOCAL+2,0:NY+1) :: rhop_old, p_old, vx_old, vy_old

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
  
  ! vx
  call MPI_SENDRECV(vx(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,vx(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
         
    
  call MPI_SENDRECV(vx(NX_LOCAL,:),number_of_values/2,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,vx(0,:),number_of_values/2, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
         

  ! windx       
  call MPI_SENDRECV(windx(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,windx(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(windx(NX_LOCAL,:),number_of_values/2,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,windx(0,:),number_of_values/2, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)

  ! vy
  call MPI_SENDRECV(vy(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,vy(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(vy(1,:),number_of_values/2,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,vy(NX_LOCAL+1,:),number_of_values/2, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  ! windy       
  call MPI_SENDRECV(windy(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,windy(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
   call MPI_SENDRECV(windy(1,:),number_of_values/2,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,windy(NX_LOCAL+1,:),number_of_values/2, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)

  ! rhop 
   call MPI_SENDRECV(rhop_old(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,rhop_old(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
   call MPI_SENDRECV(rhop_old(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,rhop_old(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)

  ! pp
   call MPI_SENDRECV(p_old(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,p_old(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
   call MPI_SENDRECV(p_old(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,p_old(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  
  ! rho0
   call MPI_SENDRECV(rho0(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,rho0(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
   call MPI_SENDRECV(rho0(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,rho0(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
        
  ! p0
   call MPI_SENDRECV(p0(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,p0(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
   call MPI_SENDRECV(p0(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,p0(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)


  do j = 1,NY
    do i = 1,NX_LOCAL

       i_global = i + offset_i

       Im1 = i-1
       Im2 = i-2
       Ip1 = i+1
       Ip2 = i+2
       Jm1 = j-1
       Jm2 = j-2
       Jp1 = j+1
       Jp2 = j+2

       ! modify index dependingon the boundary condition
       call get_index_boundarycondition(i,j, Im2, Im1,Ip1, Ip2, Jm2, Jm1,Jp1, Jp2)

        ! interpolate material parameters at the right location in the staggered grid cell
        vx_half_x = (vx(Ip1,j) + vx(i,j)) * 0.5d0
        vy_half_y = (vy(i,j) + vy(i,Jm1)) * 0.5d0
        windx_half_x = (windx(Ip1,j) + windx(i,j)) * 0.5d0
        windy_half_y = (windy(i,j) + windy(i,Jm1)) * 0.5d0

        ! derivative computations
	u_mm = vx(Im1,j); u_m = vx(i,j); u_p = vx(Ip1,j); u_pp = vx(Ip2,j)
        call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
        
        u_mm = windx(Im1,j); u_m = windx(i,j); u_p = windx(Ip1,j); u_pp = windx(Ip2,j)
        call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

        u_mm = vy(i,Jm2); u_m = vy(i,Jm1); u_p = vy(i,j); u_pp = vy(i,Jp1)
        call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      
        u_mm = windy(i,Jm2); u_m = windy(i,Jm1); u_p = windy(i,j); u_pp = windy(i,Jp1)
        call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)

        eq1_memory_dvx_dx_fw(i,j) = b_x_half(i_global) * eq1_memory_dvx_dx_fw(i,j) + a_x_half(i_global) * value_dvx_dx
        eq1_memory_dvy_dy_fw(i,j) = b_y(j) * eq1_memory_dvy_dy_fw(i,j) + a_y(j) * value_dvy_dy
        eq1_memory_dwindx_dx_fw(i,j) = c_x_half(i_global) * value_dwindx_dx
        eq1_memory_dwindy_dy_fw(i,j) = c_y(j) * value_dwindy_dy        

        ! decentered derivative in x
        if (windx_half_x >= 0) then
          u_mm = p0(Im2,j)      ; u_m = p0(Im1,j)      ; u_p = p0(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, p0(i,j), u_p, value_dp0_dx, ONE_OVER_SIX_DELTAX)
          u_mm = p_old(Im2,j); u_m = p_old(Im1,j); u_p = p_old(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, p_old(i,j), u_p, value_dpressure_dx, ONE_OVER_SIX_DELTAX)
          u_mm = rho0(Im2,j)    ; u_m = rho0(Im1,j)    ; u_p = rho0(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, rho0(i,j), u_p, value_drho0_dx, ONE_OVER_SIX_DELTAX)
          u_mm = rhop_old(Im2,j)     ; u_m = rhop_old(Im1,j)     ; u_p = rhop_old(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, rhop_old(i,j), u_p, value_drhop_dx, ONE_OVER_SIX_DELTAX)
        else 
          u_pp = p0(Ip2,j)      ; u_p = p0(Ip1,j)      ; u_m = p0(Im1,j) 
          call compute_decentered_forward_dU(u_m, p0(i,j), u_p, u_pp, value_dp0_dx, ONE_OVER_SIX_DELTAX)
          u_pp = p_old(Ip2,j); u_p = p_old(Ip1,j); u_m = p_old(Im1,j) 
          call compute_decentered_forward_dU(u_m, p_old(i,j), u_p, u_pp, value_dpressure_dx, ONE_OVER_SIX_DELTAX)
          u_pp = rho0(Ip2,j)    ; u_p = rho0(Ip1,j)    ; u_m = rho0(Im1,j) 
          call compute_decentered_forward_dU(u_m, rho0(i,j), u_p, u_pp,value_drho0_dx, ONE_OVER_SIX_DELTAX)
          u_pp = rhop_old(Ip2,j)     ; u_p = rhop_old(Ip1,j)     ; u_m = rhop_old(Im1,j) 
          call compute_decentered_forward_dU(u_m, rhop_old(i,j), u_p, u_pp, value_drhop_dx, ONE_OVER_SIX_DELTAX)
        endif
     
        eq1_memory_dpressure_dx_fw(i,j) = b_x_half(i_global)*eq1_memory_dpressure_dx_fw(i,j) + a_x_half(i_global)*value_dpressure_dx
        eq1_memory_drhop_dx_fw(i,j) = b_x_half(i_global) * eq1_memory_drhop_dx_fw(i,j) + a_x_half(i_global) * value_drhop_dx
        eq1_memory_dp0_dx_fw(i,j) = c_x_half(i_global) * value_dp0_dx
        eq1_memory_drho0_dx_fw(i,j) = c_x_half(i_global) * value_drho0_dx 
             
                
        if (windy_half_y >= 0) then
          u_mm = p0(i,Jm2)      ; u_m = p0(i,Jm1)      ; u_p = p0(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, p0(i,j), u_p, value_dp0_dy, ONE_OVER_SIX_DELTAY)
          u_mm = p_old(i,Jm2); u_m = p_old(i,Jm1); u_p = p_old(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, p_old(i,j), u_p, value_dpressure_dy, ONE_OVER_SIX_DELTAY)
          u_mm = rho0(i,Jm2)    ; u_m = rho0(i,Jm1)    ; u_p = rho0(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rho0(i,j), u_p, value_drho0_dy, ONE_OVER_SIX_DELTAY)
          u_mm = rhop_old(i,Jm2)     ; u_m = rhop_old(i,Jm1)     ; u_p = rhop_old(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rhop_old(i,j), u_p, value_drhop_dy, ONE_OVER_SIX_DELTAY)
        else 
          u_pp = p0(i,Jp2)      ; u_p = p0(i,Jp1)      ; u_m = p0(i,Jm1) 
          call compute_decentered_forward_dU(u_m, p0(i,j), u_p, u_pp, value_dp0_dy, ONE_OVER_SIX_DELTAY)
          u_pp = p_old(i,Jp2); u_p = p_old(i,Jp1); u_m = p_old(i,Jm1) 
          call compute_decentered_forward_dU(u_m, p_old(i,j), u_p, u_pp, value_dpressure_dy, ONE_OVER_SIX_DELTAY)
          u_pp = rho0(i,Jp2)    ; u_p = rho0(i,Jp1)    ; u_m = rho0(i,Jm1) 
          call compute_decentered_forward_dU(u_m, rho0(i,j), u_p, u_pp,value_drho0_dy, ONE_OVER_SIX_DELTAY)
          u_pp = rhop_old(i,Jp2)     ; u_p = rhop_old(i,Jp1)       ; u_m = rhop_old(i,Jm1) 
          call compute_decentered_forward_dU(u_m, rhop_old(i,j), u_p, u_pp, value_drhop_dy, ONE_OVER_SIX_DELTAY)
        endif
        
        eq1_memory_dpressure_dy_fw(i,j) = b_y(j) * eq1_memory_dpressure_dy_fw(i,j) + a_y(j) * value_dpressure_dy
        eq1_memory_drhop_dy_fw(i,j) = b_y(j) * eq1_memory_drhop_dy_fw(i,j) + a_y(j) * value_drhop_dy
        eq1_memory_dp0_dy_fw(i,j) = c_y(j) * value_dp0_dy
        eq1_memory_drho0_dy_fw(i,j) = c_y(j) * value_drho0_dy 
        
        value_dvx_dx = value_dvx_dx * one_over_K_x_half(i_global) + eq1_memory_dvx_dx_fw(i,j)
        value_dvy_dy = value_dvy_dy * one_over_K_y(j) + eq1_memory_dvy_dy_fw(i,j) 
        value_dwindx_dx = value_dwindx_dx * one_over_K_x_half(i_global) + eq1_memory_dwindx_dx_fw(i,j)
        value_dwindy_dy = value_dwindy_dy * one_over_K_y(j) + eq1_memory_dwindy_dy_fw(i,j) 
        
        value_dp0_dx = value_dp0_dx * one_over_K_x_half(i_global) + eq1_memory_dp0_dx_fw(i,j)
        value_drho0_dx = value_drho0_dx * one_over_K_x_half(i_global) + eq1_memory_drho0_dx_fw(i,j)
        value_dpressure_dx = value_dpressure_dx * one_over_K_x_half(i_global) + eq1_memory_dpressure_dx_fw(i,j)
        value_drhop_dx = value_drhop_dx * one_over_K_x_half(i_global) + eq1_memory_drhop_dx_fw(i,j)
        
        value_dp0_dy = value_dp0_dy * one_over_K_y(j) + eq1_memory_dp0_dy_fw(i,j) 
        value_drho0_dy = value_drho0_dy * one_over_K_y(j) + eq1_memory_drho0_dy_fw(i,j) 
        value_dpressure_dy = value_dpressure_dy * one_over_K_y(j) + eq1_memory_dpressure_dy_fw(i,j) 
        value_drhop_dy = value_drhop_dy * one_over_K_y(j) + eq1_memory_drhop_dy_fw(i,j)    
        
        
        ! intermediate computations
        value_vdp0   = (vx_half_x  * value_dp0_dx)         + (vy_half_y  * value_dp0_dy)
        value_v0dp   = (windx_half_x * value_dpressure_dx) + (windy_half_y * value_dpressure_dy)

        value_vdrho0   = (vx_half_x  * value_drho0_dx)   + (vy_half_y  * value_drho0_dy)
        value_v0drhop = (windx_half_x * value_drhop_dx)  + (windy_half_y * value_drhop_dy)
        
        ! updateT
        pressure(i,j) = pressure(i,j) - (kappa_unrelaxed(i,j) * (value_dvx_dx + value_dvy_dy)) * DELTAT
        pressure(i,j) = pressure(i,j) - value_vdp0 * DELTAT ! * gamma_chimie
        pressure(i,j) = pressure(i,j) - gamma_chimie * pressure(i,j) * (value_dwindx_dx + value_dwindy_dy) * DELTAT
        pressure(i,j) = pressure(i,j) - value_v0dp * DELTAT

        rhop(i,j) = rhop(i,j) - rho0(i,j) * (value_dvx_dx + value_dvy_dy) * DELTAT
        rhop(i,j) = rhop(i,j) - value_vdrho0 * DELTAT
        rhop(i,j) = rhop(i,j) - rhop(i,j) * (value_dwindx_dx + value_dwindy_dy) * DELTAT
        rhop(i,j) = rhop(i,j) - value_v0drhop * DELTAT



      enddo
    enddo

 ! Dircihlet conditions  
 if (USE_PML_XMIN) then
 
    if (rank == nb_procs-1) then
      rhop(NX_LOCAL:NX_LOCAL+2,:) = ZERO
      pressure(NX_LOCAL:NX_LOCAL+2,:) = ZERO
    else if (rank == 0) then
      rhop(-1:1,:) = ZERO
      pressure(-1:1,:) = ZERO    
    endif
   !else
   !  if (rank == nb_procs - 1) then
   !    pressure(NX_LOCAL,:) = pressure(0,:)
   !    rhop(NX_LOCAL,:) = rhop(0,:)
    ! endif
  endif
  
  if (USE_PML_YMIN) then
    rhop(:,1) = ZERO
    pressure(:,1) = ZERO
    
    rhop(:,NY) = ZERO
    pressure(:,NY) = ZERO
  else
    pressure(:,NY) = pressure(:,1)
    rhop(:,NY) = rhop(:,1)
  endif
  

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
  i = ISOURCE - offset_i
  if (0 <= i .and. i < NX_LOCAL) then ! TODO
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
      rhop(i,j) = rhop(i,j) + source_term * factor_ssf * DELTAT * (rho0(i,j)/gamma_chimie / p0(i,j))
    !enddo
  endif

  ! write the source
  open(unit=211,file='./OUTPUT/source_time_function_model.dat',status='unknown',position="append")
    write(211,*) (t-t0), source_term
  close(211)


!--------------------------------------------------------
! compute velocity and update memory variables for C-PML
!--------------------------------------------------------

  vx_old(:,:) = vx(:,:)
  vy_old(:,:) = vy(:,:)

 ! rhop
  call MPI_SENDRECV(rhop(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,rhop(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(rhop(1,:),number_of_values/2,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,rhop(NX_LOCAL+1,:),number_of_values/2, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  ! pressure
  call MPI_SENDRECV(pressure(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,pressure(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(pressure(1,:),number_of_values/2,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,pressure(NX_LOCAL+1,:),number_of_values/2, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)

 ! rho0
   call MPI_SENDRECV(rho0(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,rho0(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(rho0(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,rho0(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
         
 ! p0
   call MPI_SENDRECV(p0(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,p0(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(p0(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,p0(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
         
          
 ! windx
  call MPI_SENDRECV(windx(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,windx(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(windx(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,windx(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)

  ! vx
  call MPI_SENDRECV(vx_old(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,vx_old(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(vx_old(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,vx_old(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)

  ! windy
  call MPI_SENDRECV(windy(NX_LOCAL,:),number_of_values/2,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,windy(0,:),number_of_values/2, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)


 ! vy
  call MPI_SENDRECV(vy_old(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,vy_old(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
         
 ! TODO : manque le cas Im1, Jm1
         
  
  do j = 1,NY
    do i = 1,NX_LOCAL
    
      i_global = i + offset_i

       Im1 = i-1
       Im2 = i-2
       Ip1 = i+1
       Ip2 = i+2
       Jm1 = j-1
       Jm2 = j-2
       Jp1 = j+1
       Jp2 = j+2
      
      ! modify index depending on the boundary condition
      call get_index_boundarycondition(i,j, Im2, Im1,Ip1, Ip2, Jm2, Jm1,Jp1, Jp2) 

      ! some interpolations
      rho0_half_x  = 0.5d0 * (rho0(i,j)  + rho0(Im1,j))
      rhop_half_x = 0.5d0 * (rhop(i,j) + rhop(Im1,j))

      vy_half_x_half_y  = (vy(i,j)    + vy(i,Jm1)  + vy(Im1,j) + vy(Im1,Jm1))   * 0.25d0
      windy_half_x_half_y = (windy(i,j)   + windy(i,Jm1) + windy(Im1,j) + windy(Im1,Jm1)) * 0.25d0

      ! compute derivative of the pressure, density according to x
      u_mm = rho0(Im2,j); u_m = rho0(Im1,j); u_p = rho0(i,j); u_pp = rho0(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drho0_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      u_mm = pressure(Im2,j); u_m = pressure(Im1,j); u_p = pressure(i,j); u_pp = pressure(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dpressure_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      
      eq2_memory_drho0_dx_fw(i,j) = c_x(i_global) * value_drho0_dx
      eq2_memory_dpressure_dx_fw(i,j) = b_x(i_global) * eq2_memory_dpressure_dx_fw(i,j) + a_x(i_global) * value_dpressure_dx
      
      
      ! compute derivative of windx, vx according to x
      if (windx(i,j) >= 0) then
          u_mm = windx(Im2,j)      ; u_m = windx(Im1,j)      ; u_p = windx(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, windx(i,j), u_p, value_dwindx_dx, ONE_OVER_SIX_DELTAX)
          u_mm = vx_old(Im2,j); u_m = vx_old(Im1,j); u_p = vx_old(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, vx_old(i,j), u_p, value_dvx_dx, ONE_OVER_SIX_DELTAX)
     else 
          u_pp = windx(Ip2,j)    ; u_p = windx(Ip1,j)    ; u_m = windx(Im1,j) 
          call compute_decentered_forward_dU(u_m, windx(i,j), u_p, u_pp, value_dwindx_dx, ONE_OVER_SIX_DELTAX)
          u_pp = vx_old(Ip2,j)     ; u_p = vx_old(Ip1,j)     ; u_m = vx_old(Im1,j) 
          call compute_decentered_forward_dU(u_m, vx_old(i,j), u_p, u_pp, value_dvx_dx, ONE_OVER_SIX_DELTAX)
     endif
     
      eq2_memory_dwindx_dx_fw(i,j) = c_x(i_global) * value_dwindx_dx
      eq2_memory_dvx_dx_fw(i,j) = b_x(i_global) * eq2_memory_dvx_dx_fw(i,j) + a_x(i_global) * value_dvx_dx
      

     ! compute derivative of windy, vx according to y
     if (windy_half_x_half_y >= 0) then
          u_mm = windx(i,Jm2)      ; u_m = windx(i,Jm1)      ; u_p = windx(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, windx(i,j), u_p, value_dwindx_dy, ONE_OVER_SIX_DELTAY)
          u_mm = vx_old(i,Jm2); u_m = vx_old(i,Jm1); u_p = vx_old(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, vx_old(i,j), u_p, value_dvx_dy, ONE_OVER_SIX_DELTAY)
     else 
          u_pp = windx(i,Jp2)      ; u_p = windx(i,Jp1)      ; u_m = windx(i,Jm1) 
          call compute_decentered_forward_dU(u_m, windx(i,j), u_p, u_pp, value_dwindx_dy, ONE_OVER_SIX_DELTAY)
          u_pp = vx_old(i,Jp2); u_p = vx_old(i,Jp1); u_m = vx_old(i,Jm1) 
          call compute_decentered_forward_dU(u_m, vx_old(i,j), u_p, u_pp, value_dvx_dy, ONE_OVER_SIX_DELTAY)
     endif
     
      eq2_memory_dwindx_dy_fw(i,j) = c_y(j) * value_dwindx_dy
      eq2_memory_dvx_dy_fw(i,j) = b_y(j) * eq2_memory_dvx_dy_fw(i,j) + a_y(j) * value_dvx_dy
     
     
     ! compute derivative of rho0 according to y (evaluating in vx point) 
     if (windy_half_x_half_y >= 0) then
          u_mm = rho0(Im1,Jm2)      ; u_m = rho0(Im1,Jm1)      ; u_p = rho0(Im1,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rho0(Im1,j), u_p, value_drho0_dy_prec, ONE_OVER_SIX_DELTAY)
          u_mm = rho0(i,Jm2)      ; u_m = rho0(i,Jm1)      ; u_p = rho0(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rho0(i,j), u_p, value_drho0_dy_next, ONE_OVER_SIX_DELTAY)
     else 
          u_pp = rho0(Im1,Jp2)      ; u_p = rho0(Im1,Jp1)      ; u_m = rho0(Im1,Jm1) 
          call compute_decentered_forward_dU(u_m, rho0(Im1,j), u_p, u_pp, value_drho0_dy_prec, ONE_OVER_SIX_DELTAY)
          u_pp = rho0(i,Jp2); u_p = rho0(i,Jp1); u_m = rho0(i,Jm1) 
          call compute_decentered_forward_dU(u_m, rho0(i,j), u_p, u_pp, value_drho0_dy_next, ONE_OVER_SIX_DELTAY)
     endif
     value_drho0_dy = 0.5d0 * (value_drho0_dy_next + value_drho0_dy_prec)
     
     eq2_memory_drho0_dy_fw(i,j) = c_y(j) * value_drho0_dy


     ! compute windy according to y (evaluating in vx point)
     u_mm = windy(Im1,Jm2); u_m = windy(Im1,Jm1); u_p = windy(Im1,j); u_pp = windy(Im1,Jp1)
     call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
     u_mm = windy(i,Jm2); u_m = windy(i,Jm1); u_p = windy(i,j); u_pp = windy(i,Jp1)
     call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
     value_dwindy_dy = 0.5d0 * (value_dwindy_dy_next + value_dwindy_dy_prec)

     eq2_memory_dwindy_dy_fw(i,j) = c_y(j) * value_dwindy_dy



     !
     value_dpressure_dx = value_dpressure_dx * one_over_K_x(i_global) + eq2_memory_dpressure_dx_fw(i,j)
     value_drho0_dx = value_drho0_dx * one_over_K_x(i_global) + eq2_memory_drho0_dx_fw(i,j)
     
     value_dvx_dx = value_dvx_dx * one_over_K_x(i_global) + eq2_memory_dvx_dx_fw(i,j)
     value_dwindx_dx = value_dwindx_dx * one_over_K_x(i_global) + eq2_memory_dwindx_dx_fw(i,j)
     
     value_dvx_dy = value_dvx_dy * one_over_K_y(j) + eq2_memory_dvx_dy_fw(i,j)
     value_dwindx_dy = value_dwindx_dy * one_over_K_y(j) + eq2_memory_dwindx_dy_fw(i,j)
     
     value_drho0_dy = value_drho0_dy * one_over_K_y(j) + eq2_memory_drho0_dy_fw(i,j)

     value_dwindy_dy = value_dwindy_dy * one_over_K_y(j) + eq2_memory_dwindy_dy_fw(i,j)
     
     
     
     ! intermediate computations: (v0 . nabla) v0_x ; (v' . nabla) v0_x ; (v0 . nabla) v0_x
     value_v0dwindx = windx(i,j)     * value_dwindx_dx + windy_half_x_half_y * value_dwindx_dy
     value_vdwindx  = vx_old(i,j)  * value_dwindx_dx + vy_half_x_half_y  * value_dwindx_dy
     value_v0dvx  = windx(i,j)     * value_dvx_dx  + windy_half_x_half_y * value_dvx_dy

     ! compute rho0 div v0 ; v0 div rho0
     value_rho0dv = rho0_half_x * value_dwindx_dx + rho0_half_x        * value_dwindy_dy
     value_vdrho0 = windx(i,j)   * value_drho0_dx + windy_half_x_half_y * value_drho0_dy


     
     ! update
     vx(i,j) = vx(i,j) - value_dpressure_dx * DELTAT / rho0_half_x
     !  VERSION NON CONSERVATIVE
     vx(i,j) = vx(i,j) - vx_old(i,j)  * (value_vdrho0 + value_rho0dv) * DELTAT / rho0_half_x
     vx(i,j) = vx(i,j) - value_v0dvx  * DELTAT
     
     vx(i,j) = vx(i,j) - rhop_half_x * value_v0dwindx * DELTAT / rho0_half_x
     vx(i,j) = vx(i,j) - value_vdwindx  * DELTAT
   enddo
  enddo



 ! rho0
   call MPI_SENDRECV(rho0(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,rho0(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(rho0(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,rho0(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
         
 ! p0
   call MPI_SENDRECV(p0(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,p0(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(p0(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,p0(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
         
          
 ! windy
  call MPI_SENDRECV(windy(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,windy(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(windy(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,windy(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)

  ! vy
  call MPI_SENDRECV(vy_old(1:2,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,vy_old(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(vy_old(NX_LOCAL-1:NX_LOCAL,:),number_of_values,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,vy_old(-1:0,:),number_of_values, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)

  ! windx
  call MPI_SENDRECV(windx(NX_LOCAL,:),number_of_values/2,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,windx(0,:),number_of_values/2, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)


  do j = 1,NY
    do i = 1,NX_LOCAL

      i_global = offset_i + i 

       Im1 = i-1
       Im2 = i-2
       Ip1 = i+1
       Ip2 = i+2
       Jm1 = j-1
       Jm2 = j-2
       Jp1 = j+1
       Jp2 = j+2
      
       ! modify index dependingon the boundary condition
       call get_index_boundarycondition(i,j, Im2, Im1,Ip1, Ip2, Jm2, Jm1,Jp1, Jp2) 


!     interpolate density at the right location in the staggered grid cell
      rho0_half_y  = 0.5d0 * (rho0(i,j)  + rho0(i,Jp1))
      rhop_half_y = 0.5d0 * (rhop(i,j) + rhop(i,Jp1))

      vx_half_x_half_y  = (vx_old(Ip1,Jp1)  + vx_old(i,Jp1)  +vx_old(Ip1,j)  + vx_old(i,j))  * 0.25d0
      vy_half_y         = (vy_old(i,j)      + vy_old(i,Jp1))                                 * 0.5d0
      windx_half_x_half_y = (windx(Ip1,Jp1) + windx(i,Jp1) +windx(Ip1,j) + windx(i,j)) * 0.25d0
      windy_half_y        = (windy(i,j)     + windy(i,Jm1))                        * 0.5d0

      ! compute derivative of the pressure and density according to y
      u_mm = pressure(i,Jm1); u_m = pressure(i,j); u_p = pressure(i,Jp1); u_pp = pressure(i,Jp2)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dpressure_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      u_mm = rho0(i,Jm1); u_m = rho0(i,j); u_p = rho0(i,Jp1); u_pp = rho0(i,Jp2)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drho0_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)

      eq3_memory_dpressure_dy_fw(i,j) = b_y_half(j) * eq3_memory_dpressure_dy_fw(i,j) + a_y_half(j) * value_dpressure_dy
      eq3_memory_drho0_dy_fw(i,j) = c_y_half(j) * value_drho0_dy


      ! compute derivative of windy, vy according to x
      if (windx_half_x_half_y >= 0) then
          u_mm = windy(Im2,j)      ; u_m = windy(Im1,j)      ; u_p = windy(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, windy(i,j), u_p, value_dwindy_dx, ONE_OVER_SIX_DELTAX)
          u_mm = vy_old(Im2,j); u_m = vy_old(Im1,j); u_p = vy_old(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, vy_old(i,j), u_p, value_dvy_dx, ONE_OVER_SIX_DELTAX)
      else 
          u_pp = windy(Ip2,j)      ; u_p = windy(Ip1,j)      ; u_m = windy(Im1,j) 
          call compute_decentered_forward_dU(u_m, windy(i,j), u_p, u_pp, value_dwindy_dx, ONE_OVER_SIX_DELTAX)
          u_pp = vy_old(Ip2,j); u_p = vy_old(Ip1,j); u_m = vy_old(Im1,j) 
          call compute_decentered_forward_dU(u_m, vy_old(i,j), u_p, u_pp, value_dvy_dx, ONE_OVER_SIX_DELTAX)
      endif

      eq3_memory_dwindy_dx_fw(i,j) = c_x_half(i_global) * value_dwindy_dx
      

      if (windx_half_x_half_y >= 0) then
          u_mm = rho0(Im2,j)      ; u_m = rho0(Im1,j)      ; u_p = rho0(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, rho0(i,j), u_p, value_drho0_dx_prec, ONE_OVER_SIX_DELTAX)
          u_mm = rho0(Im2,Jp1); u_m = rho0(Im1,Jp1); u_p = rho0(Ip1,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rho0(i,Jp1), u_p, value_drho0_dx_next, ONE_OVER_SIX_DELTAX)
      else 
          u_pp = rho0(Ip2,j)      ; u_p = rho0(Ip1,j)      ; u_m = rho0(Im1,j) 
          call compute_decentered_forward_dU(u_m, rho0(i,j), u_p, u_pp, value_drho0_dx_prec, ONE_OVER_SIX_DELTAX)
          u_pp = rho0(Ip2,Jp1); u_p = rho0(Ip1,Jp1); u_m = rho0(Im1,Jp1) 
          call compute_decentered_forward_dU(u_m, rho0(i,Jp1), u_p, u_pp, value_drho0_dx_next, ONE_OVER_SIX_DELTAX)
      endif
      value_drho0_dx = 0.5d0 * (value_drho0_dx_next + value_drho0_dx_prec)

      eq3_memory_drho0_dx_fw(i,j) = c_x_half(i_global) * value_drho0_dx
      

      ! compute derivative of windy, vy according to y
      if (windy(i,j) >= 0) then
          u_mm = windy(i,Jm2)      ; u_m = windy(i,Jm1)      ; u_p = windy(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, windy(i,j), u_p, value_dwindy_dy, ONE_OVER_SIX_DELTAX)
          u_mm = vy_old(i,Jm2); u_m = vy_old(i,Jm1); u_p = vy_old(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, vy_old(i,j), u_p, value_dvy_dy, ONE_OVER_SIX_DELTAX)
      else 
          u_pp = windy(i,Jp2)      ; u_p = windy(i,Jp1)      ; u_m = windy(i,Jm1) 
          call compute_decentered_forward_dU(u_m, windy(i,j), u_p, u_pp, value_dwindy_dy, ONE_OVER_SIX_DELTAX)
          u_pp = vy_old(i,Jp2); u_p = vy_old(i,Jp1); u_m = vy_old(i,Jm1) 
          call compute_decentered_forward_dU(u_m, vy_old(i,j), u_p, u_pp, value_dvy_dy, ONE_OVER_SIX_DELTAX)
      endif
      
      eq3_memory_dwindy_dy_fw(i,j) = c_y_half(j) * value_dwindy_dy
      eq3_memory_dvy_dy_fw(i,j) = b_y_half(j) * eq3_memory_dvy_dy_fw(i,j) + a_y_half(j) * value_dvy_dy
            

     ! compute derivative of windx according to x
      u_mm = windx(Im1,j); u_m = windx(i,j); u_p = windx(Ip1,j); u_pp = windx(Ip2,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
      u_mm = windx(Im1,Jp1); u_m = windx(i,Jp1); u_p = windx(Ip1,Jp1); u_pp = windx(Ip2,Jp1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
      value_dwindx_dx = 0.5d0 * (value_dwindx_dx_next + value_dwindx_dx_prec)
      
      eq3_memory_dwindx_dx_fw(i,j) = c_x_half(i_global) * value_dwindx_dx
      
      !
      value_dpressure_dy = value_dpressure_dy * one_over_K_y_half(j) + eq3_memory_dpressure_dy_fw(i,j)
      value_drho0_dy = value_drho0_dy * one_over_K_y_half(j) + eq3_memory_drho0_dy_fw(i,j)
      
      value_dvy_dy = value_dvy_dy *  one_over_K_y_half(j) + eq3_memory_dvy_dy_fw(i,j)
      value_dwindy_dy = value_dwindy_dy *  one_over_K_y_half(j) + eq3_memory_dwindy_dy_fw(i,j)
      
      value_dvy_dx = value_dvy_dx *  one_over_K_x_half(i_global) + eq3_memory_dvy_dx_fw(i,j)
      value_dwindy_dx = value_dwindy_dx *  one_over_K_x_half(i_global) + eq3_memory_dwindy_dx_fw(i,j)
      
      value_drho0_dx = value_drho0_dx * one_over_K_x_half(i_global) + eq3_memory_drho0_dx_fw(i,j)
      
      value_dwindx_dx = value_dwindx_dx * one_over_K_x_half(i_global) + eq3_memory_dwindx_dx_fw(i,j)
      
      
      ! intermediate computations: (v0 . nabla) v0_y ; (v' . nabla) v0_y ; (v0 . nabla) v0_y
      value_v0dwindy = windx_half_x_half_y * value_dwindy_dx + windy(i,j)     * value_dwindy_dy
      value_vdwindy  = vx_half_x_half_y  * value_dwindy_dx + vy_old(i,j)  * value_dwindy_dy
      value_v0dvy  = windx_half_x_half_y * value_dvy_dx  + windy(i,j)     * value_dvy_dy


      ! compute rho0 div v0 ; v0 div rho0
      value_rho0dv = rho0_half_y        * value_dwindx_dx  + rho0_half_y * value_dwindy_dy
      value_vdrho0 = windx_half_x_half_y * value_drho0_dx  + windy(i,j)   * value_drho0_dy


      ! update
      vy(i,j) = vy(i,j) - value_dpressure_dy * DELTAT / rho0_half_y
        
      ! VERSION NON CONSERVATIVE
      vy(i,j) = vy(i,j) - vy_old(i,j)  * (value_vdrho0 + value_rho0dv) * DELTAT / rho0_half_y
      vy(i,j) = vy(i,j) - value_v0dvy  * DELTAT

      vy(i,j) = vy(i,j) - rhop_half_y * value_v0dwindy * DELTAT / rho0_half_y
      vy(i,j) = vy(i,j) - value_vdwindy  * DELTAT

    enddo
  enddo

  ! Dircihlet conditions
  if (USE_PML_XMIN) then
    if (rank == 0) then
      vx(-1:1,:) = ZERO
      vy(-1:1,:) = ZERO
    else if (rank == nb_procs-1) then
      vx(NX_LOCAL:NX_LOCAL+2,:) = ZERO
      vy(NX_LOCAL:NX_LOCAL+2,:) = ZERO   
    endif 
  !else
    !if (rank == nb_procs -1) then
    !  vx(NX_LOCAL,:) = vx(0,:)
   !   vy(NX_LOCAL,:) = vy(0,:)
   ! endif
  endif
  
  if (USE_PML_YMIN) then
    vx(:,1) = ZERO
    vx(:,NY) = ZERO
    
    vy(:,1) = ZERO
    vy(:,NY) = ZERO
  else
      vx(:,NY) = vx(:,1)
      vy(:,NY) = vy(:,1)
  endif
  
! store seismograms
  do irec = 1,NREC

! beware here that the two components of the velocity vector are not defined at the same point
! in a staggered grid, and thus the two components of the velocity vector are recorded at slightly different locations,
! vy is staggered by half a grid cell along X and along Y with respect to vx
  if ((ix_rec(irec) / NX_LOCAL) == rank  ) then
    sisvx(it,irec) = vx(ix_rec(irec),iy_rec(irec))
    sisvy(it,irec) = vy(ix_rec(irec),iy_rec(irec))
    sispressure(it,irec) = pressure(ix_rec(irec),iy_rec(irec))
    sisrhop(it,irec) = rhop(ix_rec(irec),iy_rec(irec))
  endif
 
 enddo


! output information
  if ((mod(it,IT_DISPLAY) == 0 .or. it == 5) .and. (save_sismos)) then
 
    call MPI_REDUCE(maxval(abs(pressure(1:NX_LOCAL,:))), &
        pressurenorm,1,MPI_DOUBLE_PRECISION,MPI_MAX,rank_cut_plane,MPI_COMM_WORLD,code)
        
    call MPI_REDUCE(maxval(sqrt(vx(1:NX_LOCAL,:)**2 + vy(1:NX_LOCAL,:)**2)), &
        velocnorm,1,MPI_DOUBLE_PRECISION,MPI_MAX,rank_cut_plane,MPI_COMM_WORLD,code)

    if (rank==0) then
    
    ! count elapsed wall-clock time
    call date_and_time(datein,timein,zone,time_values)
      ! time_values(3): day of the month
      ! time_values(5): hour of the day
      ! time_values(6): minutes of the hour
      ! time_values(7): seconds of the minute
      ! time_values(8): milliseconds of the second
      ! this fails if we cross the end of the month
    time_end = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
               60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0
     ! elapsed time since beginning of the simulation
    tCPU = time_end - time_start
    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes

    ! print maximum of pressure and of norm of velocity
    pressurenorm = maxval(abs(pressure))
    velocnorm = maxval(sqrt(vx**2 + vy**2))
    print *,'Time step # ',it,' out of ',NSTEP
    print *,'Time: ',sngl((it-1)*DELTAT),' seconds'
    print *,'Max absolute value of pressure = ',pressurenorm
    print *,'Max norm velocity vector V (m/s) = ',velocnorm
    print *, 'Elapsed time in seconds = ',tCPU
    write(*,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    print *,'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
    print *
    
    ! check stability of the code, exit if unstable
    if (pressurenorm > STABILITY_THRESHOLD .or. velocnorm > STABILITY_THRESHOLD) stop 'code became unstable and blew up'

    !call create_color_image(pressure,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
    !                     NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,field_number)
    
    ! save the part of the seismograms that has been computed so far, so that users can monitor the progress of the simulation
    call write_seismograms(sisvx,sisvy,sispressure,sisrhop,NSTEP,NREC,DELTAT,t0)

  endif
  
  call gather_and_generate_image(it)
  
  endif

  enddo   ! end of time loop


 if (rank==0) then
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
    
  endif


  endsubroutine forwardproblem
  
  

