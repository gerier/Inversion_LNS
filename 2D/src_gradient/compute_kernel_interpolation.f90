subroutine compute_kernel()

  use parameters
  implicit none

  integer :: it, time_last_frame
  integer :: ii, jj, i_start, i_end, j_start, j_end
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: &
      !vx_half_t, vy_half_t,                     &
      !vax_half_t, vay_half_t,                   &
      value_dvx_dt, value_dvy_dt

 double precision , dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) ::  &
   rhop_half_t, rhoa_half_t, p_half_t, pa_half_t,&
   rhop_old, rhoa_old, p_old, pa_old,            &
   vx_old, vy_old


  ! INITIALISATION
  call reset_kernel()

  ! init info about time of computation
  if (rank == 0 .and. method == 2) then
    print *, "[ Start computating kernel ]"
  endif
  call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! this fails if we cross the end of the month
  time_start = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
               60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0
  
 
  ! To have an exact backward wavefiled, we use checKpointing
  ! We make a first foward simulation to save all the frames
  if (rank == 0 .and. method == 2) then
    print *, "----------- Checkpointing part -----------"
  endif
  call save_frames()
  sispressure_prior(:,:) = sispressure(:,:)
  
  ! From this point of the code, we do not need anymore to save signals recorded at receivers
  save_sismos = .False.
  
  ! Start the kernel computation 
  ! Kernel is the sum over time of correlation between adjoint and backward field
  if (rank == 0 .and. method == 2) then
    print *, "----------- Kernel part -----------"
  endif
  
  do it=1,NSTEP

    ! update old values
    ! old values are useful since the leap frog scheme and the derivative in time in
    ! the density kernel expression
    rhop_old(:,:) = rhop(:,:)
    rhoa_old(:,:) = rhoa(:,:)
    p_old(:,:)  = pressure(:,:)
    pa_old(:,:) = pa(:,:)
    vx_old(:,:) = vx(:,:)
    vy_old(:,:) = vy(:,:)

    ! to compute an exact backward simulation, we use checkpointing
    if (NSTEP-it == NSTEP-1 .or. modulo(NSTEP-it,NSTEP/NFRAMES) == NSTEP/NFRAMES-1 ) then
        ! in that case, the closest frame is the global frame, but we can compute local frames to avoid all computations
        call save_local_frames(NSTEP-it)
    else if (NSTEP-it /= NSTEP) then
        ! in that case, the closest frame is the global or local frame ; we select the good frames and compute the solution to the time (NSTEP-it)*dt    
        ! select the frame that is the closest of the desired time and who is happening before the desired time
        call load_frame(NSTEP-it, time_last_frame)
        ! Launch a forward simulation from the previously selected frame
        call forwardproblem(p0_prior,rho0_prior,windx_prior,windy_prior,time_last_frame, NSTEP-it,1)
    endif
    
    ! compute adjoint field
    call compute_adjoint(it)

    if ((mod(it,IT_DISPLAY) == 0 .or. it == 5) .and. method == 2) then
      call gather_and_generate_image(vax,vay,pa,rhoa,it,4)
    endif

    ! make interpolation (on time), to have all field at the exact same time
    rhop_half_t(:,:)  = 0.5 * (rhop(:,:)     + rhop_old(:,:))
    rhoa_half_t(:,:)  = 0.5 * (rhoa(:,:)     + rhoa_old(:,:))
    p_half_t(:,:)     = 0.5 * (pressure(:,:) + p_old(:,:))
    pa_half_t(:,:)    = 0.5 * (pa(:,:)       + pa_old(:,:))

    value_dvy_dt(:,:) = (vy(:,:) - vy_old(:,:)) * ONE_OVER_DELTAT
    value_dvx_dt(:,:) = (vx(:,:) - vx_old(:,:)) * ONE_OVER_DELTAT

    ! compute for the time it the correlation between adjoint and backward wavefield
    call compute_kernel_iter(rho0_prior, p0_prior, windx_prior, windy_prior, rhop_half_t, p_half_t, vx, vy,&
             rhoa_half_t, pa_half_t, vax, vay, value_dvx_dt, value_dvy_dt, it)

    ! write information on the pogress in the kernel computation
    if ((rank == 0) .and. (mod(it,IT_DISPLAY) == 0 .or. it == 5) .and. method == 2) then
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
    
      print *,'Time step # ',it,' out of ',NSTEP
      print *,'Time: ',sngl((it-1)*DELTAT),' seconds'
      print *, 'Elapsed time in seconds = ',tCPU
      write(*,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
      print *,'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
      print *
    endif
     
  enddo

  ! Boundary conditions
  i_start = 1 +  NPOINTS_PML + 3 
  i_end = NX_LOCAL - NPOINTS_PML - 3 - 1
  j_start = 1 + NPOINTS_PML + 3 
  j_end = NY_LOCAL - NPOINTS_PML - 3 - 1
 
  !! Left boundary
  if (USE_PML_XMIN .and. i_rank == 0) then
     do jj=1,NY_LOCAL
        K_p0(1:i_start-1,jj) = K_p0(i_start,jj)
        K_rho0(1:i_start-1,jj) = K_rho0(i_start,jj)
        K_windx(1:i_start-1,jj) = K_windx(i_start,jj)
        K_windy(1:i_start-1,jj) = K_windy(i_start,jj)
     enddo
  endif
  !! Right boundary  
  if (USE_PML_XMAX .and. i_rank == NPROC_X-1) then
     do jj=1,NY_LOCAL
         K_p0(i_end+1:NX_LOCAL,jj) = K_p0(i_end,jj)
         K_rho0(i_end+1:NX_LOCAL,jj) = K_rho0(i_end,jj)
         K_windx(i_end+1:NX_LOCAL,jj) = K_windx(i_end,jj)
         K_windy(i_end+1:NX_LOCAL,jj) = K_windy(i_end,jj)
     enddo
  endif
  !! Bottom boundary
  if (USE_PML_YMIN .and. j_rank == 0) then
     do ii=1,NX_LOCAL
        K_p0(ii,1:j_start-1) = K_p0(ii,j_start)
        K_rho0(ii,1:j_start-1) = K_rho0(ii,j_start)
        K_windx(ii,1:j_start-1) = K_windx(ii,j_start)
        K_windy(ii,1:j_start-1) = K_windy(ii,j_start)
     enddo
  endif
  !! Top boundary
  if (USE_PML_YMAX .and. j_rank == NPROC_Y-1) then 
     do ii=1,NX_LOCAL 
       K_p0(ii,j_end+1:NY_LOCAL) = K_p0(ii,j_end)
       K_rho0(ii,j_end+1:NY_LOCAL) = K_rho0(ii,j_end)
       K_windx(ii,j_end+1:NY_LOCAL) = K_windx(ii,j_end)
       K_windy(ii,j_end+1:NY_LOCAL) = K_windy(ii,j_end)
     enddo
  endif 
      
  ! Corner Left/Bottom
  if (USE_PML_XMIN .and. USE_PML_YMIN .and. rank == 0) then
    K_p0(1:i_start-1,1:j_start-1) = K_p0(i_start,j_start)
    K_rho0(1:i_start-1,1:j_start-1) = K_rho0(i_start,j_start)
    K_windx(1:i_start-1,1:j_start-1) = K_windx(i_start,j_start)
    K_windy(1:i_start-1,1:j_start-1) = K_windy(i_start,j_start) 
  endif 
  ! Corner Left/Top
  if (USE_PML_XMIN .and. USE_PML_YMAX .and. j_rank == NPROC_Y-1 .and. i_rank == 0) then 
    K_p0(1:i_start-1,j_end+1:NY_LOCAL) = K_p0(i_start,j_end)
    K_rho0(1:i_start-1,j_end+1:NY_LOCAL) = K_rho0(i_start,j_end)
    K_windx(1:i_start-1,j_end+1:NY_LOCAL) = K_windx(i_start,j_end)
    K_windy(1:i_start-1,j_end+1:NY_LOCAL) = K_windy(i_start,j_end) 
  endif
  ! Corner Right/Bottom
  if (USE_PML_XMAX .and. USE_PML_YMIN .and. i_rank == NPROC_X-1 .and. j_rank == 0) then 
    K_p0(i_end+1:NX_LOCAL,1:j_start-1) = K_p0(i_end,j_start)
    K_rho0(i_end+1:NX_LOCAL,1:j_start-1) = K_rho0(i_end,j_start)
    K_windx(i_end+1:NX_LOCAL,1:j_start-1) = K_windx(i_end,j_start)
    K_windy(i_end+1:NX_LOCAL,1:j_start-1) = K_windy(i_end,j_start) 
  endif
  ! Corner Right/Top
  if (USE_PML_XMAX .and. USE_PML_YMAX .and. rank == NPROC-1) then 
    K_p0(i_end+1:NX_LOCAL,j_end+1:NY_LOCAL) = K_p0(i_end,j_end)
    K_rho0(i_end+1:NX_LOCAL,j_end+1:NY_LOCAL) = K_rho0(i_end,j_end)
    K_windx(i_end+1:NX_LOCAL,j_end+1:NY_LOCAL) = K_windx(i_end,j_end)
    K_windy(i_end+1:NX_LOCAL,j_end+1:NY_LOCAL) = K_windy(i_end,j_end) 
  endif


  ! Save kernel information
  call gather_and_generate_image(K_windx,K_windy,K_p0,K_rho0,it,3)

endsubroutine compute_kernel


subroutine compute_kernel_iter(rho0, p0, windx, windy, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
              value_dvx_dt, value_dvy_dt, it_time)


  use parameters, only :  K_rho0, K_p0, K_windx, gamma_chimie,  &
                          DELTAT, NX, NY, NSTEP,   &
                          NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX,    &
                          NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,    &
                          ONE_OVER_SIX_DELTAX, ONE_OVER_SIX_DELTAY, &
                          DELTAX, DELTAY,                           &
                          a, t0, factor, ISOURCE, JSOURCE, source_term, &
                          type_source, wavefront,                               &
                          distance2, factor_ssf, SSF_Sigma,                     &
                          NPOINTS_PML, USE_PML_XMIN, USE_PML_YMIN, USE_PML_XMAX, USE_PML_YMAX,                  &
                          NX_LOCAL, NY_LOCAL, i_global, j_global, NPROC_X, NPROC_Y, I_RANK, J_RANK, offset_i, offset_j
  implicit none

  integer :: i,j

  integer :: it_time
  double precision :: t

  double precision :: save_source_term

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: &
     rho0, p0, windx, windy,                    &
     rhop, pressure, vx, vy,                    &
     rhoa, pa, vax, vay,                        &
     value_dvx_dt, value_dvy_dt

  double precision :: &
      value_dvx_dx,     value_dvx_dy,     &
      value_dvy_dx,     value_dvy_dy,     &
      value_dvax_dx,    value_dvax_dy,    &
      value_dvay_dx,                      &
      value_dwindx_dx,  value_dwindx_dy,  &
      value_dwindy_dx,  value_dwindy_dy,  &
      value_vax_windx_dvx_dx,             &
      value_vax_windy_dvx_dy,             &
      value_vay_windx_dvy_dx,             &
      value_vay_windy_dvy_dy,             &
      value_vax_vx_dwindx_dx,             &
      value_vax_vy_dwindx_dy,             &
      value_vay_vx_dwindy_dx,             &
      value_vay_vy_dwindy_dy

  double precision ::   &
      value_drhoarhop_dx, &
      value_drhop_dx,     value_drhop_dy,     &
      value_drho0_dx,     value_drho0_dy,     &
      value_drhoa_dx,     value_drhoa_dy,     &
      value_dgammapap_dx,                     &
      value_dp_dx,                            &
      value_dpa_dx,       value_dpa_dy,       &
      value_vx_dvax_dx,                       &
      value_vy_dvay_dx,                       &
      value_vax_dvx_dx,                    &
      value_vay_dvy_dx,                    &
      value_dv,                               &
      value_dwind,                            &
      value_windx_drhop_dx,                   &
      value_windy_drhop_dy,                   &
      value_vx_drho0_dx,                      &
      value_vy_drho0_dy,                      &
      value_wind_dvax,                        &
      value_v_dvax
      

  double precision :: &
    rho0_half_x,                  &
    rhoa_half_x,                  &
    rhop_half_x,                  &
    pa_half_x,                    &
    windx_half_x,        windy_half_y,        &
    vax_half_x,          vay_half_y,          &
    vx_half_x,           vy_half_y,           &
    windy_half_x_half_y, &
    vay_half_x_half_y,   &
    vy_half_x_half_y,    &
    vax_dvx_dt_half_x, vay_dvy_dt_half_y,     &
    value_dvxvax_dx,   value_dvxvax_dy,       &
    value_dvyvay_dx,   value_dvyvay_dy,       &
    value_dvy_dy_half_x_half_y,               &
    value_dwindy_dy_half_x_half_y
    

   double precision :: & 
       value_dwindx_dy_next, value_dwindx_dy_prec,          &
       value_dwindy_dx_next, value_dwindy_dx_prec,          &
       value_dwindy_dy_next, value_dwindy_dy_prec,          &
       value_dvx_dy_next,    value_dvx_dy_prec,             &
       value_dvy_dx_next,    value_dvy_dx_prec,             &
       value_dvy_dy_next,    value_dvy_dy_prec,             &
       value_dvay_dx_next,   value_dvay_dx_prec,            &
       value_dvyvay_dx_next, value_dvyvay_dx_prec,          &
       value_dvxvax_dy_next, value_dvxvax_dy_prec,          &
       value_drho0_dy_next,  value_drho0_dy_prec,           &
       value_drhop_dy_next,  value_drhop_dy_prec           
       
  double precision, dimension(-1:NX_LOCAL+2, -1:NY_LOCAL+2) :: &
      vax_dvx_dt, vay_dvy_dt

  double precision  :: u_mm, u_m, u_p, u_pp
  integer :: Ip1,Im1, Ip2, Im2, Jp1,Jm1, Jp2, Jm2
  
  integer :: i_start, i_end, j_start, j_end
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! init the source 
  t = dble(NSTEP-it_time+1-0.5)*DELTAT
  save_source_term = - 2 * a* (t-t0) *  factor * exp(- a*(t-t0)*(t-t0))
  !!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Index for defining Boundaries 
  !!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Bottom boundary
  j_start = 1
  if (USE_PML_YMAX .and. USE_PML_YMIN .and. j_rank == 0) then
    j_start = NPOINTS_PML + 3 + 1
  else if (USE_PML_YMAX .and. .not.USE_PML_YMIN .and. j_rank == 0) then
    j_start = 3
  endif
  ! Top boundary
  j_end = NY_LOCAL
  if (USE_PML_YMAX .and.j_rank == NPROC_Y-1) then
    j_end = NY_LOCAL - NPOINTS_PML - 3  - 1
  endif
  ! Left boundary
  i_start = 1
  if (USE_PML_XMIN .and. i_rank == 0) then
    i_start = NPOINTS_PML + 1 + 3 
  endif
  ! Right boundary
  i_end = NX_LOCAL  
  if (USE_PML_XMAX .and. i_rank == NPROC_X -1) then
    i_end = NX_LOCAL - NPOINTS_PML - 3 - 1 
  endif
  


  ! Preliminary computations
  call send_receive_rightleft(vax)
  call send_receive_topbottom(vay)
  call send_receive_topbottom(value_dvy_dt)
  call send_receive_rightleft(value_dvx_dt)
  
  do j=j_start-1,j_end+1
   do i=i_start-1,i_end+1

   vax_dvx_dt(i,j) = vax(i,j) * value_dvx_dt(i,j)
   vay_dvy_dt(i,j) = vay(i,j) * value_dvy_dt(i,j)

   enddo
  enddo
 
  ! Information from and for other processus
  call send_receive_rightleft(rhoa)
  call send_receive_rightleft(pa)
  call send_receive_rightleft(vax)
  call send_receive_rightleft(vay)

  call send_receive_rightleft(rhop)
  call send_receive_rightleft(pressure)
  call send_receive_rightleft(vx)
  call send_receive_rightleft(vy)
  
  call send_receive_left(vax_dvx_dt)
  call send_receive_top(vay_dvy_dt)
  
  call send_receive_topbottom(rhoa)
  call send_receive_topbottom(pa)
  call send_receive_topbottom(vax)
  call send_receive_topbottom(vay)

  call send_receive_topbottom(rhop)
  call send_receive_topbottom(pressure)
  call send_receive_topbottom(vx)
  call send_receive_topbottom(vy)
  
  call send_receive_leftbottom(vax)
  call send_receive_lefttop(vax)
  call send_receive_leftbottom(vx)
  call send_receive_lefttop(vx)
  
  call send_receive_righttop(vay)
  call send_receive_lefttop(vay)
  call send_receive_righttop(vy)
  call send_receive_lefttop(vy)
  
  call send_receive_righttop(rhop)
  call send_receive_rightbottom(rhop)
  
  call send_receive_rightbottom(vx)
  call send_receive_rightbottom(vy)
   
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Kernel of wind x
  !!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Treatment different of the bottom boundary in case of a reflective ground 
  if (.not. USE_PML_YMIN .and. j_rank == 0) then  
     call derivative_first_raw_kernel_kvx(rho0, p0, windx, windy, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
               value_dvx_dt, value_dvy_dt, it_time, i_start, i_end)
  endif    


  do j=j_start,j_end
   do i =i_start,i_end

      ! Index
      i_global = i + offset_i
      j_global = j + offset_j
      
      Im1 = i-1
      Im2 = i-2
      Ip1 = i+1
      Ip2 = i+2
      Jm1 = j-1
      Jm2 = j-2
      Jp1 = j+1
      Jp2 = j+2

    
      ! Interpolations 
      rhoa_half_x        = 0.5d0 * (rhoa(i,j) + rhoa(Im1,j))
      rho0_half_x        = 0.5d0 * (rho0(i,j) + rho0(Im1,j))
      rhop_half_x        = 0.5d0 * (rhop(i,j) + rhop(Im1,j))
      pa_half_x          = 0.5d0 * (pa(i,j)   + pa(Im1,j))

      vay_half_x_half_y = 0.25d0 * (vay(i,j) + vay(Im1,j)+ vay(i,Jm1) + vay(Im1,Jm1))
      windy_half_x_half_y  = 0.25d0 * (windy(i,j) + windy(Im1,j)+ windy(i,Jm1) + windy(Im1,Jm1))
      vy_half_x_half_y  = 0.25d0 * (vy(i,j) + vy(Im1,j)+ vy(i,Jm1) + vy(Im1,Jm1))


      ! compute derivative of rhoprhoa, gammapap and pressure according to x
      u_mm = pressure(Im2,j); u_m = pressure(Im1,j); u_p = pressure(i,j); u_pp = pressure(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      u_mm = rhop(Im2,j); u_m = rhop(Im1,j); u_p = rhop(i,j); u_pp = rhop(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhop_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      u_mm = rho0(Im2,j); u_m = rho0(Im1,j); u_p = rho0(i,j); u_pp = rho0(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drho0_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      u_mm = rhop(Im2,j) * rhoa(Im2,j); u_m = rhop(Im1,j) * rhoa(Im1,j)
      u_p  = rhop(i,j) * rhoa(i,j); u_pp = rhop(Ip1,j) * rhoa(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhoarhop_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      u_mm = pressure(Im2,j) * pa(Im2,j) * gamma_chimie(Im2,j); u_m = pressure(Im1,j) * pa(Im1,j) * gamma_chimie(Im1,j)
      u_p  = pressure(i,j) * pa(i,j) * gamma_chimie(i,j)   ; u_pp = pressure(Ip1,j) * pa(Ip1,j) * gamma_chimie(IP1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dgammapap_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)


      ! compute derivative of vax and windx according to x
      if (windx(i,j) >= 0) then
          u_mm = windx(Im2,j)      ; u_m = windx(Im1,j)      ; u_p = windx(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, windx(i,j), u_p, value_dwindx_dx, ONE_OVER_SIX_DELTAX)
          u_mm = vax(Im2,j); u_m = vax(Im1,j); u_p = vax(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, vax(i,j), u_p, value_dvax_dx, ONE_OVER_SIX_DELTAX)
          u_mm = vx(Im2,j); u_m = vx(Im1,j); u_p = vx(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, vx(i,j), u_p, value_dvx_dx, ONE_OVER_SIX_DELTAX)
      else 
          u_pp = windx(Ip2,j)    ; u_p = windx(Ip1,j)    ; u_m = windx(Im1,j) 
          call compute_decentered_forward_dU(u_m, windx(i,j), u_p, u_pp, value_dwindx_dx, ONE_OVER_SIX_DELTAX)
          u_pp = vax(Ip2,j)     ; u_p = vax(Ip1,j)     ; u_m = vax(Im1,j) 
          call compute_decentered_forward_dU(u_m, vax(i,j), u_p, u_pp, value_dvax_dx, ONE_OVER_SIX_DELTAX)
          u_pp = vx(Ip2,j)     ; u_p = vx(Ip1,j)     ; u_m = vx(Im1,j) 
          call compute_decentered_forward_dU(u_m, vx(i,j), u_p, u_pp, value_dvx_dx, ONE_OVER_SIX_DELTAX)
      endif
       
        
      ! compute derivative of vay and windy according to x
      u_mm = vay(Im2,Jm1); u_m = vay(Im1,Jm1); u_p = vay(i,Jm1); u_pp = vay(Ip1,Jm1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvay_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      u_mm = vay(Im2,j); u_m = vay(Im1,j); u_p = vay(i,j); u_pp = vay(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvay_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      value_dvay_dx = 0.5d0 * ( value_dvay_dx_prec + value_dvay_dx_next)

      u_mm = windy(Im2,Jm1); u_m = windy(Im1,Jm1); u_p = windy(i,Jm1); u_pp = windy(Ip1,Jm1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      u_mm = windy(Im2,j); u_m = windy(Im1,j); u_p = windy(i,j); u_pp = windy(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      value_dwindy_dx = 0.5d0 * ( value_dwindy_dx_prec + value_dwindy_dx_next)
        
      u_mm = vy(Im2,Jm1); u_m = vy(Im1,Jm1); u_p = vy(i,Jm1); u_pp = vy(Ip1,Jm1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      u_mm = vy(Im2,j); u_m = vy(Im1,j); u_p = vy(i,j); u_pp = vy(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      value_dvy_dx = 0.5d0 * ( value_dvy_dx_prec + value_dvy_dx_next)
        
          
      ! compute derivative of windy, vy according to y
      u_mm = windy(Im1,Jm2); u_m = windy(Im1,Jm1); u_p = windy(Im1,j); u_pp = windy(Im1,Jp1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      u_mm = windy(i,Jm2); u_m = windy(i,Jm1); u_p = windy(i,j); u_pp = windy(i,Jp1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      value_dwindy_dy_half_x_half_y = 0.5d0 * ( value_dwindy_dy_prec + value_dwindy_dy_next)
    
      u_mm = vy(Im1,Jm2); u_m = vy(Im1,Jm1); u_p = vy(Im1,j); u_pp = vy(Im1,Jp1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      u_mm = vy(i,Jm2); u_m = vy(i,Jm1); u_p = vy(i,j); u_pp = vy(i,Jp1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      value_dvy_dy_half_x_half_y = 0.5d0 * ( value_dvy_dy_prec + value_dvy_dy_next)


      ! compute derivative of vax according to y
      if (windy(i,j) >= 0 .and. .not.(.not. USE_PML_YMIN .and. j_rank==0 .and. j == 3)) then
          u_mm = vax(i,Jm2); u_m = vax(i,Jm1); u_p = vax(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, vax(i,j), u_p, value_dvax_dy, ONE_OVER_SIX_DELTAX)
      else 
          u_pp = vax(i,Jp2)     ; u_p = vax(i,Jp1)     ; u_m = vax(i,Jm1) 
          call compute_decentered_forward_dU(u_m, vax(i,j), u_p, u_pp, value_dvax_dy, ONE_OVER_SIX_DELTAX)
      endif
    
    
      ! compute derivative of rhop, rho0 according to y
      if (windy(i,j) >= 0) then
          u_mm = rhop(Im1,Jm2); u_m = rhop(Im1,Jm1); u_p = rhop(Im1,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rhop(Im1,j), u_p, value_drhop_dy_prec, ONE_OVER_SIX_DELTAY)
          u_mm = rhop(i,Jm2); u_m = rhop(i,Jm1); u_p = rhop(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rhop(i,j), u_p, value_drhop_dy_next, ONE_OVER_SIX_DELTAY)
          
          u_mm = rho0(Im1,Jm2); u_m = rho0(Im1,Jm1); u_p = rho0(Im1,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rho0(Im1,j), u_p, value_drho0_dy_prec, ONE_OVER_SIX_DELTAY)
          u_mm = rho0(i,Jm2); u_m = rho0(i,Jm1); u_p = rho0(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rho0(i,j), u_p, value_drho0_dy_next, ONE_OVER_SIX_DELTAY)
      else 
          u_pp = rhop(Im1,Jp2)     ; u_p = rhop(Im1,Jp1)     ; u_m = rhop(Im1,Jm1) 
          call compute_decentered_forward_dU(u_m, rhop(Im1,j), u_p, u_pp, value_drhop_dy_prec, ONE_OVER_SIX_DELTAY)
          u_pp = rhop(i,Jp2)     ; u_p = rhop(i,Jp1)     ; u_m = rhop(i,Jm1) 
          call compute_decentered_forward_dU(u_m, rhop(i,j), u_p, u_pp, value_drhop_dy_next, ONE_OVER_SIX_DELTAY)
          
          u_pp = rho0(Im1,Jp2)     ; u_p = rho0(Im1,Jp1)     ; u_m = rho0(Im1,Jm1) 
          call compute_decentered_forward_dU(u_m, rho0(Im1,j), u_p, u_pp, value_drho0_dy_prec, ONE_OVER_SIX_DELTAY)
          u_pp = rho0(i,Jp2)     ; u_p = rho0(i,Jp1)     ; u_m = rho0(i,Jm1) 
          call compute_decentered_forward_dU(u_m, rho0(i,j), u_p, u_pp, value_drho0_dy_next, ONE_OVER_SIX_DELTAY)
      endif
      value_drhop_dy = 0.5d0 * (value_drhop_dy_prec+value_drhop_dy_next) 
      value_drho0_dy = 0.5d0 * (value_drho0_dy_prec+value_drho0_dy_next) 
    
    
      ! intermediate computations
      value_vx_dvax_dx = vx(i,j) * value_dvax_dx
      value_vy_dvay_dx = vy_half_x_half_y * value_dvay_dx

      value_vax_dvx_dx = vax(i,j) * value_dvx_dx
      value_vay_dvy_dx = vay_half_x_half_y * value_dvy_dx

      value_dwind = value_dwindx_dx + value_dwindy_dy_half_x_half_y
      value_dv = value_dvx_dx + value_dvy_dy_half_x_half_y
    
      value_windx_drhop_dx = windx(i,j) * value_drhop_dx
      value_vx_drho0_dx = vx(i,j) * value_drho0_dx
    
      value_windy_drhop_dy = windy_half_x_half_y * value_drhop_dy
      value_vy_drho0_dy = vy_half_x_half_y * value_drho0_dy

      value_wind_dvax = windx(i,j) * value_dvax_dx + windy_half_x_half_y * value_dvax_dy
      value_v_dvax = vx(i,j) * value_dvax_dx + vy_half_x_half_y * value_dvax_dy


      ! Kernel of x-velocity
      K_windx(i,j) = K_windx(i,j) + (rhoa_half_x     * value_drhop_dx)         * DELTAT
      K_windx(i,j) = K_windx(i,j) + rho0_half_x      * (value_vax_dvx_dx + value_vay_dvy_dx)  * DELTAT
      K_windx(i,j) = K_windx(i,j) - vax(i,j) * rho0_half_x * value_dvy_dy_half_x_half_y *DELTAT
      K_windx(i,j) = K_windx(i,j) - vy_half_x_half_y * rho0_half_x * value_dvax_dy *DELTAT
      K_windx(i,j) = K_windx(i,j) + (pa_half_x       * value_dp_dx)  * DELTAT
    
   enddo
  enddo



  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Kernel of rho0 and p0
  !!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Treatment different of the bottom boundary in case of a reflective ground 
  if (.not. USE_PML_YMIN .and. j_rank == 0) then  
    call derivative_first_raw_kernel_kp0krho0(rho0, p0, windx, windy, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
              value_dvx_dt, value_dvy_dt, it_time, i_start, i_end, save_source_term)
  endif              
    
              
  do j=j_start,j_end
   do i =i_start,i_end

      ! Index
      i_global = i + offset_i
      j_global = j + offset_j
     
      Im1 = i-1
      Im2 = i-2
      Ip1 = i+1
      Ip2 = i+2
      Jm1 = j-1
      Jm2 = j-2
      Jp1 = j+1
      Jp2 = j+2


      ! Interpolations
      vx_half_x          = 0.5d0 * (vx(i,j)  + vx(Ip1,j))
      vax_half_x         = 0.5d0 * (vax(i,j) + vax(Ip1,j))
      windx_half_x       = 0.5d0 * (windx(i,j) + windx(Ip1,j))
      vax_dvx_dt_half_x  = 0.5d0 * (vax_dvx_dt(i,j) + vax_dvx_dt(Ip1,j))

      vy_half_y          = 0.5d0 * (vy(i,j)  + vy(i,Jm1))
      vay_half_y         = 0.5d0 * (vay(i,j) + vay(i,Jm1))
      windy_half_y       = 0.5d0 * (windy(i,j) + windy(i,Jm1))
      vay_dvy_dt_half_y  = 0.5d0 * (vay_dvy_dt(i,j) + vay_dvy_dt(i,Jm1))

	
      ! Source function 
      source_term = 0.d0
      if (type_source == 1) then      
         if ((wavefront == 1 .and. i_global == ISOURCE ) .or. (wavefront == 2 .and. j_global == JSOURCE)) then
           factor_ssf = 1.d0
          source_term = save_source_term
         endif
      
      elseif (type_source ==2) then
        if (i_global == ISOURCE .and. j_global == JSOURCE) then
          factor_ssf = 1.d0
          source_term = save_source_term
        endif
      
      elseif (type_source == 3) then
        distance2 = ((i_global - Isource) * DELTAX)**2 + ((j_global - Jsource) * DELTAY)**2
        factor_ssf = exp( - distance2 / SSF_Sigma**2 )
        source_term = save_source_term
      endif
    
    
      ! compute derivative of windx, vx according to x 
      u_mm = windx(Im1,j); u_m = windx(i,j); u_p = windx(Ip1,j); u_pp = windx(Ip2,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    
      u_mm = vx(Im1,j); u_m = vx(i,j); u_p = vx(Ip1,j); u_pp = vx(Ip2,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
       
        
      ! compute derivative of windy, vy according to y
      u_mm = vy(i,Jm2); u_m = vy(i,Jm1); u_p = vy(i,j); u_pp = vy(i,Jp1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      
      u_mm = windy(i,Jm2); u_m = windy(i,Jm1); u_p = windy(i,j); u_pp = windy(i,Jp1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
       
       
      ! compute derivative of pa,rhoa according to x 
      if (windx_half_x >= 0) then
          u_mm = pa(Im2,j)      ; u_m = pa(Im1,j)      ; u_p = pa(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, pa(i,j), u_p, value_dpa_dx, ONE_OVER_SIX_DELTAX)
          u_mm = rhoa(Im2,j)      ; u_m = rhoa(Im1,j)      ; u_p = rhoa(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, rhoa(i,j), u_p, value_drhoa_dx, ONE_OVER_SIX_DELTAX)
          
      else 
          u_pp = pa(Ip2,j)      ; u_p = pa(Ip1,j)      ; u_m = pa(Im1,j) 
          call compute_decentered_forward_dU(u_m, pa(i,j), u_p, u_pp, value_dpa_dx, ONE_OVER_SIX_DELTAX)
          u_pp = rhoa(Ip2,j)      ; u_p = rhoa(Ip1,j)      ; u_m = rhoa(Im1,j) 
          call compute_decentered_forward_dU(u_m, rhoa(i,j), u_p, u_pp, value_drhoa_dx, ONE_OVER_SIX_DELTAX)
      endif
    
    
      ! compute derivative of windy, vy according to x 
      if (windx_half_x >= 0) then
          u_mm = windy(Im2,Jm1)      ; u_m = windy(Im1,Jm1)      ; u_p = windy(Ip1,Jm1) 
          call compute_decentered_backward_dU(u_mm, u_m, windy(i,Jm1), u_p, value_dwindy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_mm = windy(Im2,j)      ; u_m = windy(Im1,j)      ; u_p = windy(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, windy(i,j), u_p, value_dwindy_dx_next, ONE_OVER_SIX_DELTAX)
          
          u_mm = vy(Im2,Jm1)      ; u_m = vy(Im1,Jm1)      ; u_p = vy(Ip1,Jm1) 
          call compute_decentered_backward_dU(u_mm, u_m, vy(i,Jm1), u_p, value_dvy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_mm = vy(Im2,j)      ; u_m = vy(Im1,j)      ; u_p = vy(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, vy(i,j), u_p, value_dvy_dx_next, ONE_OVER_SIX_DELTAX)
          
      else 
          u_pp = windy(Ip2,Jm1)      ; u_p = windy(Ip1,Jm1)      ; u_m = windy(Im1,Jm1) 
          call compute_decentered_forward_dU(u_m, windy(i,Jm1), u_p, u_pp, value_dwindy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_pp = windy(Ip2,j)      ; u_p = windy(Ip1,j)      ; u_m = windy(Im1,j) 
          call compute_decentered_forward_dU(u_m, windy(i,j), u_p, u_pp, value_dwindy_dx_next, ONE_OVER_SIX_DELTAX)
          
          u_pp = vy(Ip2,Jm1)      ; u_p = vy(Ip1,Jm1)      ; u_m = vy(Im1,Jm1) 
          call compute_decentered_forward_dU(u_m, windy(i,Jm1), u_p, u_pp, value_dvy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_pp = vy(Ip2,j)      ; u_p = vy(Ip1,j)      ; u_m = vy(Im1,j) 
          call compute_decentered_forward_dU(u_m, vy(i,j), u_p, u_pp, value_dvy_dx_next, ONE_OVER_SIX_DELTAX)
      endif
      value_dwindy_dx = 0.5d0 * (value_dwindy_dx_prec + value_dwindy_dx_next)
      value_dvy_dx = 0.5d0 * (value_dvy_dx_prec + value_dvy_dx_next)
    
    
      ! compute derivative of rhoa, pa according to y
      if (windy_half_y >= 0 .and. .not.(.not. USE_PML_YMIN .and. j_rank==0 .and. j == 3)) then
          u_mm = pa(i,Jm2)      ; u_m = pa(i,Jm1)      ; u_p = pa(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, pa(i,j), u_p, value_dpa_dy, ONE_OVER_SIX_DELTAY)
          u_mm = rhoa(i,Jm2)      ; u_m = rhoa(i,Jm1)      ; u_p = rhoa(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rhoa(i,j), u_p, value_drhoa_dy, ONE_OVER_SIX_DELTAY)
      else 
          u_pp = pa(i,Jp2)      ; u_p = pa(i,Jp1)      ; u_m = pa(i,Jm1) 
          call compute_decentered_forward_dU(u_m, pa(i,j), u_p, u_pp, value_dpa_dy, ONE_OVER_SIX_DELTAY)
          u_pp = rhoa(i,Jp2)      ; u_p = rhoa(i,Jp1)      ; u_m = rhoa(i,Jm1)         
          call compute_decentered_forward_dU(u_m, rhoa(i,j), u_p, u_pp, value_drhoa_dy, ONE_OVER_SIX_DELTAY)
      endif
    
    
     ! compute derivative of windx, vx according to y
     if (windy_half_y >= 0) then
          u_mm = windx(i,Jm2)      ; u_m = windx(i,Jm1)      ; u_p = windx(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, windx(i,j), u_p, value_dwindx_dy_prec, ONE_OVER_SIX_DELTAY)
          u_mm = windx(Ip1,Jm2)      ; u_m = windx(Ip1,Jm1)      ; u_p = windx(Ip1,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, windx(Ip1,j), u_p, value_dwindx_dy_next, ONE_OVER_SIX_DELTAY)
          
          u_mm = vx(i,Jm2)      ; u_m = vx(i,Jm1)      ; u_p = vx(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, vx(i,j), u_p, value_dvx_dy_prec, ONE_OVER_SIX_DELTAY)
          u_mm = vx(Ip1,Jm2)      ; u_m = vx(Ip1,Jm1)      ; u_p = vx(Ip1,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, vx(Ip1,j), u_p, value_dvx_dy_next, ONE_OVER_SIX_DELTAY)
          
       else 
          u_pp = windx(i,Jp2)      ; u_p = windx(i,Jp1)      ; u_m = windx(i,Jm1) 
          call compute_decentered_forward_dU(u_m, windx(i,j), u_p, u_pp, value_dwindx_dy_prec, ONE_OVER_SIX_DELTAY)
          u_pp = windx(Ip1,Jp2)      ; u_p = windx(Ip1,Jp1)      ; u_m = windx(Ip1,Jm1) 
          call compute_decentered_forward_dU(u_m, windx(Ip1,j), u_p, u_pp, value_dwindx_dy_next, ONE_OVER_SIX_DELTAY)
          
          u_pp = vx(i,Jp2)      ; u_p = vx(i,Jp1)      ; u_m = vx(i,Jm1)                 ! SOLENE has changed Jm into Jp
          call compute_decentered_forward_dU(u_m, vx(i,j), u_p, u_pp, value_dvx_dy_prec, ONE_OVER_SIX_DELTAY)
          u_pp = vx(Ip1,Jp2)      ; u_p = vx(Ip1,Jp1)      ; u_m = vx(Ip1,Jm1)           ! SOLENE has changed Jm into Jp
          call compute_decentered_forward_dU(u_m, vx(Ip1,j), u_p, u_pp, value_dvx_dy_next, ONE_OVER_SIX_DELTAY)

      endif
      value_dwindx_dy = 0.5d0 * (value_dwindx_dy_prec + value_dwindx_dy_next)
      value_dvx_dy = 0.5d0 * (value_dvx_dy_prec + value_dvx_dy_next)
    
    
      ! compute derivatives of vpx * vax according to x
      u_mm = vx(Im1,j)*vax(Im1,j); u_m = vx(i,j)*vax(i,j); u_p = vx(IP1,j)*vax(Ip1,j); u_pp = vx(Ip2,j)*vax(Ip2,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvxvax_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    
    
      ! compute derivatives of vpy * vay according to y
      u_mm = vy(i,Jm2)*vay(i,Jm2); u_m = vy(i,Jm1)*vay(i,Jm1)
      u_p = vy(i,j)*vay(i,j); u_pp = vy(i,Jp1)*vay(i,Jp1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvyvay_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      
      
      !  compute derivative of vy*vay according to x
      if (windx_half_x >= 0) then
          u_mm = vy(Im2,Jm1)*vay(Im2,Jm1)      ; u_m = vy(Im1,Jm1)*vay(Im1,Jm1)
          u_p = vy(Ip1,Jm1)*vay(Ip1,Jm1) 
          call compute_decentered_backward_dU(u_mm, u_m, vy(i,Jm1)*vay(i,Jm1), u_p, value_dvyvay_dx_prec, ONE_OVER_SIX_DELTAX)
          u_mm = vy(Im2,j)*vay(Im2,j)      ; u_m = vy(Im1,j)*vay(Im1,j)
          u_p = vy(Ip1,j)*vay(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, vy(i,j)*vay(i,j), u_p, value_dvyvay_dx_next, ONE_OVER_SIX_DELTAX) 
      else 
          u_pp = vy(Ip2,Jm1)*vay(Ip2,Jm1)      ; u_p = vy(Ip1,Jm1)*vay(Ip1,Jm1)
          u_m = vy(Im1,Jm1)*vay(Im1,Jm1) 
          call compute_decentered_forward_dU(u_m, vy(i,Jm1)*vay(i,Jm1), u_p, u_pp, value_dvyvay_dx_prec, ONE_OVER_SIX_DELTAX)
          u_pp = vy(Ip2,j)*vay(Ip2,j)      ; u_p = vy(Ip1,j)*vay(Ip1,j)
          u_m = vy(Im1,j)*vay(Im1,j) 
          call compute_decentered_forward_dU(u_m, vy(i,j)*vay(i,j), u_p, u_pp, value_dvyvay_dx_next, ONE_OVER_SIX_DELTAX)
      endif
      value_dvyvay_dx = 0.5d0 * (value_dvyvay_dx_prec + value_dvyvay_dx_next)
    
    
      !  compute derivative of vx*vax according to y
      if (windy_half_y >= 0 .and. .not.(.not. USE_PML_YMIN .and. j_rank==0 .and. j == 3)) then
          u_mm = vx(i,Jm2)*vax(i,Jm2)      ; u_m = vx(i,Jm1)*vax(i,Jm1)      ; u_p = vx(i,Jp1)*vax(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, vx(i,j)*vax(i,j), u_p, value_dvxvax_dy_prec, ONE_OVER_SIX_DELTAY)
          u_mm = vx(Ip1,Jm2)*vax(Ip1,Jm2)      ; u_m = vx(Ip1,Jm1)*vax(Ip1,Jm1)      ; u_p = vx(Ip1,Jp1)*vax(Ip1,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, vx(Ip1,j)*vax(Ip1,j), u_p, value_dvxvax_dy_next, ONE_OVER_SIX_DELTAY)
      else 
          u_pp = vx(i,Jp2)*vax(i,Jp2)      ; u_p = vx(i,Jp1)*vax(i,Jp1)      ; u_m = vx(i,Jm1)*vax(i,Jm1) 
          call compute_decentered_forward_dU(u_m, vx(i,j)*vax(i,j), u_p, u_pp, value_dvxvax_dy_prec, ONE_OVER_SIX_DELTAY)
          u_pp = vx(Ip1,Jp2)*vax(Ip1,Jp2)      ; u_p = vx(Ip1,Jp1) *vax(Ip1,Jp1)      ; u_m = vx(Ip1,Jm1)*vax(Ip1,Jm1) 
          call compute_decentered_forward_dU(u_m, vx(i,j)*vax(Ip1,j), u_p, u_pp, value_dvxvax_dy_next, ONE_OVER_SIX_DELTAY)
      endif
      value_dvxvax_dy = 0.5d0 * (value_dvxvax_dy_prec + value_dvxvax_dy_next)


      ! intermediate computations
      value_vax_windx_dvx_dx = vax_half_x * windx_half_x * value_dvx_dx
      value_vax_windy_dvx_dy = vax_half_x * windy_half_y * value_dvx_dy
      value_vay_windx_dvy_dx = vay_half_y * windx_half_x * value_dvy_dx
      value_vay_windy_dvy_dy = vay_half_y * windy_half_y * value_dvy_dy
     
      value_vax_vx_dwindx_dx = vax_half_x * vx_half_x * value_dwindx_dx
      value_vax_vy_dwindx_dy = vax_half_x * vy_half_y * value_dwindx_dy
      value_vay_vx_dwindy_dx = vay_half_y * vx_half_x * value_dwindy_dx
      value_vay_vy_dwindy_dy = vay_half_y * vy_half_y * value_dwindy_dy


      ! Kernel of density
      K_rho0(i,j) = K_rho0(i,j) + (rhoa(i,j) * value_dvx_dx - vy_half_y * value_drhoa_dy) * DELTAT
      K_rho0(i,j) = K_rho0(i,j) - (vax_dvx_dt_half_x + vay_dvy_dt_half_y) * DELTAT
      K_rho0(i,j) = K_rho0(i,j) + (value_vax_windx_dvx_dx + value_vay_windx_dvy_dx) * DELTAT
      K_rho0(i,j) = K_rho0(i,j) + (value_vax_vy_dwindx_dy) * DELTAT
      K_rho0(i,j) = K_rho0(i,j) + pa(i,j) * factor_ssf * source_term * (gamma_chimie(i,j) * p0(i,j) / rho0(i,j)**2) * DELTAT


      ! Kernel of pressure
      K_p0(i,j) = K_p0(i,j) - (vy_half_y * value_dpa_dy) * DELTAT
      K_p0(i,j) = K_p0(i,j) + gamma_chimie(i,j) * pa(i,j) * (value_dvx_dx + value_dvy_dy)        * DELTAT 
      K_p0(i,j) = K_p0(i,j) - pa(i,j) * value_dvy_dy                                             * DELTAT 
      K_p0(i,j) = K_p0(i,j) - pa(i,j) * factor_ssf * source_term * gamma_chimie(i,j) / rho0(i,j) * DELTAT

   enddo
  enddo

endsubroutine compute_kernel_iter

