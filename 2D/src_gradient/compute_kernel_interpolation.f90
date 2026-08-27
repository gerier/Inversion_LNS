subroutine compute_kernel()
!==============================================================================
! Compute sensitivity kernels using the adjoint-state method.
!
! Description:
!   This routine reconstructs the forward wavefield backward in time using
!   checkpointing and computes the adjoint wavefield driven by the
!   data residuals.
!
!   At each time step, time correlations between forward and adjoint fields are
!   integrated in time to update the sensitivity kernels.
!
! Sensitivity kernels are computed for:
!
!   - background pressure       : K_p0
!   - background density        : K_rho0
!   - background wind components: K_windx
!
! The computation is parallelized using MPI domain decomposition.
! Kernel values inside PML regions are extrapolated.
!
! Main steps:
!   1. Initialize kernel arrays.
!   2. Perform a forward simulation and store wavefield checkpoints.
!   3. Compute adjoint sources.
!   4. Reconstruct forward wavefields backward in time.
!   5. Accumulate kernel contributions at each time step.
!   6. Apply boundary treatment and save results.
!
!==============================================================================
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


  ! Reset variables required for kernel computations.
  call reset_kernel()

  ! Initialize timing information
  if (rank == 0 .and. method == 2) then
    print *, "[ Start computing kernel ]"
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


  ! To reconstruct the exact backward wavefield, checkpointing is used.
  ! Perform a first forward simulation and store wavefield checkpoints
  if (rank == 0 .and. method == 2) then
    print *, "----------- Checkpointing part -----------"
  endif
  call save_frames()
  sispressure_prior(:,:) = sispressure(:,:)

  ! From this point onward, receiver signals no longer need to be stored.
  ! Signals have been saved during the checkpointing phase.
  save_sismos = .FALSE.

  ! Prepare adjoint source
  call prepare_adjoint_source()

  ! The sensitivity kernel is obtained by time integration of
  ! the correlation between forward perturbation fields and adjoint fields.
  
  ! Start the kernel computation
  if (rank == 0 .and. method == 2) then
    print *, "----------- Kernel part -----------"
  endif

  do it=1,NSTEP

    ! Store previous time-step values.
    ! Previous time levels are required by the leapfrog scheme and by temporal
    ! derivatives involved in the kernel formulation.
    rhop_old(:,:) = rhop(:,:)
    rhoa_old(:,:) = rhoa(:,:)
    p_old(:,:)  = pressure(:,:)
    pa_old(:,:) = pa(:,:)
    vx_old(:,:) = vx(:,:)
    vy_old(:,:) = vy(:,:)

    ! To compute an exact backward simulation, we use checkpointing
    if (NSTEP-it == NSTEP-1 .or. modulo(NSTEP-it,NSTEP/N_GLOB_FRAMES) == NSTEP/N_GLOB_FRAMES-1 ) then
        ! The closest checkpoint is a global frame. A local checkpoint is stored
        ! to avoid unnecessary forward recomputations.
        call save_local_frames(NSTEP-it)
    else if (NSTEP-it /= NSTEP) then
        ! Select the closest available checkpoint occurring before the target time.
        !  and compute the solution to the time (NSTEP-it)*dt
        call load_frame(NSTEP-it, time_last_frame)
        ! Launch a forward simulation from the previously selected frame
        call forwardproblem(p0_prior,rho0_prior,windx_prior,windy_prior,time_last_frame, NSTEP-it,1)
    endif

    ! Compute the adjoint wavefield at the current time step.
    call compute_adjoint(it)

    if ((mod(it,IT_DISPLAY) == 0 .or. it == 5) .and. method == 2) then
      call gather_and_generate_image(vax,vay,pa,rhoa,it,4)
    endif

    ! Temporal interpolation to bring all variables to the same physical time level.
    rhop_half_t(:,:)  = 0.5 * (rhop(:,:)     + rhop_old(:,:))
    rhoa_half_t(:,:)  = 0.5 * (rhoa(:,:)     + rhoa_old(:,:))
    p_half_t(:,:)     = 0.5 * (pressure(:,:) + p_old(:,:))
    pa_half_t(:,:)    = 0.5 * (pa(:,:)       + pa_old(:,:))

    value_dvy_dt(:,:) = (vy(:,:) - vy_old(:,:)) * ONE_OVER_DELTAT
    value_dvx_dt(:,:) = (vx(:,:) - vx_old(:,:)) * ONE_OVER_DELTAT

    ! Accumulate kernel contributions from the correlation between forward and
    ! adjoint wavefields at the current time step.
    call compute_kernel_iter(rho0_prior, p0_prior, windx_prior, windy_prior, rhop_half_t, p_half_t, vx, vy,&
             rhoa_half_t, pa_half_t, vax, vay, value_dvx_dt, value_dvy_dt, it)

    ! Write information about the progress of the kernel computation
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

  ! Extrapolate kernel values inside PML regions
  i_start = 1 +  NPOINTS_PML + 3
  i_end = NX_LOCAL - NPOINTS_PML - 3 - 1
  j_start = 1 + NPOINTS_PML + 3
  j_end = NY_LOCAL - NPOINTS_PML - 3 - 1

  ! Left PML boundary
  if (USE_PML_XMIN .and. i_rank == 0) then
     do jj=1,NY_LOCAL
        K_p0(1:i_start-1,jj) = K_p0(i_start,jj)
        K_rho0(1:i_start-1,jj) = K_rho0(i_start,jj)
        K_windx(1:i_start-1,jj) = K_windx(i_start,jj)
     enddo
  endif
  ! Right PML boundary
  if (USE_PML_XMAX .and. i_rank == NPROC_X-1) then
     do jj=1,NY_LOCAL
         K_p0(i_end+1:NX_LOCAL,jj) = K_p0(i_end,jj)
         K_rho0(i_end+1:NX_LOCAL,jj) = K_rho0(i_end,jj)
         K_windx(i_end+1:NX_LOCAL,jj) = K_windx(i_end,jj)
     enddo
  endif
  ! Bottom PML boundary
  if (USE_PML_YMIN .and. j_rank == 0) then
     do ii=1,NX_LOCAL
        K_p0(ii,1:j_start-1) = K_p0(ii,j_start)
        K_rho0(ii,1:j_start-1) = K_rho0(ii,j_start)
        K_windx(ii,1:j_start-1) = K_windx(ii,j_start)
     enddo
  endif
  ! Top PML boundary
  if (USE_PML_YMAX .and. j_rank == NPROC_Y-1) then
     do ii=1,NX_LOCAL
       K_p0(ii,j_end+1:NY_LOCAL) = K_p0(ii,j_end)
       K_rho0(ii,j_end+1:NY_LOCAL) = K_rho0(ii,j_end)
       K_windx(ii,j_end+1:NY_LOCAL) = K_windx(ii,j_end)
     enddo
  endif

  ! Bottom-Left PML corner
  if (USE_PML_XMIN .and. USE_PML_YMIN .and. rank == 0) then
    K_p0(1:i_start-1,1:j_start-1) = K_p0(i_start,j_start)
    K_rho0(1:i_start-1,1:j_start-1) = K_rho0(i_start,j_start)
    K_windx(1:i_start-1,1:j_start-1) = K_windx(i_start,j_start)
  endif
  ! Top-Left PML corner
  if (USE_PML_XMIN .and. USE_PML_YMAX .and. j_rank == NPROC_Y-1 .and. i_rank == 0) then
    K_p0(1:i_start-1,j_end+1:NY_LOCAL) = K_p0(i_start,j_end)
    K_rho0(1:i_start-1,j_end+1:NY_LOCAL) = K_rho0(i_start,j_end)
    K_windx(1:i_start-1,j_end+1:NY_LOCAL) = K_windx(i_start,j_end)
  endif
  ! Bottom-Right PML corner
  if (USE_PML_XMAX .and. USE_PML_YMIN .and. i_rank == NPROC_X-1 .and. j_rank == 0) then
    K_p0(i_end+1:NX_LOCAL,1:j_start-1) = K_p0(i_end,j_start)
    K_rho0(i_end+1:NX_LOCAL,1:j_start-1) = K_rho0(i_end,j_start)
    K_windx(i_end+1:NX_LOCAL,1:j_start-1) = K_windx(i_end,j_start)
  endif
  ! Top-Right PML corner
  if (USE_PML_XMAX .and. USE_PML_YMAX .and. rank == NPROC-1) then
    K_p0(i_end+1:NX_LOCAL,j_end+1:NY_LOCAL) = K_p0(i_end,j_end)
    K_rho0(i_end+1:NX_LOCAL,j_end+1:NY_LOCAL) = K_rho0(i_end,j_end)
    K_windx(i_end+1:NX_LOCAL,j_end+1:NY_LOCAL) = K_windx(i_end,j_end)
  endif

  ! Multiply K by ΔxΔy to account for the control volume (cell area in 2D) in the discrete formulation.
  K_p0(:,:)    = DELTAX * DELTAY * K_p0(:,:)
  K_rho0(:,:)  = DELTAX * DELTAY * K_rho0(:,:)
  K_windx(:,:) = DELTAX * DELTAY * K_windx(:,:)

  ! Save kernel information
  call gather_and_generate_image(K_windx,K_windy,K_p0,K_rho0,it,3)

endsubroutine compute_kernel


subroutine compute_kernel_iter(rho0, p0, windx, windy, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
              value_dvx_dt, value_dvy_dt, it_time)
!==============================================================================
! Compute and accumulate sensitivity kernels contributions at a given time step.
!
! Descriptions:
!   This routine evaluates the time-dependent contribution to the sensitivity
!   kernels using the correlation between reconstructed forward wavefields and
!   adjoint wavefields.
!
! Kernel contributions are computed from correlations involving:
!
!   - pressure perturbations
!   - density perturbations
!   - velocity fields
!   - adjoint variables
!   - background model parameters
!
! The kernel accumulation follows the discrete time integration:
! For each grid point:
!       K = K + dt * kernel_contribution
!
! Inputs:
!   rho0, p0
!       Background density and pressure models.
!   windx, windy=0
!       Background velocity components.
!   rhop, pressure, vx, vy
!       Forward wavefield variables.
!   rhoa, pa, vax, vay
!       Adjoint wavefield variables.
!   value_dvx_dt, value_dvy_dt
!       Time derivatives of the forward velocity fields.
!   it_time
!       Current iteration index.
!
! Outputs:
!   K_rho0
!       Sensitivity kernel contribution associated with the background density.
!   K_p0
!       Sensitivity kernel contribution associated with the background pressure.
!   K_windx
!       Sensitivity kernel contribution associated with the x-component of the background wind.
!
!==============================================================================
  use parameters, only :  K_rho0, K_p0, K_windx, gamma_chemestry,  &
                          DELTAT, NSTEP,   &
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
      value_dvx_dx,                       &
      value_dvy_dx,     value_dvy_dy,     &
      value_dvax_dy,                      &
      value_dwindx_dy,                    &
      value_vax_windx_dvx_dx,             &
      value_vay_windx_dvy_dx,             &
      value_vax_vy_dwindx_dy

  double precision ::   &
      value_drhop_dx,   &
      value_drhoa_dy,   &
      value_dp_dx,      &
      value_dpa_dy,     &
      value_vax_dvx_dx, &
      value_vay_dvy_dx


  double precision :: &
    rho0_half_x,                  &
    rhoa_half_x,                  &
    rhop_half_x,                  &
    pa_half_x,                    &
    windx_half_x,        &
    vax_half_x,          vay_half_y,          &
    vx_half_x,           vy_half_y,           &
    vay_half_x_half_y,   &
    vy_half_x_half_y,    &
    vax_dvx_dt_half_x, vay_dvy_dt_half_y,     &
    value_dvy_dy_half_x_half_y


   double precision :: &
       value_dwindx_dy_next, value_dwindx_dy_prec,          &
       value_dvy_dx_next,    value_dvy_dx_prec,             &
       value_dvy_dy_next,    value_dvy_dy_prec

  double precision, dimension(-1:NX_LOCAL+2, -1:NY_LOCAL+2) :: &
      vax_dvx_dt, vay_dvy_dt

  double precision  :: u_mm, u_m, u_p, u_pp
  integer :: Ip1,Im1, Ip2, Im2, Jp1,Jm1, Jp2, Jm2

  integer :: i_start, i_end, j_start, j_end

  ! Compute the time-dependent source term used in the forward simulation.
  t = dble(NSTEP-it_time+1-0.5)*DELTAT
  save_source_term = - 2 * a* (t-t0) *  factor * exp(- a*(t-t0)*(t-t0))


  ! Define computational boundaries excluding PML regions.

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



  ! Compute auxiliary products
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

  ! Exchange halo values with neighboring MPI processes
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


  !----------------------------------------------
  ! Contribution to the wind-x sensitivity kernel
  !----------------------------------------------

  ! Special treatment of the bottom boundary for a reflective ground condition.
  if (.not. USE_PML_YMIN .and. j_rank == 0) then
     call derivative_first_row_kernel_kvx(rho0, p0, windx, windy, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
               value_dvx_dt, value_dvy_dt, it_time, i_start, i_end)
  endif


  do j=j_start,j_end
   do i =i_start,i_end

      ! Compute global and local neighbor indices
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


      ! Interpolate variables to the staggered-grid locations required by the
      ! finite-difference formulation.
      rhoa_half_x        = 0.5d0 * (rhoa(i,j) + rhoa(Im1,j))
      rho0_half_x        = 0.5d0 * (rho0(i,j) + rho0(Im1,j))
      rhop_half_x        = 0.5d0 * (rhop(i,j) + rhop(Im1,j))
      pa_half_x          = 0.5d0 * (pa(i,j)   + pa(Im1,j))

      vay_half_x_half_y = 0.25d0 * (vay(i,j) + vay(Im1,j)+ vay(i,Jm1) + vay(Im1,Jm1))
      vy_half_x_half_y  = 0.25d0 * (vy(i,j) + vy(Im1,j)+ vy(i,Jm1) + vy(Im1,Jm1))


      ! Compute derivatives of rhoprhoa, gammapap and pressure with respect to x
      u_mm = pressure(Im2,j); u_m = pressure(Im1,j); u_p = pressure(i,j); u_pp = pressure(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      u_mm = rhop(Im2,j); u_m = rhop(Im1,j); u_p = rhop(i,j); u_pp = rhop(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhop_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

      ! Compute derivatives of vx with respect to x
      if (windx(i,j) >= 0) then
          u_mm = vx(Im2,j); u_m = vx(Im1,j); u_p = vx(Ip1,j)
          call compute_decentered_backward_dU(u_mm, u_m, vx(i,j), u_p, value_dvx_dx, ONE_OVER_SIX_DELTAX)
      else
          u_pp = vx(Ip2,j)     ; u_p = vx(Ip1,j)     ; u_m = vx(Im1,j)
          call compute_decentered_forward_dU(u_m, vx(i,j), u_p, u_pp, value_dvx_dx, ONE_OVER_SIX_DELTAX)
      endif


      ! Compute derivatives of vy with respect to x
      u_mm = vy(Im2,Jm1); u_m = vy(Im1,Jm1); u_p = vy(i,Jm1); u_pp = vy(Ip1,Jm1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      u_mm = vy(Im2,j); u_m = vy(Im1,j); u_p = vy(i,j); u_pp = vy(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      value_dvy_dx = 0.5d0 * ( value_dvy_dx_prec + value_dvy_dx_next)


      ! Compute derivatives of vy with respect to y
      u_mm = vy(Im1,Jm2); u_m = vy(Im1,Jm1); u_p = vy(Im1,j); u_pp = vy(Im1,Jp1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      u_mm = vy(i,Jm2); u_m = vy(i,Jm1); u_p = vy(i,j); u_pp = vy(i,Jp1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      value_dvy_dy_half_x_half_y = 0.5d0 * ( value_dvy_dy_prec + value_dvy_dy_next)


      ! Compute derivatives of vax with respect to y
      u_mm = vax(i,Jm2); u_m = vax(i,Jm1); u_p = vax(i,Jp1)
      call compute_decentered_backward_dU(u_mm, u_m, vax(i,j), u_p, value_dvax_dy, ONE_OVER_SIX_DELTAX)


      ! Compute intermediate terms appearing in the kernel expression
      value_vax_dvx_dx = vax(i,j) * value_dvx_dx
      value_vay_dvy_dx = vay_half_x_half_y * value_dvy_dx


      ! UPdate the horizontal wind sensitivity kernel
      K_windx(i,j) = K_windx(i,j) + (rhoa_half_x     * value_drhop_dx)                       * DELTAT
      K_windx(i,j) = K_windx(i,j) + rho0_half_x      * (value_vax_dvx_dx + value_vay_dvy_dx) * DELTAT
      K_windx(i,j) = K_windx(i,j) - vax(i,j) * rho0_half_x * value_dvy_dy_half_x_half_y      * DELTAT
      K_windx(i,j) = K_windx(i,j) - vy_half_x_half_y * rho0_half_x * value_dvax_dy           * DELTAT
      K_windx(i,j) = K_windx(i,j) + (pa_half_x       * value_dp_dx)                          * DELTAT

   enddo
  enddo



!---------------------------------------------------------------------
! Compute density and pressure kernel contributions
!---------------------------------------------------------------------

  ! Special treatment required for a reflective bottom boundary
  if (.not. USE_PML_YMIN .and. j_rank == 0) then
    call derivative_first_row_kernel_kp0krho0(rho0, p0, windx, windy, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
              value_dvx_dt, value_dvy_dt, it_time, i_start, i_end, save_source_term)
  endif


  do j=j_start,j_end
   do i =i_start,i_end

      ! Compute global and local neighbors index
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


      ! Interpolate variables to the staggered-grid locations required by the
      ! finite-difference formulation.
      vx_half_x          = 0.5d0 * (vx(i,j)  + vx(Ip1,j))
      vax_half_x         = 0.5d0 * (vax(i,j) + vax(Ip1,j))
      windx_half_x       = 0.5d0 * (windx(i,j) + windx(Ip1,j))
      vax_dvx_dt_half_x  = 0.5d0 * (vax_dvx_dt(i,j) + vax_dvx_dt(Ip1,j))

      vy_half_y          = 0.5d0 * (vy(i,j)  + vy(i,Jm1))
      vay_half_y         = 0.5d0 * (vay(i,j) + vay(i,Jm1))
      vay_dvy_dt_half_y  = 0.5d0 * (vay_dvy_dt(i,j) + vay_dvy_dt(i,Jm1))


      ! Evaluate the source function
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


      ! Compute derivatives of vx with respect to x
      u_mm = vx(Im1,j); u_m = vx(i,j); u_p = vx(Ip1,j); u_pp = vx(Ip2,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)


      ! Compute derivatives of vy with respect to y
      u_mm = vy(i,Jm2); u_m = vy(i,Jm1); u_p = vy(i,j); u_pp = vy(i,Jp1)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)



      ! Compute derivatives of vy with respect to x
      if (windx_half_x >= 0) then
          u_mm = vy(Im2,Jm1)      ; u_m = vy(Im1,Jm1)      ; u_p = vy(Ip1,Jm1)
          call compute_decentered_backward_dU(u_mm, u_m, vy(i,Jm1), u_p, value_dvy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_mm = vy(Im2,j)      ; u_m = vy(Im1,j)      ; u_p = vy(Ip1,j)
          call compute_decentered_backward_dU(u_mm, u_m, vy(i,j), u_p, value_dvy_dx_next, ONE_OVER_SIX_DELTAX)

      else
          u_pp = vy(Ip2,Jm1)      ; u_p = vy(Ip1,Jm1)      ; u_m = vy(Im1,Jm1)
          call compute_decentered_forward_dU(u_m, vy(i,Jm1), u_p, u_pp, value_dvy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_pp = vy(Ip2,j)      ; u_p = vy(Ip1,j)      ; u_m = vy(Im1,j)
          call compute_decentered_forward_dU(u_m, vy(i,j), u_p, u_pp, value_dvy_dx_next, ONE_OVER_SIX_DELTAX)
      endif
      value_dvy_dx = 0.5d0 * (value_dvy_dx_prec + value_dvy_dx_next)


      ! Compute derivatives of rhoa, pa with respect to y
      u_mm = pa(i,Jm2)      ; u_m = pa(i,Jm1)      ; u_p = pa(i,Jp1)
      call compute_decentered_backward_dU(u_mm, u_m, pa(i,j), u_p, value_dpa_dy, ONE_OVER_SIX_DELTAY)
      u_mm = rhoa(i,Jm2)      ; u_m = rhoa(i,Jm1)      ; u_p = rhoa(i,Jp1)
      call compute_decentered_backward_dU(u_mm, u_m, rhoa(i,j), u_p, value_drhoa_dy, ONE_OVER_SIX_DELTAY)



     ! Compute derivatives of windx, vx with respect to y
     u_mm = windx(i,Jm2)      ; u_m = windx(i,Jm1)      ; u_p = windx(i,Jp1)
     call compute_decentered_backward_dU(u_mm, u_m, windx(i,j), u_p, value_dwindx_dy_prec, ONE_OVER_SIX_DELTAY)
     u_mm = windx(Ip1,Jm2)      ; u_m = windx(Ip1,Jm1)      ; u_p = windx(Ip1,Jp1)
     call compute_decentered_backward_dU(u_mm, u_m, windx(Ip1,j), u_p, value_dwindx_dy_next, ONE_OVER_SIX_DELTAY)

      value_dwindx_dy = 0.5d0 * (value_dwindx_dy_prec + value_dwindx_dy_next)
      

      ! Compute intermediate terms appearing in the kernel expression
      value_vax_windx_dvx_dx = vax_half_x * windx_half_x * value_dvx_dx
      value_vay_windx_dvy_dx = vay_half_y * windx_half_x * value_dvy_dx
      value_vax_vy_dwindx_dy = vax_half_x * vy_half_y * value_dwindx_dy



      ! Update the density sensitivity kernel 
      K_rho0(i,j) = K_rho0(i,j) + (rhoa(i,j) * value_dvx_dx - vy_half_y * value_drhoa_dy) * DELTAT
      K_rho0(i,j) = K_rho0(i,j) - (vax_dvx_dt_half_x + vay_dvy_dt_half_y) * DELTAT
      K_rho0(i,j) = K_rho0(i,j) + (value_vax_windx_dvx_dx + value_vay_windx_dvy_dx) * DELTAT
      K_rho0(i,j) = K_rho0(i,j) + (value_vax_vy_dwindx_dy) * DELTAT
      K_rho0(i,j) = K_rho0(i,j) + pa(i,j) * factor_ssf * source_term * (gamma_chemestry(i,j) * p0(i,j) / rho0(i,j)**2) * DELTAT


      ! Update the pressure sensitivity kernel 
      K_p0(i,j) = K_p0(i,j) - (vy_half_y * value_dpa_dy) * DELTAT
      K_p0(i,j) = K_p0(i,j) + gamma_chemestry(i,j) * pa(i,j) * (value_dvx_dx + value_dvy_dy)        * DELTAT
      K_p0(i,j) = K_p0(i,j) - pa(i,j) * value_dvy_dy                                             * DELTAT
      K_p0(i,j) = K_p0(i,j) - pa(i,j) * factor_ssf * source_term * gamma_chemestry(i,j) / rho0(i,j) * DELTAT

   enddo
  enddo

endsubroutine compute_kernel_iter
