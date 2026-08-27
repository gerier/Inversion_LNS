subroutine compute_adjoint(it_step)
!==============================================================================
! Compute one iteration of the adjoint simulation.
!
! This routine updates the adjoint field by one time step.
! It computes the adjoint source term associated with the current iteration,
! updates the adjoint velocity and pressure variables, applies the PML
! formulation, exchanges ghost cells between MPI subdomains, and enforces
! boundary conditions.
!
! Input:
!   it_step
!       Current adjoint iteration index.
!
! Updated fields:
!   - adjoint pressure
!   - adjoint density perturbation
!   - adjoint velocity components
!   - PML recursive-convolution memory variables
!
!==============================================================================
  use parameters 
  implicit none

  integer it_step

  ! source term
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: adjoint_source_term

  ! Temporary variables used to avoid overwriting the current fields.
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: rhoa_old, pa_old, vax_old, vay_old

  ! derivatives for equation on vax and vay
  double precision ::                                 &
    value_drhoa_dx,     value_drhoa_dy,               &
    value_dgammap0pa_dx,     value_dgammap0pa_dy,     &
    value_dp0_dy,       &
    value_dwindx_dy,    &
    value_dvax_dx ,     &
    value_dvay_dx,      value_dvay_dy,      &
    value_windxdvax_dx, &
    value_windxdvay_dx

  ! Derivatives required for the rhoa and pa equations.
  double precision ::          &
    value_dpa_dx

  ! Interpolated quantities.
  double precision ::                         &
    rho0_half_x,         rho0_half_y,         &
    rhoa_half_x,         rhoa_half_y,         &
    pa_half_x,           pa_half_y,           &
    vax_half_x_half_y,   vay_half_x_half_y,   &
    vax_half_x,          vay_half_y,          &
    windx_half_x_half_y,                      &
    windx_half_x
  
  ! Intermediate variables used to interpolate derivative quantities.
  double precision :: &
    value_dwindx_dy_next, value_dwindx_dy_prec


  integer :: i,j, j_start
  integer :: Ip1,Im1, Ip2, Im2, Jp1,Jm1, Jp2, Jm2
  double precision :: u_mm, u_m, u_p, u_pp


  ! Build the adjoint source associated with the selected objective function.
  call compute_adjoint_source(adjoint_source_term, it_step)

  !---
  ! vax and vay computations
  !---

  ! Store the current fields before updating them.
  vax_old(:,:) = vax(:,:)
  vay_old(:,:) = vay(:,:)

  call send_receive_rightleft(rhoa)
  call send_receive_rightleft(pa)
  call send_receive_rightleft(vax_old)
  call send_receive_rightleft(vay_old)

  call send_receive_topbottom(rhoa)
  call send_receive_topbottom(pa)
  call send_receive_topbottom(vax_old)
  call send_receive_topbottom(vay_old)

  call send_receive_righttop(vax_old)
  call send_receive_leftbottom(vay_old)

  if (.not. USE_PML_YMIN .and. j_rank == 0) then
    ! Use a dedicated finite-difference stencil to enforce the free-surface
    ! (ground reflection) boundary condition.
    j_start = 3
    call derivative_first_row_vaxvay()
  else
    j_start = 1
  endif

  do j=j_start,NY_LOCAL
   do i=1,NX_LOCAL

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


    ! Compute quantities interpolated at staggered-grid locations.
    rhoa_half_x = 0.5d0 * (rhoa(i,j) + rhoa(Im1,j))
    rho0_half_x  = 0.5d0 * (rho0_prior(i,j)  + rho0_prior(Im1,j))
    pa_half_x   = 0.5d0 * (pa(i,j)   + pa(Im1,j))
    vay_half_x_half_y = 0.25d0 * (vay_old(i,j) + vay_old(Im1,j) + vay_old(i,Jm1) + vay_old(Im1,Jm1))

    rhoa_half_y = 0.5d0 * (rhoa(i,j) + rhoa(i,Jp1))
    rho0_half_y  = 0.5d0 * (rho0_prior(i,j)  + rho0_prior(i,Jp1))
    pa_half_y   = 0.5d0 * (pa(i,j)   + pa(i,Jp1))
    vax_half_x_half_y = 0.25d0 * (vax_old(i,j) + vax_old(Ip1,j) + vax_old(i,Jp1) + vax_old(Ip1,Jp1))
    windx_half_x_half_y = 0.25d0 * (windx_prior(i,j)     + windx_prior(Ip1,j)     + windx_prior(i,Jp1)     + windx_prior(Ip1,Jp1))


    ! Compute derivatives of rho0rhoa, gammap0pa with respect to x
    u_mm = rhoa(Im2,j); u_m = rhoa(Im1,j); u_p = rhoa(i,j); u_pp = rhoa(Ip1,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhoa_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

    u_mm = p0_prior(Im2,j) * pa(Im2,j) * gamma_chemestry(Im2,j); u_m = p0_prior(Im1,j) * pa(Im1,j) * gamma_chemestry(Im1,j)
    u_p  = p0_prior(i,j) * pa(i,j) * gamma_chemestry(i,j)   ; u_pp = p0_prior(Ip1,j) * pa(Ip1,j) * gamma_chemestry(Ip1,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dgammap0pa_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

    eq1_memory_drhoa_dx_adj(i,j)      = b_x(i_global) * eq1_memory_drhoa_dx_adj(i,j) + a_x(i_global) * value_drhoa_dx
    eq1_memory_dgammap0pa_dx_adj(i,j) = b_x(i_global) * eq1_memory_dgammap0pa_dx_adj(i,j) &
                                        + a_x(i_global) * value_dgammap0pa_dx


    ! Compute derivatives of windx and vax with respect to x
    if (windx_prior(i,j) <= 0) then 
    ! In the adjoint problem, the upwind discretization is reversed
    ! with respect to the forward problem.
          u_mm = vax_old(Im2,j); u_m = vax_old(Im1,j); u_p = vax_old(Ip1,j)
          call compute_decentered_backward_dU(u_mm, u_m, vax_old(i,j), u_p, value_dvax_dx, ONE_OVER_SIX_DELTAX)
     else
          u_pp = vax_old(Ip2,j)     ; u_p = vax_old(Ip1,j)     ; u_m = vax_old(Im1,j)
          call compute_decentered_forward_dU(u_m, vax_old(i,j), u_p, u_pp, value_dvax_dx, ONE_OVER_SIX_DELTAX)
     endif
     
    eq1_memory_dvax_dx_adj(i,j) = b_x(i_global) * eq1_memory_dvax_dx_adj(i,j) + a_x(i_global) * value_dvax_dx



    ! Compute derivatives of rho0, p0, rho0rhoa, gammap0pa with respect to y
    u_mm = rhoa(i,Jm1); u_m = rhoa(i,j); u_p = rhoa(i,Jp1); u_pp = rhoa(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhoa_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = p0_prior(i,Jm1); u_m = p0_prior(i,j); u_p = p0_prior(i,Jp1); u_pp = p0_prior(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp0_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)

    u_mm = p0_prior(i,Jm1) * pa(i,Jm1) * gamma_chemestry(i,Jm1); u_m = p0_prior(i,j) * pa(i,j) * gamma_chemestry(i,j)
    u_p  = p0_prior(i,Jp1) * pa(i,Jp1) * gamma_chemestry(i,Jp1); u_pp = p0_prior(i,Jp2) * pa(i,Jp2) * gamma_chemestry(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dgammap0pa_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)

    eq1_memory_drhoa_dy_adj(i,j)     = b_y_half(j_global)*eq1_memory_drhoa_dy_adj(i,j)     + a_y_half(j_global)*value_drhoa_dy
    ! eq1_memory_dp0_dy_adj(i,j)       = b_y_half(j_global)*eq1_memory_dp0_dy_adj(i,j)       + a_y_half(j_global)*value_dp0_dy
    eq1_memory_dgammap0pa_dy_adj(i,j)     = b_y_half(j_global)*eq1_memory_dgammap0pa_dy_adj(i,j)     &
                                      + a_y_half(j_global)*value_dgammap0pa_dy


    ! Compute derivatives of vay with respect to x
    if (windx_half_x_half_y <= 0) then 
    ! In the adjoint problem, the upwind discretization is reversed
    ! with respect to the forward problem.
          u_mm = vay_old(Im2,j); u_m = vay_old(Im1,j); u_p = vay_old(Ip1,j)
          call compute_decentered_backward_dU(u_mm, u_m, vay_old(i,j), u_p, value_dvay_dx, ONE_OVER_SIX_DELTAX)
    else
          u_pp = vay_old(Ip2,j); u_p = vay_old(Ip1,j); u_m = vay_old(Im1,j)
          call compute_decentered_forward_dU(u_m, vay_old(i,j), u_p, u_pp, value_dvay_dx, ONE_OVER_SIX_DELTAX)
    endif

    eq1_memory_dvay_dx_adj(i,j) = b_x_half(i_global) * eq1_memory_dvay_dx_adj(i,j) + a_x_half(i_global) * value_dvay_dx


    ! Compute derivatives of windx with respect to y, at the vay location
    u_mm = windx_prior(i,Jm1); u_m = windx_prior(i,j); u_p = windx_prior(i,Jp1); u_pp = windx_prior(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = windx_prior(Ip1,Jm1); u_m = windx_prior(Ip1,j); u_p = windx_prior(Ip1,Jp1); u_pp = windx_prior(Ip1,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    value_dwindx_dy = 0.5d0 * ( value_dwindx_dy_prec + value_dwindx_dy_next)


    ! Apply the PML correction to the computed derivatives.
    value_drhoa_dx     = value_drhoa_dx     * one_over_K_x(i_global)      + eq1_memory_drhoa_dx_adj(i,j)
    value_dgammap0pa_dx= value_dgammap0pa_dx* one_over_K_x(i_global)      + eq1_memory_dgammap0pa_dx_adj(i,j)
    value_dvax_dx      = value_dvax_dx      * one_over_K_x(i_global)      + eq1_memory_dvax_dx_adj(i,j)
    value_drhoa_dy     = value_drhoa_dy     * one_over_K_y_half(j_global) + eq1_memory_drhoa_dy_adj(i,j)
    value_dgammap0pa_dy= value_dgammap0pa_dy* one_over_K_y_half(j_global) + eq1_memory_dgammap0pa_dy_adj(i,j)
    value_dvay_dx      = value_dvay_dx      * one_over_K_y_half(j_global) + eq1_memory_dvay_dx_adj(i,j)
    value_dp0_dy       = value_dp0_dy       * one_over_Kdalpha_y_half(j_global)
    value_dwindx_dy    = value_dwindx_dy    * one_over_Kdalpha_y_half(j_global)
    

    ! Compute intermediate terms appearing in the adjoint expression
    value_windxdvax_dx = windx_prior(i,j) * value_dvax_dx
    value_windxdvay_dx = windx_half_x_half_y * value_dvay_dx

    ! Update the x-component of the adjoint velocity.
    vax(i,j) = vax(i,j) +  value_drhoa_dx * DELTAT
    vax(i,j) = vax(i,j) + value_windxdvax_dx * DELTAT
    vax(i,j) = vax(i,j) + value_dgammap0pa_dx * DELTAT / rho0_half_x

    ! Update the y-component of the adjoint velocity.
    vay(i,j) = vay(i,j) +  value_drhoa_dy * DELTAT
    vay(i,j) = vay(i,j) + value_windxdvay_dx * DELTAT
    vay(i,j) = vay(i,j) - (vax_half_x_half_y * value_dwindx_dy) * DELTAT
    vay(i,j) = vay(i,j) + (value_dgammap0pa_dy - pa_half_y * value_dp0_dy) * DELTAT / rho0_half_y

  enddo
  enddo


  ! Apply Dirichlet boundary conditions.
  !! Left boundary
  if (USE_PML_XMIN .and. i_rank == 0) then
      vax(-1:1,:) = ZERO
      vay(-1:1,:) = ZERO
  endif
  !! Right boundary
  if (USE_PML_XMAX .and. i_rank == NPROC_X-1) then
      vax(NX_LOCAL:NX_LOCAL+2,:) = ZERO
      vay(NX_LOCAL:NX_LOCAL+2,:) = ZERO
  endif
  !! Bottom boundary
  if (USE_PML_YMIN .and. j_rank == 0) then
      vax(:,-1:1) = ZERO
      vay(:,-1:1) = ZERO
  else if (j_rank == 0) then
     do j=-1,0
       vay(:,j) = -vay(:,2-j)
     enddo
     vay(:,1) = ZERO
  endif
  !! Top boundary
  if (USE_PML_YMAX .and. j_rank == NPROC_Y -1) then
      vax(:,NY_LOCAL:NY_LOCAL+2) = ZERO
      vay(:,NY_LOCAL:NY_LOCAL+2) = ZERO
  endif

  !---
  ! pa and rhoa computations
  !---

  ! Store the current fields before updating them.
  rhoa_old(:,:) = rhoa(:,:)
  pa_old(:,:) = pa(:,:)

  call send_receive_rightleft(vax)
  call send_receive_topbottom(vay)

  call send_receive_rightleft(pa)
  call send_receive_rightleft(rhoa)
  call send_receive_topbottom(pa)
  call send_receive_topbottom(rhoa)

  if ((.not. USE_PML_YMIN) .and. j_rank == 0) then
   ! Use a dedicated finite-difference stencil to enforce the free-surface
   ! (ground reflection) boundary condition.
   j_start = 3
   call derivative_first_row_parhoa(adjoint_source_term)
  else
    j_start = 1
  endif

  do j=j_start,NY_LOCAL
   do i=1,NX_LOCAL

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


    ! Compute staggered-grid interpolated quantities.
    vax_half_x = 0.5d0 * (vax(i,j) + vax(Ip1,j))
    vay_half_y = 0.5d0 * (vay(i,j) + vay(i,Jm1))
    windx_half_x = 0.5d0 * (windx_prior(i,j) + windx_prior(Ip1,j))

    ! Compute derivatives of rhoa and pa with respect to x
    if (windx_half_x <= 0) then 
    ! In the adjoint problem, the upwind discretization is reversed
    ! with respect to the forward problem.
          u_mm = pa_old(Im2,j)      ; u_m = pa_old(Im1,j)      ; u_p = pa_old(Ip1,j)
          call compute_decentered_backward_dU(u_mm, u_m, pa_old(i,j), u_p, value_dpa_dx, ONE_OVER_SIX_DELTAX)
          u_mm = rhoa_old(Im2,j); u_m = rhoa_old(Im1,j); u_p = rhoa_old(Ip1,j)
          call compute_decentered_backward_dU(u_mm, u_m, rhoa_old(i,j), u_p, value_drhoa_dx, ONE_OVER_SIX_DELTAX)
    else
          u_pp = pa_old(Ip2,j)      ; u_p = pa_old(Ip1,j)      ; u_m = pa_old(Im1,j)
          call compute_decentered_forward_dU(u_m, pa_old(i,j), u_p, u_pp, value_dpa_dx, ONE_OVER_SIX_DELTAX)
          u_pp = rhoa_old(Ip2,j); u_p = rhoa_old(Ip1,j); u_m = rhoa_old(Im1,j)
          call compute_decentered_forward_dU(u_m, rhoa_old(i,j), u_p, u_pp, value_drhoa_dx, ONE_OVER_SIX_DELTAX)
    endif

    eq2_memory_dpa_dx_adj(i,j)   = b_x_half(i_global)*eq2_memory_dpa_dx_adj(i,j)   &
                                        + a_x_half(i_global)*value_dpa_dx
    eq2_memory_drhoa_dx_adj(i,j) = b_x_half(i_global)*eq2_memory_drhoa_dx_adj(i,j) &
                                        + a_x_half(i_global)*value_drhoa_dx


    ! Compute derivatives of vax with respect to x
    u_mm = vax(Im1,j); u_m = vax(i,j); u_p = vax(Ip1,j); u_pp = vax(Ip2,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvax_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

    eq2_memory_dvax_dx_adj(i,j) = b_x_half(i_global) * eq2_memory_dvax_dx_adj(i,j) + a_x_half(i_global) * value_dvax_dx


    ! Compute derivatives of vay with respect to y
    u_mm = vay(i,Jm2); u_m = vay(i,Jm1); u_p = vay(i,j); u_pp = vay(i,Jp1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvay_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)

    eq2_memory_dvay_dy_adj(i,j) = b_y(j_global) * eq2_memory_dvay_dy_adj(i,j) + a_y(j_global) * value_dvay_dy


    ! PML treatment of the derivatives
    value_dpa_dx   = value_dpa_dx   * one_over_K_x_half(i_global) + eq2_memory_dpa_dx_adj(i,j)
    value_drhoa_dx = value_drhoa_dx * one_over_K_x_half(i_global) + eq2_memory_drhoa_dx_adj(i,j)
    value_dvax_dx  = value_dvax_dx  * one_over_K_x_half(i_global) + eq2_memory_dvax_dx_adj(i,j)
    value_dvay_dy  = value_dvay_dy  * one_over_K_y(j_global)      + eq2_memory_dvay_dy_adj(i,j)


    ! Update the adjoint density perturbation and pressure.
    rhoa(i,j) = rhoa(i,j) + windx_half_x * value_drhoa_dx            * DELTAT
    rhoa(i,j) = rhoa(i,j) - vay_half_y * g(i,j)                      * DELTAT

    pa(i,j) = pa(i,j) + windx_half_x * value_dpa_dx                  * DELTAT
    pa(i,j) = pa(i,j) + (value_dvax_dx + value_dvay_dy)              * DELTAT
    pa(i,j) = pa(i,j) + adjoint_source_term(i,j)                     * DELTAT

   enddo
  enddo

  ! Exchange ghost-cell values with neighboring MPI subdomains.
  call send_receive_rightleft(rhoa)
  call send_receive_rightleft(pa)
  call send_receive_rightleft(vax)
  call send_receive_rightleft(vay)

  call send_receive_topbottom(rhoa)
  call send_receive_topbottom(pa)
  call send_receive_topbottom(vax)
  call send_receive_topbottom(vay)


  ! Dirichlet conditions
  !! Left boundary
  if (USE_PML_XMIN .and. i_rank == 0) then
      pa(-1:1,:) = ZERO
      rhoa(-1:1,:) = ZERO
  endif
  !! Right boundary
  if (USE_PML_XMAX .and. i_rank == NPROC_X-1) then
      pa(NX_LOCAL:NX_LOCAL+2,:) = ZERO
      rhoa(NX_LOCAL:NX_LOCAL+2,:) = ZERO
  endif
  !! Bottom boundary
  if (USE_PML_YMIN .and. j_rank == 0) then
      pa(:,-1:1) = ZERO
      rhoa(:,-1:1) = ZERO
  endif
  !! Top boundary
  if (USE_PML_YMAX .and. j_rank == NPROC_Y-1) then
      pa(:,NY_LOCAL:NY_LOCAL+2) = ZERO
      rhoa(:,NY_LOCAL:NY_LOCAL+2) = ZERO
  endif


endsubroutine compute_adjoint



subroutine prepare_adjoint_source()
!==============================================================================
! Prepare the adjoint source time functions.
!
! This routine constructs the adjoint source from the difference between
! observed and simulated pressure signals (or from travel-time residuals,
! depending on the selected objective function).
!
! The source can optionally be:
!   - normalized receiver by receiver,
!   - tapered during the first iterations,
!   - windowed around selected arrivals.
!
! The resulting adjoint source is stored for subsequent adjoint simulations.
!
!==============================================================================
use parameters, only : sispressure_true, sispressure_prior, normsq_pressure_true_per_rec,&
                      NSTEP, NREC, t0, DELTAT, PI, observation, adjoint_source, TINYVAL, &
                      wr, window_waveform, REC_wr, save_adjoint_source
implicit none
integer :: irec, it_step
double precision, dimension(NSTEP) :: coef_damping
integer :: it_t0, it, i_tmin, i_delta_tmin
integer :: i,j

  if (observation == 0) then
    coef_damping(:) = 1.0d0
    do it_step=1,NSTEP
      ! Apply an initial cosine taper to reduce high-frequency artifacts.
      if (it_step*DELTAT < t0) then
        coef_damping(NSTEP-it_step+1) = 0.5 - 0.5 * cos(2*PI*it_step*DELTAT/(2*t0))
      endif
    enddo

    adjoint_source(:,:) = (sispressure_true(:,:) - sispressure_prior(:,:))
    do irec=1,NREC

      if (window_waveform == 1) then
         wr(:) = 0.d0
         i_tmin = REC_wr(irec,1) / DELTAT
         i_delta_tmin = REC_wr(irec,2) / DELTAT
         wr(i_tmin: i_tmin+i_delta_tmin) = 1.0d0
         ! Smoothly taper both ends of the correlation window with a Hann window.
         ! This ensures that the signal vanishes at the beginning and end of the window.
         it_t0 = t0/DELTAT
         do it=1,it_t0
            wr(i_tmin - it) = 0.5 + 0.5 * cos(PI*it*DELTAT/t0)
            wr(it + i_tmin + i_delta_tmin) = 0.5 + 0.5 * cos(PI*it*DELTAT/t0)
         enddo
      else
        wr(:) = 1.0
      endif

      if (normsq_pressure_true_per_rec(irec) > TINYVAL) then
        adjoint_source(:,irec) = adjoint_source(:,irec) / normsq_pressure_true_per_rec(irec)
        adjoint_source(:,irec) =  adjoint_source(:,irec) * coef_damping(:) * wr(:)
      endif
    enddo

  else

   call get_Tr_adjoint_source()

  endif

  ! Save the adjoint source for post-processing.
  if (save_adjoint_source) then
    open(unit=10, file='./OUTPUT/source_adjointe.txt', status='replace', action='write')
    do i = 1, NSTEP
       write(10, '(50(ES16.8,1X))') (adjoint_source(i,j), j=1,NREC)
    end do
    close(10)
  endif

endsubroutine prepare_adjoint_source


subroutine compute_adjoint_source(adjoint_source_term, it_step)
!==============================================================================
! Compute the spatial adjoint source term.
!
! This routine injects the adjoint source corresponding to the current
! iteration at the receiver locations owned by the current MPI subdomain.
!
! Inputs:
!   it_step
!       Current adjoint iteration index.
!
! Output:
!   adjoint_source_term
!       Spatial distribution of the adjoint source term.
!
!==============================================================================
 use parameters, only : ix_rec, iy_rec, NREC, NSTEP, adjoint_source, &
                        i_rank, j_rank,NX_LOCAL,NY_LOCAL, offset_i, offset_j

 double precision, dimension(-1:NX_LOCAL+2, -1:NY_LOCAL+2) :: adjoint_source_term
 integer :: it_step, irec
 integer :: i,j

  adjoint_source_term(:,:) = 0.0d0
  do irec=1,NREC

     if (i_rank == (ix_rec(irec)-1)/NX_LOCAL .and. j_rank == (iy_rec(irec)-1)/NY_LOCAL) then
       i = ix_rec(irec) - offset_i
       j = iy_rec(irec) - offset_j

       adjoint_source_term(i,j) = adjoint_source(NSTEP-it_step+1,irec)
     endif

  enddo

endsubroutine compute_adjoint_source
