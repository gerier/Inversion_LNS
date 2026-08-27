!!======================================================================
!! Adjoint and sensitivity-kernel routines for the first grid row
!! (reflective ground boundary condition)
!!======================================================================


!======================================================================
! Adjoint routines
!======================================================================

subroutine derivative_first_row_vaxvay()
!======================================================================
!> Update the adjoint velocity components (vax, vay) on the first grid row.
!!
!! This routine computes the spatial derivatives required by the adjoint
!! momentum equations on the reflective-ground boundary. Spatial
!! derivatives are evaluated using centered or reversed-upwind finite
!! differences depending on the local flow direction.
!!======================================================================
  use parameters
  implicit none

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: vax_old, vay_old

  ! derivatives for equation on vax and vay
  double precision ::                       &
    value_drhoa_dx,     value_drhoa_dy,     &
    value_dgammap0pa_dx,     value_dgammap0pa_dy,     &
    value_dp0_dy,                           &
    value_dwindx_dy,                        &
    value_dvax_dx,                          &
    value_dvay_dx,                          &
    value_windxdvax_dx, value_windxdvay_dx

  ! Interpolated variables
  double precision ::                         &
    rho0_half_x,         rho0_half_y,         &
    rhoa_half_x,         rhoa_half_y,         &
    pa_half_x,           pa_half_y,           &
    vax_half_x_half_y,   vay_half_x_half_y,   &
    windx_half_x_half_y

  ! Intermediate variables used to interpolate derivative quantities.
  double precision :: &
    value_dwindx_dy_next, value_dwindx_dy_prec

  ! Compute global and local neighbor indices
  integer :: i
  integer :: Ip1,Im1, Ip2, Im2
  double precision :: u_mm, u_m, u_p, u_pp


  !!!!!!!!!!!!!!!!!!!!!!!
  ! Update adjoint velocity components
  !!!!!!!!!!!!!!!!!!!!!!!

  ! Store the current fields before the in-place update.
  vax_old(:,:) = vax(:,:)
  vay_old(:,:) = vay(:,:)

  do i=1,NX_LOCAL

    ! Compute global and local neighbor indices
    i_global = i + offset_i

    Im1 = i-1
    Im2 = i-2
    Ip1 = i+1
    Ip2 = i+2


    ! Compute quantities interpolated at staggered-grid locations.
    rhoa_half_x = 0.5d0 * (rhoa(i,2) + rhoa(Im1,2))
    rho0_half_x  = 0.5d0 * (rho0_prior(i,2)  + rho0_prior(Im1,2))
    pa_half_x   = 0.5d0 * (pa(i,2)   + pa(Im1,2))
    vay_half_x_half_y = 0.25d0 * (vay_old(i,2) + vay_old(Im1,2) + vay_old(i,1) + vay_old(Im1,1))
    
    rhoa_half_y = 0.5d0 * (rhoa(i,2) + rhoa(i,3))
    rho0_half_y  = 0.5d0 * (rho0_prior(i,2)  + rho0_prior(i,3))
    pa_half_y   = 0.5d0 * (pa(i,2)   + pa(i,3))
    vax_half_x_half_y = 0.25d0 * (vax_old(i,2) + vax_old(Ip1,2) + vax_old(i,3) + vax_old(Ip1,3))
    windx_half_x_half_y = 0.25d0 * (windx_prior(i,2)     + windx_prior(Ip1,2)     + windx_prior(i,3)     + windx_prior(Ip1,3))


    ! Compute derivatives of rho0, p0, rho0rhoa, gammap0pa with respect to x
    u_mm = rhoa(Im2,2); u_m = rhoa(Im1,2); u_p = rhoa(i,2); u_pp = rhoa(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhoa_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    
    u_mm = p0_prior(Im2,2) * pa(Im2,2) * gamma_chemestry(Im2,2)
    u_m  = p0_prior(Im1,2) * pa(Im1,2) * gamma_chemestry(Im1,2)
    u_p  = p0_prior(i,2)   * pa(i,2)   * gamma_chemestry(i,2)  
    u_pp = p0_prior(Ip1,2) * pa(Ip1,2) * gamma_chemestry(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dgammap0pa_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

    eq1_memory_drhoa_dx_adj(i,2)      = b_x(i_global) * eq1_memory_drhoa_dx_adj(i,2)     + a_x(i_global) * value_drhoa_dx
    eq1_memory_dgammap0pa_dx_adj(i,2) = b_x(i_global) * eq1_memory_dgammap0pa_dx_adj(i,2)    &
                                      + a_x(i_global) * value_dgammap0pa_dx


    ! Compute derivatives of windx and vax with respect to x
    if (windx_prior(i,2) <= 0) then 
    ! In the adjoint problem, the upwind discretization is reversed
    ! with respect to the forward problem.
        u_mm = vax_old(Im2,2); u_m = vax_old(Im1,2); u_p = vax_old(Ip1,2)
        call compute_decentered_backward_dU(u_mm, u_m, vax_old(i,2), u_p, value_dvax_dx, ONE_OVER_SIX_DELTAX)
     else
        u_pp = vax_old(Ip2,2)     ; u_p = vax_old(Ip1,2)     ; u_m = vax_old(Im1,2)
        call compute_decentered_forward_dU(u_m, vax_old(i,2), u_p, u_pp, value_dvax_dx, ONE_OVER_SIX_DELTAX)
     endif

    eq1_memory_dvax_dx_adj(i,2) = b_x(i_global) * eq1_memory_dvax_dx_adj(i,2) + a_x(i_global) * value_dvax_dx


    ! Compute derivatives of rho0, p0, rho0rhoa, gammap0pa with respect to y
    value_drhoa_dy = (rhoa(i,3) - rhoa(i,2)) *  ONE_OVER_DELTAY
    value_dgammap0pa_dy = (p0_prior(i,3)*pa(i,3)*gamma_chemestry(i,3) - p0_prior(i,2)*pa(i,2)*gamma_chemestry(i,2)) &
                       *  ONE_OVER_DELTAY

    u_mm = p0_prior(i,1); u_m = p0_prior(i,2); u_p = p0_prior(i,3); u_pp = p0_prior(i,4)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp0_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)

    eq1_memory_drhoa_dy_adj(i,2)     = b_y_half(2)*eq1_memory_drhoa_dy_adj(i,2)     + a_y_half(2)*value_drhoa_dy
    !eq1_memory_dp0_dy_adj(i,2)       = b_y_half(2)*eq1_memory_dp0_dy_adj(i,2)       + a_y_half(2)*value_dp0_dy
    eq1_memory_dgammap0pa_dy_adj(i,2)= b_y_half(2)*eq1_memory_dgammap0pa_dy_adj(i,2)+ a_y_half(2)*value_dgammap0pa_dy


    ! Compute derivatives of vay with respect to x
    if (windx_half_x_half_y <= 0) then 
    ! In the adjoint problem, the upwind discretization is reversed
    ! with respect to the forward problem.
        u_mm = vay_old(Im2,2); u_m = vay_old(Im1,2); u_p = vay_old(Ip1,2)
        call compute_decentered_backward_dU(u_mm, u_m, vay_old(i,2), u_p, value_dvay_dx, ONE_OVER_SIX_DELTAX)
    else
        u_pp = vay_old(Ip2,2); u_p = vay_old(Ip1,2); u_m = vay_old(Im1,2)
        call compute_decentered_forward_dU(u_m, vay_old(i,2), u_p, u_pp, value_dvay_dx, ONE_OVER_SIX_DELTAX)
    endif

    eq1_memory_dvay_dx_adj(i,2) = b_x_half(i_global) * eq1_memory_dvay_dx_adj(i,2) + a_x_half(i_global) * value_dvay_dx


    ! Compute derivatives of windx with respect to y, at vay location
    u_mm = windx_prior(i,1); u_m = windx_prior(i,2); u_p = windx_prior(i,3); u_pp = windx_prior(i,4)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = windx_prior(Ip1,1); u_m = windx_prior(Ip1,2); u_p = windx_prior(Ip1,3); u_pp = windx_prior(Ip1,4)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    value_dwindx_dy = 0.5d0 * ( value_dwindx_dy_prec + value_dwindx_dy_next)
    !eq1_memory_dwindx_dy_adj(i,2) = c_y_half(2) * value_dwindx_dy


    ! Apply the PML correction to the spatial derivatives
    value_drhoa_dx     = value_drhoa_dx     * one_over_K_x(i_global) + eq1_memory_drhoa_dx_adj(i,2)
    value_dgammap0pa_dx= value_dgammap0pa_dx* one_over_K_x(i_global) + eq1_memory_dgammap0pa_dx_adj(i,2)
    value_dvax_dx      = value_dvax_dx      * one_over_K_x(i_global) + eq1_memory_dvax_dx_adj(i,2)
    value_drhoa_dy     = value_drhoa_dy     * one_over_K_y_half(2)   + eq1_memory_drhoa_dy_adj(i,2)
    value_dgammap0pa_dy= value_dgammap0pa_dy* one_over_K_y_half(2)   + eq1_memory_dgammap0pa_dy_adj(i,2)
    value_dvay_dx      = value_dvay_dx      * one_over_K_y_half(2)   + eq1_memory_dvay_dx_adj(i,2)
    value_dp0_dy       = value_dp0_dy       * one_over_Kdalpha_y_half(2)   !+ eq1_memory_dp0_dy_adj(i,2)
    value_dwindx_dy    = value_dwindx_dy    * one_over_Kdalpha_y_half(2)   !+ eq1_memory_dwindx_dy_adj(i,2)


    ! Compute intermediate terms appearing in the adjoint expression
    value_windxdvax_dx = windx_prior(i,2) * value_dvax_dx
    value_windxdvay_dx = windx_half_x_half_y * value_dvay_dx
  

    ! Update the x-component of adjoint velocity
    vax(i,2) = vax(i,2) + value_drhoa_dx * DELTAT
    vax(i,2) = vax(i,2) + value_windxdvax_dx * DELTAT
    vax(i,2) = vax(i,2) + value_dgammap0pa_dx * DELTAT / rho0_half_x

    ! Update the y-component of adjoint velocity
    vay(i,2) = vay(i,2) + value_drhoa_dy * DELTAT
    vay(i,2) = vay(i,2) + value_windxdvay_dx * DELTAT
    vay(i,2) = vay(i,2) - (vax_half_x_half_y * value_dwindx_dy) * DELTAT
    vay(i,2) = vay(i,2) + (value_dgammap0pa_dy - pa_half_y * value_dp0_dy) * DELTAT / rho0_half_y

  enddo

endsubroutine derivative_first_row_vaxvay





subroutine derivative_first_row_parhoa(adjoint_source_term)
!======================================================================
!> Update the adjoint density and pressure on the first grid row.
!!
!! This routine advances the adjoint continuity and pressure equations on
!! the reflective-ground boundary. Spatial derivatives are computed with
!! centered or reversed-upwind finite differences.
!!
!! Arguments:
!!   adjoint_source_term : Adjoint source contribution.
!!======================================================================
use parameters

implicit none

  ! Source term
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: adjoint_source_term

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: rhoa_old, pa_old

  ! Derivatives for equation on vax and vay
  double precision ::       &
    value_drhoa_dx,         &
    value_dvax_dx,      value_dvay_dy

  ! Derivative for equation on rhoa and pa
  double precision ::     &
    value_dpa_dx

  ! Interpolated quantities
  double precision ::                         &
    vax_half_x,          vay_half_y,          &
    windx_half_x
  
  ! Compute global and local neighbor indices
  integer :: i
  integer :: Ip1,Im1, Ip2, Im2
  double precision :: u_mm, u_m, u_p, u_pp


  !!!!!!!!!!!!!!!!!!!!!!!
  ! Update adjoint pressure and density
  !!!!!!!!!!!!!!!!!!!!!!!

  ! Store the current fields before the in-place update.
  pa_old(:,:) = pa(:,:)
  rhoa_old(:,:) = rhoa(:,:)


  do i=1,NX_LOCAL

    ! Compute global and local neighbor indices
    i_global = i + offset_i

    Im1 = i-1
    Im2 = i-2
    Ip1 = i+1
    Ip2 = i+2


    ! Compute quantities interpolated at staggered-grid locations.
    vax_half_x = 0.5d0 * (vax(i,2) + vax(Ip1,2))
    vay_half_y = 0.5d0 * (vay(i,2) + vay(i,1))
    windx_half_x = 0.5d0 * (windx_prior(i,2) + windx_prior(Ip1,2))

    ! Compute derivatives of rhoa and pa with respect to x
    if (windx_half_x <= 0) then 
    ! In the adjoint problem, the upwind discretization is reversed
    ! with respect to the forward problem.
         u_mm = pa_old(Im2,2)      ; u_m = pa_old(Im1,2)      ; u_p = pa_old(Ip1,2)
         call compute_decentered_backward_dU(u_mm, u_m, pa_old(i,2), u_p, value_dpa_dx, ONE_OVER_SIX_DELTAX)
         u_mm = rhoa_old(Im2,2); u_m = rhoa_old(Im1,2); u_p = rhoa_old(Ip1,2)
         call compute_decentered_backward_dU(u_mm, u_m, rhoa_old(i,2), u_p, value_drhoa_dx, ONE_OVER_SIX_DELTAX)
    else
         u_pp = pa_old(Ip2,2)      ; u_p = pa_old(Ip1,2)      ; u_m = pa_old(Im1,2)
         call compute_decentered_forward_dU(u_m, pa_old(i,2), u_p, u_pp, value_dpa_dx, ONE_OVER_SIX_DELTAX)
         u_pp = rhoa_old(Ip2,2); u_p = rhoa_old(Ip1,2); u_m = rhoa_old(Im1,2)
         call compute_decentered_forward_dU(u_m, rhoa_old(i,2), u_p, u_pp, value_drhoa_dx, ONE_OVER_SIX_DELTAX)
    endif

    eq2_memory_dpa_dx_adj(i,2)   = b_x_half(i_global)*eq2_memory_dpa_dx_adj(i,2)   &
                                       + a_x_half(i_global)*value_dpa_dx
    eq2_memory_drhoa_dx_adj(i,2) = b_x_half(i_global)*eq2_memory_drhoa_dx_adj(i,2) &
                                       + a_x_half(i_global)*value_drhoa_dx


     ! Compute derivatives of vax with respect to x
     u_mm = vax(Im1,2); u_m = vax(i,2); u_p = vax(Ip1,2); u_pp = vax(Ip2,2)
     call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvax_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

     eq2_memory_dvax_dx_adj(i,2)   = b_x_half(i_global) * eq2_memory_dvax_dx_adj(i,2)   + a_x_half(i_global) * value_dvax_dx


     ! Compute derivatives of vay with respect to y
     u_mm = vay(i,0); u_m = vay(i,1); u_p = vay(i,2); u_pp = vay(i,3)
     call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvay_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)

     eq2_memory_dvay_dy_adj(i,2)   = b_y(2) * eq2_memory_dvay_dy_adj(i,2)   + a_y(2) * value_dvay_dy


     ! Apply PML correction to spatial derivatives
     value_dpa_dx   = value_dpa_dx    * one_over_K_x_half(i_global) + eq2_memory_dpa_dx_adj(i,2)
     value_drhoa_dx = value_drhoa_dx  * one_over_K_x_half(i_global) + eq2_memory_drhoa_dx_adj(i,2)
     
     value_dvax_dx  = value_dvax_dx * one_over_K_x_half(i_global) + eq2_memory_dvax_dx_adj(i,2)
     value_dvay_dy  = value_dvay_dy * one_over_K_y(2)             + eq2_memory_dvay_dy_adj(i,2)

   
     ! Update adjoint density
     rhoa(i,2) = rhoa(i,2) + (windx_half_x * value_drhoa_dx)  * DELTAT
     rhoa(i,2) = rhoa(i,2) - vay_half_y * g(i,2)              * DELTAT


     ! Update adjoint pressure
     pa(i,2) = pa(i,2) + (windx_half_x * value_dpa_dx)     * DELTAT
     pa(i,2) = pa(i,2) + (value_dvax_dx + value_dvay_dy)   * DELTAT
     pa(i,2) = pa(i,2) + adjoint_source_term(i,2)          * DELTAT


  enddo

endsubroutine derivative_first_row_parhoa




!======================================================================
! Sensitivity-kernel routines
!======================================================================

subroutine derivative_first_row_kernel_kvx(rho0, p0, windx, windy, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
              value_dvx_dt, value_dvy_dt, it_time, i_start, i_end)
!======================================================================
!> Compute the wind-x sensitivity kernel on the first grid row.
!!
!! This routine evaluates the contributions to the horizontal wind kernel
!! associated with the reflective-ground boundary. The required forward
!! and adjoint spatial derivatives are evaluated on the staggered grid
!! and accumulated into K_windx.
!!
!! The kernel is accumulated over the current time step.
!!======================================================================
   use parameters, only :  K_windx,                    &
                          DELTAT, NX_LOCAL, NY_LOCAL,  &
                          NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX,    &
                          NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,    &
                          ONE_OVER_SIX_DELTAX,     &
                          ONE_OVER_DELTAY,         &
                          i_global, offset_i
  implicit none

  integer :: i
  integer :: it_time
  integer :: i_start, i_end


  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: &
     rho0, p0, windx, windy,                    &
     rhop, pressure, vx, vy,                    &
     rhoa, pa, vax, vay,                        &
     value_dvx_dt, value_dvy_dt

  double precision :: &
      value_dvx_dx,                       &
      value_dvy_dx,                       &
      value_dvax_dy

  double precision ::   &
      value_drhop_dx,                         &
      value_dp_dx,                            &
      value_vax_dvx_dx,                       &
      value_vay_dvy_dx

  double precision :: &
    rho0_half_x,                              &
    rhoa_half_x,                              &
    rhop_half_x,                              &
    pa_half_x,                                &
    vay_half_x_half_y,   &
    vy_half_x_half_y,    &
    value_dvy_dy_half_x_half_y

  double precision :: &
    value_dvy_dy_next,    value_dvy_dy_prec,       &
    value_dvy_dx_next, value_dvy_dx_prec


  double precision  :: u_mm, u_m, u_p, u_pp
  integer :: Ip1,Im1, Ip2, Im2


  do i =i_start,i_end

    ! Compute global and local neighbor indices
    i_global = i + offset_i

    Im1 = i-1
    Im2 = i-2
    Ip1 = i+1
    Ip2 = i+2


    ! Interpolate variables to the staggered-grid locations required by the
    ! finite-difference formulation.
    rhoa_half_x        = 0.5d0 * (rhoa(i,2) + rhoa(Im1,2))
    rho0_half_x        = 0.5d0 * (rho0(i,2) + rho0(Im1,2))
    rhop_half_x        = 0.5d0 * (rhop(i,2) + rhop(Im1,2))
    pa_half_x          = 0.5d0 * (pa(i,2)   + pa(Im1,2))

    vay_half_x_half_y = 0.25d0 * (vay(i,2) + vay(Im1,2)+ vay(i,1) + vay(Im1,1))
    vy_half_x_half_y  = 0.25d0 * (vy(i,2) + vy(Im1,2)+ vy(i,1) + vy(Im1,1))


    ! Compute derivatives of rhoprhoa, gammapap and pressure with respect to x
    u_mm = pressure(Im2,2); u_m = pressure(Im1,2); u_p = pressure(i,2); u_pp = pressure(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    
    u_mm = rhop(Im2,2); u_m = rhop(Im1,2); u_p = rhop(i,2); u_pp = rhop(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhop_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    

    ! Compute derivatives of vx with respect to x
    if (windx(i,2) >= 0) then
        u_mm = vx(Im2,2); u_m = vx(Im1,2); u_p = vx(Ip1,2)
        call compute_decentered_backward_dU(u_mm, u_m, vx(i,2), u_p, value_dvx_dx, ONE_OVER_SIX_DELTAX)
    else
        u_pp = vx(Ip2,2)     ; u_p = vx(Ip1,2)     ; u_m = vx(Im1,2)
        call compute_decentered_forward_dU(u_m, vx(i,2), u_p, u_pp, value_dvx_dx, ONE_OVER_SIX_DELTAX)
    endif


    ! Compute derivatives of vy with respect to x
    u_mm = vy(Im2,1); u_m = vy(Im1,1); u_p = vy(i,1); u_pp = vy(Ip1,1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = vy(Im2,2); u_m = vy(Im1,2); u_p = vy(i,2); u_pp = vy(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    value_dvy_dx = 0.5d0 * ( value_dvy_dx_prec + value_dvy_dx_next)


    ! Compute derivatives of vy with respect to y
    u_mm = vy(Im1,0); u_m = vy(Im1,1); u_p = vy(Im1,2); u_pp = vy(Im1,3)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = vy(i,0); u_m = vy(i,1); u_p = vy(i,2); u_pp = vy(i,3)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    value_dvy_dy_half_x_half_y = 0.5d0 * ( value_dvy_dy_prec + value_dvy_dy_next)


    ! Compute derivatives of vax with respect to y
    u_pp = vax(i,4)     ; u_p = vax(i,3)
    call compute_decentered_forward_dU_o2(vax(i,2), u_p, u_pp, value_dvax_dy, ONE_OVER_DELTAY)


    ! Compute intermediate terms appearing in the adjoint expression
    value_vax_dvx_dx = vax(i,2) * value_dvx_dx
    value_vay_dvy_dx = vay_half_x_half_y * value_dvy_dx

  
    ! Update the horizontal wind sensitivity kernel
    K_windx(i,2) = K_windx(i,2) + (rhoa_half_x     * value_drhop_dx)                        * DELTAT
    K_windx(i,2) = K_windx(i,2) + rho0_half_x      * (value_vax_dvx_dx + value_vay_dvy_dx)  * DELTAT
    K_windx(i,2) = K_windx(i,2) - vax(i,2) * rho0_half_x * value_dvy_dy_half_x_half_y       * DELTAT
    K_windx(i,2) = K_windx(i,2) - vy_half_x_half_y * rho0_half_x * value_dvax_dy            * DELTAT
    K_windx(i,2) = K_windx(i,2) + (pa_half_x       * value_dp_dx)                           * DELTAT

enddo

endsubroutine derivative_first_row_kernel_kvx


subroutine derivative_first_row_kernel_kp0krho0(rho0, p0, windx, windy, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
              value_dvx_dt, value_dvy_dt, it_time, i_start, i_end, save_source_term)
!======================================================================
!> Compute the density and pressure sensitivity kernels on the first grid row.
!!
!! This routine evaluates the contributions to the background density
!! (K_rho0) and pressure (K_p0) sensitivity kernels at the reflective
!! boundary. Source-term contributions are reconstructed consistently
!! with the forward simulation before being accumulated into the kernels.
!!
!! Arguments:
!!   save_source_term : Source amplitude stored from the forward run.
!!======================================================================
 use parameters, only :  K_rho0, K_p0, gamma_chemestry,             &
                          NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX,    &
                          NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,    &
                          ONE_OVER_SIX_DELTAX, ONE_OVER_SIX_DELTAY, &
                          ONE_OVER_DELTAY,                          &
                          DELTAT, DELTAX, DELTAY,                   &
                          ISOURCE, JSOURCE, source_term,            &
                          type_source, wavefront,                   &
                          distance2, factor_ssf, SSF_Sigma,         &
                          NX_LOCAL, NY_LOCAL, i_global, offset_i
  implicit none

  integer :: i
  integer :: it_time
  integer :: i_start, i_end

  double precision :: save_source_term

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: &
     rho0, p0, windx, windy,                    &
     rhop, pressure, vx, vy,                    &
     rhoa, pa, vax, vay,                        &
     value_dvx_dt, value_dvy_dt

  double precision :: &
      value_dvx_dx,                       &
      value_dvy_dx,     value_dvy_dy,     &
      value_dwindx_dy,                    &
      value_vax_windx_dvx_dx,             &
      value_vay_windx_dvy_dx,             &
      value_vax_vy_dwindx_dy

  double precision ::   &
      value_drhoa_dy,     &
      value_dpa_dy

  double precision :: &
    windx_half_x,                             &
    vax_half_x,          vay_half_y,          &
    vx_half_x,           vy_half_y,           &
    vax_dvx_dt_half_x, vay_dvy_dt_half_y

   double precision :: &
       value_dwindx_dy_next, value_dwindx_dy_prec,          &
       value_dvy_dx_next,    value_dvy_dx_prec

  double precision  :: u_mm, u_m, u_p, u_pp
  integer :: Ip1,Im1, Ip2, Im2


  do i =i_start,i_end

    ! Compute global and local neighbor indices
    i_global = i + offset_i

    Im1 = i-1
    Im2 = i-2
    Ip1 = i+1
    Ip2 = i+2


    ! Interpolate variables to the staggered-grid locations required by the
    ! finite-difference formulation.
    vx_half_x          = 0.5d0 * (vx(i,2)  + vx(Ip1,2))
    vax_half_x         = 0.5d0 * (vax(i,2) + vax(Ip1,2))
    windx_half_x       = 0.5d0 * (windx(i,2) + windx(Ip1,2))
    vax_dvx_dt_half_x  = 0.5d0 * (vax(i,2)*value_dvx_dt(i,2) + vax(Ip1,2)*value_dvx_dt(Ip1,2))

    vy_half_y          = 0.5d0 * (vy(i,2)  + vy(i,1))
    vay_half_y         = 0.5d0 * (vay(i,2) + vay(i,1))
    vay_dvy_dt_half_y  = 0.5d0 * (vay(i,2)*value_dvy_dt(i,2) + vay(i,1)*value_dvy_dt(i,1))


    ! Evaluate the source term
    source_term = 0.d0
    if (type_source == 1) then
      if ((wavefront == 1 .and. i_global == ISOURCE ) .or. (wavefront == 2 .and. 2 == JSOURCE)) then
        factor_ssf = 1.d0
        source_term = save_source_term
      endif

    elseif (type_source ==2) then
      if (i_global == ISOURCE .and. 2 == JSOURCE) then
        factor_ssf = 1.d0
        source_term = save_source_term
      endif

    elseif (type_source == 3) then
        distance2 = ((i_global - Isource) * DELTAX)**2 + ((2 - Jsource) * DELTAY)**2
        factor_ssf = exp( - distance2 / SSF_Sigma**2 )
        source_term = save_source_term
    endif


    ! Compute derivatives of vx with respect to x
    u_mm = vx(Im1,2); u_m = vx(i,2); u_p = vx(Ip1,2); u_pp = vx(Ip2,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)


    ! Compute derivatives of vy with respect to y
    u_mm = vy(i,0); u_m = vy(i,1); u_p = vy(i,2); u_pp = vy(i,3)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)


    ! Compute derivatives of vy with respect to x
    if (windx_half_x >= 0) then
          u_mm = vy(Im2,1)      ; u_m = vy(Im1,1)      ; u_p = vy(Ip1,1)
          call compute_decentered_backward_dU(u_mm, u_m, vy(i,1), u_p, value_dvy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_mm = vy(Im2,2)      ; u_m = vy(Im1,2)      ; u_p = vy(Ip1,2)
          call compute_decentered_backward_dU(u_mm, u_m, vy(i,2), u_p, value_dvy_dx_next, ONE_OVER_SIX_DELTAX)

    else
          u_pp = vy(Ip2,1)      ; u_p = vy(Ip1,1)      ; u_m = vy(Im1,1)
          call compute_decentered_forward_dU(u_m, vy(i,1), u_p, u_pp, value_dvy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_pp = vy(Ip2,2)      ; u_p = vy(Ip1,2)      ; u_m = vy(Im1,2)
          call compute_decentered_forward_dU(u_m, vy(i,2), u_p, u_pp, value_dvy_dx_next, ONE_OVER_SIX_DELTAX)
    endif
    value_dvy_dx = 0.5d0 * (value_dvy_dx_prec + value_dvy_dx_next)


    ! Compute derivatives of rhoa, pa with respect to y
    u_pp = pa(i,4)      ; u_p = pa(i,3)
    call compute_decentered_forward_dU_o2(pa(i,2), u_p, u_pp, value_dpa_dy, ONE_OVER_DELTAY)

    u_pp = rhoa(i,4)      ; u_p = rhoa(i,3)
    call compute_decentered_forward_dU_o2(rhoa(i,2), u_p, u_pp, value_drhoa_dy, ONE_OVER_DELTAY)


    ! Compute derivatives of windx, vx with respect to y
    u_mm = windx(i,0)      ; u_m = windx(i,1)      ; u_p = windx(i,3)
    call compute_decentered_backward_dU(u_mm, u_m, windx(i,2), u_p, value_dwindx_dy_prec, ONE_OVER_SIX_DELTAY)
    u_mm = windx(Ip1,0)      ; u_m = windx(Ip1,1)      ; u_p = windx(Ip1,3)
    call compute_decentered_backward_dU(u_mm, u_m, windx(Ip1,2), u_p, value_dwindx_dy_next, ONE_OVER_SIX_DELTAY)

    value_dwindx_dy = 0.5d0 * (value_dwindx_dy_prec + value_dwindx_dy_next)


     ! Compute intermediate terms appearing in the kernel expression
    value_vax_windx_dvx_dx = vax_half_x * windx_half_x * value_dvx_dx
    value_vay_windx_dvy_dx = vay_half_y * windx_half_x * value_dvy_dx
    value_vax_vy_dwindx_dy = vax_half_x * vy_half_y * value_dwindx_dy


    ! Update the density sensitivity kernel
    K_rho0(i,2) = K_rho0(i,2) + (rhoa(i,2) * value_dvx_dx - vy_half_y * value_drhoa_dy) * DELTAT
    K_rho0(i,2) = K_rho0(i,2) - (vax_dvx_dt_half_x + vay_dvy_dt_half_y) * DELTAT
    K_rho0(i,2) = K_rho0(i,2) + (value_vax_windx_dvx_dx + value_vay_windx_dvy_dx) * DELTAT
    K_rho0(i,2) = K_rho0(i,2) + (value_vax_vy_dwindx_dy) * DELTAT
    K_rho0(i,2) = K_rho0(i,2) + pa(i,2) * factor_ssf * source_term * (gamma_chemestry(i,2) * p0(i,2) / rho0(i,2)**2) * DELTAT


    ! Update the pressure sensitivity kernel
    K_p0(i,2) = K_p0(i,2) - (vy_half_y * value_dpa_dy)                                         * DELTAT
    K_p0(i,2) = K_p0(i,2) + gamma_chemestry(i,2) * pa(i,2) * (value_dvx_dx + value_dvy_dy)        * DELTAT
    K_p0(i,2) = K_p0(i,2) - pa(i,2) * value_dvy_dy                                             * DELTAT
    K_p0(i,2) = K_p0(i,2) - pa(i,2) * factor_ssf * source_term * gamma_chemestry(i,2) / rho0(i,2) * DELTAT

  enddo

endsubroutine derivative_first_row_kernel_kp0krho0
