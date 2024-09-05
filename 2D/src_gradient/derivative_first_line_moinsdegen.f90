!!!!!!!!!!!!!!!!!!!!!!!
!!  Specific functions of adjoint and sensitivity kernel computation for the first row to model reflective ground
!!!!!!!!!!!!!!!!!!!!!!!


!---------------
! ADJOINT 
!---------------
subroutine derivative_first_line_vaxvay()

  use parameters
  implicit none

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: vax_old, vay_old
   
  ! derivative for equation on vax and vay
  double precision ::                       &
    value_drhoa_dx,     value_drhoa_dy,     &
    value_dgammap0pa_dx,     value_dgammap0pa_dy,     &
    value_dp0_dx,       value_dp0_dy,       &
    value_dwindx_dx,    value_dwindx_dy,    &
    value_dwindy_dx,    value_dwindy_dy,     &
    value_dvax_dx,      value_dvax_dy,      &
    value_dvay_dx,      value_dvay_dy,      &
    value_windxdvax_dx, value_windydvax_dy, &
    value_windxdvay_dx, value_windydvay_dy

  ! interpolated variable
  double precision ::                         &
    rho0_half_x,         rho0_half_y,         &
    rhoa_half_x,         rhoa_half_y,         &
    pa_half_x,           pa_half_y,           &
    vax_half_x_half_y,   vay_half_x_half_y,   &
    windx_half_x_half_y, windy_half_x_half_y
    
  ! intermediate variable to create interpolation on derivative 
  double precision :: &
    value_dwindx_dy_next, value_dwindx_dy_prec,&
    value_dwindy_dx_next, value_dwindy_dx_prec 

  ! index
  integer :: i
  integer :: Ip1,Im1, Ip2, Im2
  double precision :: u_mm, u_m, u_p, u_pp
  
  
  !!!!!!!!!!!!!!!!!!!!!!!
  ! vax and vay computations
  !!!!!!!!!!!!!!!!!!!!!!!
     
  ! save value of previous iterations in an auxiliary variable    
  vax_old(:,:) = vax(:,:)
  vay_old(:,:) = vay(:,:) 
  
  do i=1,NX_LOCAL

    ! Index
    i_global = i + offset_i
    
    Im1 = i-1
    Im2 = i-2
    Ip1 = i+1
    Ip2 = i+2
 
 
    ! useful interpolation
    rhoa_half_x = 0.5d0 * (rhoa(i,2) + rhoa(Im1,2))
    rho0_half_x  = 0.5d0 * (rho0_prior(i,2)  + rho0_prior(Im1,2))
    pa_half_x   = 0.5d0 * (pa(i,2)   + pa(Im1,2))
    vay_half_x_half_y = 0.25d0 * (vay_old(i,2) + vay_old(Im1,2) + vay_old(i,1) + vay_old(Im1,1))
    windy_half_x_half_y = 0.25d0 * (windy_prior(i,2)     + windy_prior(Im1,2)     + windy_prior(i,1)     + windy_prior(Im1,1))

    rhoa_half_y = 0.5d0 * (rhoa(i,2) + rhoa(i,3))
    rho0_half_y  = 0.5d0 * (rho0_prior(i,2)  + rho0_prior(i,3))
    pa_half_y   = 0.5d0 * (pa(i,2)   + pa(i,3))
    vax_half_x_half_y = 0.25d0 * (vax_old(i,2) + vax_old(Ip1,2) + vax_old(i,3) + vax_old(Ip1,3))
    windx_half_x_half_y = 0.25d0 * (windx_prior(i,2)     + windx_prior(Ip1,2)     + windx_prior(i,3)     + windx_prior(Ip1,3))


    ! compute derivative of rho0, p0, rho0rhoa, gammap0pa accroding to x
    u_mm = rhoa(Im2,2); u_m = rhoa(Im1,2); u_p = rhoa(i,2); u_pp = rhoa(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhoa_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = p0_prior(Im2,2); u_m = p0_prior(Im1,2); u_p = p0_prior(i,2); u_pp = p0_prior(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp0_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = p0_prior(Im2,2) * pa(Im2,2) * gamma_chimie(Im2,2); u_m = p0_prior(Im1,2) * pa(Im1,2) * gamma_chimie(Im1,2)
    u_p  = p0_prior(i,2) * pa(i,2) * gamma_chimie(i,2)   ; u_pp = p0_prior(Ip1,2) * pa(Ip1,2) * gamma_chimie(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dgammap0pa_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

    eq1_memory_drhoa_dx_adj(i,2)     = b_x(i_global) * eq1_memory_drhoa_dx_adj(i,2)     + a_x(i_global) * value_drhoa_dx
    eq1_memory_dp0_dx_adj(i,2)       = b_x(i_global) * eq1_memory_dp0_dx_adj(i,2)       + a_x(i_global) * value_dp0_dx
    eq1_memory_dgammap0pa_dx_adj(i,2)     = b_x(i_global) * eq1_memory_dgammap0pa_dx_adj(i,2)    &
                                      + a_x(i_global) * value_dgammap0pa_dx

    ! compute derivative of  vax according to y
    u_pp = vax_old(i,4); u_p = vax_old(i,3)
    call compute_decentered_forward_dU_o2(vax_old(i,2), u_p, u_pp, value_dvax_dy, ONE_OVER_DELTAY)

    eq1_memory_dvax_dy_adj(i,2) = b_y(2) * eq1_memory_dvax_dy_adj(i,2) + a_y(2) * value_dvax_dy


    ! derivative of windx and vax according to x
    if (windx_prior(i,2) <= 0) then !! Solene a change a cause du signe
        u_mm = windx_prior(Im2,2)      ; u_m = windx_prior(Im1,2)      ; u_p = windx_prior(Ip1,2)
        call compute_decentered_backward_dU(u_mm, u_m, windx_prior(i,2), u_p, value_dwindx_dx, ONE_OVER_SIX_DELTAX)
        u_mm = vax_old(Im2,2); u_m = vax_old(Im1,2); u_p = vax_old(Ip1,2)
        call compute_decentered_backward_dU(u_mm, u_m, vax_old(i,2), u_p, value_dvax_dx, ONE_OVER_SIX_DELTAX)
     else
        u_pp = windx_prior(Ip2,2)    ; u_p = windx_prior(Ip1,2)    ; u_m = windx_prior(Im1,2)
        call compute_decentered_forward_dU(u_m, windx_prior(i,2), u_p, u_pp, value_dwindx_dx, ONE_OVER_SIX_DELTAX)
        u_pp = vax_old(Ip2,2)     ; u_p = vax_old(Ip1,2)     ; u_m = vax_old(Im1,2)
        call compute_decentered_forward_dU(u_m, vax_old(i,2), u_p, u_pp, value_dvax_dx, ONE_OVER_SIX_DELTAX)
     endif

    eq1_memory_dwindx_dx_adj(i,2) = c_x(i_global) * value_dwindx_dx
    eq1_memory_dvax_dx_adj(i,2) = b_x(i_global) * eq1_memory_dvax_dx_adj(i,2) + a_x(i_global) * value_dvax_dx


    ! compute derivative of windy according to x (evaluating in vax)
    u_mm = windy_prior(Im2,1); u_m = windy_prior(Im1,1); u_p = windy_prior(i,1); u_pp = windy_prior(Ip1,1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = windy_prior(Im2,2); u_m = windy_prior(Im1,2); u_p = windy_prior(i,2); u_pp = windy_prior(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    value_dwindy_dx = 0.5d0 * ( value_dwindy_dx_prec + value_dwindy_dx_next)

    eq1_memory_dwindy_dx_adj(i,2) = c_x(i_global) * value_dwindy_dx


    ! compute derivative of rho0, p0, rho0rhoa, gammap0pa accroding to y
    value_drhoa_dy = (rhoa(i,3) - rhoa(i,2)) *  ONE_OVER_DELTAY
    value_dgammap0pa_dy = (p0_prior(i,3)*pa(i,3)*gamma_chimie(i,3) - p0_prior(i,2)*pa(i,2)*gamma_chimie(i,2)) *  ONE_OVER_DELTAY
  
    u_mm = p0_prior(i,1); u_m = p0_prior(i,2); u_p = p0_prior(i,3); u_pp = p0_prior(i,4)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp0_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
  
    eq1_memory_drhoa_dy_adj(i,2)     = b_y_half(2)*eq1_memory_drhoa_dy_adj(i,2)     + a_y_half(2)*value_drhoa_dy
    eq1_memory_dp0_dy_adj(i,2)       = b_y_half(2)*eq1_memory_dp0_dy_adj(i,2)       + a_y_half(2)*value_dp0_dy
    eq1_memory_dgammap0pa_dy_adj(i,2)= b_y_half(2)*eq1_memory_dgammap0pa_dy_adj(i,2)+ a_y_half(2)*value_dgammap0pa_dy


    ! compute derivative of vay according to x
    if (windx_half_x_half_y <= 0) then !! Solene a change a cause du signe
        u_mm = vay_old(Im2,2); u_m = vay_old(Im1,2); u_p = vay_old(Ip1,2)
        call compute_decentered_backward_dU(u_mm, u_m, vay_old(i,2), u_p, value_dvay_dx, ONE_OVER_SIX_DELTAX)
    else
        u_pp = vay_old(Ip2,2); u_p = vay_old(Ip1,2); u_m = vay_old(Im1,2)
        call compute_decentered_forward_dU(u_m, vay_old(i,2), u_p, u_pp, value_dvay_dx, ONE_OVER_SIX_DELTAX)
    endif

    eq1_memory_dvay_dx_adj(i,2) = b_x_half(i_global) * eq1_memory_dvay_dx_adj(i,2) + a_x_half(i_global) * value_dvay_dx


    ! compute derivative of windy, vay according to y
    if (windy_prior(i,2) <= 0) then  !! Solene a change a cause du signe
          u_mm = windy_prior(i,0)      ; u_m = windy_prior(i,1)      ; u_p = windy_prior(i,3) 
          call compute_decentered_backward_dU(u_mm, u_m, windy_prior(i,2), u_p, value_dwindy_dy, ONE_OVER_SIX_DELTAY) 
    else 
          u_pp = windy_prior(i,4)      ; u_p = windy_prior(i,3)      ; u_m = windy_prior(i,1) 
          call compute_decentered_forward_dU(u_m, windy_prior(i,2), u_p, u_pp, value_dwindy_dy, ONE_OVER_SIX_DELTAY)
    endif

    u_pp = vay_old(i,4); u_p = vay_old(i,3)
    call compute_decentered_forward_dU_o2(vay_old(i,2), u_p, u_pp, value_dvay_dy, ONE_OVER_DELTAY)

    eq1_memory_dwindy_dy_adj(i,2) = c_y_half(2) * value_dwindy_dy
    eq1_memory_dvay_dy_adj(i,2) = b_y_half(2) * eq1_memory_dvay_dy_adj(i,2) + a_y_half(2) * value_dvay_dy


    ! computing derivative of windx according to y, evaluating in vay
    u_mm = windx_prior(i,1); u_m = windx_prior(i,2); u_p = windx_prior(i,3); u_pp = windx_prior(i,4)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = windx_prior(Ip1,1); u_m = windx_prior(Ip1,2); u_p = windx_prior(Ip1,3); u_pp = windx_prior(Ip1,4)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    value_dwindx_dy = 0.5d0 * ( value_dwindx_dy_prec + value_dwindx_dy_next)

    eq1_memory_dwindx_dy_adj(i,2) = c_y_half(2) * value_dwindx_dy


    ! Treatment of derivatives for PML 
    value_drhoa_dx     = value_drhoa_dx     * one_over_K_x(i_global) + eq1_memory_drhoa_dx_adj(i,2)
    value_dp0_dx       =  value_dp0_dx      * one_over_K_x(i_global) + eq1_memory_dp0_dx_adj(i,2)
    value_dgammap0pa_dx= value_dgammap0pa_dx* one_over_K_x(i_global) + eq1_memory_dgammap0pa_dx_adj(i,2)

    value_dvax_dy      = value_dvax_dy      * one_over_K_x(i_global) + eq1_memory_dvax_dy_adj(i,2)

    value_dwindx_dx    = value_dwindx_dx    * one_over_K_x(i_global) + eq1_memory_dwindx_dx_adj(i,2)
    value_dvax_dx      = value_dvax_dx      * one_over_K_x(i_global) + eq1_memory_dvax_dx_adj(i,2)

    value_dwindy_dx    = value_dwindy_dx    * one_over_K_x(i_global) + eq1_memory_dwindy_dx_adj(i,2)

    value_drhoa_dy     = value_drhoa_dy     * one_over_K_y_half(2) + eq1_memory_drhoa_dy_adj(i,2)
    value_dp0_dy       = value_dp0_dy       * one_over_K_y_half(2) + eq1_memory_dp0_dy_adj(i,2)
    value_dgammap0pa_dy     = value_dgammap0pa_dy     * one_over_K_y_half(2) + eq1_memory_dgammap0pa_dy_adj(i,2)

    value_dvay_dx      = value_dvay_dx      * one_over_K_y_half(2) + eq1_memory_dvay_dx_adj(i,2)

    value_dwindy_dy    = value_dwindy_dy    * one_over_K_y_half(2) + eq1_memory_dwindy_dy_adj(i,2)
    value_dvay_dy      = value_dvay_dy      * one_over_K_y_half(2) + eq1_memory_dvay_dy_adj(i,2)

    value_dwindx_dy    = value_dwindx_dy    * one_over_K_y_half(2) + eq1_memory_dwindx_dy_adj(i,2)


    ! intermediate computations
    value_windxdvax_dx = windx_prior(i,2) * value_dvax_dx
    value_windydvax_dy = windy_half_x_half_y * value_dvax_dy

    value_windxdvay_dx = windx_half_x_half_y * value_dvay_dx
    value_windydvay_dy = windy_prior(i,2) * value_dvay_dy


    ! x component
    vax(i,2) = vax(i,2) +  value_drhoa_dx * DELTAT
    vax(i,2) = vax(i,2) + (value_windxdvax_dx + value_windydvax_dy) * DELTAT
    vax(i,2) = vax(i,2) - (vax_old(i,2) * value_dwindx_dx + vay_half_x_half_y * value_dwindy_dx) * DELTAT
    vax(i,2) = vax(i,2) + (value_dgammap0pa_dx - pa_half_x * value_dp0_dx) * DELTAT / rho0_half_x

    ! y component
    vay(i,2) = vay(i,2) +  value_drhoa_dy * DELTAT
    vay(i,2) = vay(i,2) + (value_windxdvay_dx + value_windydvay_dy) * DELTAT
    vay(i,2) = vay(i,2) - (vax_half_x_half_y * value_dwindx_dy + vay_old(i,2) * value_dwindy_dy) * DELTAT
    vay(i,2) = vay(i,2) + (value_dgammap0pa_dy - pa_half_y * value_dp0_dy) * DELTAT / rho0_half_y
  
  enddo
  
endsubroutine derivative_first_line_vaxvay





subroutine derivative_first_line_parhoa(adjoint_source_term)

use parameters

implicit none

  ! source term
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: adjoint_source_term
  
 double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: rhoa_old, pa_old
   
  ! derivative for equation on vax and vay
  double precision ::                       &
    value_drhoa_dx,     value_drhoa_dy,     &
    value_dwindx_dx,    value_dwindx_dy,    &
    value_dwindy_dx,    value_dwindy_dy,     &
    value_dvax_dx,      value_dvay_dy
    
  ! derivatie for equation on rhoa and pa
  double precision ::          &
    value_vax_windx_dwindx_dx, & 
    value_vax_windy_dwindx_dy, &
    value_vay_windx_dwindy_dx, &
    value_vay_windy_dwindy_dy, &
    value_dpa_dx,         &
    value_dpa_dy

  ! interpolated variable
  double precision ::                         &
    vax_half_x,          vay_half_y,          &
    windx_half_x,        windy_half_y
  ! intermediate variable to create interpolation on derivative 
  double precision :: &
    value_dwindx_dy_next, value_dwindx_dy_prec,&
    value_dwindy_dx_next, value_dwindy_dx_prec 

  ! index
  integer :: i
  integer :: Ip1,Im1, Ip2, Im2
  double precision :: u_mm, u_m, u_p, u_pp
  
  
  !!!!!!!!!!!!!!!!!!!!!!!
  ! pa and rhoa computations
  !!!!!!!!!!!!!!!!!!!!!!!
  
  ! save value of previous iterations in an auxiliary variable    
  pa_old(:,:) = pa(:,:)
  rhoa_old(:,:) = rhoa(:,:)
  
  
  do i=1,NX_LOCAL

    ! Index
    i_global = i + offset_i
        
    Im1 = i-1
    Im2 = i-2
    Ip1 = i+1
    Ip2 = i+2
   
   
    ! Interpolations   
    vax_half_x = 0.5d0 * (vax(i,2) + vax(Ip1,2))
    vay_half_y = 0.5d0 * (vay(i,2) + vay(i,1))
    windx_half_x = 0.5d0 * (windx_prior(i,2) + windx_prior(Ip1,2))
    windy_half_y = 0.5d0 * (windy_prior(i,2) + windy_prior(i,1))

  
    ! compute derivative of rhoa and pa according to x
    if (windx_half_x <= 0) then !! Solene a change a cause du signe
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
   
   
    ! compute derivative of rhoa and p according to y
    u_pp = pa_old(i,4)      ; u_p = pa_old(i,3)  
    call compute_decentered_forward_dU_o2(pa_old(i,2), u_p, u_pp, value_dpa_dy, ONE_OVER_DELTAY)
   
    if (windy_half_y <= 0) then  !! Solene a change a cause du signe
          u_mm = rhoa_old(i,0); u_m = rhoa_old(i,1); u_p = rhoa_old(i,3) 
          call compute_decentered_backward_dU(u_mm, u_m, rhoa_old(i,2), u_p, value_drhoa_dy, ONE_OVER_SIX_DELTAY)
     else 
          
          u_pp = rhoa_old(i,4); u_p = rhoa_old(i,3); u_m = rhoa_old(i,1) 
          call compute_decentered_forward_dU(u_m, rhoa_old(i,2), u_p, u_pp, value_drhoa_dy, ONE_OVER_SIX_DELTAY)
     endif
 

     eq2_memory_dpa_dy_adj(i,2)   = b_y(2) * eq2_memory_dpa_dy_adj(i,2) + a_y(2) * value_dpa_dy
     eq2_memory_drhoa_dy_adj(i,2) = b_y(2) * eq2_memory_drhoa_dy_adj(i,2) + a_y(2) * value_drhoa_dy


     ! compute derivative of windx and vax according to x
     u_mm = vax(Im1,2); u_m = vax(i,2); u_p = vax(Ip1,2); u_pp = vax(Ip2,2)
     call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvax_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
       
     u_mm = windx_prior(Im1,2); u_m = windx_prior(i,2); u_p = windx_prior(Ip1,2); u_pp = windx_prior(Ip2,2)
     call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

     eq2_memory_dwindx_dx_adj(i,2) = c_x_half(i_global) * value_dwindx_dx
     eq2_memory_dvax_dx_adj(i,2)   = b_x_half(i_global) * eq2_memory_dvax_dx_adj(i,2)   + a_x_half(i_global) * value_dvax_dx
   
   
     ! compute derivative of windy and vay according to y
     u_mm = vay(i,0); u_m = vay(i,1); u_p = vay(i,2); u_pp = vay(i,3)
     call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvay_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
     u_mm = windy_prior(i,0); u_m = windy_prior(i,1); u_p = windy_prior(i,2); u_pp = windy_prior(i,3)
     call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
   
     eq2_memory_dwindy_dy_adj(i,2) = c_y(2) * value_dwindy_dy
     eq2_memory_dvay_dy_adj(i,2)   = b_y(2) * eq2_memory_dvay_dy_adj(i,2)   + a_y(2) * value_dvay_dy


     ! compute derivative of windy according to x 
     if (windx_half_x <= 0) then !! Solene a change a cause du signe
         u_mm = windy_prior(Im2,1)      ; u_m = windy_prior(Im1,1)      ; u_p = windy_prior(Ip1,1) 
         call compute_decentered_backward_dU(u_mm, u_m, windy_prior(i,1), u_p, value_dwindy_dx_prec, ONE_OVER_SIX_DELTAX)
         u_mm = windy_prior(Im2,2)      ; u_m = windy_prior(Im1,2)      ; u_p = windy_prior(Ip1,2) 
         call compute_decentered_backward_dU(u_mm, u_m, windy_prior(i,2), u_p, value_dwindy_dx_next, ONE_OVER_SIX_DELTAX)
     else 
         u_pp = windy_prior(Ip2,1)      ; u_p = windy_prior(Ip1,1)      ; u_m = windy_prior(Im1,1) 
         call compute_decentered_forward_dU(u_m, windy_prior(i,1), u_p, u_pp, value_dwindy_dx_prec, ONE_OVER_SIX_DELTAX)
         u_pp = windy_prior(Ip2,2)      ; u_p = windy_prior(Ip1,2)      ; u_m = windy_prior(Im1,2) 
         call compute_decentered_forward_dU(u_m, windy_prior(i,2), u_p, u_pp, value_dwindy_dx_next, ONE_OVER_SIX_DELTAX)
     endif
     value_dwindy_dx = 0.5d0 * (value_dwindy_dx_prec + value_dwindy_dx_next)
   
     eq2_memory_dwindy_dx_adj(i,2) = b_x_half(i_global) * eq2_memory_dwindy_dx_adj(i,2) + a_x_half(i_global) * value_dwindy_dx 
   
   
     ! compute derivative of windx according to y
     if (windy_half_y <= 0) then !! Solene a change a cause du signe
         u_mm = windx_prior(i,0)      ; u_m = windx_prior(i,1)      ; u_p = windx_prior(i,3) 
         call compute_decentered_backward_dU(u_mm, u_m, windx_prior(i,2), u_p, value_dwindx_dy_prec, ONE_OVER_SIX_DELTAY)
         u_mm = windx_prior(Ip1,0)      ; u_m = windx_prior(Ip1,1)      ; u_p = windx_prior(Ip1,3) 
         call compute_decentered_backward_dU(u_mm, u_m, windx_prior(Ip1,2), u_p, value_dwindx_dy_next, ONE_OVER_SIX_DELTAY)
     else 
         u_pp = windx_prior(i,4)      ; u_p = windx_prior(i,3)      ; u_m = windx_prior(i,1) 
         call compute_decentered_forward_dU(u_m, windx_prior(i,2), u_p, u_pp, value_dwindx_dy_prec, ONE_OVER_SIX_DELTAY)
         u_pp = windx_prior(Ip1,4)      ; u_p = windx_prior(Ip1,3)      ; u_m = windx_prior(Ip1,1) 
         call compute_decentered_forward_dU(u_m, windx_prior(Ip1,2), u_p, u_pp, value_dwindx_dy_next, ONE_OVER_SIX_DELTAY)
     endif
     value_dwindx_dy = 0.5d0 * (value_dwindx_dy_prec + value_dwindx_dy_next)

     eq2_memory_dwindx_dy_adj(i,2) = b_y(2) * eq2_memory_dwindx_dy_adj(i,2) + a_y(2) * value_dwindx_dy
   
   
     ! Treatment of derivatives for PML   
     value_dpa_dx   = value_dpa_dx   * one_over_K_x_half(i_global) + eq2_memory_dpa_dx_adj(i,2) 
     value_dpa_dy   = value_dpa_dy   * one_over_K_y(2)      + eq2_memory_dpa_dy_adj(i,2)
   
     value_drhoa_dx   = value_drhoa_dx   * one_over_K_x_half(i_global) + eq2_memory_drhoa_dx_adj(i,2) 
     value_drhoa_dy   = value_drhoa_dy   * one_over_K_y(2)      + eq2_memory_drhoa_dy_adj(i,2)
    
     value_dwindx_dx     = value_dwindx_dx * one_over_K_x_half(i_global) + eq2_memory_dwindx_dx_adj(i,2)
     value_dvax_dx       = value_dvax_dx   * one_over_K_x_half(i_global) + eq2_memory_dvax_dx_adj(i,2)
     value_dwindy_dy     = value_dwindy_dy * one_over_K_y(2)      + eq2_memory_dwindy_dy_adj(i,2)
     value_dvay_dy       = value_dvay_dy   * one_over_K_y(2)      + eq2_memory_dvay_dy_adj(i,2)
    
     value_dwindy_dx     = value_dwindy_dx * one_over_K_x_half(i_global) + eq2_memory_dwindy_dx_adj(i,2)
     value_dwindx_dy     = value_dwindx_dy * one_over_K_y(2)      + eq2_memory_dwindx_dy_adj(i,2)
   

     ! intermediate computation
     value_vax_windx_dwindx_dx = vax_half_x * windx_half_x * value_dwindx_dx
     value_vax_windy_dwindx_dy = vax_half_x * windy_half_y * value_dwindx_dy
     value_vay_windx_dwindy_dx = vay_half_y * windx_half_x * value_dwindy_dx
     value_vay_windy_dwindy_dy = vay_half_y * windy_half_y * value_dwindy_dy


     ! compute rhoa
     rhoa(i,2) = rhoa(i,2) + (windx_half_x * value_drhoa_dx + windy_half_y * value_drhoa_dy)          * DELTAT
     rhoa(i,2) = rhoa(i,2) - (value_vax_windx_dwindx_dx + value_vax_windy_dwindx_dy)                  * DELTAT
     rhoa(i,2) = rhoa(i,2) - (value_vay_windx_dwindy_dx + value_vay_windy_dwindy_dy)                  * DELTAT
     rhoa(i,2) = rhoa(i,2) - vay_half_y * g(i,2)                                                      * DELTAT


     ! compute pa
     pa(i,2) = pa(i,2) -  (gamma_chimie(i,2) - 1.0d0) * pa_old(i,2) * (value_dwindx_dx + value_dwindy_dy) * DELTAT
     pa(i,2) = pa(i,2) + (windx_half_x * value_dpa_dx + windy_half_y * value_dpa_dy)                      * DELTAT
     pa(i,2) = pa(i,2) + (value_dvax_dx + value_dvay_dy)                                                  * DELTAT
     pa(i,2) = pa(i,2) + adjoint_source_term(i,2) * DELTAT

   
  enddo

endsubroutine derivative_first_line_parhoa




!---------------
! KERNELS 
!---------------

subroutine derivative_first_raw_kernel_kvx(rho0, p0, windx, windy, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
              value_dvx_dt, value_dvy_dt, it_time, i_start, i_end)

   use parameters, only :  K_windx, gamma_chimie,  &
                          DELTAT, NX, NY,   &
                          NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX,    &
                          NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,    &
                          ONE_OVER_SIX_DELTAX, ONE_OVER_SIX_DELTAY, &
                          ONE_OVER_DELTAY,         &
                          NX_LOCAL, NY_LOCAL, i_global, offset_i
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
      value_dvax_dx,    value_dvax_dy,    &
      value_dvay_dx,                      &
      value_dwindx_dx,                    &
      value_dwindy_dx

  double precision ::   &
      value_drhoarhop_dx, &
      value_drhop_dx,     value_drhop_dy,     &
      value_drho0_dx,     value_drho0_dy,     &
      value_dgammapap_dx,            &
      value_dp_dx,              &
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
    rho0_half_x,                              &
    rhoa_half_x,                              &
    rhop_half_x,                              &
    pa_half_x,                                &
    windy_half_x_half_y, &
    vay_half_x_half_y,   &
    vy_half_x_half_y,    &
    value_dvy_dy_half_x_half_y,               &
    value_dwindy_dy_half_x_half_y
    
   double precision :: & 
       value_dwindy_dx_next, value_dwindy_dx_prec,          &
       value_dwindy_dy_next, value_dwindy_dy_prec,          &
       value_dvy_dy_next,    value_dvy_dy_prec,             &
       value_dvy_dx_next, value_dvy_dx_prec,          &
       value_dvay_dx_next,   value_dvay_dx_prec,            &
       value_drho0_dy_next,  value_drho0_dy_prec,           &
       value_drhop_dy_next,  value_drhop_dy_prec           
       

  double precision  :: u_mm, u_m, u_p, u_pp
  integer :: Ip1,Im1, Ip2, Im2
  
  
  do i =i_start,i_end

    !! Index
    i_global = i + offset_i
    
    Im1 = i-1
    Im2 = i-2
    Ip1 = i+1
    Ip2 = i+2


    ! Interpolation
    rhoa_half_x        = 0.5d0 * (rhoa(i,2) + rhoa(Im1,2))
    rho0_half_x        = 0.5d0 * (rho0(i,2) + rho0(Im1,2))
    rhop_half_x        = 0.5d0 * (rhop(i,2) + rhop(Im1,2))
    pa_half_x          = 0.5d0 * (pa(i,2)   + pa(Im1,2))

    vay_half_x_half_y = 0.25d0 * (vay(i,2) + vay(Im1,2)+ vay(i,1) + vay(Im1,1))
    windy_half_x_half_y  = 0.25d0 * (windy(i,2) + windy(Im1,2)+ windy(i,1) + windy(Im1,1))
    vy_half_x_half_y  = 0.25d0 * (vy(i,2) + vy(Im1,2)+ vy(i,1) + vy(Im1,1))


    ! compute derivative of rhoprhoa, gammapap and pressure according to x
    u_mm = pressure(Im2,2); u_m = pressure(Im1,2); u_p = pressure(i,2); u_pp = pressure(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = rhop(Im2,2); u_m = rhop(Im1,2); u_p = rhop(i,2); u_pp = rhop(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhop_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = rho0(Im2,2); u_m = rho0(Im1,2); u_p = rho0(i,2); u_pp = rho0(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drho0_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = rhop(Im2,2) * rhoa(Im2,2); u_m = rhop(Im1,2) * rhoa(Im1,2)
    u_p  = rhop(i,2) * rhoa(i,2); u_pp = rhop(Ip1,2) * rhoa(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhoarhop_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = pressure(Im2,2) * pa(Im2,2) * gamma_chimie(Im2,2); u_m = pressure(Im1,2) * pa(Im1,2) * gamma_chimie(Im1,2)
    u_p  = pressure(i,2) * pa(i,2) * gamma_chimie(i,2)   ; u_pp = pressure(Ip1,2) * pa(Ip1,2) * gamma_chimie(IP1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dgammapap_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

    
    ! compute derivative of vax and windx according to x
    if (windx(i,2) >= 0) then
        u_mm = windx(Im2,2)      ; u_m = windx(Im1,2)      ; u_p = windx(Ip1,2) 
        call compute_decentered_backward_dU(u_mm, u_m, windx(i,2), u_p, value_dwindx_dx, ONE_OVER_SIX_DELTAX)
        u_mm = vax(Im2,2); u_m = vax(Im1,2); u_p = vax(Ip1,2) 
        call compute_decentered_backward_dU(u_mm, u_m, vax(i,2), u_p, value_dvax_dx, ONE_OVER_SIX_DELTAX)
        u_mm = vx(Im2,2); u_m = vx(Im1,2); u_p = vx(Ip1,2) 
        call compute_decentered_backward_dU(u_mm, u_m, vx(i,2), u_p, value_dvx_dx, ONE_OVER_SIX_DELTAX)
    else 
        u_pp = windx(Ip2,2)    ; u_p = windx(Ip1,2)    ; u_m = windx(Im1,2) 
        call compute_decentered_forward_dU(u_m, windx(i,2), u_p, u_pp, value_dwindx_dx, ONE_OVER_SIX_DELTAX)
        u_pp = vax(Ip2,2)     ; u_p = vax(Ip1,2)     ; u_m = vax(Im1,2) 
        call compute_decentered_forward_dU(u_m, vax(i,2), u_p, u_pp, value_dvax_dx, ONE_OVER_SIX_DELTAX)
        u_pp = vx(Ip2,2)     ; u_p = vx(Ip1,2)     ; u_m = vx(Im1,2) 
        call compute_decentered_forward_dU(u_m, vx(i,2), u_p, u_pp, value_dvx_dx, ONE_OVER_SIX_DELTAX)
    endif
    
      
    ! compute derivative of vay and windy according to x
    u_mm = vay(Im2,1); u_m = vay(Im1,1); u_p = vay(i,1); u_pp = vay(Ip1,1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvay_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = vay(Im2,2); u_m = vay(Im1,2); u_p = vay(i,2); u_pp = vay(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvay_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    value_dvay_dx = 0.5d0 * ( value_dvay_dx_prec + value_dvay_dx_next)

    u_mm = windy(Im2,1); u_m = windy(Im1,1); u_p = windy(i,1); u_pp = windy(Ip1,1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = windy(Im2,2); u_m = windy(Im1,2); u_p = windy(i,2); u_pp = windy(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    value_dwindy_dx = 0.5d0 * ( value_dwindy_dx_prec + value_dwindy_dx_next)

    u_mm = vy(Im2,1); u_m = vy(Im1,1); u_p = vy(i,1); u_pp = vy(Ip1,1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = vy(Im2,2); u_m = vy(Im1,2); u_p = vy(i,2); u_pp = vy(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    value_dvy_dx = 0.5d0 * ( value_dvy_dx_prec + value_dvy_dx_next)
        
          
    ! compute derivative of windy, vy according to y
    u_mm = windy(Im1,0); u_m = windy(Im1,1); u_p = windy(Im1,2); u_pp = windy(Im1,3)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = windy(i,0); u_m = windy(i,1); u_p = windy(i,2); u_pp = windy(i,3)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    value_dwindy_dy_half_x_half_y = 0.5d0 * ( value_dwindy_dy_prec + value_dwindy_dy_next)
  
    u_mm = vy(Im1,0); u_m = vy(Im1,1); u_p = vy(Im1,2); u_pp = vy(Im1,3)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = vy(i,0); u_m = vy(i,1); u_p = vy(i,2); u_pp = vy(i,3)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    value_dvy_dy_half_x_half_y = 0.5d0 * ( value_dvy_dy_prec + value_dvy_dy_next)


    ! compute derivative of vax according to y
    u_pp = vax(i,4)     ; u_p = vax(i,3)  
    call compute_decentered_forward_dU_o2(vax(i,2), u_p, u_pp, value_dvax_dy, ONE_OVER_DELTAY)
  
  
    ! compute derivative of rhop, rho0 according to y
      if (windy(i,2) >= 0) then
        u_mm = rhop(Im1,0); u_m = rhop(Im1,1); u_p = rhop(Im1,3) 
        call compute_decentered_backward_dU(u_mm, u_m, rhop(Im1,2), u_p, value_drhop_dy_prec, ONE_OVER_SIX_DELTAY)
        u_mm = rhop(i,0); u_m = rhop(i,1); u_p = rhop(i,3) 
        call compute_decentered_backward_dU(u_mm, u_m, rhop(i,2), u_p, value_drhop_dy_next, ONE_OVER_SIX_DELTAY)
        
        u_mm = rho0(Im1,0); u_m = rho0(Im1,1); u_p = rho0(Im1,3) 
        call compute_decentered_backward_dU(u_mm, u_m, rho0(Im1,2), u_p, value_drho0_dy_prec, ONE_OVER_SIX_DELTAY)
        u_mm = rho0(i,0); u_m = rho0(i,1); u_p = rho0(i,3) 
        call compute_decentered_backward_dU(u_mm, u_m, rho0(i,2), u_p, value_drho0_dy_next, ONE_OVER_SIX_DELTAY)
    else 
        u_pp = rhop(Im1,4)     ; u_p = rhop(Im1,3)     ; u_m = rhop(Im1,1) 
        call compute_decentered_forward_dU(u_m, rhop(Im1,2), u_p, u_pp, value_drhop_dy_prec, ONE_OVER_SIX_DELTAY)
        u_pp = rhop(i,4)     ; u_p = rhop(i,3)     ; u_m = rhop(i,1) 
        call compute_decentered_forward_dU(u_m, rhop(i,2), u_p, u_pp, value_drhop_dy_next, ONE_OVER_SIX_DELTAY)
        
        u_pp = rho0(Im1,4)     ; u_p = rho0(Im1,3)     ; u_m = rho0(Im1,1) 
        call compute_decentered_forward_dU(u_m, rho0(Im1,2), u_p, u_pp, value_drho0_dy_prec, ONE_OVER_SIX_DELTAY)
        u_pp = rho0(i,4)     ; u_p = rho0(i,3)     ; u_m = rho0(i,1) 
        call compute_decentered_forward_dU(u_m, rho0(i,2), u_p, u_pp, value_drho0_dy_next, ONE_OVER_SIX_DELTAY)
    endif
    value_drhop_dy = 0.5d0 * (value_drhop_dy_prec+value_drhop_dy_next) 
    value_drho0_dy = 0.5d0 * (value_drho0_dy_prec+value_drho0_dy_next) 
 
   
    ! intermediate computations
    value_vx_dvax_dx = vx(i,2) * value_dvax_dx
    value_vy_dvay_dx = vy_half_x_half_y * value_dvay_dx

    value_vax_dvx_dx = vax(i,2) * value_dvx_dx
    value_vay_dvy_dx = vay_half_x_half_y * value_dvy_dx

    value_dwind = value_dwindx_dx + value_dwindy_dy_half_x_half_y
    value_dv = value_dvx_dx + value_dvy_dy_half_x_half_y
  
    value_windx_drhop_dx = windx(i,2) * value_drhop_dx
    value_vx_drho0_dx = vx(i,2) * value_drho0_dx
  
    value_windy_drhop_dy = windy_half_x_half_y * value_drhop_dy
    value_vy_drho0_dy = vy_half_x_half_y * value_drho0_dy
  
    value_wind_dvax = windx(i,2) * value_dvax_dx + windy_half_x_half_y * value_dvax_dy
    value_v_dvax = vx(i,2) * value_dvax_dx + vy_half_x_half_y * value_dvax_dy


    ! Kernel of x-velocity
    K_windx(i,2) = K_windx(i,2) + (rhoa_half_x     * value_drhop_dx)         * DELTAT
    K_windx(i,2) = K_windx(i,2) + rho0_half_x      * (value_vax_dvx_dx + value_vay_dvy_dx)  * DELTAT
    K_windx(i,2) = K_windx(i,2) - vax(i,2) * rho0_half_x * value_dvy_dy_half_x_half_y *DELTAT
    K_windx(i,2) = K_windx(i,2) - vy_half_x_half_y * rho0_half_x * value_dvax_dy *DELTAT
    K_windx(i,2) = K_windx(i,2) + (pa_half_x       * value_dp_dx)  * DELTAT

enddo

endsubroutine derivative_first_raw_kernel_kvx




subroutine derivative_first_raw_kernel_kvy(rho0, p0, windx, windy, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
              value_dvx_dt, value_dvy_dt, it_time, i_start, i_end)

  use parameters, only :  K_windy, gamma_chimie,  &
                          DELTAT, NX, NY,   &
                          NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX,    &
                          NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,    &
                          ONE_OVER_SIX_DELTAX, ONE_OVER_SIX_DELTAY, &
                          ONE_OVER_DELTAY,         &
                          NX_LOCAL, NY_LOCAL, i_global, offset_i
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
      value_dvy_dy,     &
      value_dvax_dy,    &
      value_dvay_dx,    value_dvay_dy,    &
      value_dwindx_dy,  &
      value_dwindy_dy

  double precision ::   &
      value_drhoarhop_dy, &
      value_drhop_dx,     value_drhop_dy,     &
      value_drho0_dx,     value_drho0_dy,     &
      value_dgammapap_dy,      &
      value_dp_dy,        &
      value_vx_dvax_dy,                       &
      value_vy_dvay_dy,                       &
      value_vax_dwindx_dy,                    &
      value_vay_dwindy_dy,                    &
      value_dv,                               &
      value_dwind,                            &
      value_windx_drhop_dx,                   &
      value_windy_drhop_dy,                   &
      value_vx_drho0_dx,                      &
      value_vy_drho0_dy,                      &
      value_wind_dvay,                        &
      value_v_dvay 
      
  double precision :: &
    rho0_half_y,         &
    rhoa_half_y,         &
    rhop_half_y,         &
    pa_half_y,           &
    windx_half_x_half_y, &
    vax_half_x_half_y,   &
    vx_half_x_half_y,    &
    value_dvx_dx_half_x_half_y,               &
    value_dwindx_dx_half_x_half_y

   double precision :: & 
       value_dwindx_dy_next, value_dwindx_dy_prec,          &
       value_dwindx_dx_next, value_dwindx_dx_prec,          &
       value_dvx_dx_next,    value_dvx_dx_prec,             &
       value_dvax_dy_next,   value_dvax_dy_prec,            &
       value_drho0_dx_next,  value_drho0_dx_prec,           &
       value_drhop_dx_next,  value_drhop_dx_prec          

  double precision  :: u_mm, u_m, u_p, u_pp
  integer :: Ip1,Im1, Ip2, Im2
  
 
  do i =i_start,i_end

    ! Index
    i_global = i + offset_i
    Im1 = i-1
    Im2 = i-2
    Ip1 = i+1
    Ip2 = i+2

   
    ! Interpolations
    rhoa_half_y        = 0.5d0 * (rhoa(i,2) + rhoa(i,3))
    rho0_half_y        = 0.5d0 * (rho0(i,2) + rho0(i,3))
    rhop_half_y        = 0.5d0 * (rhop(i,2) + rhop(i,3))
    pa_half_y          = 0.5d0 * (pa(i,2)   + pa(i,3))

    vax_half_x_half_y = 0.25d0 * (vax(i,2) + vax(Ip1,2)+ vax(i,3) + vax(Ip1,3))
    windx_half_x_half_y = 0.25d0 * (windx(i,2) + windx(Ip1,2)+ windx(i,3) + windx(Ip1,3))
    vx_half_x_half_y  = 0.25d0 * (vx(i,2) + vx(Ip1,2)+ vx(i,3) + vx(Ip1,3))


    ! compute derivative of rho0, rhop and pressure according to y
    u_mm = pressure(i,1); u_m = pressure(i,2); u_p = pressure(i,3); u_pp = pressure(i,4)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = rhop(i,1); u_m = rhop(i,2); u_p = rhop(i,3); u_pp = rhop(i,4)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhop_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = rho0(i,1); u_m = rho0(i,2); u_p = rho0(i,3); u_pp = rho0(i,4)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drho0_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)


    ! compute derivative of rhoprhoa, gammapap according to y
    u_m = rhop(i,2) * rhoa(i,2); u_p  = rhop(i,3) * rhoa(i,3)
    value_drhoarhop_dy = (u_p - u_m) * ONE_OVER_DELTAY      
    u_m = pressure(i,2) * pa(i,2) *  gamma_chimie(i,2); u_p  = pressure(i,3) * pa(i,3) * gamma_chimie(i,3)
    value_dgammapap_dy = (u_p - u_m) * ONE_OVER_DELTAY


    !compute derivative of vay, vy and windy according to y
    if (windy(i,2) >= 0) then
        u_mm = windy(i,0)      ; u_m = windy(i,1)      ; u_p = windy(i,3) 
        call compute_decentered_backward_dU(u_mm, u_m, windy(i,2), u_p, value_dwindy_dy, ONE_OVER_SIX_DELTAY)
        u_mm = vy(i,0); u_m = vy(i,1); u_p = vy(i,3) 
        call compute_decentered_backward_dU(u_mm, u_m, vy(i,2), u_p, value_dvy_dy, ONE_OVER_SIX_DELTAY)
        u_mm = vay(i,0); u_m = vay(i,1); u_p = vay(i,3) 
        call compute_decentered_backward_dU(u_mm, u_m, vay(i,2), u_p, value_dvay_dy, ONE_OVER_SIX_DELTAY)         
    else 
        u_pp = windy(i,4)      ; u_p = windy(i,3)      ; u_m = windy(i,1) 
        call compute_decentered_forward_dU(u_m, windy(i,2), u_p, u_pp, value_dwindy_dy, ONE_OVER_SIX_DELTAY)
        u_pp = vy(i,4); u_p = vy(i,3); u_m = vy(i,1) 
        call compute_decentered_forward_dU(u_m, vy(i,2), u_p, u_pp, value_dvy_dy, ONE_OVER_SIX_DELTAY)
        u_pp = vay(i,4); u_p = vay(i,3); u_m = vay(i,1) 
        call compute_decentered_forward_dU(u_m, vay(i,2), u_p, u_pp, value_dvay_dy, ONE_OVER_SIX_DELTAY)      
    endif


    !  compute derivative of windx according to y
    u_mm = windx(i,1); u_m = windx(i,2); u_p = windx(i,3); u_pp = windx(i,4)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = windx(Ip1,1); u_m = windx(Ip1,2); u_p = windx(Ip1,3); u_pp = windx(Ip1,4)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    value_dwindx_dy = 0.5d0 * ( value_dwindx_dy_prec + value_dwindx_dy_next)


    !  compute derivative of vax according to y
    u_m = vax(i,2); u_p = vax(i,3)
    value_dvax_dy_prec = (u_p - u_m) * ONE_OVER_DELTAY
    u_m = vax(Ip1,2); u_p = vax(Ip1,3)
    value_dvax_dy_next = (u_p - u_m) * ONE_OVER_DELTAY
    value_dvax_dy = 0.5d0 * ( value_dvax_dy_prec + value_dvax_dy_next)


    ! compute derivative of windx, vx according to x
    u_mm = windx(Im2,2); u_m = windx(Im1,2); u_p = windx(i,2); u_pp = windx(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = windx(Im2,3); u_m = windx(Im1,3); u_p = windx(i,3); u_pp = windx(Ip1,3)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    value_dwindx_dx_half_x_half_y = 0.5d0 * ( value_dwindx_dx_prec + value_dwindx_dx_next)

    u_mm = vx(Im2,2); u_m = vx(Im1,2); u_p = vx(i,2); u_pp = vx(Ip1,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvx_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = vx(Im2,3); u_m = vx(Im1,3); u_p = vx(i,3); u_pp = vx(Ip1,3)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvx_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    value_dvx_dx_half_x_half_y = 0.5d0 * ( value_dvx_dx_prec + value_dvx_dx_next)


    ! compute derivative of vay according to x
    if (windx(i,2) >= 0) then
        u_mm = vay(Im2,2)      ; u_m = vay(Im1,2)      ; u_p = vay(Ip1,2) 
        call compute_decentered_backward_dU(u_mm, u_m, vay(i,2), u_p, value_dvay_dx, ONE_OVER_SIX_DELTAY)
    else 
        u_pp = vay(Ip2,2)      ; u_p = vay(Ip1,2)      ; u_m = vay(Im1,2) 
        call compute_decentered_forward_dU(u_m, vay(i,2), u_p, u_pp, value_dvay_dx, ONE_OVER_SIX_DELTAY)
    endif

    ! compute derivative of rhop, rho0 according to x
    if (windx(i,2) >= 0) then
        u_mm = rhop(Im2,2)      ; u_m = rhop(Im1,2)      ; u_p = rhop(Ip1,2) 
        call compute_decentered_backward_dU(u_mm, u_m, rhop(i,2), u_p, value_drhop_dx_prec, ONE_OVER_SIX_DELTAX)
        u_mm = rhop(Im2,3)      ; u_m = rhop(Im1,3)      ; u_p = rhop(Ip1,3) 
        call compute_decentered_backward_dU(u_mm, u_m, rhop(i,3), u_p, value_drhop_dx_next, ONE_OVER_SIX_DELTAX)
                
        u_mm = rho0(Im2,2)      ; u_m = rho0(Im1,2)      ; u_p = rho0(Ip1,2) 
        call compute_decentered_backward_dU(u_mm, u_m, rho0(i,2), u_p, value_drho0_dx_prec, ONE_OVER_SIX_DELTAX)
        u_mm = rho0(Im2,3)      ; u_m = rho0(Im1,3)      ; u_p = rho0(Ip1,3) 
        call compute_decentered_backward_dU(u_mm, u_m, rho0(i,3), u_p, value_drho0_dx_next, ONE_OVER_SIX_DELTAX)
    else 
        u_pp = rhop(Ip2,2)      ; u_p = rhop(Ip1,2)      ; u_m = rhop(Im1,2) 
        call compute_decentered_forward_dU(u_m, rhop(i,2), u_p, u_pp, value_drhop_dx_prec, ONE_OVER_SIX_DELTAX)
        u_pp = rhop(Ip2,3)      ; u_p = rhop(Ip1,3)      ; u_m = rhop(Im1,3) 
        call compute_decentered_forward_dU(u_m, rhop(i,3), u_p, u_pp, value_drhop_dx_next, ONE_OVER_SIX_DELTAX)
     
        u_pp = rho0(Ip2,2)      ; u_p = rho0(Ip1,2)      ; u_m = rho0(Im1,2) 
        call compute_decentered_forward_dU(u_m, rho0(i,2), u_p, u_pp, value_drho0_dx_prec, ONE_OVER_SIX_DELTAX)
        u_pp = rho0(Ip2,3)      ; u_p = rho0(Ip1,3)      ; u_m = rho0(Im1,3) 
        call compute_decentered_forward_dU(u_m, rho0(i,3), u_p, u_pp, value_drho0_dx_next, ONE_OVER_SIX_DELTAX)
    endif
    value_drhop_dx = 0.50d0 * (value_drhop_dx_prec+value_drhop_dx_next)
    value_drho0_dx = 0.50d0 * (value_drho0_dx_prec+value_drho0_dx_next)
 
 
    ! intermediate computations
    value_vx_dvax_dy = vx_half_x_half_y * value_dvax_dy
    value_vy_dvay_dy = vy(i,2) * value_dvay_dy

    value_vax_dwindx_dy = vax_half_x_half_y * value_dwindx_dy
    value_vay_dwindy_dy = vay(i,2) * value_dwindy_dy
    
    value_dwind = value_dwindx_dx_half_x_half_y + value_dwindy_dy
    value_dv = value_dvx_dx_half_x_half_y + value_dvy_dy

    value_windx_drhop_dx = windx_half_x_half_y * value_drhop_dx
    value_vx_drho0_dx = vx_half_x_half_y * value_drho0_dx

    value_windy_drhop_dy = windy(i,2) * value_drhop_dy
    value_vy_drho0_dy = vy(i,2) * value_drho0_dy

    value_wind_dvay = windx_half_x_half_y * value_dvay_dx + windy(i,2) * value_dvay_dy
    value_v_dvay = vx_half_x_half_y * value_dvay_dx + vy(i,2) * value_dvay_dy


    ! Kernel of y-velocity
    K_windy(i,2) = K_windy(i,2) + (rhoa_half_y     * value_drhop_dy - value_drhoarhop_dy)         * DELTAT
    K_windy(i,2) = K_windy(i,2) - rho0_half_y      * (value_vx_dvax_dy + value_vy_dvay_dy)        * DELTAT !
    K_windy(i,2) = K_windy(i,2) + rhop_half_y      * (value_vax_dwindx_dy + value_vay_dwindy_dy)  * DELTAT

    K_windy(i,2) = K_windy(i,2) - vay(i,2) * (rhop_half_y * value_dwind + rho0_half_y * value_dv) * DELTAT
    K_windy(i,2) = K_windy(i,2) - vay(i,2) * (value_windx_drhop_dx + value_windy_drhop_dy)        * DELTAT 
    K_windy(i,2) = K_windy(i,2) - vay(i,2) * (value_vx_drho0_dx + value_vy_drho0_dy)              * DELTAT
    K_windy(i,2) = K_windy(i,2) - (rhop_half_y * value_wind_dvay + rho0_half_y * value_v_dvay)    * DELTAT

    K_windy(i,2) = K_windy(i,2) + (pa_half_y       * value_dp_dy - value_dgammapap_dy)  * DELTAT !

  enddo

endsubroutine derivative_first_raw_kernel_kvy




subroutine derivative_first_raw_kernel_kp0krho0(rho0, p0, windx, windy, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
              value_dvx_dt, value_dvy_dt, it_time, i_start, i_end, save_source_term)

 use parameters, only :  K_rho0, K_p0, gamma_chimie,  &
                          DELTAT,  &
                          NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX,    &
                          NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,    &
                          ONE_OVER_SIX_DELTAX, ONE_OVER_SIX_DELTAY, &
                          ONE_OVER_DELTAY,         &
                          DELTAX, DELTAY,                           &
                          ISOURCE, JSOURCE, source_term, &
                          type_source, wavefront,                               &
                          distance2, factor_ssf, SSF_Sigma,                     &
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
      value_dvx_dx,     value_dvx_dy,     &
      value_dvy_dx,     value_dvy_dy,     &
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
      value_drhoa_dx,     value_drhoa_dy,     &
      value_dpa_dx,       value_dpa_dy
      
  double precision :: &
    windx_half_x,        windy_half_y,        &
    vax_half_x,          vay_half_y,          &
    vx_half_x,           vy_half_y,           &
    vax_dvx_dt_half_x, vay_dvy_dt_half_y,     &
    value_dvxvax_dx,   value_dvxvax_dy,       &
    value_dvyvay_dx,   value_dvyvay_dy
    
   double precision :: & 
       value_dwindx_dy_next, value_dwindx_dy_prec,          &
       value_dwindy_dx_next, value_dwindy_dx_prec,          &
       value_dvx_dy_next,    value_dvx_dy_prec,             &
       value_dvy_dx_next,    value_dvy_dx_prec,             &
       value_dvyvay_dx_next, value_dvyvay_dx_prec,          &
       value_dvxvax_dy_next, value_dvxvax_dy_prec

  double precision  :: u_mm, u_m, u_p, u_pp
  integer :: Ip1,Im1, Ip2, Im2
  
 
  do i =i_start,i_end

    ! Index
    i_global = i + offset_i
    
    Im1 = i-1
    Im2 = i-2
    Ip1 = i+1
    Ip2 = i+2


    ! Interpolations
    vx_half_x          = 0.5d0 * (vx(i,2)  + vx(Ip1,2))
    vax_half_x         = 0.5d0 * (vax(i,2) + vax(Ip1,2))
    windx_half_x       = 0.5d0 * (windx(i,2) + windx(Ip1,2))
    vax_dvx_dt_half_x  = 0.5d0 * (vax(i,2)*value_dvx_dt(i,2) + vax(Ip1,2)*value_dvx_dt(Ip1,2))

    vy_half_y          = 0.5d0 * (vy(i,2)  + vy(i,1))
    vay_half_y         = 0.5d0 * (vay(i,2) + vay(i,1))
    windy_half_y       = 0.5d0 * (windy(i,2) + windy(i,1))
    vay_dvy_dt_half_y  = 0.5d0 * (vay(i,2)*value_dvy_dt(i,2) + vay(i,1)*value_dvy_dt(i,1))


    ! Source term estimation
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
    
    
    ! compute derivative of windx, vx according to x 
    u_mm = windx(Im1,2); u_m = windx(i,2); u_p = windx(Ip1,2); u_pp = windx(Ip2,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    
    u_mm = vx(Im1,2); u_m = vx(i,2); u_p = vx(Ip1,2); u_pp = vx(Ip2,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
        
        
    ! compute derivative of windy, vy according to y
    u_mm = vy(i,0); u_m = vy(i,1); u_p = vy(i,2); u_pp = vy(i,3)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      
    u_mm = windy(i,0); u_m = windy(i,1); u_p = windy(i,2); u_pp = windy(i,3)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      
      
    ! compute derivative of pa,rhoa according to x 
    if (windx_half_x >= 0) then
          u_mm = pa(Im2,2)      ; u_m = pa(Im1,2)      ; u_p = pa(Ip1,2) 
          call compute_decentered_backward_dU(u_mm, u_m, pa(i,2), u_p, value_dpa_dx, ONE_OVER_SIX_DELTAX)
          u_mm = rhoa(Im2,2)      ; u_m = rhoa(Im1,2)      ; u_p = rhoa(Ip1,2) 
          call compute_decentered_backward_dU(u_mm, u_m, rhoa(i,2), u_p, value_drhoa_dx, ONE_OVER_SIX_DELTAX)
          
    else 
          u_pp = pa(Ip2,2)      ; u_p = pa(Ip1,2)      ; u_m = pa(Im1,2) 
          call compute_decentered_forward_dU(u_m, pa(i,2), u_p, u_pp, value_dpa_dx, ONE_OVER_SIX_DELTAX)
          u_pp = rhoa(Ip2,2)      ; u_p = rhoa(Ip1,2)      ; u_m = rhoa(Im1,2) 
          call compute_decentered_forward_dU(u_m, rhoa(i,2), u_p, u_pp, value_drhoa_dx, ONE_OVER_SIX_DELTAX)
    endif
    
    
    ! compute derivative of windy, vy according to x 
    if (windx_half_x >= 0) then
          u_mm = windy(Im2,1)      ; u_m = windy(Im1,1)      ; u_p = windy(Ip1,1) 
          call compute_decentered_backward_dU(u_mm, u_m, windy(i,1), u_p, value_dwindy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_mm = windy(Im2,2)      ; u_m = windy(Im1,2)      ; u_p = windy(Ip1,2) 
          call compute_decentered_backward_dU(u_mm, u_m, windy(i,2), u_p, value_dwindy_dx_next, ONE_OVER_SIX_DELTAX)
          
          u_mm = vy(Im2,1)      ; u_m = vy(Im1,1)      ; u_p = vy(Ip1,1) 
          call compute_decentered_backward_dU(u_mm, u_m, vy(i,1), u_p, value_dvy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_mm = vy(Im2,2)      ; u_m = vy(Im1,2)      ; u_p = vy(Ip1,2) 
          call compute_decentered_backward_dU(u_mm, u_m, vy(i,2), u_p, value_dvy_dx_next, ONE_OVER_SIX_DELTAX)
          
    else 
          u_pp = windy(Ip2,1)      ; u_p = windy(Ip1,1)      ; u_m = windy(Im1,1) 
          call compute_decentered_forward_dU(u_m, windy(i,1), u_p, u_pp, value_dwindy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_pp = windy(Ip2,2)      ; u_p = windy(Ip1,2)      ; u_m = windy(Im1,2) 
          call compute_decentered_forward_dU(u_m, windy(i,2), u_p, u_pp, value_dwindy_dx_next, ONE_OVER_SIX_DELTAX)
          
          u_pp = vy(Ip2,1)      ; u_p = vy(Ip1,1)      ; u_m = vy(Im1,1) 
          call compute_decentered_forward_dU(u_m, windy(i,1), u_p, u_pp, value_dvy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_pp = vy(Ip2,2)      ; u_p = vy(Ip1,2)      ; u_m = vy(Im1,2) 
          call compute_decentered_forward_dU(u_m, vy(i,2), u_p, u_pp, value_dvy_dx_next, ONE_OVER_SIX_DELTAX)
    endif
    value_dwindy_dx = 0.5d0 * (value_dwindy_dx_prec + value_dwindy_dx_next)
    value_dvy_dx = 0.5d0 * (value_dvy_dx_prec + value_dvy_dx_next)
    
    
    ! compute derivative of rhoa, pa according to y
    u_pp = pa(i,4)      ; u_p = pa(i,3)   
    call compute_decentered_forward_dU_o2(pa(i,2), u_p, u_pp, value_dpa_dy, ONE_OVER_DELTAY)
          
    u_pp = rhoa(i,4)      ; u_p = rhoa(i,3)              
    call compute_decentered_forward_dU_o2(rhoa(i,2), u_p, u_pp, value_drhoa_dy, ONE_OVER_DELTAY)
    
    
    ! compute derivative of windx, vx according to y
    if (windy_half_y >= 0) then
          u_mm = windx(i,0)      ; u_m = windx(i,1)      ; u_p = windx(i,3) 
          call compute_decentered_backward_dU(u_mm, u_m, windx(i,2), u_p, value_dwindx_dy_prec, ONE_OVER_SIX_DELTAY)
          u_mm = windx(Ip1,0)      ; u_m = windx(Ip1,1)      ; u_p = windx(Ip1,3) 
          call compute_decentered_backward_dU(u_mm, u_m, windx(Ip1,2), u_p, value_dwindx_dy_next, ONE_OVER_SIX_DELTAY)
          
          u_mm = vx(i,0)      ; u_m = vx(i,1)      ; u_p = vx(i,3) 
          call compute_decentered_backward_dU(u_mm, u_m, vx(i,2), u_p, value_dvx_dy_prec, ONE_OVER_SIX_DELTAY)
          u_mm = vx(Ip1,0)      ; u_m = vx(Ip1,1)      ; u_p = vx(Ip1,3) 
          call compute_decentered_backward_dU(u_mm, u_m, vx(Ip1,2), u_p, value_dvx_dy_next, ONE_OVER_SIX_DELTAY)
          
     else 
          u_pp = windx(i,4)      ; u_p = windx(i,3)      ; u_m = windx(i,1) 
          call compute_decentered_forward_dU(u_m, windx(i,2), u_p, u_pp, value_dwindx_dy_prec, ONE_OVER_SIX_DELTAY)
          u_pp = windx(Ip1,4)      ; u_p = windx(Ip1,3)      ; u_m = windx(Ip1,1) 
          call compute_decentered_forward_dU(u_m, windx(Ip1,2), u_p, u_pp, value_dwindx_dy_next, ONE_OVER_SIX_DELTAY)
          
          u_pp = vx(i,4)      ; u_p = vx(i,3)      ; u_m = vx(i,1)                 ! SOLENE has changed Jm into Jp
          call compute_decentered_forward_dU(u_m, vx(i,2), u_p, u_pp, value_dvx_dy_prec, ONE_OVER_SIX_DELTAY)
          u_pp = vx(Ip1,4)      ; u_p = vx(Ip1,3)      ; u_m = vx(Ip1,1)           ! SOLENE has changed Jm into Jp
          call compute_decentered_forward_dU(u_m, vx(Ip1,2), u_p, u_pp, value_dvx_dy_next, ONE_OVER_SIX_DELTAY)

    endif
    value_dwindx_dy = 0.5d0 * (value_dwindx_dy_prec + value_dwindx_dy_next)
    value_dvx_dy = 0.5d0 * (value_dvx_dy_prec + value_dvx_dy_next)
    
    
    ! compute derivatives of vpx * vax according to x
    u_mm = vx(Im1,2)*vax(Im1,2); u_m = vx(i,2)*vax(i,2); u_p = vx(Ip1,2)*vax(Ip1,2); u_pp = vx(Ip2,2)*vax(Ip2,2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvxvax_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    
    
    ! compute derivatives of vpy * vay according to y
    u_mm = vy(i,0)*vay(i,0); u_m = vy(i,1)*vay(i,1)
    u_p = vy(i,2)*vay(i,2); u_pp = vy(i,3)*vay(i,3)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvyvay_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      
      
    ! compute derivatives of vpy * vay according to x  
    if (windx_half_x >= 0) then
          u_mm = vy(Im2,1)*vay(Im2,1)      ; u_m = vy(Im1,1)*vay(Im1,1)
          u_p = vy(Ip1,1)*vay(Ip1,1) 
          call compute_decentered_backward_dU(u_mm, u_m, vy(i,1)*vay(i,1), u_p, value_dvyvay_dx_prec, ONE_OVER_SIX_DELTAX)
          u_mm = vy(Im2,2)*vay(Im2,2)      ; u_m = vy(Im1,2)*vay(Im1,2)
          u_p = vy(Ip1,2)*vay(Ip1,2) 
          call compute_decentered_backward_dU(u_mm, u_m, vy(i,2)*vay(i,2), u_p, value_dvyvay_dx_next, ONE_OVER_SIX_DELTAX) 
    else 
          u_pp = vy(Ip2,1)*vay(Ip2,1)      ; u_p = vy(Ip1,1)*vay(Ip1,1)
          u_m = vy(Im1,1)*vay(Im1,1) 
          call compute_decentered_forward_dU(u_m, vy(i,1)*vay(i,1), u_p, u_pp, value_dvyvay_dx_prec, ONE_OVER_SIX_DELTAX)
          u_pp = vy(Ip2,2)*vay(Ip2,2)      ; u_p = vy(Ip1,2)*vay(Ip1,2)
          u_m = vy(Im1,2)*vay(Im1,2) 
          call compute_decentered_forward_dU(u_m, vy(i,2)*vay(i,2), u_p, u_pp, value_dvyvay_dx_next, ONE_OVER_SIX_DELTAX)
    endif
    value_dvyvay_dx = 0.5d0 * (value_dvyvay_dx_prec + value_dvyvay_dx_next)
    
    
    ! compute derivatives of vpx * vax according to y
    u_pp = vx(i,4)*vax(i,4)      ; u_p = vx(i,3)*vax(i,3)  
    call compute_decentered_forward_dU_o2(vx(i,2)*vax(i,2), u_p, u_pp, value_dvxvax_dy_prec, ONE_OVER_DELTAY)
    u_pp = vx(Ip1,4)*vax(Ip1,4)      ; u_p = vx(Ip1,3) *vax(Ip1,3)
    call compute_decentered_forward_dU_o2(vx(i,2)*vax(Ip1,2), u_p, u_pp, value_dvxvax_dy_next, ONE_OVER_DELTAY)
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
    K_rho0(i,2) = K_rho0(i,2) + (rhoa(i,2) * value_dvx_dx - vy_half_y * value_drhoa_dy) * DELTAT
    K_rho0(i,2) = K_rho0(i,2) - (vax_dvx_dt_half_x + vay_dvy_dt_half_y) * DELTAT
    K_rho0(i,2) = K_rho0(i,2) + (value_vax_windx_dvx_dx + value_vay_windx_dvy_dx) * DELTAT
    K_rho0(i,2) = K_rho0(i,2) + (value_vax_vy_dwindx_dy) * DELTAT
    K_rho0(i,2) = K_rho0(i,2) + pa(i,2) * factor_ssf * source_term * (gamma_chimie(i,2) * p0(i,2) / rho0(i,2)**2) * DELTAT



    ! Kernel of pressure
    K_p0(i,2) = K_p0(i,2) - (vy_half_y * value_dpa_dy)                                         * DELTAT
    K_p0(i,2) = K_p0(i,2) + gamma_chimie(i,2) * pa(i,2) * (value_dvx_dx + value_dvy_dy)        * DELTAT 
    K_p0(i,2) = K_p0(i,2) - pa(i,2) * value_dvy_dy                                             * DELTAT 
    K_p0(i,2) = K_p0(i,2) - pa(i,2) * factor_ssf * source_term * gamma_chimie(i,2) / rho0(i,2) * DELTAT

  enddo

endsubroutine derivative_first_raw_kernel_kp0krho0
