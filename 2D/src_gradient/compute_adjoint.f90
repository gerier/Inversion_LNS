subroutine compute_adjoint(it_step)

  use parameters !, only : pa, rhoa, vax, vay &
                 !      p0_prior, rho0_prior, windx_prior, windy_prior,&
                 !      NX, NY, NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX,       &
                 !      NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,ONE_OVER_SIX_DELTAX,ONE_OVER_SIX_DELTAY, &
                 !      DELTAT, NSTEP,   &
                 !      ZERO, gamma_chimie
  implicit none

  integer it_step

  ! source term
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: adjoint_source_term
  
  ! variable to avoid overwritting data
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: rhoa_old, pa_old, vax_old, vay_old
   
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
    rho0_half_x,         rho0_half_y,         &
    rhoa_half_x,         rhoa_half_y,         &
    pa_half_x,           pa_half_y,           &
    vax_half_x_half_y,   vay_half_x_half_y,   &
    vax_half_x,          vay_half_y,          &
    windx_half_x_half_y, windy_half_x_half_y, &
    windx_half_x,        windy_half_y
  ! intermediate variable to create interpolation on derivative 
  double precision :: &
    value_dwindx_dy_next, value_dwindx_dy_prec,&
    value_dwindy_dx_next, value_dwindy_dx_prec 

  integer :: i,j, j_start
  integer :: Ip1,Im1, Ip2, Im2, Jp1,Jm1, Jp2, Jm2
  double precision :: u_mm, u_m, u_p, u_pp


  ! adjoint source
  call compute_adjoint_source(adjoint_source_term, it_step)
   
  !---
  ! vax and vay computations
  !--- 
  
  ! initialisation  
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
    ! for boundary condition (ground reflection), different numerical scheme 
    j_start = 3
    call derivative_first_line_vaxvay()
  else
    j_start = 1
  endif
    
  do j=j_start,NY_LOCAL
   do i=1,NX_LOCAL

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
      
      
    ! Useful interpolations
    rhoa_half_x = 0.5d0 * (rhoa(i,j) + rhoa(Im1,j))
    rho0_half_x  = 0.5d0 * (rho0_prior(i,j)  + rho0_prior(Im1,j))
    pa_half_x   = 0.5d0 * (pa(i,j)   + pa(Im1,j))
    vay_half_x_half_y = 0.25d0 * (vay_old(i,j) + vay_old(Im1,j) + vay_old(i,Jm1) + vay_old(Im1,Jm1))
    windy_half_x_half_y = 0.25d0 * (windy_prior(i,j)     + windy_prior(Im1,j)     + windy_prior(i,Jm1)     + windy_prior(Im1,Jm1))

    rhoa_half_y = 0.5d0 * (rhoa(i,j) + rhoa(i,Jp1))
    rho0_half_y  = 0.5d0 * (rho0_prior(i,j)  + rho0_prior(i,Jp1))
    pa_half_y   = 0.5d0 * (pa(i,j)   + pa(i,Jp1))
    vax_half_x_half_y = 0.25d0 * (vax_old(i,j) + vax_old(Ip1,j) + vax_old(i,Jp1) + vax_old(Ip1,Jp1))
    windx_half_x_half_y = 0.25d0 * (windx_prior(i,j)     + windx_prior(Ip1,j)     + windx_prior(i,Jp1)     + windx_prior(Ip1,Jp1))


    ! compute derivative of rho0, p0, rho0rhoa, gammap0pa accroding to x
    u_mm = rhoa(Im2,j); u_m = rhoa(Im1,j); u_p = rhoa(i,j); u_pp = rhoa(Ip1,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhoa_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = p0_prior(Im2,j); u_m = p0_prior(Im1,j); u_p = p0_prior(i,j); u_pp = p0_prior(Ip1,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp0_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = p0_prior(Im2,j) * pa(Im2,j) * gamma_chimie(Im2,j); u_m = p0_prior(Im1,j) * pa(Im1,j) * gamma_chimie(Im1,j)
    u_p  = p0_prior(i,j) * pa(i,j) * gamma_chimie(i,j)   ; u_pp = p0_prior(Ip1,j) * pa(Ip1,j) * gamma_chimie(Ip1,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dgammap0pa_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

    eq1_memory_drhoa_dx_adj(i,j)     = b_x(i_global) * eq1_memory_drhoa_dx_adj(i,j)     + a_x(i_global) * value_drhoa_dx
    eq1_memory_dp0_dx_adj(i,j)       = b_x(i_global) * eq1_memory_dp0_dx_adj(i,j)       + a_x(i_global) * value_dp0_dx
    eq1_memory_dgammap0pa_dx_adj(i,j)     = b_x(i_global) * eq1_memory_dgammap0pa_dx_adj(i,j)    &
                                        + a_x(i_global) * value_dgammap0pa_dx 


    ! compute derivative of vax according to y
    if (windy_half_x_half_y >= 0 .and. .not.(.not. USE_PML_YMIN .and. j_rank==0 .and. j == 3)) then
          u_mm = vax_old(i,Jm2); u_m = vax_old(i,Jm1); u_p = vax_old(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, vax_old(i,j), u_p, value_dvax_dy, ONE_OVER_SIX_DELTAY)
    else 
          u_pp = vax_old(i,Jp2); u_p = vax_old(i,Jp1); u_m = vax_old(i,Jm1) 
          call compute_decentered_forward_dU(u_m, vax_old(i,j), u_p, u_pp, value_dvax_dy, ONE_OVER_SIX_DELTAY)
    endif

    eq1_memory_dvax_dy_adj(i,j) = b_y(j_global) * eq1_memory_dvax_dy_adj(i,j) + a_y(j_global) * value_dvax_dy


    ! derivative of windx and vax according to x
    if (windx_prior(i,j) <= 0) then !! Solene a change a cause du signe
          u_mm = windx_prior(Im2,j)      ; u_m = windx_prior(Im1,j)      ; u_p = windx_prior(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, windx_prior(i,j), u_p, value_dwindx_dx, ONE_OVER_SIX_DELTAX)
          u_mm = vax_old(Im2,j); u_m = vax_old(Im1,j); u_p = vax_old(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, vax_old(i,j), u_p, value_dvax_dx, ONE_OVER_SIX_DELTAX)
     else 
          u_pp = windx_prior(Ip2,j)    ; u_p = windx_prior(Ip1,j)    ; u_m = windx_prior(Im1,j) 
          call compute_decentered_forward_dU(u_m, windx_prior(i,j), u_p, u_pp, value_dwindx_dx, ONE_OVER_SIX_DELTAX)
          u_pp = vax_old(Ip2,j)     ; u_p = vax_old(Ip1,j)     ; u_m = vax_old(Im1,j) 
          call compute_decentered_forward_dU(u_m, vax_old(i,j), u_p, u_pp, value_dvax_dx, ONE_OVER_SIX_DELTAX)
     endif

    eq1_memory_dwindx_dx_adj(i,j) = c_x(i_global) * value_dwindx_dx
    eq1_memory_dvax_dx_adj(i,j) = b_x(i_global) * eq1_memory_dvax_dx_adj(i,j) + a_x(i_global) * value_dvax_dx


    ! compute derivative of windy according to x (evaluating in vax)
    u_mm = windy_prior(Im2,Jm1); u_m = windy_prior(Im1,Jm1); u_p = windy_prior(i,Jm1); u_pp = windy_prior(Ip1,Jm1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = windy_prior(Im2,j); u_m = windy_prior(Im1,j); u_p = windy_prior(i,j); u_pp = windy_prior(Ip1,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    value_dwindy_dx = 0.5d0 * ( value_dwindy_dx_prec + value_dwindy_dx_next)

    eq1_memory_dwindy_dx_adj(i,j) = c_x(i_global) * value_dwindy_dx


    ! compute derivative of rho0, p0, rho0rhoa, gammap0pa accroding to y
    u_mm = rhoa(i,Jm1); u_m = rhoa(i,j); u_p = rhoa(i,Jp1); u_pp = rhoa(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhoa_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)      
    u_mm = p0_prior(i,Jm1); u_m = p0_prior(i,j); u_p = p0_prior(i,Jp1); u_pp = p0_prior(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp0_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
   
    u_mm = p0_prior(i,Jm1) * pa(i,Jm1) * gamma_chimie(i,Jm1); u_m = p0_prior(i,j) * pa(i,j) * gamma_chimie(i,j)
    u_p  = p0_prior(i,Jp1) * pa(i,Jp1) * gamma_chimie(i,Jp1); u_pp = p0_prior(i,Jp2) * pa(i,Jp2) * gamma_chimie(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dgammap0pa_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
 
    eq1_memory_drhoa_dy_adj(i,j)     = b_y_half(j_global)*eq1_memory_drhoa_dy_adj(i,j)     + a_y_half(j_global)*value_drhoa_dy
    eq1_memory_dp0_dy_adj(i,j)       = b_y_half(j_global)*eq1_memory_dp0_dy_adj(i,j)       + a_y_half(j_global)*value_dp0_dy
    eq1_memory_dgammap0pa_dy_adj(i,j)     = b_y_half(j_global)*eq1_memory_dgammap0pa_dy_adj(i,j)     &
                                      + a_y_half(j_global)*value_dgammap0pa_dy
 
 
    ! compute derivative of vay according to x
    if (windx_half_x_half_y <= 0) then !! Solene a change a cause du signe
          u_mm = vay_old(Im2,j); u_m = vay_old(Im1,j); u_p = vay_old(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, vay_old(i,j), u_p, value_dvay_dx, ONE_OVER_SIX_DELTAX)
    else 
          u_pp = vay_old(Ip2,j); u_p = vay_old(Ip1,j); u_m = vay_old(Im1,j) 
          call compute_decentered_forward_dU(u_m, vay_old(i,j), u_p, u_pp, value_dvay_dx, ONE_OVER_SIX_DELTAX)
    endif

    eq1_memory_dvay_dx_adj(i,j) = b_x_half(i_global) * eq1_memory_dvay_dx_adj(i,j) + a_x_half(i_global) * value_dvay_dx


    ! compute derivative of windy, vay according to y
    if (windy_prior(i,j) <= 0) then  !! Solene a change a cause du signe
          u_mm = windy_prior(i,Jm2)      ; u_m = windy_prior(i,Jm1)      ; u_p = windy_prior(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, windy_prior(i,j), u_p, value_dwindy_dy, ONE_OVER_SIX_DELTAY)
          u_mm = vay_old(i,Jm2); u_m = vay_old(i,Jm1); u_p = vay_old(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, vay_old(i,j), u_p, value_dvay_dy, ONE_OVER_SIX_DELTAY)
    else 
          u_pp = windy_prior(i,Jp2)      ; u_p = windy_prior(i,Jp1)      ; u_m = windy_prior(i,Jm1) 
          call compute_decentered_forward_dU(u_m, windy_prior(i,j), u_p, u_pp, value_dwindy_dy, ONE_OVER_SIX_DELTAY)
          u_pp = vay_old(i,Jp2); u_p = vay_old(i,Jp1); u_m = vay_old(i,Jm1) 
          call compute_decentered_forward_dU(u_m, vay_old(i,j), u_p, u_pp, value_dvay_dy, ONE_OVER_SIX_DELTAY)
    endif

    eq1_memory_dwindy_dy_adj(i,j) = c_y_half(j_global) * value_dwindy_dy
    eq1_memory_dvay_dy_adj(i,j) = b_y_half(j_global) * eq1_memory_dvay_dy_adj(i,j) + a_y_half(j_global) * value_dvay_dy
        
        
    ! computing derivative of windx according to y, evaluating in vay
    u_mm = windx_prior(i,Jm1); u_m = windx_prior(i,j); u_p = windx_prior(i,Jp1); u_pp = windx_prior(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = windx_prior(Ip1,Jm1); u_m = windx_prior(Ip1,j); u_p = windx_prior(Ip1,Jp1); u_pp = windx_prior(Ip1,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    value_dwindx_dy = 0.5d0 * ( value_dwindx_dy_prec + value_dwindx_dy_next)

    eq1_memory_dwindx_dy_adj(i,j) = c_y_half(j_global) * value_dwindx_dy


    ! PML treatment of derivatives
    value_drhoa_dx     = value_drhoa_dx     * one_over_K_x(i_global) + eq1_memory_drhoa_dx_adj(i,j)  
    value_dp0_dx       =  value_dp0_dx      * one_over_K_x(i_global) + eq1_memory_dp0_dx_adj(i,j)
    value_dgammap0pa_dx= value_dgammap0pa_dx* one_over_K_x(i_global) + eq1_memory_dgammap0pa_dx_adj(i,j) 
    
    value_dvax_dy      = value_dvax_dy      * one_over_K_x(i_global) + eq1_memory_dvax_dy_adj(i,j)

    value_dwindx_dx    = value_dwindx_dx    * one_over_K_x(i_global) + eq1_memory_dwindx_dx_adj(i,j)
    value_dvax_dx      = value_dvax_dx      * one_over_K_x(i_global) + eq1_memory_dvax_dx_adj(i,j) 

    value_dwindy_dx    = value_dwindy_dx    * one_over_K_x(i_global) + eq1_memory_dwindy_dx_adj(i,j)
    
    value_drhoa_dy     = value_drhoa_dy     * one_over_K_y_half(j_global) + eq1_memory_drhoa_dy_adj(i,j)
    value_dp0_dy       = value_dp0_dy       * one_over_K_y_half(j_global) + eq1_memory_dp0_dy_adj(i,j)  
    value_dgammap0pa_dy     = value_dgammap0pa_dy     * one_over_K_y_half(j_global) + eq1_memory_dgammap0pa_dy_adj(i,j)  

    value_dvay_dx      = value_dvay_dx      * one_over_K_y_half(j_global) + eq1_memory_dvay_dx_adj(i,j) 

    value_dwindy_dy    = value_dwindy_dy    * one_over_K_y_half(j_global) + eq1_memory_dwindy_dy_adj(i,j) 
    value_dvay_dy      = value_dvay_dy      * one_over_K_y_half(j_global) + eq1_memory_dvay_dy_adj(i,j)
       
    value_dwindx_dy    = value_dwindx_dy    * one_over_K_y_half(j_global) + eq1_memory_dwindx_dy_adj(i,j)
    
    
    ! intermediate computations    
    value_windxdvax_dx = windx_prior(i,j) * value_dvax_dx
    value_windydvax_dy = windy_half_x_half_y * value_dvax_dy

    value_windxdvay_dx = windx_half_x_half_y * value_dvay_dx
    value_windydvay_dy = windy_prior(i,j) * value_dvay_dy


    ! x component
    vax(i,j) = vax(i,j) +  value_drhoa_dx * DELTAT 
    vax(i,j) = vax(i,j) + (value_windxdvax_dx + value_windydvax_dy) * DELTAT 
    vax(i,j) = vax(i,j) - (vax_old(i,j) * value_dwindx_dx + vay_half_x_half_y * value_dwindy_dx) * DELTAT
    vax(i,j) = vax(i,j) + (value_dgammap0pa_dx - pa_half_x * value_dp0_dx) * DELTAT / rho0_half_x

    ! y component
    vay(i,j) = vay(i,j) +  value_drhoa_dy * DELTAT
    vay(i,j) = vay(i,j) + (value_windxdvay_dx + value_windydvay_dy) * DELTAT
    vay(i,j) = vay(i,j) - (vax_half_x_half_y * value_dwindx_dy + vay_old(i,j) * value_dwindy_dy) * DELTAT
    vay(i,j) = vay(i,j) + (value_dgammap0pa_dy - pa_half_y * value_dp0_dy) * DELTAT / rho0_half_y

  enddo
  enddo


  ! Dircihlet conditions
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

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !---
  ! pa and rhoa computations
  !--- 
  
  rhoa_old(:,:) = rhoa(:,:)
  pa_old(:,:) = pa(:,:)
  
  call send_receive_rightleft(vax)  
  call send_receive_topbottom(vay)
  
  call send_receive_rightleft(pa) 
  call send_receive_rightleft(rhoa)   
  call send_receive_topbottom(pa)
  call send_receive_topbottom(rhoa)
  
  if ((.not. USE_PML_YMIN) .and. j_rank == 0) then
   ! for boundary condition (ground reflection), different numerical scheme
   j_start = 3
   call derivative_first_line_parhoa(adjoint_source_term)
  else
    j_start = 1
  endif
  
  do j=j_start,NY_LOCAL
   do i=1,NX_LOCAL

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
    vax_half_x = 0.5d0 * (vax(i,j) + vax(Ip1,j))
    vay_half_y = 0.5d0 * (vay(i,j) + vay(i,Jm1))
    windx_half_x = 0.5d0 * (windx_prior(i,j) + windx_prior(Ip1,j))
    windy_half_y = 0.5d0 * (windy_prior(i,j) + windy_prior(i,Jm1))

    
    ! compute derivative of rhoa and pa according to x
    if (windx_half_x <= 0) then !! Solene a change a cause du signe
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
    
    
    ! compute derivative of rhoa and p according to y
    if (windy_half_y <= 0 .and. .not.(.not. USE_PML_YMIN .and. j_rank==0 .and. j==3)) then  !! Solene a change a cause du signe
          u_mm = pa_old(i,Jm2)      ; u_m = pa_old(i,Jm1)      ; u_p = pa_old(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, pa_old(i,j), u_p, value_dpa_dy, ONE_OVER_SIX_DELTAY)
          u_mm = rhoa_old(i,Jm2); u_m = rhoa_old(i,Jm1); u_p = rhoa_old(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rhoa_old(i,j), u_p, value_drhoa_dy, ONE_OVER_SIX_DELTAY)
    else 
          u_pp = pa_old(i,Jp2)      ; u_p = pa_old(i,Jp1)      ; u_m = pa_old(i,Jm1) 
          call compute_decentered_forward_dU(u_m, pa_old(i,j), u_p, u_pp, value_dpa_dy, ONE_OVER_SIX_DELTAY)
          u_pp = rhoa_old(i,Jp2); u_p = rhoa_old(i,Jp1); u_m = rhoa_old(i,Jm1) 
          call compute_decentered_forward_dU(u_m, rhoa_old(i,j), u_p, u_pp, value_drhoa_dy, ONE_OVER_SIX_DELTAY)
    endif

    eq2_memory_dpa_dy_adj(i,j)   = b_y(j_global) * eq2_memory_dpa_dy_adj(i,j) + a_y(j_global) * value_dpa_dy
    eq2_memory_drhoa_dy_adj(i,j) = b_y(j_global) * eq2_memory_drhoa_dy_adj(i,j) + a_y(j_global) * value_drhoa_dy


    ! compute derivative of windx and vax according to x
    u_mm = vax(Im1,j); u_m = vax(i,j); u_p = vax(Ip1,j); u_pp = vax(Ip2,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvax_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
        
    u_mm = windx_prior(Im1,j); u_m = windx_prior(i,j); u_p = windx_prior(Ip1,j); u_pp = windx_prior(Ip2,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

    eq2_memory_dwindx_dx_adj(i,j) = c_x_half(i_global) * value_dwindx_dx
                !b_x_half(i_global) * eq2_memory_dwindx_dx_adj(i,j) + a_x_half(i_global) * value_dwindx_dx 
    eq2_memory_dvax_dx_adj(i,j)   = b_x_half(i_global) * eq2_memory_dvax_dx_adj(i,j)   + a_x_half(i_global) * value_dvax_dx
    
    
    ! compute derivative of windy and vay according to y
    u_mm = vay(i,Jm2); u_m = vay(i,Jm1); u_p = vay(i,j); u_pp = vay(i,Jp1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvay_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      
    u_mm = windy_prior(i,Jm2); u_m = windy_prior(i,Jm1); u_p = windy_prior(i,j); u_pp = windy_prior(i,Jp1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)

    eq2_memory_dwindy_dy_adj(i,j) = c_y(j_global) * value_dwindy_dy
!b_y(j_global) * eq2_memory_dwindy_dy_adj(i,j) + a_y(j_global) * value_dwindy_dy
    eq2_memory_dvay_dy_adj(i,j)   = b_y(j_global) * eq2_memory_dvay_dy_adj(i,j)   + a_y(j_global) * value_dvay_dy


    ! compute derivative of windy according to x 
    if (windx_half_x <= 0) then !! Solene a change a cause du signe
          u_mm = windy_prior(Im2,Jm1)      ; u_m = windy_prior(Im1,Jm1)      ; u_p = windy_prior(Ip1,Jm1) 
          call compute_decentered_backward_dU(u_mm, u_m, windy_prior(i,Jm1), u_p, value_dwindy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_mm = windy_prior(Im2,j)      ; u_m = windy_prior(Im1,j)      ; u_p = windy_prior(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, windy_prior(i,j), u_p, value_dwindy_dx_next, ONE_OVER_SIX_DELTAX)
    else 
          u_pp = windy_prior(Ip2,Jm1)      ; u_p = windy_prior(Ip1,Jm1)      ; u_m = windy_prior(Im1,Jm1) 
          call compute_decentered_forward_dU(u_m, windy_prior(i,Jm1), u_p, u_pp, value_dwindy_dx_prec, ONE_OVER_SIX_DELTAX)
          u_pp = windy_prior(Ip2,j)      ; u_p = windy_prior(Ip1,j)      ; u_m = windy_prior(Im1,j) 
          call compute_decentered_forward_dU(u_m, windy_prior(i,j), u_p, u_pp, value_dwindy_dx_next, ONE_OVER_SIX_DELTAX)
    endif
    value_dwindy_dx = 0.5d0 * (value_dwindy_dx_prec + value_dwindy_dx_next)
    
    eq2_memory_dwindy_dx_adj(i,j) = b_x_half(i_global) * eq2_memory_dwindy_dx_adj(i,j) + a_x_half(i_global) * value_dwindy_dx 
    
    
    ! compute derivative of windx according to y
    if (windy_half_y <= 0) then !! Solene a change a cause du signe
          u_mm = windx_prior(i,Jm2)      ; u_m = windx_prior(i,Jm1)      ; u_p = windx_prior(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, windx_prior(i,j), u_p, value_dwindx_dy_prec, ONE_OVER_SIX_DELTAY)
          u_mm = windx_prior(Ip1,Jm2)      ; u_m = windx_prior(Ip1,Jm1)      ; u_p = windx_prior(Ip1,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, windx_prior(Ip1,j), u_p, value_dwindx_dy_next, ONE_OVER_SIX_DELTAY)
     else 
          u_pp = windx_prior(i,Jp2)      ; u_p = windx_prior(i,Jp1)      ; u_m = windx_prior(i,Jm1) 
          call compute_decentered_forward_dU(u_m, windx_prior(i,j), u_p, u_pp, value_dwindx_dy_prec, ONE_OVER_SIX_DELTAY)
          u_pp = windx_prior(Ip1,Jp2)      ; u_p = windx_prior(Ip1,Jp1)      ; u_m = windx_prior(Ip1,Jm1) 
          call compute_decentered_forward_dU(u_m, windx_prior(Ip1,j), u_p, u_pp, value_dwindx_dy_next, ONE_OVER_SIX_DELTAY)

    endif
    value_dwindx_dy = 0.5d0 * (value_dwindx_dy_prec + value_dwindx_dy_next)

    eq2_memory_dwindx_dy_adj(i,j) = b_y(j_global) * eq2_memory_dwindx_dy_adj(i,j) + a_y(j_global) * value_dwindx_dy
    
    
    ! PML treatment of the derivatives   
    value_dpa_dx   = value_dpa_dx   * one_over_K_x_half(i_global) + eq2_memory_dpa_dx_adj(i,j) 
    value_dpa_dy   = value_dpa_dy   * one_over_K_y(j_global)      + eq2_memory_dpa_dy_adj(i,j)
    
    value_drhoa_dx   = value_drhoa_dx   * one_over_K_x_half(i_global) + eq2_memory_drhoa_dx_adj(i,j) 
    value_drhoa_dy   = value_drhoa_dy   * one_over_K_y(j_global)      + eq2_memory_drhoa_dy_adj(i,j)
    
    value_dwindx_dx     = value_dwindx_dx * one_over_K_x_half(i_global) + eq2_memory_dwindx_dx_adj(i,j)
    value_dvax_dx       = value_dvax_dx   * one_over_K_x_half(i_global) + eq2_memory_dvax_dx_adj(i,j)
    value_dwindy_dy     = value_dwindy_dy * one_over_K_y(j_global)      + eq2_memory_dwindy_dy_adj(i,j)
    value_dvay_dy       = value_dvay_dy   * one_over_K_y(j_global)      + eq2_memory_dvay_dy_adj(i,j)
    
    value_dwindy_dx     = value_dwindy_dx * one_over_K_x_half(i_global) + eq2_memory_dwindy_dx_adj(i,j)
    value_dwindx_dy     = value_dwindx_dy * one_over_K_y(j_global)      + eq2_memory_dwindx_dy_adj(i,j)
    

    ! intermediate computation
    value_vax_windx_dwindx_dx = vax_half_x * windx_half_x * value_dwindx_dx
    value_vax_windy_dwindx_dy = vax_half_x * windy_half_y * value_dwindx_dy
    value_vay_windx_dwindy_dx = vay_half_y * windx_half_x * value_dwindy_dx
    value_vay_windy_dwindy_dy = vay_half_y * windy_half_y * value_dwindy_dy


    ! compute rhoa and pa
    rhoa(i,j) = rhoa(i,j) + (windx_half_x * value_drhoa_dx + windy_half_y * value_drhoa_dy)          * DELTAT
    rhoa(i,j) = rhoa(i,j) - (value_vax_windx_dwindx_dx + value_vax_windy_dwindx_dy)                  * DELTAT
    rhoa(i,j) = rhoa(i,j) - (value_vay_windx_dwindy_dx + value_vay_windy_dwindy_dy)                  * DELTAT
    rhoa(i,j) = rhoa(i,j) - vay_half_y * g(i,j)                                                      * DELTAT
 
    pa(i,j) = pa(i,j) -  (gamma_chimie(i,j) - 1.0d0) * pa_old(i,j) * (value_dwindx_dx + value_dwindy_dy) * DELTAT
    pa(i,j) = pa(i,j) + (windx_half_x * value_dpa_dx + windy_half_y * value_dpa_dy)                      * DELTAT
    pa(i,j) = pa(i,j) + (value_dvax_dx + value_dvay_dy)                                                  * DELTAT
    pa(i,j) = pa(i,j) + adjoint_source_term(i,j) * DELTAT
    
   enddo
  enddo


  call send_receive_rightleft(rhoa)
  call send_receive_rightleft(pa)
  call send_receive_rightleft(vax)
  call send_receive_rightleft(vay)
  
  call send_receive_topbottom(rhoa)
  call send_receive_topbottom(pa)
  call send_receive_topbottom(vax)
  call send_receive_topbottom(vay)
 
  
  ! Dircihlet conditions
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



subroutine compute_adjoint_source(adjoint_source_term, it_step)

 use parameters, only : sispressure_true, sispressure_prior, normsq_pressure_true_per_rec,t0,DELTAT,PI, &
                        ix_rec, iy_rec, NREC, NSTEP, TINYVAL,&
                        i_rank, j_rank,NX_LOCAL,NY_LOCAL, offset_i, offset_j
 
 double precision, dimension(-1:NX_LOCAL+2, -1:NY_LOCAL+2) :: adjoint_source_term
  double precision :: diff
 integer :: it_step, irec
 integer :: i,j
 double precision :: coef_damping = 1.0d0
  adjoint_source_term(:,:) = 0.0d0

  ! add a damping in the first iteration to avoid exciting high frequencies
  if (it_step*DELTAT < t0) then
    coef_damping = 0.5 - 0.5 * cos(2*PI*it_step*DELTAT/(2*t0))
  endif 
  
  do irec=1,NREC
   
     if (i_rank == (ix_rec(irec)-1)/NX_LOCAL .and. j_rank == (iy_rec(irec)-1)/NY_LOCAL) then
     
       i = ix_rec(irec) - offset_i 
       j = iy_rec(irec) - offset_j
   
       diff = (sispressure_true(NSTEP-it_step+1,irec) - sispressure_prior(NSTEP-it_step+1,irec))
       
       if (abs(diff)  > TINYVAL) then
         adjoint_source_term(i,j) = coef_damping * diff
       endif
     
       if (normsq_pressure_true_per_rec(irec) > TINYVAL .and. abs(diff) > TINYVAL) then
         adjoint_source_term(i,j) = adjoint_source_term(i,j) / normsq_pressure_true_per_rec(irec)    
       endif
     endif
       
  enddo
 
endsubroutine compute_adjoint_source

