subroutine compute_kernel()
   
  use parameters
  implicit none
                   
  integer :: it, time_last_frame
  integer :: j, ii, jj
  double precision, dimension(0:NX+1,0:NY+1) :: &
      vx_half_t, vy_half_t,                     &
      vax_half_t, vay_half_t,                   &
      value_dvx_dt, value_dvy_dt
 
 double precision , dimension(0:NX+1,0:NY+1) :: &
   rhop_half_t, rhoa_half_t, p_half_t, pa_half_t,&
   rhop_old, rhoa_old, p_old, pa_old

double precision, dimension(0:NX+1, 0:NY+1) :: adjoint_source_term

  use_checkpoint = .False.
  call save_frames()
  calcpressure(:,:) = sispressure(:,:)

  use_checkpoint = .True.


  do it=1,NSTEP
    rhop_old(:,:) = rhop(:,:)
    rhoa_old(:,:) = rhoa(:,:)
    p_old(:,:)  = pressure(:,:)
    pa_old(:,:) = pa(:,:)
    vx_old(:,:) = vx(:,:)
    vy_old(:,:) = vy(:,:)
    
    
    call load_last_frame(NSTEP-it, time_last_frame, pressure, rhop, vx, vy)
    call forwardproblem(pressure, rhop, vx, vy, p0, rho, v0x, v0y, kappa_unrelaxed, time_last_frame, NSTEP-it,3)
    call compute_adjoint(pa, rhoa, vax, vay, it)
    
    if ((mod(it,IT_DISPLAY) == 0 .or. it == 5)) then
        call create_color_image(pa,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
                         0,.FALSE.,.FALSE.,.FALSE.,.FALSE.,10)
    endif
    
    rhop_half_t(:,:)  = 0.5 * (rhop(:,:)     + rhop_old(:,:))
    rhoa_half_t(:,:)  = 0.5 * (rhoa(:,:)     + rhoa_old(:,:))
    p_half_t(:,:)     = 0.5 * (pressure(:,:) + p_old(:,:))
    pa_half_t(:,:)    = 0.5 * (pa(:,:)       + pa_old(:,:))
    
    value_dvy_dt(:,:) = (vy(:,:) - vy_old(:,:)) * ONE_OVER_DELTAT
    value_dvx_dt(:,:) = (vx(:,:) - vx_old(:,:)) * ONE_OVER_DELTAT

    call compute_kernel_iter(rho, p0, v0x, v0y, rhop_half_t, p_half_t, vx, vy,&
              rhoa_half_t, pa_half_t, vax, vay, value_dvx_dt, value_dvy_dt, NSTEP-it)
      
  enddo

  call create_color_image(Krho,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              0,.FALSE.,.FALSE.,.FALSE.,.FALSE.,4)
  call create_color_image(Kp,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              0,.FALSE.,.FALSE.,.FALSE.,.FALSE.,5)
  call create_color_image(Kvy,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              0,.FALSE.,.FALSE.,.FALSE.,.FALSE.,7)
  call create_color_image(Kvx,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              0,.FALSE.,.FALSE.,.FALSE.,.FALSE.,6)
  !OPEN(UNIT=124, FILE="aoutput.txt", ACTION="write", STATUS="replace")
  !DO j=1,NX
  !  WRITE(124,*) Krho(j,100), Kp(j,100), Kvx(j,100), Kvy(j,100) 
  !END DO
  
 OPEN(UNIT=12, FILE="OUTPUT/Kv0x.txt", ACTION="write")
  DO ii=1,NX
    WRITE(12,*) (Kvx(ii,jj), jj=1,NY)
  END DO
  CLOSE(12)
  
  OPEN(UNIT=12, FILE="OUTPUT/Kv0y.txt", ACTION="write")
  DO ii=1,NX
    WRITE(12,*) (Kvy(ii,jj), jj=1,NY)
  END DO
  CLOSE(12)
  
  OPEN(UNIT=12, FILE="OUTPUT/Krho0.txt", ACTION="write")
  DO ii=1,NX
    WRITE(12,*) (Krho(ii,jj), jj=1,NY)
  END DO
  CLOSE(12)
  
  OPEN(UNIT=12, FILE="OUTPUT/Kp0.txt", ACTION="write")
  DO ii=1,NX
    WRITE(12,*) (Kp(ii,jj), jj=1,NY)
  END DO
  CLOSE(12)

endsubroutine compute_kernel


subroutine compute_kernel_iter(rho, p0, v0x, v0y, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
              value_dvx_dt, value_dvy_dt, it_time)


  use parameters, only :  Krho, Kp, Kvx, Kvy, gamma_chimie,         &
                          ONE_OVER_DELTAT, DELTAT, NX, NY,          &
                          NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX,    &
                          NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,    &
                          ONE_OVER_SIX_DELTAX, ONE_OVER_SIX_DELTAY, &
                          DELTAX, DELTAY,                           &
                          a, f0, t0, pi, factor, ISOURCE, JSOURCE, source_term, &
                          distance2, factor_ssf, SSF_Sigma
                          
  implicit none

  integer :: Im1, Ip1, Jp1, Jm1
  integer :: i,j
  
  integer :: it_time
  double precision :: t
  
  double precision, dimension(0:NX+1,0:NY+1) :: &
     rho, p0, v0x, v0y,      &
     rhop, pressure, vx, vy, &
     rhoa, pa, vax, vay,     &
     value_dvx_dt, value_dvy_dt
  
  double precision :: & 
      value_dvx_dx,     &
      value_dvy_dy,     &
      value_dvy_dx,     &
      value_dvx_dy,     &
      value_dvax_dx,    &
      value_dvay_dy,    &
      value_dvay_dx,    &
      value_dvax_dy,    &
      value_dv0x_dx,    &
      value_dv0x_dy,    &
      value_dv0y_dx,    &
      value_dv0y_dy,    &
      value_drhoavx_dx, &
      value_drhoavy_dy, &
      value_vax_v0x_dvx_dx,  &
      value_vax_v0y_dvx_dy,  &
      value_vay_v0x_dvy_dx,  &
      value_vay_v0y_dvy_dy,  &
      value_vax_vx_dv0x_dx,  &
      value_vax_vy_dv0x_dy,  &
      value_vay_vx_dv0y_dx,  &
      value_vay_vy_dv0y_dy,  &
      value_dpavx_dx,        &
      value_dpavy_dy,        &
      value_dsumxvax_dx,     &
      value_dsumxvay_dx,     &
      value_dsumyvax_dy,     &
      value_dsumyvay_dy
      
  double precision ::   &
      value_drho_dx,      &
      value_drho_dy,      &
      value_drhoarhop_dx, &
      value_drhop_dx,     &
      value_vx_dvax_dx,   &
      value_vy_dvay_dx,   &
      value_vax_dv0x_dx,   &
      value_vay_dv0y_dx,   &
      value_dsumva_dx,    &
      value_dpap_dx,      &
      value_dp_dx,        &
      value_drhoarhop_dy, &
      value_drhop_dy,     &
      value_vx_dvax_dy,   &
      value_vy_dvay_dy,   &
      value_vax_dv0x_dy,   &
      value_vay_dv0y_dy,   &
      value_dsumva_dy,    &
      value_dpap_dy,      &
      value_dp_dy,        &
      value_dsumvax,      &
      value_dsumvay
      
  double precision :: &
    rho_half_x,         &
    rho_half_y,         &   
    rhoa_half_x,        &
    rhoa_half_y,        &  
    rhop_half_x,        &
    rhop_half_y,        &     
    pa_half_x,          &
    pa_half_y,          &
    v0x_half_x,         &
    vax_half_x,         &
    vx_half_x,          &
    v0x_half_x_half_y,  &
    v0y_half_x_half_y,  &
    vax_half_x_half_y,  &
    vay_half_x_half_y,  &
    v0y_half_y,         &
    vay_half_y,         &
    vy_half_y,          &
    vx_half_x_half_y,   &
    vy_half_x_half_y,   &
    vax_dvx_dt_half_x,  &
    vay_dvy_dt_half_y
    
   double precision :: value_dp0_dy, value_dp0pa_dy, &
       value_dv0x_dy_next, value_dv0x_dy_prec,       &
       value_dv0y_dx_next, value_dv0y_dx_prec,       &
       value_dvx_dy_next, value_dvx_dy_prec,         &
       value_dvy_dx_next, value_dvy_dx_prec,         &
       value_dvax_dy_next, value_dvax_dy_prec,         &
       value_dvay_dx_next, value_dvay_dx_prec
    
  double precision, dimension(0:NX+1, 0:NY+1) :: pap, rhoprhoa, &
      sumx_vax, sumx_vay, sumy_vax, sumy_vay, &
      pavx, pavy, rhoavx, rhoavy, &
      vax_dvx_dt, vay_dvy_dt
    
    
  !!!!!!!!!!!!!!!!!!
  do j=1,NY
   do i=1,NX
    !! Boundary condition
    Im1 = i-1
    Ip1 = i+1
    Jm1 = j-1
    Jp1 = j+1
    if (i == NX) then
      Ip1 = 2
    elseif (i == 1) then
      Im1 = NX-1
    endif
    if (j == NY) then
      Jp1 = 2
    elseif (j == 1) then
      Jm1 = NY-1
    endif   
 
    rho_half_x  = 0.5d0 * (rho(i,j)  + rho(Im1,j))
    rhop_half_x = 0.5d0 * (rhop(i,j) + rhop(Im1,j))
    vx_half_x   = 0.5d0 * (vx(i,j)   + vx(Ip1,j))
    
    rho_half_y  = 0.5d0 * (rho(i,j)  + rho(i,Jp1)) 
    rhop_half_y = 0.5d0 * (rhop(i,j) + rhop(i,Jp1)) 
    vy_half_y   = 0.5d0 * (vy(i,j)   + vy(i,Jm1))
    
    vx_half_x_half_y  = 0.25d0 * (vx(i,j)  + vx(Ip1,j) + vx(i,Jp1)  + vx(Ip1,Jp1))
    vy_half_x_half_y  = 0.25d0 * (vy(i,j)  + vy(Im1,j) + vy(i,Jm1)  + vy(Im1,Jm1))
    v0x_half_x_half_y = 0.25d0 * (v0x(i,j) + v0x(Ip1,j)+ v0x(i,Jp1) + v0x(Ip1,Jp1))
    v0y_half_x_half_y = 0.25d0 * (v0y(i,j) + v0y(Im1,j)+ v0y(i,Jm1) + v0y(Im1,Jm1))
    
   !rhoavx = 
   !rhoavy =    
   
   rhoprhoa(i,j) = rhop(i,j) * rhoa(i,j)
   pap(i,j)      = pa(i,j)   * pressure(i,j)
   
   sumx_vax(i,j) = (rhop_half_x*v0x(i,j)          + rho_half_x*vx(i,j))          * vax(i,j)
   sumx_vay(i,j) = (rhop_half_y*v0x_half_x_half_y + rho_half_y*vx_half_x_half_y) * vay(i,j)
   sumy_vax(i,j) = (rhop_half_x*v0y_half_x_half_y + rho_half_x*vy_half_x_half_y) * vax(i,j)
   sumy_vay(i,j) = (rhop_half_y*v0y(i,j)          + rho_half_y*vy(i,j))          * vay(i,j) 
   
   pavx(i,j) = pa(i,j) * vx_half_x
   pavy(i,j) = pa(i,j) * vy_half_y
   
   vax_dvx_dt(i,j) = vax(i,j) * value_dvx_dt(i,j)
   vay_dvy_dt(i,j) = vay(i,j) * value_dvy_dt(i,j)
   
  enddo
 enddo 
  !!!!!!!!!!!!!!!!!!  
  
  do j=1,NY
   do i =1,NX

    !! Boundary condition
    Im1 = i-1
    Ip1 = i+1
    Jm1 = j-1
    Jp1 = j+1
    if (i == NX) then
      Ip1 = 2
    elseif (i == 1) then
      Im1 = NX-1
    endif
    if (j == NY) then
      Jp1 = 2
    elseif (j == 1) then
      Jm1 = NY-1
    endif    
    
    rhoa_half_x        = 0.5d0 * (rhoa(i,j) + rhoa(Im1,j))
    pa_half_x          = 0.5d0 * (pa(i,j)   + pa(Im1,j))
    rhoa_half_y        = 0.5d0 * (rhoa(i,j) + rhoa(i,Jp1))
    pa_half_y          = 0.5d0 * (pa(i,j)   + pa(i,Jp1))
    
    vax_half_x_half_y = 0.25d0 * (vax(i,j) + vax(Ip1,j)+ vax(i,Jp1) + vax(Ip1,Jp1))
    vay_half_x_half_y = 0.25d0 * (vay(i,j) + vay(Im1,j)+ vay(i,Jm1) + vay(Im1,Jm1))
    vx_half_x_half_y  = 0.25d0 * (vx(i,j) + vx(Ip1,j)+ vx(i,Jp1) + vx(Ip1,Jp1))
    vy_half_x_half_y  = 0.25d0 * (vy(i,j) + vy(Im1,j)+ vy(i,Jm1) + vy(Im1,Jm1))
    
 
    call compute_centered_dU_dx_in_i(rhoprhoa, Im1, j, value_drhoarhop_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dx_in_i(rhop, Im1, j, value_drho_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dx_in_i(pap, Im1, j, value_dpap_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dx_in_i(pressure, Im1, j, value_dp_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    
    call compute_decentered_dU_dx_in_i(vax, i, j, value_dvax_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x(i,j))
    call compute_centered_dU_dx_in_i(vay, Im1, Jm1, value_dvay_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dx_in_i(vay, Im1, j, value_dvay_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    value_dvay_dx = 0.5d0 * ( value_dvay_dx_prec + value_dvay_dx_next) 
    
    value_vx_dvax_dx = vx(i,j) * value_dvax_dx
    value_vy_dvay_dx = vy_half_x_half_y * value_dvay_dx
    
    call compute_decentered_dU_dx_in_i(v0x, i, j, value_dv0x_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x(i,j))
    call compute_centered_dU_dx_in_i(v0y, Im1, Jm1, value_dv0y_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dx_in_i(v0y, Im1, j, value_dv0y_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    value_dv0y_dx = 0.5d0 * ( value_dv0y_dx_prec + value_dv0y_dx_next) 
    
    value_vax_dv0x_dx = vax(i,j) * value_dv0x_dx
    value_vay_dv0y_dx = vay_half_x_half_y * value_dv0y_dx
    
    call compute_decentered_dU_dx_in_i(sumx_vax, i, j, value_dsumxvax_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x(i,j))
    call compute_decentered_dU_dz_in_i(sumy_vax, i, j, value_dsumyvax_dy, ONE_OVER_SIX_DELTAX,NX,NY,v0y_half_x_half_y)
    value_dsumvax = value_dsumxvax_dx + value_dsumyvax_dy

    call compute_centered_dU_dz_in_i(rhoprhoa, i, j, value_drhoarhop_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    call compute_centered_dU_dz_in_i(rhop, i, j, value_drhop_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    call compute_centered_dU_dz_in_i(pap, i, j, value_dpap_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    call compute_centered_dU_dz_in_i(pressure, i, j, value_dp_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    
    call compute_decentered_dU_dz_in_i(vay, i, j, value_dvay_dy, ONE_OVER_SIX_DELTAY,NX,NY,v0y(i,j))
    call compute_centered_dU_dz_in_i(vax, i, j, value_dvax_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    call compute_centered_dU_dz_in_i(vax, Ip1, j, value_dvax_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    value_dvax_dy = 0.5d0 * ( value_dvax_dy_prec + value_dvax_dy_next) 
    
    value_vx_dvax_dy = vx_half_x_half_y * value_dvax_dy
    value_vy_dvay_dy = vy(i,j) * value_dvay_dy
    
    call compute_decentered_dU_dz_in_i(v0y, i, j, value_dv0y_dy, ONE_OVER_SIX_DELTAY,NX,NY,v0y(i,j))
    call compute_centered_dU_dz_in_i(v0x, i, j, value_dv0x_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    call compute_centered_dU_dz_in_i(v0x, Ip1, j, value_dv0x_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    value_dv0x_dy = 0.5d0 * ( value_dv0y_dx_prec + value_dv0y_dx_next) 
    
    value_vax_dv0x_dy = vax_half_x_half_y * value_dv0x_dy
    value_vay_dv0y_dy = vay(i,j) * value_dv0y_dy
    
    call compute_decentered_dU_dx_in_i(sumx_vay, i, j, value_dsumxvay_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x(i,j))
    call compute_decentered_dU_dz_in_i(sumy_vay, i, j, value_dsumyvay_dy, ONE_OVER_SIX_DELTAX,NX,NY,v0y(i,j))
    value_dsumvay = value_dsumxvay_dx + value_dsumyvay_dy


    ! Kernel of x-velocity
    !Kvx(i,j) = Kvx(i,j) + (rhoa_half_x     * value_drhop_dx - value_drhoarhop_dx)         * DELTAT 
    !Kvx(i,j) = Kvx(i,j) - rho_half_x       * (value_vx_dvax_dx + value_vy_dvay_dx)        * DELTAT !
    !Kvx(i,j) = Kvx(i,j) + rhop_half_x      * (value_vax_dv0x_dx + value_vay_dv0y_dx)      * DELTAT 
    !Kvx(i,j) = Kvx(i,j) - value_dsumvax                                                   * DELTAT !
    !Kvx(i,j) = Kvx(i,j) + (pa_half_x       * value_dp_dx - gamma_chimie * value_dpap_dx)  * DELTAT  
   
    ! Kernel of y-velocity
    !Kvy(i,j) = Kvy(i,j) + (rhoa_half_y     * value_drhop_dy - value_drhoarhop_dy)         * DELTAT 
    !Kvy(i,j) = Kvy(i,j) - rho_half_y       * (value_vx_dvax_dy + value_vy_dvay_dy)        * DELTAT !
    !Kvy(i,j) = Kvy(i,j) + rhop_half_y      * (value_vax_dv0x_dy + value_vay_dv0y_dy)      * DELTAT 
    !Kvy(i,j) = Kvy(i,j) - value_dsumvay                                                   * DELTAT !
    !Kvy(i,j) = Kvy(i,j) + (pa_half_y       * value_dp_dy - gamma_chimie * value_dpap_dy)  * DELTAT ! 
  

    vx_half_x          = 0.5d0 * (vx(i,j)  + vx(Ip1,j))
    vax_half_x         = 0.5d0 * (vax(i,j) + vax(Ip1,j))
    v0x_half_x         = 0.5d0 * (v0x(i,j) + v0x(Ip1,j))
    vax_dvx_dt_half_x  = 0.5d0 * (vax_dvx_dt(i,j) + vax_dvx_dt(Ip1,j))

    vy_half_y          = 0.5d0 * (vy(i,j)  + vy(i,Jm1))  
    vay_half_y         = 0.5d0 * (vay(i,j) + vay(i,Jm1))
    v0y_half_y         = 0.5d0 * (v0y(i,j) + v0y(i,Jm1))
    vay_dvy_dt_half_y  = 0.5d0 * (vay_dvy_dt(i,j) + vay_dvy_dt(i,Jm1))

    call compute_centered_dU_dx_in_i(v0x, i, j, value_dv0x_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dz_in_i(v0y, i, Jm1, value_dv0y_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    
    call compute_centered_dU_dx_in_i(vx, i, j, value_dvx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dz_in_i(vy,i,Jm1, value_dvy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    
    call compute_decentered_dU_dx_in_i(rhoavx, i, j, value_drhoavx_dx,ONE_OVER_SIX_DELTAX,NX,NY,v0x(i,j))
    call compute_decentered_dU_dz_in_i(rhoavy, i, j, value_drhoavy_dy,ONE_OVER_SIX_DELTAY,NX,NY,v0y(i,j))

    call compute_decentered_dU_dx_in_i(vy, i, Jm1, value_dvy_dx_prec,ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
    call compute_decentered_dU_dz_in_i(vx, i, j, value_dvx_dy_prec,ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y)
    call compute_decentered_dU_dx_in_i(vy, i, j, value_dvy_dx_next,ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
    call compute_decentered_dU_dz_in_i(vx, Ip1, j, value_dvx_dy_next,ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y)
    value_dvy_dx = 0.5d0 * (value_dvy_dx_prec + value_dvy_dx_next)  
    value_dvx_dy = 0.5d0 * (value_dvx_dy_prec + value_dvx_dy_next)
    
    value_vax_v0x_dvx_dx = vax_half_x * v0x_half_x * value_dvx_dx
    value_vax_v0y_dvx_dy = vax_half_x * v0y_half_y * value_dvx_dy
    value_vay_v0x_dvy_dx = vay_half_y * v0x_half_x * value_dvy_dx 
    value_vay_v0y_dvy_dy = vay_half_y * v0y_half_y * value_dvy_dy
    
    call compute_decentered_dU_dx_in_i(v0y, i, Jm1, value_dv0y_dx_prec,ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
    call compute_decentered_dU_dz_in_i(v0x, i, j, value_dv0x_dy_prec,ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y)
    call compute_decentered_dU_dx_in_i(v0y, i, j, value_dv0y_dx_next,ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
    call compute_decentered_dU_dz_in_i(v0x, Ip1, j, value_dv0x_dy_next,ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y)
    value_dv0y_dx = 0.5d0 * (value_dv0y_dx_prec + value_dv0y_dx_next)  
    value_dv0x_dy = 0.5d0 * (value_dv0x_dy_prec + value_dv0x_dy_next)
    
    value_vax_vx_dv0x_dx = vax_half_x * vx_half_x * value_dv0x_dx
    value_vax_vy_dv0x_dy = vax_half_x * vy_half_y * value_dv0x_dy
    value_vay_vx_dv0y_dx = vay_half_y * vx_half_x * value_dv0y_dx 
    value_vay_vy_dv0y_dy = vay_half_y * vy_half_y * value_dv0y_dy
        
    call compute_decentered_dU_dx_in_i(pavx, i, j, value_dpavx_dx,ONE_OVER_SIX_DELTAX,NX,NY,v0x(i,j))
    call compute_decentered_dU_dz_in_i(pavy, i, j, value_dpavy_dy,ONE_OVER_SIX_DELTAY,NX,NY,v0y(i,j))

    
    
    a = pi*pi*f0*f0
    t = dble(it_time-1)*DELTAT
    source_term = -8 * a* (t-t0) *  factor * exp(- 4 * a*(t-t0)*(t-t0))
    distance2 = ((i - Isource) * DELTAX)**2 + ((j - Jsource) * DELTAY)**2
    factor_ssf = exp( - distance2 / SSF_Sigma**2 )
      
    ! Kernel of density
    Krho(i,j) = Krho(i,j) + rhoa(i,j) * (value_dvx_dx + value_dvy_dy) * DELTAT
    Krho(i,j) = Krho(i,j) - (value_drhoavx_dx + value_drhoavy_dy) * DELTAT
    Krho(i,j) = Krho(i,j) - (vax_dvx_dt_half_x + vay_dvy_dt_half_y) * DELTAT
    !Krho(i,j) = Krho(i,j) + (value_vax_v0x_dvx_dx + value_vax_v0y_dvx_dy) * DELTAT
    !Krho(i,j) = Krho(i,j) + (value_vay_v0x_dvy_dx + value_vay_v0y_dvy_dy) * DELTAT
    !Krho(i,j) = Krho(i,j) + (value_vax_vx_dv0x_dx + value_vax_vy_dv0x_dy) * DELTAT
    !Krho(i,j) = Krho(i,j) + (value_vay_vx_dv0y_dx + value_vay_vy_dv0y_dy) * DELTAT
    Krho(i,j) = Krho(i,j) + factor_ssf * source_term / (gamma_chimie * p0(i,j)) * DELTAT * rhoa(i,j)
    
    ! Kernel of pressure
    Kp(i,j) = Kp(i,j) - (value_dpavx_dx + value_dpavy_dy) * DELTAT 
    Kp(i,j) = Kp(i,j) + gamma_chimie * pa(i,j) * (value_dvx_dx + value_dvy_dy) * DELTAT ! 
    Kp(i,j) = Kp(i,j) - factor_ssf * source_term * rho(i,j) / (gamma_chimie * p0(i,j)**2) * rhoa(i,j) * DELTAT
    
   enddo
  enddo
  

endsubroutine compute_kernel_iter
