subroutine compute_kernel()

  use parameters
  implicit none

  integer :: it, time_last_frame
  integer :: ii, jj, i_start, i_end, j_start, j_end,i,j
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: &
      !vx_half_t, vy_half_t,                     &
      !vax_half_t, vay_half_t,                   &
      value_dvx_dt, value_dvy_dt

 double precision , dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) ::  &
   rhop_half_t, rhoa_half_t, p_half_t, pa_half_t,&
   rhop_old, rhoa_old, p_old, pa_old,            &
   vx_old, vy_old

 character(len=100) :: file_name

  ! INITIALISATION
  call reset_kernel()
 
  ! To have an exact backward wavefiled, we use checKpointing
  ! We make a first foward simulation to save all the frames
  call save_frames()
  sispressure_prior(:,:) = sispressure(:,:)
    write(file_name, "('./OUTPUT/p_test_',i6.6,'.txt')") rank
    OPEN(UNIT=12, FILE=file_name, ACTION="write")
    DO ii=1,NX_LOCAL
      WRITE(12,*) (pressure(ii,jj), jj=1,NY_LOCAL)
    END DO
    CLOSE(12)
    
  !OPEN(UNIT=12, FILE="OUTPUT/vx.txt", ACTION="write")
  !DO ii=1,NX_LOCAL
  !  WRITE(12,*) (vx(ii,jj), jj=1,NY)
  !END DO
  !CLOSE(12)
  
  ! From this point of the code, we do not need anymore to save signals recorded at receivers
  save_sismos = .False.


  ! Start the kernel computation 
  ! Kernel is the sum over time of correlation between adjoint and backward field
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
    ! select the frame that is the closest of the desired time and who is happening before the desired time
    call load_frame(NSTEP-it, time_last_frame)
    ! Launch a forward simulation from the previously selected frame
    call forwardproblem(p0_prior,rho0_prior,windx_prior,windy_prior,kappa_unrelaxed_prior,time_last_frame, NSTEP-it,3)
    ! compute adjoint field
    call compute_adjoint(it)

    ! draw image of the adjoint field
    !if ((mod(it,IT_DISPLAY) == 0 .or. it == 5)) then
    !    call create_color_image(pa,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
    !                     NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,3333)
    !endif
    if ((mod(it,IT_DISPLAY) == 0 .or. it == 5)) then
      call gather_and_generate_image(vax,vay,pa,rhoa,it,4)
    endif

    ! make interpolation (on time), to hhave all field at the exact same time
    rhop_half_t(:,:)  = 0.5 * (rhop(:,:)     + rhop_old(:,:))
    rhoa_half_t(:,:)  = 0.5 * (rhoa(:,:)     + rhoa_old(:,:))
    p_half_t(:,:)     = 0.5 * (pressure(:,:) + p_old(:,:))
    pa_half_t(:,:)    = 0.5 * (pa(:,:)       + pa_old(:,:))

    value_dvy_dt(:,:) = (vy(:,:) - vy_old(:,:)) * ONE_OVER_DELTAT
    value_dvx_dt(:,:) = (vx(:,:) - vx_old(:,:)) * ONE_OVER_DELTAT

    ! compute for the time it the correlation between adjoint and backward wavefield
    call compute_kernel_iter(rho0_prior, p0_prior, windx_prior, windy_prior, rhop_half_t, p_half_t, vx, vy,&
             rhoa_half_t, pa_half_t, vax, vay, value_dvx_dt, value_dvy_dt, it)

  enddo

 i_start = 1
 i_end = NX_LOCAL 
 j_start = 1
 j_end = NY_LOCAL
 
 if (USE_PML_XMIN) then
    if (i_rank == 0) then
       i_start = NPOINTS_PML + 1
    endif
    if (i_rank == NPROC_X-1) then
       i_end = NX_LOCAL - NPOINTS_PML
     endif 
   
    do jj=j_start,j_end
     K_p0(1:i_start-1,jj) = K_p0(i_start,jj)
     K_rho0(1:i_start-1,jj) = K_rho0(i_start,jj)
     K_windx(1:i_start-1,jj) = K_windx(i_start,jj)
     K_windy(1:i_start-1,jj) = K_windy(i_start,jj)
    
     K_p0(i_end+1:NX_LOCAL,jj) = K_p0(i_end,jj)
     K_rho0(i_end+1:NX_LOCAL,jj) = K_rho0(i_end,jj)
     K_windx(i_end+1:NX_LOCAL,jj) = K_windx(i_end,jj)
     K_windy(i_end+1:NX_LOCAL,jj) = K_windy(i_end,jj)
    enddo
  
  endif
  
  if (USE_PML_YMIN) then
   if (j_rank == 0) then
    j_start = NPOINTS_PML + 1
   endif
   if (j_rank == NPROC_Y-1) then 
    j_end = NY_LOCAL - NPOINTS_PML
   endif 
   
   do ii=i_start,i_end
    K_p0(ii,1:j_start-1) = K_p0(ii,j_start)
    K_rho0(ii,1:j_start-1) = K_rho0(ii,j_start)
    K_windx(ii,1:j_start-1) = K_windx(ii,j_start)
    K_windy(ii,1:j_start-1) = K_windy(ii,j_start)
    
    K_p0(ii,j_end+1:NY_LOCAL) = K_p0(ii,j_end)
    K_rho0(ii,j_end+1:NY_LOCAL) = K_rho0(ii,j_end)
    K_windx(ii,j_end+1:NY_LOCAL) = K_windx(ii,j_end)
    K_windy(ii,j_end+1:NY_LOCAL) = K_windy(ii,j_end)
   enddo
  endif

  if (USE_PML_XMIN .and. USE_PML_YMIN) then
    K_p0(1:i_start-1,1:j_start-1) = K_p0(i_start,j_start)
    K_rho0(1:i_start-1,1:j_start-1) = K_rho0(i_start,j_start)
    K_windx(1:i_start-1,1:j_start-1) = K_windx(i_start,j_start)
    K_windy(1:i_start-1,1:j_start-1) = K_windy(i_start,j_start) 
   
    K_p0(1:i_start-1,j_end+1:NY_LOCAL) = K_p0(i_start,j_end)
    K_rho0(1:i_start-1,j_end+1:NY_LOCAL) = K_rho0(i_start,j_end)
    K_windx(1:i_start-1,j_end+1:NY_LOCAL) = K_windx(i_start,j_end)
    K_windy(1:i_start-1,j_end+1:NY_LOCAL) = K_windy(i_start,j_end) 
   
    K_p0(i_end+1:NX_LOCAL,1:j_start-1) = K_p0(i_end,j_start)
    K_rho0(i_end+1:NX_LOCAL,1:j_start-1) = K_rho0(i_end,j_start)
    K_windx(i_end+1:NX_LOCAL,1:j_start-1) = K_windx(i_end,j_start)
    K_windy(i_end+1:NX_LOCAL,1:j_start-1) = K_windy(i_end,j_start) 
   
    K_p0(i_end+1:NX_LOCAL,j_end+1:NY_LOCAL) = K_p0(i_end,j_end)
    K_rho0(i_end+1:NX_LOCAL,j_end+1:NY_LOCAL) = K_rho0(i_end,j_end)
    K_windx(i_end+1:NX_LOCAL,j_end+1:NY_LOCAL) = K_windx(i_end,j_end)
    K_windy(i_end+1:NX_LOCAL,j_end+1:NY_LOCAL) = K_windy(i_end,j_end) 
  endif

  ! Save kernel information
  !call create_color_image(K_rho0,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
  !            NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,444)
  !call create_color_image(K_p0,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
  !            NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,111)
  !call create_color_image(K_windy,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
  !            NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,222)
  !call create_color_image(K_windx,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
  !            NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,333)
  call gather_and_generate_image(K_windx,K_windy,K_p0,K_rho0,it,3)

  !OPEN(UNIT=12, FILE="OUTPUT/K_windx.txt", ACTION="write")
  !DO ii=1,NX
  !  WRITE(12,*) (K_windx(ii,jj), jj=1,NY)
  !END DO
  !CLOSE(12)
    write(file_name, "('./OUTPUT/Kwindx_true_',i6.6,'.txt')") rank
    OPEN(UNIT=12, FILE=file_name, ACTION="write")
    DO ii=1,NX_LOCAL
      WRITE(12,*) (K_windx(ii,jj), jj=1,NY_LOCAL)
    END DO
    CLOSE(12)

  !OPEN(UNIT=12, FILE="OUTPUT/K_windy.txt", ACTION="write")
  !DO ii=1,NX
  !  WRITE(12,*) (K_windy(ii,jj), jj=1,NY)
  !END DO
  !CLOSE(12)
    write(file_name, "('./OUTPUT/Kwindy_true_',i6.6,'.txt')") rank
    OPEN(UNIT=12, FILE=file_name, ACTION="write")
    DO ii=1,NX_LOCAL
      WRITE(12,*) (K_windy(ii,jj), jj=1,NY_LOCAL)
    END DO
    CLOSE(12)
    
  !OPEN(UNIT=12, FILE="OUTPUT/K_rho0.txt", ACTION="write")
  !DO ii=1,NX
  !  WRITE(12,*) (K_rho0(ii,jj), jj=1,NY)
  !END DO
  !CLOSE(12)
    write(file_name, "('./OUTPUT/Krho0_true_',i6.6,'.txt')") rank
    OPEN(UNIT=12, FILE=file_name, ACTION="write")
    DO ii=1,NX_LOCAL
      WRITE(12,*) (K_rho0(ii,jj), jj=1,NY_LOCAL)
    END DO
    CLOSE(12)
    
  !OPEN(UNIT=12, FILE="OUTPUT/K_p0.txt", ACTION="write")
  !DO ii=1,NX
  !  WRITE(12,*) (K_p0(ii,jj), jj=1,NY)
  !END DO
  !CLOSE(12)
    write(file_name, "('./OUTPUT/Kp0_true_',i6.6,'.txt')") rank
    OPEN(UNIT=12, FILE=file_name, ACTION="write")
    DO ii=1,NX_LOCAL
      WRITE(12,*) (K_p0(ii,jj), jj=1,NY_LOCAL)
    END DO
    CLOSE(12)

  !OPEN(UNIT=12, FILE="OUTPUT/va.txt", ACTION="write")
  !DO ii=1,NX
  !  WRITE(12,*) (vax(ii,jj), jj=1,NY)
  !END DO
  !CLOSE(12)

  !  OPEN(UNIT=12, FILE="OUTPUT/pa.txt", ACTION="write")
  !DO ii=1,NX
  !  WRITE(12,*) (pa(ii,jj), jj=1,NY)
  !END DO
  !CLOSE(12)
      write(file_name, "('./OUTPUT/pa_true_',i6.6,'.txt')") rank
    OPEN(UNIT=12, FILE=file_name, ACTION="write")
    DO ii=1,NX_LOCAL
      WRITE(12,*) (pa(ii,jj), jj=1,NY_LOCAL)
    END DO
    CLOSE(12)
  
  !    OPEN(UNIT=12, FILE="OUTPUT/rhoa.txt", ACTION="write")
  !DO ii=1,NX
  !  WRITE(12,*) (rhoa(ii,jj), jj=1,NY)
  !END DO
  !CLOSE(12)
  
  
endsubroutine compute_kernel


subroutine compute_kernel_iter(rho0, p0, windx, windy, rhop, pressure, vx, vy, rhoa, pa, vax, vay, &
              value_dvx_dt, value_dvy_dt, it_time)


  use parameters, only :  K_rho0, K_p0, K_windx, K_windy, gamma_chimie,         &
                          ONE_OVER_DELTAT, DELTAT, NX, NY,          &
                          NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX,    &
                          NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,    &
                          ONE_OVER_SIX_DELTAX, ONE_OVER_SIX_DELTAY, &
                          DELTAX, DELTAY,                           &
                          a, f0, t0, pi, factor, ISOURCE, JSOURCE, source_term, &
                          distance2, factor_ssf, SSF_Sigma,                     &
                          NPOINTS_PML, USE_PML_XMIN, USE_PML_YMIN,              &
                          NX_LOCAL, NY_LOCAL, i_global, j_global, NPROC_X, NPROC_Y, I_RANK, J_RANK, offset_i, offset_j

  implicit none


  integer :: i,j

  integer :: it_time
  double precision :: t

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: &
     rho0, p0, windx, windy,                    &
     rhop, pressure, vx, vy,                    &
     rhoa, pa, vax, vay,                        &
     value_dvx_dt, value_dvy_dt

  double precision :: &
      value_dvx_dx,     value_dvx_dy,     &
      value_dvy_dx,     value_dvy_dy,     &
      value_dvax_dx,    value_dvax_dy,    &
      value_dvay_dx,    value_dvay_dy,    &
      value_dwindx_dx,  value_dwindx_dy,  &
      value_dwindy_dx,  value_dwindy_dy,  &
      value_drhoavx_dx, value_drhoavy_dy, &
      value_dpavx_dx,   value_dpavy_dy,   &
      value_vax_windx_dvx_dx,             &
      value_vax_windy_dvx_dy,             &
      value_vay_windx_dvy_dx,             &
      value_vay_windy_dvy_dy,             &
      value_vax_vx_dwindx_dx,             &
      value_vax_vy_dwindx_dy,             &
      value_vay_vx_dwindy_dx,             &
      value_vay_vy_dwindy_dy,             &
      value_dsumxvax_dx,                  &
      value_dsumxvay_dx,                  &
      value_dsumyvax_dy,                  &
      value_dsumyvay_dy

  double precision ::   &
      value_drhoarhop_dx, value_drhoarhop_dy, &
      value_drhop_dx,     value_drhop_dy,     &
      value_dpap_dx,      value_dpap_dy,      &
      value_dp_dx,        value_dp_dy,        &
      value_vx_dvax_dx,                       &
      value_vy_dvay_dx,                       &
      value_vx_dvax_dy,                       &
      value_vy_dvay_dy,                       &
      value_vax_dwindx_dx,                    &
      value_vay_dwindy_dx,                    &
      value_vax_dwindx_dy,                    &
      value_vay_dwindy_dy,                    &
      value_dsumvax, value_dsumvay     
      

  double precision :: &
    rho0_half_x,         rho0_half_y,         &
    rhoa_half_x,         rhoa_half_y,         &
    rhop_half_x,         rhop_half_y,         &
    pa_half_x,           pa_half_y,           &
    windx_half_x,        windy_half_y,        &
    vax_half_x,          vay_half_y,          &
    vx_half_x,           vy_half_y,           &
    windx_half_x_half_y, windy_half_x_half_y, &
    vax_half_x_half_y,   vay_half_x_half_y,   &
    vx_half_x_half_y,    vy_half_x_half_y,    &
    vax_dvx_dt_half_x, vay_dvy_dt_half_y
    

   double precision :: & 
       value_dwindx_dy_next, value_dwindx_dy_prec,          &
       value_dwindy_dx_next, value_dwindy_dx_prec,          &
       value_dvx_dy_next,    value_dvx_dy_prec,             &
       value_dvy_dx_next,    value_dvy_dx_prec,             &
       value_dvax_dy_next,   value_dvax_dy_prec,            &
       value_dvay_dx_next,   value_dvay_dx_prec

  double precision, dimension(-1:NX_LOCAL+2, -1:NY_LOCAL+2) :: &
      pap, rhoprhoa, &
      sumx_vax, sumx_vay, sumy_vax, sumy_vay, &
      pavx, pavy, rhoavx, rhoavy, &
      vax_dvx_dt, vay_dvy_dt


  double precision  :: u_mm, u_m, u_p, u_pp
  integer :: Ip1,Im1, Ip2, Im2, Jp1,Jm1, Jp2, Jm2
  
  
  integer :: i_start, i_end, j_start, j_end
  
  
  j_start = 1
  j_end = NY_LOCAL
  
  if (USE_PML_YMIN) then
  if (j_rank == 0) then
    j_start = NPOINTS_PML + 1
  endif
  if (j_rank == NPROC_Y-1) then
    j_end = NY_LOCAL - NPOINTS_PML
  endif
  endif
  
  i_start = 1
  i_end = NX_LOCAL
  if (USE_PML_XMIN) then
  if (i_rank == 0) then
    i_start = NPOINTS_PML + 1
  endif
  if (i_rank == NPROC_X -1) then
    i_end = NX_LOCAL - NPOINTS_PML
  endif
  endif
  

  !!!!!!!!!!!!!!!!!!
  
  call send_receive_left(vx)
  call send_receive_right(vy)  
  call send_receive_right(rhop)
  
  call send_receive_top(vy)
  call send_receive_bottom(vx)
  call send_receive_bottom(rhop)
  
  call send_receive_righttop(vy)
  call send_receive_leftbottom(vx)
  
  do j=j_start,j_end
   do i=i_start,i_end
   
    !! Boundary condition
       Im1 = i-1
       Im2 = i-2
       Ip1 = i+1
       Ip2 = i+2
       Jm1 = j-1
       Jm2 = j-2
       Jp1 = j+1
       Jp2 = j+2
      
      
    ! modify index dependingon the boundary condition
    ! note, only useful when periodic conditions
    call get_index_boundarycondition(i,j, Im2, Im1,Ip1, Ip2, Jm2, Jm1,Jp1, Jp2)

    rho0_half_x  = 0.5d0 * (rho0(i,j)  + rho0(Im1,j))
    rhop_half_x = 0.5d0 * (rhop(i,j) + rhop(Im1,j))
    vx_half_x   = 0.5d0 * (vx(i,j)   + vx(Ip1,j))

    rho0_half_y  = 0.5d0 * (rho0(i,j)  + rho0(i,Jp1))
    rhop_half_y = 0.5d0 * (rhop(i,j) + rhop(i,Jp1))
    vy_half_y   = 0.5d0 * (vy(i,j)   + vy(i,Jm1))

    vx_half_x_half_y  = 0.25d0 * (vx(i,j)  + vx(Ip1,j) + vx(i,Jp1)  + vx(Ip1,Jp1))
    vy_half_x_half_y  = 0.25d0 * (vy(i,j)  + vy(Im1,j) + vy(i,Jm1)  + vy(Im1,Jm1))
    windx_half_x_half_y = 0.25d0 * (windx(i,j) + windx(Ip1,j)+ windx(i,Jp1) + windx(Ip1,Jp1))
    windy_half_x_half_y = 0.25d0 * (windy(i,j) + windy(Im1,j)+ windy(i,Jm1) + windy(Im1,Jm1))

   sumx_vax(i,j) = (rhop_half_x*windx(i,j)          + rho0_half_x*vx(i,j))          * vax(i,j)
   sumx_vay(i,j) = (rhop_half_y*windx_half_x_half_y + rho0_half_y*vx_half_x_half_y) * vay(i,j)
   sumy_vax(i,j) = (rhop_half_x*windy_half_x_half_y + rho0_half_x*vy_half_x_half_y) * vax(i,j)
   sumy_vay(i,j) = (rhop_half_y*windy(i,j)          + rho0_half_y*vy(i,j))          * vay(i,j)

   pavx(i,j) = pa(i,j) * vx_half_x
   pavy(i,j) = pa(i,j) * vy_half_y

   vax_dvx_dt(i,j) = vax(i,j) * value_dvx_dt(i,j)
   vay_dvy_dt(i,j) = vay(i,j) * value_dvy_dt(i,j)

  enddo
 enddo
  !!!!!!!!!!!!!!!!!!

  call send_receive_rightleft(rhoa)
  call send_receive_rightleft(pa)
  call send_receive_rightleft(vax)
  call send_receive_rightleft(vay)
  call send_receive_rightleft(rhop)
  call send_receive_rightleft(pressure)
  call send_receive_rightleft(vx)
  call send_receive_rightleft(vy)
  call send_receive_rightleft(sumx_vax)
  call send_receive_rightleft(sumx_vay)
  call send_receive_rightleft(rhoavx)
  call send_receive_rightleft(pavx)
  
  call send_receive_left(vax_dvx_dt)
  call send_receive_top(vay_dvy_dt)
  
  call send_receive_topbottom(rhoa)
  call send_receive_topbottom(pa)
  call send_receive_topbottom(vax)
  call send_receive_topbottom(vay)
  call send_receive_topbottom(rhop)
  call send_receive_topbottom(pressure)
  call send_receive_topbottom(vy)
  call send_receive_topbottom(sumy_vax)
  call send_receive_topbottom(sumy_vay)
  call send_receive_topbottom(rhoavy)
  call send_receive_topbottom(pavy)
  
  call send_receive_leftbottom(vax)
  call send_receive_lefttop(vax)
  call send_receive_leftbottom(vx)
  call send_receive_lefttop(vx)
  
  call send_receive_righttop(vay)
  call send_receive_lefttop(vay)
  call send_receive_righttop(vy)
  call send_receive_lefttop(vy)
  

  do j=j_start,j_end
   do i =i_start,i_end

    i_global = i + offset_i
    j_global = j + offset_j
    

    !! Boundary condition
       Im1 = i-1
       Im2 = i-2
       Ip1 = i+1
       Ip2 = i+2
       Jm1 = j-1
       Jm2 = j-2
       Jp1 = j+1
       Jp2 = j+2
     
    call get_index_boundarycondition(i,j, Im2, Im1,Ip1, Ip2, Jm2, Jm1,Jp1, Jp2)

    rhoa_half_x        = 0.5d0 * (rhoa(i,j) + rhoa(Im1,j))
    pa_half_x          = 0.5d0 * (pa(i,j)   + pa(Im1,j))
    rhoa_half_y        = 0.5d0 * (rhoa(i,j) + rhoa(i,Jp1))
    pa_half_y          = 0.5d0 * (pa(i,j)   + pa(i,Jp1))

    vax_half_x_half_y = 0.25d0 * (vax(i,j) + vax(Ip1,j)+ vax(i,Jp1) + vax(Ip1,Jp1))
    vay_half_x_half_y = 0.25d0 * (vay(i,j) + vay(Im1,j)+ vay(i,Jm1) + vay(Im1,Jm1))
    vx_half_x_half_y  = 0.25d0 * (vx(i,j) + vx(Ip1,j)+ vx(i,Jp1) + vx(Ip1,Jp1))
    vy_half_x_half_y  = 0.25d0 * (vy(i,j) + vy(Im1,j)+ vy(i,Jm1) + vy(Im1,Jm1))


    ! compute derivative of rhoprhoa, pap and pressure according to x
    u_mm = pressure(Im2,j); u_m = pressure(Im1,j); u_p = pressure(i,j); u_pp = pressure(Ip1,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = rhop(Im2,j); u_m = rhop(Im1,j); u_p = rhop(i,j); u_pp = rhop(Ip1,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhop_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = rhop(Im2,j) * rhoa(Im2,j); u_m = rhop(Im1,j) * rhoa(Im1,j)
    u_p  = rhop(i,j) * rhoa(i,j); u_pp = rhop(Ip1,j) * rhoa(Ip1,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhoarhop_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = pressure(Im2,j) * pa(Im2,j); u_m = pressure(Im1,j) * pa(Im1,j)
    u_p  = pressure(i,j) * pa(i,j)    ; u_pp = pressure(Ip1,j) * pa(Ip1,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dpap_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)

    ! compute derivative of vax and windx according to x
    if (windx(i,j) >= 0) then
          u_mm = windx(Im2,j)      ; u_m = windx(Im1,j)      ; u_p = windx(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, windx(i,j), u_p, value_dwindx_dx, ONE_OVER_SIX_DELTAX)
          u_mm = vax(Im2,j); u_m = vax(Im1,j); u_p = vax(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, vax(i,j), u_p, value_dvax_dx, ONE_OVER_SIX_DELTAX)
    else 
          u_pp = windx(Ip2,j)    ; u_p = windx(Ip1,j)    ; u_m = windx(Im1,j) 
          call compute_decentered_forward_dU(u_m, windx(i,j), u_p, u_pp, value_dwindx_dx, ONE_OVER_SIX_DELTAX)
          u_pp = vax(Ip2,j)     ; u_p = vax(Ip1,j)     ; u_m = vax(Im1,j) 
          call compute_decentered_forward_dU(u_m, vax(i,j), u_p, u_pp, value_dvax_dx, ONE_OVER_SIX_DELTAX)
    endif
        
    ! compute derivative of vay and windy according to x
    u_mm = vay(Im2,Jm1); u_m = vay(Im1,Jm1); u_p = vay(i,Jm1); u_pp = vay(Ip1,Jm1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvay_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = vay(Im2,j); u_m = vay(Im1,j); u_p = vay(i,j); u_pp = vay(Ip1,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvay_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    value_dvay_dx = 0.5d0 * ( value_dvay_dx_prec + value_dvay_dx_next)

    u_mm = windy(Im2,Jm1); u_m = windy(Im1,Jm1); u_p = windy(i,Jm1); u_pp = windy(Ip1,Jm1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    u_mm = windy(Im2,j); u_m = windy(Im1,j); u_p = windy(i,j); u_pp = windy(Ip1,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dx_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    value_dwindy_dx = 0.5d0 * ( value_dwindy_dx_prec + value_dwindy_dx_next)
        
    ! intermediate computations
    value_vx_dvax_dx = vx(i,j) * value_dvax_dx
    value_vy_dvay_dx = vy_half_x_half_y * value_dvay_dx

    value_vax_dwindx_dx = vax(i,j) * value_dwindx_dx
    value_vay_dwindy_dx = vay_half_x_half_y * value_dwindy_dx

    ! compute derivative of (rhop v0 x va_x + rho0vp x va_x)
    if (windx(i,j) >= 0) then
          u_mm = sumx_vax(Im2,j)      ; u_m = sumx_vax(Im1,j)      ; u_p = sumx_vax(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, sumx_vax(i,j), u_p, value_dsumxvax_dx, ONE_OVER_SIX_DELTAX)
    else 
          u_pp = sumx_vax(Ip2,j)      ; u_p = sumx_vax(Ip1,j)      ; u_m = sumx_vax(Im1,j) 
          call compute_decentered_forward_dU(u_m, sumx_vax(i,j), u_p, u_pp, value_dsumxvax_dx, ONE_OVER_SIX_DELTAX)
    endif
     
    if (windy_half_x_half_y >= 0) then
          u_mm = sumy_vax(i,Jm2)      ; u_m = sumy_vax(i,Jm1)      ; u_p = sumy_vax(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, sumy_vax(i,j), u_p, value_dsumyvax_dy, ONE_OVER_SIX_DELTAX)
    else 
          u_pp = sumy_vax(i,Jp2)      ; u_p = sumy_vax(i,Jp1)      ; u_m = sumy_vax(i,Jm1) 
          call compute_decentered_forward_dU(u_m, sumy_vax(i,j), u_p, u_pp, value_dsumyvax_dy, ONE_OVER_SIX_DELTAX)
    endif
    value_dsumvax = value_dsumxvax_dx + value_dsumyvax_dy
          
    ! compute derivative of rhoprhoa, pap and pressure according to y
    u_mm = pressure(i,Jm1); u_m = pressure(i,j); u_p = pressure(i,Jp1); u_pp = pressure(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dp_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
     u_mm = rhop(i,Jm1); u_m = rhop(i,j); u_p = rhop(i,Jp1); u_pp = rhop(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhop_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = rhop(i,Jm1) * rhoa(i,Jm1); u_m = rhop(i,j) * rhoa(i,j)
    u_p  = rhop(i,Jp1) * rhoa(i,Jp1); u_pp = rhop(i,Jp2) * rhoa(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drhoarhop_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)      
    u_mm = pressure(i,Jm1) * pa(i,Jm1); u_m = pressure(i,j) * pa(i,j)
    u_p  = pressure(i,Jp1) * pa(i,Jp1); u_pp = pressure(i,Jp2) * pa(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dpap_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    

    !compute derivative of vay and windy according to y
    if (windy(i,j) >= 0) then
          u_mm = windy(i,Jm2)      ; u_m = windy(i,Jm1)      ; u_p = windy(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, windy(i,j), u_p, value_dwindy_dy, ONE_OVER_SIX_DELTAY)
          u_mm = vay(i,Jm2); u_m = vay(i,Jm1); u_p = vay(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, vay(i,j), u_p, value_dvay_dy, ONE_OVER_SIX_DELTAY)
    else 
          u_pp = windy(i,Jp2)      ; u_p = windy(i,Jp1)      ; u_m = windy(i,Jm1) 
          call compute_decentered_forward_dU(u_m, windy(i,j), u_p, u_pp, value_dwindy_dy, ONE_OVER_SIX_DELTAY)
          u_pp = vay(i,Jp2); u_p = vay(i,Jp1); u_m = vay(i,Jm1) 
          call compute_decentered_forward_dU(u_m, vay(i,j), u_p, u_pp, value_dvay_dy, ONE_OVER_SIX_DELTAY)
    endif
    
    !  compute derivative of vax and windx according to y
    u_mm = windx(i,Jm1); u_m = windx(i,j); u_p = windx(i,Jp1); u_pp = windx(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = windx(Ip1,Jm1); u_m = windx(Ip1,j); u_p = windx(Ip1,Jp1); u_pp = windx(Ip1,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    value_dwindx_dy = 0.5d0 * ( value_dwindx_dy_prec + value_dwindx_dy_next)
    
    u_mm = vax(i,Jm1); u_m = vax(i,j); u_p = vax(i,Jp1); u_pp = vax(i,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvax_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    u_mm = vax(Ip1,Jm1); u_m = vax(Ip1,j); u_p = vax(Ip1,Jp1); u_pp = vax(Ip1,Jp2)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvax_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
    value_dvax_dy = 0.5d0 * ( value_dvax_dy_prec + value_dvax_dy_next)
    
    ! intermediate computations
    value_vx_dvax_dy = vx_half_x_half_y * value_dvax_dy
    value_vy_dvay_dy = vy(i,j) * value_dvay_dy

    value_vax_dwindx_dy = vax_half_x_half_y * value_dwindx_dy
    value_vay_dwindy_dy = vay(i,j) * value_dwindy_dy

    ! compute derivative of (rhop v0 x va_y + rho0vp x va_y)
    if (windx_half_x_half_y >= 0) then
          u_mm = sumx_vay(Im2,j)      ; u_m = sumx_vay(Im1,j)      ; u_p = sumx_vay(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, sumx_vay(i,j), u_p, value_dsumxvay_dx, ONE_OVER_SIX_DELTAX)
    else 
          u_pp = sumx_vay(Ip2,j)      ; u_p = sumx_vay(Ip1,j)      ; u_m = sumx_vay(Im1,j) 
          call compute_decentered_forward_dU(u_m, sumx_vay(i,j), u_p, u_pp, value_dsumxvay_dx, ONE_OVER_SIX_DELTAX)
    endif
   
    if (windy(i,j)>= 0) then
          u_mm = sumy_vay(i,Jm2)      ; u_m = sumy_vay(i,Jm1)      ; u_p = sumy_vay(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, sumy_vay(i,j), u_p, value_dsumyvay_dy, ONE_OVER_SIX_DELTAX)
    else 
          u_pp = sumy_vay(i,Jp2)      ; u_p = sumy_vay(i,Jp1)      ; u_m = sumy_vay(i,Jm1) 
          call compute_decentered_forward_dU(u_m, sumy_vay(i,j), u_p, u_pp, value_dsumyvay_dy, ONE_OVER_SIX_DELTAX)
    endif
    value_dsumvay = value_dsumxvay_dx + value_dsumyvay_dy


    ! Kernel of x-velocity
    K_windx(i,j) = K_windx(i,j) + (rhoa_half_x     * value_drhop_dx - value_drhoarhop_dx)         * DELTAT
    K_windx(i,j) = K_windx(i,j) - rho0_half_x      * (value_vx_dvax_dx + value_vy_dvay_dx)        * DELTAT !
    K_windx(i,j) = K_windx(i,j) + rhop_half_x      * (value_vax_dwindx_dx + value_vay_dwindy_dx)  * DELTAT
    K_windx(i,j) = K_windx(i,j) - value_dsumvax                                                   * DELTAT !
    K_windx(i,j) = K_windx(i,j) + (pa_half_x       * value_dp_dx - gamma_chimie * value_dpap_dx)  * DELTAT

    ! Kernel of y-velocity
    K_windy(i,j) = K_windy(i,j) + (rhoa_half_y     * value_drhop_dy - value_drhoarhop_dy)         * DELTAT
    K_windy(i,j) = K_windy(i,j) - rho0_half_y      * (value_vx_dvax_dy + value_vy_dvay_dy)        * DELTAT !
    K_windy(i,j) = K_windy(i,j) + rhop_half_y      * (value_vax_dwindx_dy + value_vay_dwindy_dy)  * DELTAT
    K_windy(i,j) = K_windy(i,j) - value_dsumvay                                                   * DELTAT !
    K_windy(i,j) = K_windy(i,j) + (pa_half_y       * value_dp_dy - gamma_chimie * value_dpap_dy)  * DELTAT !


    vx_half_x          = 0.5d0 * (vx(i,j)  + vx(Ip1,j))
    vax_half_x         = 0.5d0 * (vax(i,j) + vax(Ip1,j))
    windx_half_x       = 0.5d0 * (windx(i,j) + windx(Ip1,j))
    vax_dvx_dt_half_x  = 0.5d0 * (vax_dvx_dt(i,j) + vax_dvx_dt(Ip1,j))

    vy_half_y          = 0.5d0 * (vy(i,j)  + vy(i,Jm1))
    vay_half_y         = 0.5d0 * (vay(i,j) + vay(i,Jm1))
    windy_half_y       = 0.5d0 * (windy(i,j) + windy(i,Jm1))
    vay_dvy_dt_half_y  = 0.5d0 * (vay_dvy_dt(i,j) + vay_dvy_dt(i,Jm1))

    ! compute derivative of windx, vx according to x 
    u_mm = windx(Im1,j); u_m = windx(i,j); u_p = windx(Ip1,j); u_pp = windx(Ip2,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
    
    u_mm = vx(Im1,j); u_m = vx(i,j); u_p = vx(Ip1,j); u_pp = vx(Ip2,j)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
        
    ! compute derivative of rhoavx, pavx according to x
    if (windx_half_x >= 0) then
          u_mm = rhoavx(Im2,j); u_m = rhoavx(Im1,j); u_p = rhoavx(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, rhoavx(i,j), u_p, value_drhoavx_dx, ONE_OVER_SIX_DELTAX)
          u_mm = pavx(Im2,j); u_m = pavx(Im1,j); u_p = pavx(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, pavx(i,j), u_p, value_dpavx_dx, ONE_OVER_SIX_DELTAX)
    else 
          u_pp = rhoavx(Ip2,j); u_p = rhoavx(Ip1,j); u_m = rhoavx(Im1,j) 
          call compute_decentered_forward_dU(u_m, rhoavx(i,j), u_p, u_pp, value_drhoavx_dx, ONE_OVER_SIX_DELTAX)
          u_pp = pavx(Ip2,j); u_p = pavx(Ip1,j); u_m = pavx(Im1,j) 
          call compute_decentered_forward_dU(u_m, pavx(i,j), u_p, u_pp, value_dpavx_dx, ONE_OVER_SIX_DELTAX)
    endif
    
    ! compute derivative of windy, vy according to y
    u_mm = vy(i,Jm2); u_m = vy(i,Jm1); u_p = vy(i,j); u_pp = vy(i,Jp1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      
    u_mm = windy(i,Jm2); u_m = windy(i,Jm1); u_p = windy(i,j); u_pp = windy(i,Jp1)
    call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dwindy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      
    ! compute derivative of rhoavy, pavy according to y 
    if (windy_half_y >= 0) then
          u_mm = pavy(i,Jm2)      ; u_m = pavy(i,Jm1)      ; u_p = pavy(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, pavy(i,j), u_p, value_dpavy_dy, ONE_OVER_SIX_DELTAY)
          u_mm = rhoavy(i,Jm2); u_m = rhoavy(i,Jm1); u_p = rhoavy(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rhoavy(i,j), u_p, value_drhoavy_dy, ONE_OVER_SIX_DELTAY)
    else 
          u_pp = pavy(i,Jp2)      ; u_p = pavy(i,Jp1)      ; u_m = pavy(i,Jm1) 
          call compute_decentered_forward_dU(u_m, pavy(i,j), u_p, u_pp, value_dpavy_dy, ONE_OVER_SIX_DELTAY)
          u_pp = rhoavy(i,Jp2); u_p = rhoavy(i,Jp1); u_m = rhoavy(i,Jm1) 
          call compute_decentered_forward_dU(u_m, rhoavy(i,j), u_p, u_pp, value_drhoavy_dy, ONE_OVER_SIX_DELTAY)
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
    

    ! intermediate computations
    value_vax_windx_dvx_dx = vax_half_x * windx_half_x * value_dvx_dx
    value_vax_windy_dvx_dy = vax_half_x * windy_half_y * value_dvx_dy
    value_vay_windx_dvy_dx = vay_half_y * windx_half_x * value_dvy_dx
    value_vay_windy_dvy_dy = vay_half_y * windy_half_y * value_dvy_dy
    
    value_vax_vx_dwindx_dx = vax_half_x * vx_half_x * value_dwindx_dx
    value_vax_vy_dwindx_dy = vax_half_x * vy_half_y * value_dwindx_dy
    value_vay_vx_dwindy_dx = vay_half_y * vx_half_x * value_dwindy_dx
    value_vay_vy_dwindy_dy = vay_half_y * vy_half_y * value_dwindy_dy


    if (i_global == ISOURCE .and. j_global == JSOURCE) then
      a = pi*pi*f0*f0
      t = dble(it_time-1)*DELTAT
      source_term = -8 * a* (t-t0) *  factor * exp(- 4 * a*(t-t0)*(t-t0))
      !distance2 = ((i - Isource) * DELTAX)**2 + ((j - Jsource) * DELTAY)**2
      factor_ssf = 1!exp( - distance2 / SSF_Sigma**2 )
    else
      source_term = 0
    endif

    ! Kernel of density
    K_rho0(i,j) = K_rho0(i,j) + rhoa(i,j) * (value_dvx_dx + value_dvy_dy) * DELTAT
    K_rho0(i,j) = K_rho0(i,j) - (value_drhoavx_dx + value_drhoavy_dy) * DELTAT
    K_rho0(i,j) = K_rho0(i,j) - (vax_dvx_dt_half_x + vay_dvy_dt_half_y) * DELTAT
    K_rho0(i,j) = K_rho0(i,j) + (value_vax_windx_dvx_dx + value_vax_windy_dvx_dy) * DELTAT
    K_rho0(i,j) = K_rho0(i,j) + (value_vay_windx_dvy_dx + value_vay_windy_dvy_dy) * DELTAT
    K_rho0(i,j) = K_rho0(i,j) + (value_vax_vx_dwindx_dx + value_vax_vy_dwindx_dy) * DELTAT
    K_rho0(i,j) = K_rho0(i,j) + (value_vay_vx_dwindy_dx + value_vay_vy_dwindy_dy) * DELTAT
    K_rho0(i,j) = K_rho0(i,j) + rhoa(i,j) * factor_ssf * source_term / (gamma_chimie * p0(i,j)) * DELTAT

    ! Kernel of pressure
    K_p0(i,j) = K_p0(i,j) - (value_dpavx_dx + value_dpavy_dy) * DELTAT
    K_p0(i,j) = K_p0(i,j) + gamma_chimie * pa(i,j) * (value_dvx_dx + value_dvy_dy) * DELTAT !
    K_p0(i,j) = K_p0(i,j) - pa(i,j) * factor_ssf * source_term * rho0(i,j) / (gamma_chimie * p0(i,j)**2)   * DELTAT

   enddo
  enddo


endsubroutine compute_kernel_iter
