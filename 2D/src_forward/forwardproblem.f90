subroutine forwardproblem(p0,rho0,windx,windy, it_start, it_end, type_number)
 !

use parameters
implicit none

  integer :: it_start, it_end,&
             ix_irec, iy_irec 

  integer :: type_number

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: &
    p0, rho0, windx, windy

  double precision :: &
      value_dvx_dx, &
      value_dvy_dy, &
      value_dpressure_dx, &
      value_dpressure_dy, &
      value_dp0_dy, &
      value_drho0_dy, &
      value_drhop_dx, &
      value_drhop_dy, &
      value_vdrho0,    &
      value_v0drhop, &
      value_vdp0 , &
      value_v0dp

  double precision :: &
    value_dvy_dx, &
    value_dvx_dy, &
    value_drho0_dy_prec, &
    value_drho0_dy_next

  double precision :: &
      vx_half_x, &
      vy_half_y,&
      windx_half_x, &
      windx_half_x_half_y, &
      vy_half_x_half_y, &
      vx_half_x_half_y

   double precision :: &
      value_vdwindx, &
      value_v0dvy, &
      value_v0dvx, &
      value_dwindx_dy

  double precision :: maxvalue
        
  double precision  :: u_mm, u_m, u_p, u_pp
  
  integer :: Ip1,Im1, Ip2, Im2, Jp1,Jm1, Jp2, Jm2
  integer :: i,j,it,irec

  double precision :: velocnorm,pressurenorm

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: rhop_old, p_old, vx_old, vy_old

 ! to interpolate material parameters or velocity at the right location in the staggered grid cell
  double precision :: rho0_half_x,rho0_half_y,rhop_half_x,rhop_half_y,g_half_y

!---
!--- program starts here
!---

  do it = it_start, it_end

!-----------------------------------------------------------------------
! compute pressure, density and update memory variables for C-PML
!-----------------------------------------------------------------------

  ! store pressure and density values in axiliary variable
  p_old(:,:) = pressure(:,:)
  rhop_old(:,:) = rhop(:,:)
  
  ! send information to neighbooring processus
  call send_receive_rightleft(vx)
  call send_receive_rightleft(vy)
  call send_receive_topbottom(vy)
  call send_receive_rightleft(rhop_old)
  call send_receive_topbottom(rhop_old)
  call send_receive_rightleft(p_old)
  call send_receive_topbottom(p_old)

  do j = 1,NY_LOCAL
    do i = 1,NX_LOCAL

       ! index definition
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


       ! interpolate material parameters at the right location in the staggered grid cell
       vx_half_x = (vx(Ip1,j) + vx(i,j)) * 0.5d0
       vy_half_y = (vy(i,j) + vy(i,Jm1)) * 0.5d0
       windx_half_x = (windx(Ip1,j) + windx(i,j)) * 0.5d0


       ! derivative computations : vx, windx according to x
       u_mm = vx(Im1,j); u_m = vx(i,j); u_p = vx(Ip1,j); u_pp = vx(Ip2,j)
       call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvx_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
        
       eq1_memory_dvx_dx_fw(i,j) = b_x_half(i_global) * eq1_memory_dvx_dx_fw(i,j) + a_x_half(i_global) * value_dvx_dx
	
	
       ! derivative computations : vy, windy according to y
       u_mm = vy(i,Jm2); u_m = vy(i,Jm1); u_p = vy(i,j); u_pp = vy(i,Jp1)
       call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dvy_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      
       eq1_memory_dvy_dy_fw(i,j) = b_y(j_global) * eq1_memory_dvy_dy_fw(i,j) + a_y(j_global) * value_dvy_dy       


       ! decentered derivative of p0, p, rho0 and rhop in x
       if (windx_half_x >= 0) then
          u_mm = p_old(Im2,j); u_m = p_old(Im1,j); u_p = p_old(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, p_old(i,j), u_p, value_dpressure_dx, ONE_OVER_SIX_DELTAX)
         
          u_mm = rhop_old(Im2,j)     ; u_m = rhop_old(Im1,j)     ; u_p = rhop_old(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, rhop_old(i,j), u_p, value_drhop_dx, ONE_OVER_SIX_DELTAX)
       else 
          u_pp = p_old(Ip2,j); u_p = p_old(Ip1,j); u_m = p_old(Im1,j) 
          call compute_decentered_forward_dU(u_m, p_old(i,j), u_p, u_pp, value_dpressure_dx, ONE_OVER_SIX_DELTAX)
         
          u_pp = rhop_old(Ip2,j)     ; u_p = rhop_old(Ip1,j)     ; u_m = rhop_old(Im1,j) 
          call compute_decentered_forward_dU(u_m, rhop_old(i,j), u_p, u_pp, value_drhop_dx, ONE_OVER_SIX_DELTAX)
       endif
     
       eq1_memory_dpressure_dx_fw(i,j) = b_x_half(i_global) * eq1_memory_dpressure_dx_fw(i,j) &
                 + a_x_half(i_global) * value_dpressure_dx
       eq1_memory_drhop_dx_fw(i,j) = b_x_half(i_global) * eq1_memory_drhop_dx_fw(i,j) + a_x_half(i_global) * value_drhop_dx
             
              
       ! decentered derivative of p0, p, rho0 and rhop in y  
          u_mm = p0(i,Jm2)      ; u_m = p0(i,Jm1)      ; u_p = p0(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, p0(i,j), u_p, value_dp0_dy, ONE_OVER_SIX_DELTAY)
          u_mm = p_old(i,Jm2); u_m = p_old(i,Jm1); u_p = p_old(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, p_old(i,j), u_p, value_dpressure_dy, ONE_OVER_SIX_DELTAY)
          u_mm = rho0(i,Jm2)    ; u_m = rho0(i,Jm1)    ; u_p = rho0(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rho0(i,j), u_p, value_drho0_dy, ONE_OVER_SIX_DELTAY)
          u_mm = rhop_old(i,Jm2)     ; u_m = rhop_old(i,Jm1)     ; u_p = rhop_old(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rhop_old(i,j), u_p, value_drhop_dy, ONE_OVER_SIX_DELTAY)
    
              
       eq1_memory_dpressure_dy_fw(i,j) = b_y(j_global) * eq1_memory_dpressure_dy_fw(i,j) + a_y(j_global) * value_dpressure_dy
       eq1_memory_drhop_dy_fw(i,j) = b_y(j_global) * eq1_memory_drhop_dy_fw(i,j) + a_y(j_global) * value_drhop_dy
       eq1_memory_dp0_dy_fw(i,j) = c_y(j_global) * value_dp0_dy
       eq1_memory_drho0_dy_fw(i,j) = c_y(j_global) * value_drho0_dy
 
                      
       ! derivative updates for PML
       value_dvx_dx = value_dvx_dx * one_over_K_x_half(i_global) + eq1_memory_dvx_dx_fw(i,j)
       value_dvy_dy = value_dvy_dy * one_over_K_y(j_global) + eq1_memory_dvy_dy_fw(i,j) 
    
       value_dpressure_dx = value_dpressure_dx * one_over_K_x_half(i_global) + eq1_memory_dpressure_dx_fw(i,j)
       value_drhop_dx = value_drhop_dx * one_over_K_x_half(i_global) + eq1_memory_drhop_dx_fw(i,j)
        
       value_dp0_dy = value_dp0_dy * one_over_Kdalpha_y(j_global) 
       value_drho0_dy = value_drho0_dy * one_over_Kdalpha_y(j_global) 
       value_dpressure_dy = value_dpressure_dy * one_over_K_y(j_global) + eq1_memory_dpressure_dy_fw(i,j) 
       value_drhop_dy = value_drhop_dy * one_over_K_y(j_global) + eq1_memory_drhop_dy_fw(i,j)        
        
        
       ! intermediate computations
       value_vdp0   = (vy_half_y  * value_dp0_dy)
       value_v0dp   = (windx_half_x * value_dpressure_dx)
       value_vdrho0   =  (vy_half_y  * value_drho0_dy)
       value_v0drhop = (windx_half_x * value_drhop_dx)
        
        
       ! updateT
       pressure(i,j) = pressure(i,j) - (gamma_chimie(i,j) * p0(i,j) * (value_dvx_dx + value_dvy_dy)) * DELTAT
       pressure(i,j) = pressure(i,j) - value_vdp0 * DELTAT 
       pressure(i,j) = pressure(i,j) - value_v0dp * DELTAT

       rhop(i,j) = rhop(i,j) - rho0(i,j) * (value_dvx_dx + value_dvy_dy) * DELTAT
       rhop(i,j) = rhop(i,j) - value_vdrho0 * DELTAT
       rhop(i,j) = rhop(i,j) - value_v0drhop * DELTAT

      enddo
    enddo

    ! Dircihlet conditions  
    !! Left conditions
    if (USE_PML_XMIN .and. i_rank == 0) then
     if (type_source == 1 .and. wavefront == 2 .and. USE_PML_XMIN .and. USE_PML_XMAX) then
      do j=1,NY_LOCAL
        rhop(-1:NPOINTS_PML,j) = rhop(NPOINTS_PML,j)
        pressure(-1:NPOINTS_PML,j) = pressure(NPOINTS_PML,j)
      enddo
     else
      rhop(-1:1,:) = ZERO
      pressure(-1:1,:) = ZERO 
     endif
    endif
    !! Right conditions
    if (USE_PML_XMAX .and. i_rank == NPROC_X-1) then
     if (type_source == 1 .and. wavefront == 2 .and. USE_PML_XMIN .and. USE_PML_XMAX) then
      do j=1,NY_LOCAL
        rhop(NX_LOCAL:NX_LOCAL+2,j) = rhop(NX_LOCAL,j)
        pressure(NX_LOCAL:NX_LOCAL+2,j) = pressure(NX_LOCAL,j)
      enddo
     else
      rhop(NX_LOCAL:NX_LOCAL+2,:) = ZERO
      pressure(NX_LOCAL:NX_LOCAL+2,:) = ZERO
     endif
    endif
    !! Bottom conditions  
    if (USE_PML_YMIN .and. j_rank == 0) then
      rhop(:,-1:1) = ZERO
      pressure(:,-1:1) = ZERO
    else if (j_rank == 0) then
       do j=-1,1
         rhop(:,j) = rhop(:,3-j)
        pressure(:,j) = pressure(:,3-j)
       enddo 
    endif
    !! Top conditions
    if (USE_PML_YMAX .and. j_rank == NPROC_Y-1) then
      rhop(:,NY_LOCAL:NY_LOCAL+2) = ZERO
      pressure(:,NY_LOCAL:NY_LOCAL+2) = ZERO
    endif



  ! add source
  ! add the source (pressure located at a given grid point)
  t = dble(it-1)*DELTAT

  ! Gaussian
  !  source_term = factor * exp(- 4 * a*(t-t0)*(t-t0))  
  ! derivative of guassian
  source_term = - 2 * a* (t-t0) *  factor * exp(- a*(t-t0)*(t-t0))   
  ! define location of the source  
  if (type_source == 1) then ! plane wave case 
    if (wavefront == 1 .and. i_rank == ISOURCE / NX_LOCAL) then
      i = ISOURCE - offset_i
      do j = 1,NY_LOCAL
        pressure(i,j) = pressure(i,j) + source_term * DELTAT * gamma_chimie(i,j) * p0(i,j) / rho0(i,j)
        rhop(i,j) = rhop(i,j) + source_term * DELTAT
    enddo
    elseif (wavefront ==2 .and. j_rank == JSOURCE / NY_LOCAL) then
      j = JSOURCE - offset_j
      do i = 1,NX_LOCAL
        pressure(i,j) = pressure(i,j) + source_term * DELTAT * gamma_chimie(i,j) * p0(i,j) / rho0(i,j)
        rhop(i,j) = rhop(i,j) + source_term * DELTAT 
      enddo
    endif
    
  elseif (type_source == 2 ) then ! Point Source
    ! point source case 
    j = JSOURCE - offset_j
    i = ISOURCE - offset_i
    if (i_rank == ISOURCE / NX_LOCAL .and. j_rank == JSOURCE / NY_LOCAL) then   
      pressure(i,j) = pressure(i,j) + source_term * DELTAT * gamma_chimie(i,j) * p0(i,j) / rho0(i,j)
      rhop(i,j) = rhop(i,j) + source_term * DELTAT
    endif
     
   else if (type_source ==3) then  ! Point source case with SPREAD_SSF
      do j = 1,NY_LOCAL
       do i = 1,NX_LOCAL
 
          distance2 = ((i + offset_i - Isource) * DELTAX)**2 + ((j + offset_j - Jsource) * DELTAY)**2
          factor_ssf = exp( - distance2 / SSF_Sigma**2 )
         
          pressure(i,j) = pressure(i,j) + source_term * factor_ssf * DELTAT * gamma_chimie(i,j) * p0(i,j) / rho0(i,j)
          rhop(i,j) = rhop(i,j) + source_term * factor_ssf * DELTAT 
       enddo
      enddo

  endif  ! type_source
    
  ! write the source
  if (rank == 0) then
   open(unit=211,file='./OUTPUT/source_time_function_model.dat',status='unknown',position="append",action='write')
    write(211,*) (t-t0), source_term
   close(211)
  endif
  


!--------------------------------------------------------
! compute velocity and update memory variables for C-PML
!--------------------------------------------------------

  ! store pressure and density values in axiliary variable
  vx_old(:,:) = vx(:,:)
  vy_old(:,:) = vy(:,:)



  !  send information to neighbooring processus
  call send_receive_rightleft(rhop)
  call send_receive_rightleft(pressure)
  call send_receive_rightleft(vx_old)
  call send_receive_topbottom(vx_old)
  call send_receive_topbottom(vy_old)
  
  call MPI_SENDRECV(vy_old(NX_LOCAL-1:NX_LOCAL,:),number_of_values_x,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,vy_old(-1:0,:),number_of_values_x, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
   
  call MPI_SENDRECV(vy_old(NX_LOCAL-1:NX_LOCAL,NY_LOCAL-1:NY_LOCAL),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_right_top_shift,message_tag,vy_old(-1:0,-1:0),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_right_top_shift,message_tag,MPI_COMM_WORLD,message_status,code)
    

  do j = 1,NY_LOCAL
    do i = 1,NX_LOCAL
    
      ! index definition
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
      
      ! some interpolations
      rho0_half_x  = 0.5d0 * (rho0(i,j)  + rho0(Im1,j))
      rhop_half_x = 0.5d0 * (rhop(i,j) + rhop(Im1,j))
      vy_half_x_half_y  = (vy_old(i,j)    + vy_old(i,Jm1)+ vy_old(Im1,j)+ vy_old(Im1,Jm1))* 0.25d0

      ! compute derivative of the pressure, density according to x
      u_mm = pressure(Im2,j); u_m = pressure(Im1,j); u_p = pressure(i,j); u_pp = pressure(Ip1,j)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dpressure_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX)
      
      eq2_memory_dpressure_dx_fw(i,j) = b_x(i_global) * eq2_memory_dpressure_dx_fw(i,j) + a_x(i_global) * value_dpressure_dx
      
      
      ! compute derivative of windx, vx according to x
      if (windx(i,j) >= 0) then
          u_mm = vx_old(Im2,j); u_m = vx_old(Im1,j); u_p = vx_old(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, vx_old(i,j), u_p, value_dvx_dx, ONE_OVER_SIX_DELTAX)
     else 
          u_pp = vx_old(Ip2,j)     ; u_p = vx_old(Ip1,j)     ; u_m = vx_old(Im1,j) 
          call compute_decentered_forward_dU(u_m, vx_old(i,j), u_p, u_pp, value_dvx_dx, ONE_OVER_SIX_DELTAX)
     endif
     
      eq2_memory_dvx_dx_fw(i,j) = b_x(i_global) * eq2_memory_dvx_dx_fw(i,j) + a_x(i_global) * value_dvx_dx
      

     ! compute derivative of windx, vx according to y
          u_mm = windx(i,Jm2)      ; u_m = windx(i,Jm1)      ; u_p = windx(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, windx(i,j), u_p, value_dwindx_dy, ONE_OVER_SIX_DELTAY)
          u_mm = vx_old(i,Jm2); u_m = vx_old(i,Jm1); u_p = vx_old(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, vx_old(i,j), u_p, value_dvx_dy, ONE_OVER_SIX_DELTAY)
 
      eq2_memory_dwindx_dy_fw(i,j) = c_y(j_global) * value_dwindx_dy
      eq2_memory_dvx_dy_fw(i,j) = b_y(j_global) * eq2_memory_dvx_dy_fw(i,j) + a_y(j_global) * value_dvx_dy
     
     
     ! compute derivative of rho0 according to y (evaluating in vx point) 
          u_mm = rho0(Im1,Jm2)      ; u_m = rho0(Im1,Jm1)      ; u_p = rho0(Im1,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rho0(Im1,j), u_p, value_drho0_dy_prec, ONE_OVER_SIX_DELTAY)
          u_mm = rho0(i,Jm2)      ; u_m = rho0(i,Jm1)      ; u_p = rho0(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, rho0(i,j), u_p, value_drho0_dy_next, ONE_OVER_SIX_DELTAY)
     value_drho0_dy = 0.5d0 * (value_drho0_dy_next + value_drho0_dy_prec)
     
     eq2_memory_drho0_dy_fw(i,j) = c_y(j_global) * value_drho0_dy


     ! derivative updates for PML
     value_dpressure_dx = value_dpressure_dx * one_over_K_x(i_global) + eq2_memory_dpressure_dx_fw(i,j)
     value_dvx_dx = value_dvx_dx * one_over_K_x(i_global) + eq2_memory_dvx_dx_fw(i,j)
     value_dvx_dy = value_dvx_dy * one_over_K_y(j_global) + eq2_memory_dvx_dy_fw(i,j)
     value_dwindx_dy = value_dwindx_dy * one_over_Kdalpha_y(j_global) !+ eq2_memory_dwindx_dy_fw(i,j)
     value_drho0_dy = value_drho0_dy * one_over_Kdalpha_y(j_global)! + eq2_memory_drho0_dy_fw(i,j)

     
     ! intermediate computations: (v0 . nabla) v0_x ; (v' . nabla) v0_x ; (v0 . nabla) v0_x
     value_vdwindx  =  vy_half_x_half_y  * value_dwindx_dy
     value_v0dvx  = windx(i,j)     * value_dvx_dx  


     ! update
     vx(i,j) = vx(i,j) - value_dpressure_dx * DELTAT / rho0_half_x
     vx(i,j) = vx(i,j) - value_v0dvx  * DELTAT
     vx(i,j) = vx(i,j) - value_vdwindx  * DELTAT
     
   enddo
  enddo


  ! send information to neighbooring processus
  call send_receive_topbottom(rhop)
  call send_receive_topbottom(pressure)
  call send_receive_rightleft(vy_old)
  call send_receive_topbottom(vy_old)
  call send_receive_topbottom(vx_old)
  
  call MPI_SENDRECV(vy_old(NX_LOCAL-1:NX_LOCAL,:),number_of_values_x,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,vy_old(-1:0,:),number_of_values_x, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
   
  call MPI_SENDRECV(vx_old(1:2,1:2),number_of_values_corner,MPI_DOUBLE_PRECISION, &
        receiver_bottom_left_shift,message_tag,vx_old(NX_LOCAL+1:NX_LOCAL+2,NY_LOCAL+1:NY_LOCAL+2),number_of_values_corner, &
        MPI_DOUBLE_PRECISION,sender_bottom_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)

   
  do j = 1,NY_LOCAL
    do i = 1,NX_LOCAL

      ! index definition
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


      ! interpolate density at the right location in the staggered grid cell
      rho0_half_y  = 0.5d0 * (rho0(i,j)  + rho0(i,Jp1))
      rhop_half_y = 0.5d0 * (rhop(i,j) + rhop(i,Jp1))
      g_half_y = 0.5d0 * (g(i,j) + g(i,Jp1))

      vx_half_x_half_y  = (vx_old(Ip1,Jp1)  + vx_old(i,Jp1)  +vx_old(Ip1,j)  + vx_old(i,j))  * 0.25d0 
      windx_half_x_half_y = (windx(Ip1,Jp1) + windx(i,Jp1) +windx(Ip1,j) + windx(i,j)) * 0.25d0


      ! compute derivative of the pressure and density according to y
      u_mm = pressure(i,Jm1); u_m = pressure(i,j); u_p = pressure(i,Jp1); u_pp = pressure(i,Jp2)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_dpressure_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)
      u_mm = rho0(i,Jm1); u_m = rho0(i,j); u_p = rho0(i,Jp1); u_pp = rho0(i,Jp2)
      call compute_centered_dU(u_mm, u_m, u_p, u_pp, value_drho0_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY)

      eq3_memory_dpressure_dy_fw(i,j) = b_y_half(j_global) * eq3_memory_dpressure_dy_fw(i,j) &
                + a_y_half(j_global) * value_dpressure_dy
      eq3_memory_drho0_dy_fw(i,j) = c_y_half(j_global) * value_drho0_dy


      ! compute derivative of windy, vy according to x
      if (windx_half_x_half_y >= 0) then

          u_mm = vy_old(Im2,j); u_m = vy_old(Im1,j); u_p = vy_old(Ip1,j) 
          call compute_decentered_backward_dU(u_mm, u_m, vy_old(i,j), u_p, value_dvy_dx, ONE_OVER_SIX_DELTAX)
      else 

          u_pp = vy_old(Ip2,j); u_p = vy_old(Ip1,j); u_m = vy_old(Im1,j) 
          call compute_decentered_forward_dU(u_m, vy_old(i,j), u_p, u_pp, value_dvy_dx, ONE_OVER_SIX_DELTAX)
      endif
   
      eq3_memory_dvy_dx_fw(i,j) = b_y_half(j_global) * eq3_memory_dvy_dx_fw(i,j) &
                + a_y_half(j_global) * value_dvy_dx
   
      
  
      ! compute derivative of windy, vy according to y
          u_mm = vy_old(i,Jm2); u_m = vy_old(i,Jm1); u_p = vy_old(i,Jp1) 
          call compute_decentered_backward_dU(u_mm, u_m, vy_old(i,j), u_p, value_dvy_dy, ONE_OVER_SIX_DELTAX)
   
      eq3_memory_dvy_dy_fw(i,j) = b_y_half(j_global) * eq3_memory_dvy_dy_fw(i,j) + a_y_half(j_global) * value_dvy_dy
            

      ! derivative updates for PML
      value_dpressure_dy = value_dpressure_dy * one_over_K_y_half(j_global) + eq3_memory_dpressure_dy_fw(i,j)
      value_drho0_dy = value_drho0_dy * one_over_Kdalpha_y_half(j_global) !+ eq3_memory_drho0_dy_fw(i,j)
      
      value_dvy_dy = value_dvy_dy *  one_over_K_y_half(j_global) + eq3_memory_dvy_dy_fw(i,j)
    
      value_dvy_dx = value_dvy_dx *  one_over_K_x_half(i_global) + eq3_memory_dvy_dx_fw(i,j)
   
      
      
      ! intermediate computations: (v0 . nabla) v0_y ; (v' . nabla) v0_y ; (v0 . nabla) v0_y
      value_v0dvy  = windx_half_x_half_y * value_dvy_dx 

   
      ! update
      vy(i,j) = vy(i,j) - value_dpressure_dy * DELTAT / rho0_half_y
      vy(i,j) = vy(i,j) - value_v0dvy  * DELTAT
      vy(i,j) = vy(i,j) - rhop_half_y * g_half_y * DELTAT / rho0_half_y
 
    enddo
  enddo

  ! Dircihlet conditions
  !! Left conditions
  if (USE_PML_XMIN .and. i_rank == 0) then
    if (type_source == 1 .and. wavefront == 2 .and. USE_PML_XMIN .and. USE_PML_XMAX) then
      do j=1,NY_LOCAL
        vx(-1:NPOINTS_PML,j) = vx(NPOINTS_PML,j)
        vy(-1:NPOINTS_PML,j) = vy(NPOINTS_PML,j)
      enddo
    else 
      vx(-1:1,:) = ZERO
      vy(-1:1,:) = ZERO
    endif
  endif
  !! Right conditions
  if (USE_PML_XMAX .and. i_rank == NPROC_X-1) then
    if (type_source == 1 .and. wavefront == 2 .and. USE_PML_XMIN .and. USE_PML_XMAX) then
      do j=1,NY_LOCAL
        vx(NX_LOCAL-NPOINTS_PML:NX_LOCAL+2,j) = vx(NX_LOCAL-NPOINTS_PML,j)
        vy(NX_LOCAL-NPOINTS_PML:NX_LOCAL+2,j) = vy(NX_LOCAL-NPOINTS_PML,j)
      enddo
    else
      vx(NX_LOCAL:NX_LOCAL+2,:) = ZERO
      vy(NX_LOCAL:NX_LOCAL+2,:) = ZERO   
    endif
  endif 
  !! Bottom conditions
  if (USE_PML_YMIN .and. j_rank == 0) then ! PML
      vx(:,-1:1) = ZERO
      vy(:,-1:1) = ZERO
  elseif (j_rank == 0) then ! Reflection on the ground
     do j=-1,1
        vx(:,j) = -vx(:,3-j)
        vy(:,j) = -vy(:,2-j)
     enddo
     vy(:,1) = ZERO
  endif
  !! Top conditions
  if (USE_PML_YMAX .and. j_rank == NPROC_Y-1) then
      vx(:,NY_LOCAL:NY_LOCAL+2) = ZERO
      vy(:,NY_LOCAL:NY_LOCAL+2) = ZERO
  endif


  ! store seismograms
  do irec = 1,NREC
     ! beware here that the two components of the velocity vector are not defined at the same point
     ! in a staggered grid, and thus the two components of the velocity vector are recorded at slightly different locations,
     ! vy is staggered by half a grid cell along X and along Y with respect to vx
     if (((ix_rec(irec)-1) / NX_LOCAL) == i_rank .and. (iy_rec(irec)-1) / NY_LOCAL == j_rank ) then
  
        ix_irec = ix_rec(irec) - offset_i
        iy_irec = iy_rec(irec) - offset_j
   
        sisvx(it,irec) = vx(ix_irec,iy_irec)
        sisvy(it,irec) = vy(ix_irec,iy_irec)
        sispressure(it,irec) = pressure(ix_irec,iy_irec)
        sisrhop(it,irec) = rhop(ix_irec,iy_irec)
      endif 
  enddo

   ! update the maxvalue to create normalised image
   if (save_normimage_overtime == 1) then
     call MPI_ALLREDUCE(maxval(abs(pressure(1:NX_LOCAL,1:NY_LOCAL))), &
         maxvalue,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,code)
     if (maxval_image_p < maxvalue) then
       maxval_image_p = maxvalue
     endif

     call MPI_ALLREDUCE(maxval(abs(rhop(1:NX_LOCAL,1:NY_LOCAL))), &
         maxvalue,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,code)
     if (maxval_image_rho < maxvalue) then
       maxval_image_rho = maxvalue
     endif
 
     call MPI_ALLREDUCE(maxval(abs(vx(1:NX_LOCAL,1:NY_LOCAL))), &
         maxvalue,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,code)
     if (maxval_image_vx < maxvalue) then
       maxval_image_vx = maxvalue
     endif

     call MPI_ALLREDUCE(maxval(abs(vy(1:NX_LOCAL,1:NY_LOCAL))), &
         maxvalue,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,code)
     if (maxval_image_vy < maxvalue) then
       maxval_image_vy = maxvalue
     endif
   endif


   ! output information
   if ((mod(it,IT_DISPLAY) == 0 .or. it == 5) .and. (save_sismos)) then
 
     call MPI_ALLREDUCE(maxval(abs(pressure(1:NX_LOCAL,1:NY_LOCAL))), &
        pressurenorm,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,code)   
     call MPI_REDUCE(maxval(sqrt(vx(1:NX_LOCAL,1:NY_LOCAL)**2 + vy(1:NX_LOCAL,1:NY_LOCAL)**2)), &
        velocnorm,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,code)
     call MPI_BARRIER(MPI_COMM_WORLD, code)
      
     ! display informations
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

     endif
     
     ! create snapshots of pressure, velocity and density and save seismograms  
     call gather_and_generate_image(vx,vy,pressure,rhop,it,type_number)
     call write_seismograms(sisvx,sisvy,sispressure,sisrhop,NSTEP,NREC, &
             DELTAT,t0,type_number)
  
   endif
   
 enddo   ! end of time loop


 ! save seismograms
 if (save_sismos) then
     call write_seismograms(sisvx,sisvy,sispressure,sisrhop,NSTEP,NREC,DELTAT,t0,type_number)
 endif
   
 if (rank==0 .and. method <= 2) then 
   if (it == NSTEP) then
       print *
       print *,'End of the simulation'
       print *
   endif
 endif

endsubroutine forwardproblem
