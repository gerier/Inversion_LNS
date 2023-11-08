
subroutine load_frame(it_time, it_time_frame)
  use parameters ! , only : NX, NY, NSTEP, NFRAMES, FRAMES, ZERO
  implicit none
  
  integer :: i_frame, it_time_frame, it_time, step_between_frames
  
  ! compute the number of step between each frames
  step_between_frames = NSTEP / NFRAMES

  ! evaluate which frame is the closest and happening before it_time  
  i_frame = it_time / step_between_frames 
   
  if (i_frame == 0 ) then
  ! if it_time < step_between_frames, then the frame is the initial condition 
      call reset_forward()
      
      ! the frame is the initial condition (so iframe = i_start)
      it_time_frame = 1  
   else  
   ! if not, take the wavefield of the frame 
      pressure(:,:) = FRAMES(:,:,1,i_frame)
      rhop(:,:)     = FRAMES(:,:,2,i_frame)
      vx(:,:)       = FRAMES(:,:,3,i_frame)
      vy(:,:)       = FRAMES(:,:,4,i_frame)
      
      ! PML variables
      eq1_memory_dp0_dx_fw(:,:)       = eq1_memory_dp0_dx_fw_fr(:,:,i_frame)       
      eq1_memory_dp0_dy_fw(:,:)       = eq1_memory_dp0_dy_fw_fr(:,:,i_frame)       
      eq1_memory_drho0_dx_fw(:,:)     = eq1_memory_drho0_dx_fw_fr(:,:,i_frame)     
      eq1_memory_drho0_dy_fw(:,:)     = eq1_memory_drho0_dy_fw_fr(:,:,i_frame)     
      eq1_memory_dpressure_dx_fw(:,:) = eq1_memory_dpressure_dx_fw_fr(:,:,i_frame) 
      eq1_memory_dpressure_dy_fw(:,:) = eq1_memory_dpressure_dy_fw_fr(:,:,i_frame) 
      eq1_memory_drhop_dx_fw(:,:)     = eq1_memory_drhop_dx_fw_fr(:,:,i_frame)     
      eq1_memory_drhop_dy_fw(:,:)     = eq1_memory_drhop_dy_fw_fr(:,:,i_frame)     
      eq1_memory_dvx_dx_fw(:,:)       = eq1_memory_dvx_dx_fw_fr(:,:,i_frame)       
      eq1_memory_dvy_dy_fw(:,:)       = eq1_memory_dvy_dy_fw_fr(:,:,i_frame)       
      eq1_memory_dwindx_dx_fw(:,:)    = eq1_memory_dwindx_dx_fw_fr(:,:,i_frame)    
      eq1_memory_dwindy_dy_fw(:,:)    = eq1_memory_dwindy_dy_fw_fr(:,:,i_frame)    

      eq2_memory_dpressure_dx_fw(:,:) = eq2_memory_dpressure_dx_fw_fr(:,:,i_frame) 
      eq2_memory_drho0_dx_fw(:,:)     = eq2_memory_drho0_dx_fw_fr(:,:,i_frame)     
      eq2_memory_drho0_dy_fw(:,:)     = eq2_memory_drho0_dy_fw_fr(:,:,i_frame)     
      eq2_memory_dvx_dx_fw(:,:)       = eq2_memory_dvx_dx_fw_fr(:,:,i_frame)       
      eq2_memory_dvx_dy_fw(:,:)       = eq2_memory_dvx_dy_fw_fr(:,:,i_frame)       
      eq2_memory_dwindx_dx_fw(:,:)    = eq2_memory_dwindx_dx_fw_fr(:,:,i_frame)    
      eq2_memory_dwindx_dy_fw(:,:)    = eq2_memory_dwindx_dy_fw_fr(:,:,i_frame)    
      eq2_memory_dwindy_dy_fw(:,:)    = eq2_memory_dwindy_dy_fw_fr(:,:,i_frame)    

      eq3_memory_dpressure_dy_fw(:,:) = eq3_memory_dpressure_dy_fw_fr(:,:,i_frame) 
      eq3_memory_drho0_dy_fw(:,:)     = eq3_memory_drho0_dy_fw_fr(:,:,i_frame)     
      eq3_memory_drho0_dx_fw(:,:)     = eq3_memory_drho0_dx_fw_fr(:,:,i_frame)     
      eq3_memory_dvy_dy_fw(:,:)       = eq3_memory_dvy_dy_fw_fr(:,:,i_frame)       
      eq3_memory_dvy_dx_fw(:,:)       = eq3_memory_dvy_dx_fw_fr(:,:,i_frame)       
      eq3_memory_dwindy_dy_fw(:,:)    = eq3_memory_dwindy_dy_fw_fr(:,:,i_frame)    
      eq3_memory_dwindy_dx_fw(:,:)    = eq3_memory_dwindy_dx_fw_fr(:,:,i_frame)    
      eq3_memory_dwindx_dx_fw(:,:)    = eq3_memory_dwindx_dx_fw_fr(:,:,i_frame)  

      ! compute at which iteration the frame appears
      it_time_frame = i_frame * step_between_frames + 1
    endif
  
endsubroutine load_frame



subroutine save_frames()
  
 use parameters
 implicit none
 
 integer :: it_frame
 integer :: it_start, it_end
 integer :: step_between_frames
 
 ! compute the number of step between each frames
 step_between_frames = NSTEP / NFRAMES - 1
 
 ! initialisation
 it_end = 0
 
 ! make NFRAMES simulations (forward problem)
 ! the starting field (of the simulation) is the last frame saved
 ! the ending field (of the simulation) is the newt frame to save
 do it_frame=1,NFRAMES
   ! initialisation of starting and ending point
   it_start = it_end + 1
   it_end   = it_start + step_between_frames
   ! forward simulation
   call forwardproblem(p0_prior,rho0_prior,windx_prior,windy_prior,kappa_unrelaxed_prior, it_start, it_end,3) 
   ! save the new frame 
   FRAMES(:,:,1,it_frame) = pressure(:,:)
   FRAMES(:,:,2,it_frame) = rhop(:,:)
   FRAMES(:,:,3,it_frame) = vx(:,:)
   FRAMES(:,:,4,it_frame) = vy(:,:)

   eq1_memory_dp0_dx_fw_fr(:,:,it_frame)       = eq1_memory_dp0_dx_fw(:,:)
   eq1_memory_dp0_dy_fw_fr(:,:,it_frame)       = eq1_memory_dp0_dy_fw(:,:)
   eq1_memory_drho0_dx_fw_fr(:,:,it_frame)     = eq1_memory_drho0_dx_fw(:,:)
   eq1_memory_drho0_dy_fw_fr(:,:,it_frame)     = eq1_memory_drho0_dy_fw(:,:)
   eq1_memory_dpressure_dx_fw_fr(:,:,it_frame) = eq1_memory_dpressure_dx_fw(:,:)
   eq1_memory_dpressure_dy_fw_fr(:,:,it_frame) = eq1_memory_dpressure_dy_fw(:,:)
   eq1_memory_drhop_dx_fw_fr(:,:,it_frame)     = eq1_memory_drhop_dx_fw(:,:)
   eq1_memory_drhop_dy_fw_fr(:,:,it_frame)     = eq1_memory_drhop_dy_fw(:,:)
   eq1_memory_dvx_dx_fw_fr(:,:,it_frame)       = eq1_memory_dvx_dx_fw(:,:)
   eq1_memory_dvy_dy_fw_fr(:,:,it_frame)       = eq1_memory_dvy_dy_fw(:,:)
   eq1_memory_dwindx_dx_fw_fr(:,:,it_frame)    = eq1_memory_dwindx_dx_fw(:,:)
   eq1_memory_dwindy_dy_fw_fr(:,:,it_frame)    = eq1_memory_dwindy_dy_fw(:,:)
     
   eq2_memory_dpressure_dx_fw_fr(:,:,it_frame) = eq2_memory_dpressure_dx_fw(:,:)
   eq2_memory_drho0_dx_fw_fr(:,:,it_frame)     = eq2_memory_drho0_dx_fw(:,:)
   eq2_memory_drho0_dy_fw_fr(:,:,it_frame)     = eq2_memory_drho0_dy_fw(:,:)
   eq2_memory_dvx_dx_fw_fr(:,:,it_frame)       = eq2_memory_dvx_dx_fw(:,:)
   eq2_memory_dvx_dy_fw_fr(:,:,it_frame)       = eq2_memory_dvx_dy_fw(:,:)
   eq2_memory_dwindx_dx_fw_fr(:,:,it_frame)    = eq2_memory_dwindx_dx_fw(:,:)
   eq2_memory_dwindx_dy_fw_fr(:,:,it_frame)    = eq2_memory_dwindx_dy_fw(:,:)
   eq2_memory_dwindy_dy_fw_fr(:,:,it_frame)    = eq2_memory_dwindy_dy_fw(:,:)
      
   eq3_memory_dpressure_dy_fw_fr(:,:,it_frame) = eq3_memory_dpressure_dy_fw(:,:)
   eq3_memory_drho0_dy_fw_fr(:,:,it_frame)     = eq3_memory_drho0_dy_fw(:,:)
   eq3_memory_drho0_dx_fw_fr(:,:,it_frame)     = eq3_memory_drho0_dx_fw(:,:)
   eq3_memory_dvy_dy_fw_fr(:,:,it_frame)       = eq3_memory_dvy_dy_fw(:,:)
   eq3_memory_dvy_dx_fw_fr(:,:,it_frame)       = eq3_memory_dvy_dx_fw(:,:)
   eq3_memory_dwindy_dy_fw_fr(:,:,it_frame)    = eq3_memory_dwindy_dy_fw(:,:)
   eq3_memory_dwindy_dx_fw_fr(:,:,it_frame)    = eq3_memory_dwindy_dx_fw(:,:)
   eq3_memory_dwindx_dx_fw_fr(:,:,it_frame)    = eq3_memory_dwindx_dx_fw(:,:)
   
 enddo
 
endsubroutine save_frames

