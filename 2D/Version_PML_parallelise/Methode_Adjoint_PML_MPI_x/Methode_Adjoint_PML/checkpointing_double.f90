
subroutine load_frame(it_time, it_time_frame)
  use parameters ! , only : NX, NY, NSTEP, NFRAMES, FRAMES, ZERO
  implicit none
  
  integer :: i_frame, it_time_frame, it_time, step_frames, step_local_frames
 
  step_frames = NSTEP / NFRAMES 
  step_local_frames = step_frames / N_LOC_FRAMES 
  
  if (modulo(it_time, step_frames) < step_local_frames) then
      i_frame = it_time / step_frames
      
      it_time_frame = i_frame * step_frames + 1
      
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
  else 
      i_frame = (modulo(it_time, step_frames)) / step_local_frames
 
      it_time_frame = (it_time / (step_frames)) * step_frames + i_frame * step_local_frames + 1
      
      pressure(:,:) = LOC_FRAMES(:,:,1,i_frame)
      rhop(:,:)     = LOC_FRAMES(:,:,2,i_frame)
      vx(:,:)       = LOC_FRAMES(:,:,3,i_frame)
      vy(:,:)       = LOC_FRAMES(:,:,4,i_frame)
      
      ! PML variables
      eq1_memory_dp0_dx_fw(:,:)       = eq1_memory_dp0_dx_fw_loc_fr(:,:,i_frame)       
      eq1_memory_dp0_dy_fw(:,:)       = eq1_memory_dp0_dy_fw_loc_fr(:,:,i_frame)       
      eq1_memory_drho0_dx_fw(:,:)     = eq1_memory_drho0_dx_fw_loc_fr(:,:,i_frame)     
      eq1_memory_drho0_dy_fw(:,:)     = eq1_memory_drho0_dy_fw_loc_fr(:,:,i_frame)     
      eq1_memory_dpressure_dx_fw(:,:) = eq1_memory_dpressure_dx_fw_loc_fr(:,:,i_frame) 
      eq1_memory_dpressure_dy_fw(:,:) = eq1_memory_dpressure_dy_fw_loc_fr(:,:,i_frame) 
      eq1_memory_drhop_dx_fw(:,:)     = eq1_memory_drhop_dx_fw_loc_fr(:,:,i_frame)     
      eq1_memory_drhop_dy_fw(:,:)     = eq1_memory_drhop_dy_fw_loc_fr(:,:,i_frame)     
      eq1_memory_dvx_dx_fw(:,:)       = eq1_memory_dvx_dx_fw_loc_fr(:,:,i_frame)       
      eq1_memory_dvy_dy_fw(:,:)       = eq1_memory_dvy_dy_fw_loc_fr(:,:,i_frame)       
      eq1_memory_dwindx_dx_fw(:,:)    = eq1_memory_dwindx_dx_fw_loc_fr(:,:,i_frame)    
      eq1_memory_dwindy_dy_fw(:,:)    = eq1_memory_dwindy_dy_fw_loc_fr(:,:,i_frame)    

      eq2_memory_dpressure_dx_fw(:,:) = eq2_memory_dpressure_dx_fw_loc_fr(:,:,i_frame) 
      eq2_memory_drho0_dx_fw(:,:)     = eq2_memory_drho0_dx_fw_loc_fr(:,:,i_frame)     
      eq2_memory_drho0_dy_fw(:,:)     = eq2_memory_drho0_dy_fw_loc_fr(:,:,i_frame)     
      eq2_memory_dvx_dx_fw(:,:)       = eq2_memory_dvx_dx_fw_loc_fr(:,:,i_frame)       
      eq2_memory_dvx_dy_fw(:,:)       = eq2_memory_dvx_dy_fw_loc_fr(:,:,i_frame)       
      eq2_memory_dwindx_dx_fw(:,:)    = eq2_memory_dwindx_dx_fw_loc_fr(:,:,i_frame)    
      eq2_memory_dwindx_dy_fw(:,:)    = eq2_memory_dwindx_dy_fw_loc_fr(:,:,i_frame)    
      eq2_memory_dwindy_dy_fw(:,:)    = eq2_memory_dwindy_dy_fw_loc_fr(:,:,i_frame)    

      eq3_memory_dpressure_dy_fw(:,:) = eq3_memory_dpressure_dy_fw_loc_fr(:,:,i_frame) 
      eq3_memory_drho0_dy_fw(:,:)     = eq3_memory_drho0_dy_fw_loc_fr(:,:,i_frame)     
      eq3_memory_drho0_dx_fw(:,:)     = eq3_memory_drho0_dx_fw_loc_fr(:,:,i_frame)     
      eq3_memory_dvy_dy_fw(:,:)       = eq3_memory_dvy_dy_fw_loc_fr(:,:,i_frame)       
      eq3_memory_dvy_dx_fw(:,:)       = eq3_memory_dvy_dx_fw_loc_fr(:,:,i_frame)       
      eq3_memory_dwindy_dy_fw(:,:)    = eq3_memory_dwindy_dy_fw_loc_fr(:,:,i_frame)    
      eq3_memory_dwindy_dx_fw(:,:)    = eq3_memory_dwindy_dx_fw_loc_fr(:,:,i_frame)    
      eq3_memory_dwindx_dx_fw(:,:)    = eq3_memory_dwindx_dx_fw_loc_fr(:,:,i_frame) 
      
  endif
  
  
endsubroutine load_frame



subroutine save_frames()
  
 use parameters
 implicit none
 
 integer :: it_frame
 integer :: it_start, it_end
 integer :: step_frames
 
 ! compute the number of step between each frames
 step_frames = NSTEP / NFRAMES - 1
 
 ! initialisation
 it_end = 0
 
 ! make NFRAMES simulations (forward problem)
 ! the starting field (of the simulation) is the last frame saved
 ! the ending field (of the simulation) is the newt frame to save
 do it_frame=1,NFRAMES
   ! initialisation of starting and ending point
   it_start = it_end + 1
   it_end   = it_start + step_frames
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



subroutine save_local_frames(it_time)
  
 use parameters
 implicit none
 
 integer :: it_frame
 integer :: it_time
 integer :: it_start, it_end
 integer :: step_local_frames, step_frames
 
 ! compute the number of step between each frames
 step_frames = NSTEP / NFRAMES 
 step_local_frames = step_frames / N_LOC_FRAMES - 1
 
 ! Load last frame as initialisation
 call load_frame(it_time/step_frames*step_frames+1, it_end)
 it_end = it_end - 1
 
 ! make N_LOC_FRAMES simulations (forward problem)
 ! the starting field (of the simulation) is the last frame saved
 ! the ending field (of the simulation) is the newt frame to save
 do it_frame=1,N_LOC_FRAMES
   ! initialisation of starting and ending point
   it_start = it_end + 1
   it_end   = min(it_time, it_start + step_local_frames)
   
   ! forward simulation
   call forwardproblem(p0_prior,rho0_prior,windx_prior,windy_prior,kappa_unrelaxed_prior, it_start, it_end,3) 
   ! save the new frame 
      ! save the new frame 
   LOC_FRAMES(:,:,1,it_frame) = pressure(:,:)
   LOC_FRAMES(:,:,2,it_frame) = rhop(:,:)
   LOC_FRAMES(:,:,3,it_frame) = vx(:,:)
   LOC_FRAMES(:,:,4,it_frame) = vy(:,:)

   eq1_memory_dp0_dx_fw_loc_fr(:,:,it_frame)       = eq1_memory_dp0_dx_fw(:,:)
   eq1_memory_dp0_dy_fw_loc_fr(:,:,it_frame)       = eq1_memory_dp0_dy_fw(:,:)
   eq1_memory_drho0_dx_fw_loc_fr(:,:,it_frame)     = eq1_memory_drho0_dx_fw(:,:)
   eq1_memory_drho0_dy_fw_loc_fr(:,:,it_frame)     = eq1_memory_drho0_dy_fw(:,:)
   eq1_memory_dpressure_dx_fw_loc_fr(:,:,it_frame) = eq1_memory_dpressure_dx_fw(:,:)
   eq1_memory_dpressure_dy_fw_loc_fr(:,:,it_frame) = eq1_memory_dpressure_dy_fw(:,:)
   eq1_memory_drhop_dx_fw_loc_fr(:,:,it_frame)     = eq1_memory_drhop_dx_fw(:,:)
   eq1_memory_drhop_dy_fw_loc_fr(:,:,it_frame)     = eq1_memory_drhop_dy_fw(:,:)
   eq1_memory_dvx_dx_fw_loc_fr(:,:,it_frame)       = eq1_memory_dvx_dx_fw(:,:)
   eq1_memory_dvy_dy_fw_loc_fr(:,:,it_frame)       = eq1_memory_dvy_dy_fw(:,:)
   eq1_memory_dwindx_dx_fw_loc_fr(:,:,it_frame)    = eq1_memory_dwindx_dx_fw(:,:)
   eq1_memory_dwindy_dy_fw_loc_fr(:,:,it_frame)    = eq1_memory_dwindy_dy_fw(:,:)
     
   eq2_memory_dpressure_dx_fw_loc_fr(:,:,it_frame) = eq2_memory_dpressure_dx_fw(:,:)
   eq2_memory_drho0_dx_fw_loc_fr(:,:,it_frame)     = eq2_memory_drho0_dx_fw(:,:)
   eq2_memory_drho0_dy_fw_loc_fr(:,:,it_frame)     = eq2_memory_drho0_dy_fw(:,:)
   eq2_memory_dvx_dx_fw_loc_fr(:,:,it_frame)       = eq2_memory_dvx_dx_fw(:,:)
   eq2_memory_dvx_dy_fw_loc_fr(:,:,it_frame)       = eq2_memory_dvx_dy_fw(:,:)
   eq2_memory_dwindx_dx_fw_loc_fr(:,:,it_frame)    = eq2_memory_dwindx_dx_fw(:,:)
   eq2_memory_dwindx_dy_fw_loc_fr(:,:,it_frame)    = eq2_memory_dwindx_dy_fw(:,:)
   eq2_memory_dwindy_dy_fw_loc_fr(:,:,it_frame)    = eq2_memory_dwindy_dy_fw(:,:)
      
   eq3_memory_dpressure_dy_fw_loc_fr(:,:,it_frame) = eq3_memory_dpressure_dy_fw(:,:)
   eq3_memory_drho0_dy_fw_loc_fr(:,:,it_frame)     = eq3_memory_drho0_dy_fw(:,:)
   eq3_memory_drho0_dx_fw_loc_fr(:,:,it_frame)     = eq3_memory_drho0_dx_fw(:,:)
   eq3_memory_dvy_dy_fw_loc_fr(:,:,it_frame)       = eq3_memory_dvy_dy_fw(:,:)
   eq3_memory_dvy_dx_fw_loc_fr(:,:,it_frame)       = eq3_memory_dvy_dx_fw(:,:)
   eq3_memory_dwindy_dy_fw_loc_fr(:,:,it_frame)    = eq3_memory_dwindy_dy_fw(:,:)
   eq3_memory_dwindy_dx_fw_loc_fr(:,:,it_frame)    = eq3_memory_dwindy_dx_fw(:,:)
   eq3_memory_dwindx_dx_fw_loc_fr(:,:,it_frame)    = eq3_memory_dwindx_dx_fw(:,:)
   
 enddo
endsubroutine save_local_frames

