subroutine load_frame(it_time, it_time_frame)
!==============================================================================
! Load a saved simulation state.
!
! This routine restores the complete simulation state stored in the closest
! available frame before or at the requested simulation iteration.
!
! Two types of saved frames are considered:
!   - global frames stored at regular time intervals
!   - local frames stored between two consecutive global frames to provide
!     a finer temporal resolution
!
! Restored quantities:
!   - pressure field
!   - density perturbation field
!   - velocity components
!   - PML recursive convolution memory variables
!
! Inputs:
!   it_time
!       Requested simulation iteration index.
!
! Outputs:
!   it_time_frame
!       Iteration index corresponding to the restored frame.
!
! The restored state can be directly used as the starting point of a
! subsequent forward simulation.
!
!==============================================================================
  use parameters
  implicit none
  integer, intent(in) :: it_time
  integer, intent(out) :: it_time_frame
  integer :: i_frame, step_frames, step_local_frames
 
  ! Compute the number of iterations between two global frames.
  step_frames = NSTEP / N_GLOB_FRAMES 
  
  ! Compute the number of iterations between two local frames.
  step_local_frames = step_frames / N_LOC_FRAMES 
  
  if (modulo(it_time, step_frames) < step_local_frames) then 
  ! Determine whether the requested iteration belongs to a global frame
  ! or to a local intermediate frame.
  
      ! Index of the selected frame.
      i_frame = it_time / step_frames
  
      if (i_frame == 0) then
          ! If the requested iteration precedes the first saved frame,
          ! Restore the initial state corresponding to the beginning of the simulation
          ! (iteration index 1, representing the interval from t = 0 to t = Δt).
          call reset_forward()
          it_time_frame = 1 
      else 
        
          ! Iteration index corresponding to the selected frame.
          it_time_frame = i_frame * step_frames + 1
      
          ! Load the global frame corresponding to index i_frame.
          pressure(:,:) = GLOB_FRAMES(:,:,1,i_frame)
          rhop(:,:)     = GLOB_FRAMES(:,:,2,i_frame)
          vx(:,:)       = GLOB_FRAMES(:,:,3,i_frame)
          vy(:,:)       = GLOB_FRAMES(:,:,4,i_frame)
      
          ! Restore PML memory variables associated with the selected frame.
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
      endif  
  
  else  ! The requested state corresponds to a local frame.
  
      ! Index of the selected frame
      i_frame = (modulo(it_time, step_frames)) / step_local_frames
      ! Time corresponding to the selected frame
      it_time_frame = (it_time / step_frames) * step_frames + i_frame * step_local_frames + 1
      
      ! Load the local frame corresponding to index i_frame.
      pressure(:,:) = LOC_FRAMES(:,:,1,i_frame)
      rhop(:,:)     = LOC_FRAMES(:,:,2,i_frame)
      vx(:,:)       = LOC_FRAMES(:,:,3,i_frame)
      vy(:,:)       = LOC_FRAMES(:,:,4,i_frame)
      
      ! Restore PML memory variables associated with the selected frame.
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
!==============================================================================
! Save global simulation frames.
!
! This routine performs successive forward simulations and stores the complete
! simulation state at predefined regular time intervals.
!
! Each saved frame contains:
!
!   - pressure field
!   - density perturbation field
!   - velocity components
!   - PML recursive convolution memory variables
!
! The stored frames are later used to reconstruct the wavefield evolution
! without recomputing the complete forward simulation.
!
! Starting from the initial state, step_frames_minus_1 forward iterations
! are performed between two consecutive saved frames.
!
!==============================================================================
 use parameters
 implicit none
 
 integer :: it_frame
 integer :: it_start, it_end
 integer :: step_frames_minus_1
 
! Compute the number of time steps between two (global) saved frames.
! The initial state is already available, therefore only the remaining
! time steps between two saved frames have to be computed.
 step_frames_minus_1 = NSTEP / N_GLOB_FRAMES - 1
 
 ! Initialization
 it_end = 0
 
 ! Perform N_GLOB_FRAMES forward simulations.  
 ! The initial state corresponds to the previously saved frame
 ! (or to the initial simulation state for the first frame).
 ! Starting from this state, step_frames_minus_1 + 1 forward
 ! iterations (i_start through i_end, inclusive) are performed
 ! to reach the next saved frame.
 do it_frame=1,N_GLOB_FRAMES
   ! Define the starting and ending time steps.
   it_start = it_end + 1
   it_end   = it_start + step_frames_minus_1
   
   ! Forward simulation
   call forwardproblem(p0_prior,rho0_prior,windx_prior,windy_prior, it_start, it_end,1) 
   
   !! Store the current state (and the memory variables for PML) as a new saved frame.
   GLOB_FRAMES(:,:,1,it_frame) = pressure(:,:)
   GLOB_FRAMES(:,:,2,it_frame) = rhop(:,:)
   GLOB_FRAMES(:,:,3,it_frame) = vx(:,:)
   GLOB_FRAMES(:,:,4,it_frame) = vy(:,:)

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
!==============================================================================
! Save local intermediate simulation frames.
!
! This routine generates additional intermediate frames between two consecutive
! global frames. The previous global frame is first restored and used as the
! initial state.
!
! Several shorter forward simulations are then performed to generate local
! frames with a finer temporal resolution.
!
! Each local frame contains:
!
!   - pressure field
!   - density perturbation field
!   - velocity components
!   - PML recursive convolution memory variables
!
! Local frames improve temporal resolution during reverse-time reconstruction
! while limiting the additional memory requirements.
!
!==============================================================================
 use parameters
 implicit none
 integer, intent(in) :: it_time 
 integer :: it_frame
 integer :: it_start, it_end
 integer :: step_local_frames_minus_1, step_frames
 
! Compute the number of time steps between two saved frames.
 step_frames = NSTEP / N_GLOB_FRAMES 
! The initial state is already available, therefore only the remaining
! time steps between two saved frames have to be computed.
 step_local_frames_minus_1 = step_frames / N_LOC_FRAMES - 1
 
! Load the previous saved frame as the initial state.
 call load_frame(it_time/step_frames*step_frames+1, it_end)
 it_end = it_end - 1
 
 ! The simulation starts from the previous global frame.
 ! Perform N_LOC_FRAMES forward simulations.
 ! The initial state corresponds to the previously saved frame.
 ! Starting from this state, step_local_frames_minus_1 + 1 forward iterations
 ! are performed (i.e., i_end - it_start + 1 iterations)
 ! to reach the next local saved frame.
 do it_frame=1,N_LOC_FRAMES
   ! Define the starting and ending time steps.
   it_start = it_end + 1
   it_end   = min(it_time, it_start + step_local_frames_minus_1)
   
   ! Forward simulation
   call forwardproblem(p0_prior,rho0_prior,windx_prior,windy_prior, it_start, it_end,1) 
   
   ! Store the current state and PML memory variables as a new saved frame.  
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

