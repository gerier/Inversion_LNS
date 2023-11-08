subroutine reset_forward()

 use parameters
 implicit none

  ! sismometers
  sisvx(:,:) = ZERO
  sisvy(:,:) = ZERO
  sispressure(:,:) = ZERO
  sisrhop(:,:) = ZERO
  
  ! variable of the forward problem
  pressure(:,:)    = ZERO
  rhop(:,:)        = ZERO
  vx(:,:)          = ZERO
  vy(:,:)          = ZERO
  
  ! PML
  eq1_memory_dp0_dx_fw(:,:)       = ZERO
  eq1_memory_dp0_dy_fw(:,:)       = ZERO
  eq1_memory_drho0_dx_fw(:,:)     = ZERO
  eq1_memory_drho0_dy_fw(:,:)     = ZERO
  eq1_memory_dpressure_dx_fw(:,:) = ZERO
  eq1_memory_dpressure_dy_fw(:,:) = ZERO
  eq1_memory_drhop_dx_fw(:,:)     = ZERO
  eq1_memory_drhop_dy_fw(:,:)     = ZERO
  eq1_memory_dvx_dx_fw(:,:)       = ZERO
  eq1_memory_dvy_dy_fw(:,:)       = ZERO
  eq1_memory_dwindx_dx_fw(:,:)    = ZERO
  eq1_memory_dwindy_dy_fw(:,:)    = ZERO
  ! for equation on vx
  eq2_memory_dpressure_dx_fw(:,:) = ZERO
  eq2_memory_drho0_dx_fw(:,:)     = ZERO
  eq2_memory_drho0_dy_fw(:,:)     = ZERO
  eq2_memory_dvx_dx_fw(:,:)       = ZERO
  eq2_memory_dvx_dy_fw(:,:)       = ZERO
  eq2_memory_dwindx_dx_fw(:,:)    = ZERO
  eq2_memory_dwindx_dy_fw(:,:)    = ZERO
  eq2_memory_dwindy_dy_fw(:,:)    = ZERO
  ! for equation on vy
  eq3_memory_dpressure_dy_fw(:,:) = ZERO
  eq3_memory_drho0_dy_fw(:,:)     = ZERO
  eq3_memory_drho0_dx_fw(:,:)     = ZERO
  eq3_memory_dvy_dy_fw(:,:)       = ZERO
  eq3_memory_dvy_dx_fw(:,:)       = ZERO
  eq3_memory_dwindy_dy_fw(:,:)    = ZERO
  eq3_memory_dwindy_dx_fw(:,:)    = ZERO
  eq3_memory_dwindx_dx_fw(:,:)    = ZERO
  
endsubroutine reset_forward



