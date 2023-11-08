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


subroutine reset_kernel()
 
 use parameters
 implicit none

  
  call reset_forward()
  
  ! variables of the adjoint problem
  pa(:,:)    = ZERO
  rhoa(:,:)  = ZERO
  vax(:,:)   = ZERO
  vay(:,:)   = ZERO

  ! PML variables fo the adjoint problem
  eq1_memory_drho0_dx_adj(:,:)     = ZERO
  eq1_memory_dp0_dx_adj(:,:)       = ZERO
  eq1_memory_drhoarho0_dx_adj(:,:) = ZERO
  eq1_memory_dp0pa_dx_adj(:,:)     = ZERO
  eq1_memory_dp0pa_dy_adj(:,:)     = ZERO
  eq1_memory_dwindx_dx_adj(:,:)    = ZERO
  eq1_memory_dwindy_dx_adj(:,:)    = ZERO
  eq1_memory_dvax_dx_adj(:,:)      = ZERO
  eq1_memory_dvax_dy_adj(:,:)      = ZERO
  eq1_memory_drho0_dy_adj(:,:)     = ZERO
  eq1_memory_dp0_dy_adj(:,:)       = ZERO
  eq1_memory_drhoarho0_dy_adj(:,:) = ZERO
  eq1_memory_dp0pa_dy_adj(:,:)     = ZERO
  eq1_memory_dwindy_dy_adj(:,:)    = ZERO
  eq1_memory_dwindx_dy_adj(:,:)    = ZERO
  eq1_memory_dvay_dx_adj(:,:)      = ZERO
  eq1_memory_dvay_dy_adj(:,:)      = ZERO

  eq2_memory_dpawindx_dx_adj(:,:)  = ZERO
  eq2_memory_drhoawindx_dx_adj(:,:)= ZERO
  eq2_memory_dpawindy_dy_adj(:,:)  = ZERO
  eq2_memory_drhoawindy_dy_adj(:,:)= ZERO
  eq2_memory_dwindx_dx_adj(:,:)    = ZERO
  eq2_memory_dwindy_dy_adj(:,:)    = ZERO
  eq2_memory_dvax_dx_adj(:,:)      = ZERO
  eq2_memory_dvay_dy_adj(:,:)      = ZERO
  eq2_memory_dwindy_dx_adj(:,:)    = ZERO
  eq2_memory_dwindx_dy_adj(:,:)    = ZERO
       
  ! kernel variables  
  K_p0(:,:)    = ZERO
  K_rho0(:,:)  = ZERO
  K_windx(:,:)   = ZERO
  K_windy(:,:)   = ZERO
  
endsubroutine reset_kernel

