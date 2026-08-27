
!==========================================================
! Initialization routines for forward and adjoint simulations
!==========================================================

subroutine reset_forward()
!==========================================================
! Reset forward simulation variables
!
! Initialize all forward wavefield variables, seismograms
! and PML memory variables to zero before starting a new
! forward simulation.
!==========================================================
 use parameters
 implicit none

  ! Seismograms
  sisvx(:,:)       = ZERO
  sisvy(:,:)       = ZERO
  sispressure(:,:) = ZERO
  sisrhop(:,:)     = ZERO
  
  ! Forward problem variables
  pressure(:,:)    = ZERO
  rhop(:,:)        = ZERO
  vx(:,:)          = ZERO
  vy(:,:)          = ZERO
  
  ! Forward PML memory variables
  ! PML memory variables for the p and rho equations
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
  ! PML memory variables for the vx equation
  eq2_memory_dpressure_dx_fw(:,:) = ZERO
  eq2_memory_drho0_dx_fw(:,:)     = ZERO
  eq2_memory_drho0_dy_fw(:,:)     = ZERO
  eq2_memory_dvx_dx_fw(:,:)       = ZERO
  eq2_memory_dvx_dy_fw(:,:)       = ZERO
  eq2_memory_dwindx_dx_fw(:,:)    = ZERO
  eq2_memory_dwindx_dy_fw(:,:)    = ZERO
  eq2_memory_dwindy_dy_fw(:,:)    = ZERO
  ! PML memory variables for the vy equation
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
!==========================================================
! Reset kernel computation variables
!
! Initialize the forward and adjoint wavefields, PML memory
! variables and sensitivity kernels before computing
! adjoint-based sensitivity kernels or validation kernels.
!==========================================================
 use parameters
 implicit none

  ! Forward problem variables  
  call reset_forward()
  
  ! Adjoint problem variables
  pa(:,:)    = ZERO
  rhoa(:,:)  = ZERO
  vax(:,:)   = ZERO
  vay(:,:)   = ZERO

  ! PML memory variables for the adjoint problem
  eq1_memory_drhoa_dx_adj(:,:)     = ZERO
  eq1_memory_dp0_dx_adj(:,:)       = ZERO
  eq1_memory_dgammap0pa_dx_adj(:,:)= ZERO
  eq1_memory_dgammap0pa_dy_adj(:,:)= ZERO
  eq1_memory_dwindx_dx_adj(:,:)    = ZERO
  eq1_memory_dwindy_dx_adj(:,:)    = ZERO
  eq1_memory_dvax_dx_adj(:,:)      = ZERO
  eq1_memory_dvax_dy_adj(:,:)      = ZERO
  eq1_memory_drhoa_dy_adj(:,:)     = ZERO
  eq1_memory_dp0_dy_adj(:,:)       = ZERO
  eq1_memory_dwindy_dy_adj(:,:)    = ZERO
  eq1_memory_dwindx_dy_adj(:,:)    = ZERO
  eq1_memory_dvay_dx_adj(:,:)      = ZERO
  eq1_memory_dvay_dy_adj(:,:)      = ZERO

  eq2_memory_dpa_dx_adj(:,:)       = ZERO
  eq2_memory_drhoa_dx_adj(:,:)     = ZERO
  eq2_memory_dpa_dy_adj(:,:)       = ZERO
  eq2_memory_drhoa_dy_adj(:,:)     = ZERO
  eq2_memory_dwindx_dx_adj(:,:)    = ZERO
  eq2_memory_dwindy_dy_adj(:,:)    = ZERO
  eq2_memory_dvax_dx_adj(:,:)      = ZERO
  eq2_memory_dvay_dy_adj(:,:)      = ZERO
  eq2_memory_dwindy_dx_adj(:,:)    = ZERO
  eq2_memory_dwindx_dy_adj(:,:)    = ZERO
       
  ! Sensitivity kernels (adjoint method)
  K_p0(:,:)    = ZERO
  K_rho0(:,:)  = ZERO
  K_windx(:,:) = ZERO
  K_windy(:,:) = ZERO
  
  !  Sensitivity kernels (small perturbation method - validation) 
  K(:,:)       = ZERO
endsubroutine reset_kernel

