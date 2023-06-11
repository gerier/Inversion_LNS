module parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Paramaters 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! total number of grid points in each direction of the grid
  integer, parameter :: NY = 201
  integer, parameter :: NX = 201

! size of a grid cell
  double precision, parameter :: DELTAX = 100.d0
  double precision, parameter :: DELTAY = DELTAX

! P-velocity and density
! the unrelaxed value is the value at frequency = 0 (the relaxed value would be the value at frequency = +infinity)
  double precision, parameter :: cp_unrelaxed = 347.763977d0 
  double precision, parameter :: density = 1.13837624d0 
  double precision, parameter :: gamma_chimie = 1.4d0

! Time step in seconds.
! The CFL stability number for the O(2,2) algorithm is 1 / sqrt(2) = 0.707
! i.e. one must choose  cp * deltat / deltax < 0.707.
! For the O(2,4) algorithm used here it is a bit more restrictive,
! it is cp * deltat / deltax < 0.606  (see Levander 1988 eq (7)).
! However this only ensures that the scheme is stable. To have a scheme that is both stable and accurate,
! for O(2,4) some numerical tests show that one needs to take about half of that,
! i.e. choose deltat so that cp * deltat / deltax is equal to about 0.30 or so. (or any value below; but not above).
! Since the time scheme is only second order, this also depends on how many time steps are performed in total
! (i.e. what the value of NSTEP below is); for large values of NSTEP, of course numerical errors will start to accumulate.
  double precision, parameter :: DELTAT = 0.05d0

! total number of time steps
  integer, parameter :: NSTEP = 1200

! parameters for the source
  double precision, parameter :: f0 = 0.1d0
  double precision, parameter :: t0 = 1.20d0 / f0
  double precision, parameter :: factor = 1.d0

! source (in pressure, thus at a gridpoint rather than half a grid cell away)
  double precision, parameter :: xsource = 4000.d0
  double precision, parameter :: ysource = 10000.d0
  integer, parameter :: ISOURCE = xsource / DELTAX + 1
  integer, parameter :: JSOURCE = ysource / DELTAY + 1
  double precision, parameter :: SSF_Sigma = 500.d0
  ! spread the source spatial function
  double precision :: distance2, factor_ssf
 
  double precision, parameter :: obstacle_xstart = -1000.d0
  double precision, parameter :: obstacle_xend   = 20000.d0
  double precision, parameter :: obstacle_ystart = -1000.d0
  double precision, parameter :: obstacle_yend   = 21000.d0
  integer, parameter :: IObs_start = obstacle_xstart / DELTAX + 1
  integer, parameter :: IObs_end   = obstacle_xend   / DELTAX + 1
  integer, parameter :: JObs_start = obstacle_ystart / DELTAY + 1
  integer, parameter :: JObs_end   = obstacle_yend   / DELTAY + 1
  double precision, parameter :: obstacle_factor_rho = 1.0d0
  double precision, parameter :: obstacle_factor_c2 = 1.1d0
  logical, parameter :: add_wind_profile = .False.
  
! receivers
  integer, parameter :: NREC = 1 !201
!! DK DK I use 2301 here instead of 2300 in order to fall exactly on a grid point
  double precision, parameter :: xdeb = 12500.d0   ! first receiver x in meters
  double precision, parameter :: ydeb = 10000.d0   ! first receiver y in meters
  double precision, parameter :: xfin = 12500.d0   ! last receiver x in meters
  double precision, parameter :: yfin = 10000.d0   ! last receiver y in meters

! display information on the screen from time to time
  integer, parameter :: IT_DISPLAY = 200

! compute some constants once and for all for the fourth-order spatial scheme
! These coefficients are given for instance by Levander, Geophysics, vol. 53(11), p. 1436, equation (A-2)
  double precision, parameter :: NINE_OVER_8_DELTAX = 9.d0 / (8.d0*DELTAX)
  double precision, parameter :: NINE_OVER_8_DELTAY = 9.d0 / (8.d0*DELTAY)
  double precision, parameter :: ONE_OVER_24_DELTAX = 1.d0 / (24.d0*DELTAX)
  double precision, parameter :: ONE_OVER_24_DELTAY = 1.d0 / (24.d0*DELTAY)
  ! compute some constants once and for all for the second-order spatial scheme
  double precision, parameter :: ONE_OVER_SIX_DELTAX = 1.d0 / (6.d0*DELTAX)
  double precision, parameter :: ONE_OVER_SIX_DELTAY = 1.d0 / (6.d0*DELTAY)
  ! Some constants for the derivation in time
  double precision, parameter :: ONE_OVER_DELTAT = 1.d0 / (DELTAT)

! value of PI
  double precision, parameter :: PI = 3.141592653589793238462643d0
! zero
  double precision, parameter :: ZERO = 0.d0
! large value for maximum and minimum
  double precision, parameter :: HUGEVAL = 1.d+30
  double precision :: TINYVAL = 1e-16
! threshold above which we consider that the code became unstable
  double precision, parameter :: STABILITY_THRESHOLD = 1.d+25

! main arrays
! in order to be able to use a fourth-order spatial operator on the edges of the model
! here we define the arrays with size (0:NX+1,0:NY+1) instead of size (NX,NY) as in the second-order case
  double precision, dimension(0:NX+1,0:NY+1) :: & 
    pressure,rhop,vx,vy,                        &
    rhoa,pa,vax,vay,                            &
    kappa_unrelaxed_prior,rho0_prior,p0_prior,windx_prior,windy_prior,&
    kappa_unrelaxed_true,rho0_true, p0_true, windx_true, windy_true                               


  ! for the source
  double precision :: a,t,source_term
  
  ! for receivers
  double precision xspacerec,yspacerec,distval,dist
  integer, dimension(NREC) :: ix_rec,iy_rec
  double precision, dimension(NREC) :: xrec,yrec
  integer :: myNREC
  logical :: save_sismos ! to save or not the seismograms
  
! for seismograms
  double precision, dimension(NSTEP,NREC) :: sisvx,sisvy,sispressure,sisrhop

  double precision, dimension(NSTEP,NREC) :: sispressure_true, sispressure_prior
  double precision norm_pressure_true

  ! choose the desired kernel 
  integer ::  which_kernel = 1
  double precision, dimension(0:NX+1, 0:NY+1) :: K

  ! checkpointing
  integer, parameter :: NFRAMES = 50
  double precision, dimension(0:NX+1,0:NY+1,1:4,1:NFRAMES) :: FRAMES
  
  ! PML parameters
  ! flags to add PML layers to the edges of the grid
  logical, parameter :: USE_PML_XMIN = .true.
  logical, parameter :: USE_PML_XMAX = .true.
  logical, parameter :: USE_PML_YMIN = .true.
  logical, parameter :: USE_PML_YMAX = .true.
  ! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 20
  ! power to compute d0 profile
  double precision, parameter :: NPOWER = 2.d0
  ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
  double precision, parameter :: K_MAX_PML = 1.d0
  double precision, parameter :: ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0) ! from Festa and Vilotte
  ! 1D arrays for the damping profiles
  double precision, dimension(NX) :: d_x,K_x,alpha_x,a_x,b_x,d_x_half,K_x_half,alpha_x_half,a_x_half,b_x_half,c_x,c_x_half, &
                                     one_over_K_x,one_over_K_x_half
  double precision, dimension(NY) :: d_y,K_y,alpha_y,a_y,b_y,d_y_half,K_y_half,alpha_y_half,a_y_half,b_y_half,c_y,c_y_half, &
                                     one_over_K_y,one_over_K_y_half

  double precision :: thickness_PML_x,thickness_PML_y,xoriginleft,xoriginright,yoriginbottom,yorigintop
  double precision :: Rcoef,d0_x,d0_y,xval,yval,abscissa_in_PML,abscissa_normalized
  
  ! PML memory variables
    double precision, dimension(NX,NY) ::        &
  ! for equation on rhop and pressure
      eq1_memory_dp0_dx_fw, eq1_memory_dp0_dy_fw,            &
      eq1_memory_drho0_dx_fw, eq1_memory_drho0_dy_fw,        &
      eq1_memory_dpressure_dx_fw, eq1_memory_dpressure_dy_fw,&
      eq1_memory_drhop_dx_fw, eq1_memory_drhop_dy_fw,        &
      eq1_memory_dvx_dx_fw, eq1_memory_dvy_dy_fw,            &
      eq1_memory_dwindx_dx_fw, eq1_memory_dwindy_dy_fw,      &
  ! for equation on vx
      eq2_memory_dpressure_dx_fw,                            &
      eq2_memory_drho0_dx_fw, eq2_memory_drho0_dy_fw,        &
      eq2_memory_dvx_dx_fw, eq2_memory_dvx_dy_fw,            &
      eq2_memory_dwindx_dx_fw, eq2_memory_dwindx_dy_fw,      &
      eq2_memory_dwindy_dy_fw,                               &
  ! for equation on vy
      eq3_memory_dpressure_dy_fw,                            &
      eq3_memory_drho0_dy_fw, eq3_memory_drho0_dx_fw,        &
      eq3_memory_dvy_dy_fw, eq3_memory_dvy_dx_fw,            &
      eq3_memory_dwindy_dy_fw, eq3_memory_dwindy_dx_fw,      &
      eq3_memory_dwindx_dx_fw
  double precision, dimension(NX,NY) ::                         &
      eq1_memory_drho0_dx_adj,                                  &
      eq1_memory_dp0_dx_adj,                                    &
      eq1_memory_drhoarho0_dx_adj, eq1_memory_drhoarho0_dy_adj, &
      eq1_memory_dp0pa_dx_adj, eq1_memory_dp0pa_dy_adj,         &
      eq1_memory_dwindx_dx_adj, eq1_memory_dwindy_dx_adj,       &
      eq1_memory_dvax_dx_adj, eq1_memory_dvax_dy_adj,           &
      eq1_memory_drho0_dy_adj,                                  &
      eq1_memory_dp0_dy_adj,                                    &
      eq1_memory_dwindy_dy_adj, eq1_memory_dwindx_dy_adj,       &
      eq1_memory_dvay_dx_adj, eq1_memory_dvay_dy_adj,           &
      !
      eq2_memory_dpawindx_dx_adj,                               &
      eq2_memory_drhoawindx_dx_adj,                             &
      eq2_memory_dpawindy_dy_adj,                               &
      eq2_memory_drhoawindy_dy_adj,                             &
      eq2_memory_dwindx_dx_adj, eq2_memory_dwindy_dy_adj,       &
      eq2_memory_dvax_dx_adj, eq2_memory_dvay_dy_adj,           &
      eq2_memory_dwindy_dx_adj, eq2_memory_dwindx_dy_adj
  ! PML variables for checkpointing
      double precision, dimension(NX,NY,NFRAMES) ::                &
  ! for equation on rhop and pressure
      eq1_memory_dp0_dx_fw_fr, eq1_memory_dp0_dy_fw_fr,            &
      eq1_memory_drho0_dx_fw_fr, eq1_memory_drho0_dy_fw_fr,        &
      eq1_memory_dpressure_dx_fw_fr, eq1_memory_dpressure_dy_fw_fr,&
      eq1_memory_drhop_dx_fw_fr, eq1_memory_drhop_dy_fw_fr,        &
      eq1_memory_dvx_dx_fw_fr, eq1_memory_dvy_dy_fw_fr,            &
      eq1_memory_dwindx_dx_fw_fr, eq1_memory_dwindy_dy_fw_fr,      &
  ! for equation on vx
      eq2_memory_dpressure_dx_fw_fr,                               &
      eq2_memory_drho0_dx_fw_fr, eq2_memory_drho0_dy_fw_fr,        &
      eq2_memory_dvx_dx_fw_fr, eq2_memory_dvx_dy_fw_fr,            &
      eq2_memory_dwindx_dx_fw_fr, eq2_memory_dwindx_dy_fw_fr,      &
      eq2_memory_dwindy_dy_fw_fr,                                  &
  ! for equation on vy
      eq3_memory_dpressure_dy_fw_fr,                               &
      eq3_memory_drho0_dy_fw_fr, eq3_memory_drho0_dx_fw_fr,        &
      eq3_memory_dvy_dy_fw_fr, eq3_memory_dvy_dx_fw_fr,            &
      eq3_memory_dwindy_dy_fw_fr, eq3_memory_dwindy_dx_fw_fr,      &
      eq3_memory_dwindx_dx_fw_fr

end module parameters 
