module parameters

  use MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Paramaters 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! method
! 1: forward, 2: kernel, 3: inversion
 integer, parameter :: method = 2
 
 
! total number of grid points in each direction of the grid
  integer, parameter :: NY = 1507
  integer, parameter :: NX = 4608

  integer, parameter :: NPROC_X = 36 !! 20
  integer, parameter :: NPROC_Y = 11 !! 20
  integer, parameter :: NPROC = NPROC_X * NPROC_Y !! 20
  
  integer, parameter :: NX_LOCAL = NX / NPROC_X
  integer, parameter :: NY_LOCAL = NY / NPROC_Y 
  
  
! size of a grid cell
  double precision, parameter :: DELTAX = 100.d0
  double precision, parameter :: DELTAY = DELTAX

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
  integer, parameter :: NSTEP = 29000

! parameters for the source
  double precision, parameter :: f0 = 0.1d0
  double precision, parameter :: t0 = 1.20d0 / f0
  double precision, parameter :: factor = 1.0d0
  double precision, parameter :: PI = 3.141592653589793238462643d0
  double precision, parameter :: a = pi*pi*f0*f0

! source (in pressure, thus at a gridpoint rather than half a grid cell away)
  integer, parameter :: type_source = 2! 1. Plane wave, 2. Point source, 3. POint source with SPREAD_SSF
  ! if type_source == 1 
  integer, parameter :: wavefront = 2 ! 1. Wavefront in x direction, 2. Wavefront in y direction
  ! if type_source == 1,2 or 3
  double precision, parameter :: xsource = 100000.d0
  double precision, parameter :: ysource = 100.d0
  integer, parameter :: ISOURCE = xsource / DELTAX + 1
  integer, parameter :: JSOURCE = ysource / DELTAY + 1
  ! if type_source == 3 
  !! Standard deviation to spread the source
  !! SSF_Sigma should be around 1 wave length
  !! Attention : the source must have a greater edge distance than SSF_Sigma
  double precision, parameter :: SSF_Sigma = 500.d0
  ! spread the source spatial function
  double precision :: distance2, factor_ssf
 
 
 ! read from file or use a simulation
 ! 0: observation are created by a resolution of a simulation
 ! 1: observation are read from file
 integer :: observation_from_file = 1
 ! If observation_from_file == 1, read from a file of the form :
 character(len=200) :: path_obs_file = &
     "./Data/real_pressure_file_obs_"
        

 integer, parameter :: NREC_SET = 1
 integer, dimension(NREC_SET), parameter :: NREC_PER_SET = (/6/) 
 double precision, dimension(NREC_SET,4), parameter :: REC_SET_INFO = transpose(reshape( &
      (/420000,	100,	430000,	100 /), (/4,NREC_SET/)))

integer, dimension(100) :: list_irec_alt
integer :: dim_list_irec        

integer, parameter :: NREC = sum(NREC_PER_SET)



  logical, parameter :: atmospheric_model_file = .true.
   character(len=200) :: atmospheric_file_name_true = &
     "./atmospheric_model_Hukkakero_grav.dat"
   character(len=200) :: atmospheric_file_name_prior = & 
     "./atmospheric_model_Hukkakero.dat"

  ! P-velocity and density
  ! the unrelaxed value is the value at frequency = 0 (the relaxed value would be the value at frequency = +infinity)
    double precision, parameter :: gamma_chimie_value = 1.4d0
    double precision, parameter :: cp_unrelaxed_true = 347.763977d0 
  double precision, parameter :: density_true = 1.13837624d0 
  double precision, parameter :: windx_value_true = 0.0d0
  double precision, parameter :: cp_unrelaxed_prior = 347.763977d0 
  double precision, parameter :: density_prior = 1.13837624d0 
  double precision, parameter :: windx_value_prior = 0.0d0
 
  
  ! add perturbation into the model
  !! perturbation correspond to a square or a circle, in which one of the parameter (c,rho) are multiplied by a factor
  !! define the number of perturbation
  integer, parameter :: NPERTURB_MODEL = 0
  !! define the type of perturbation
  !!! If perturbation has the shape of a square
  !!!!   Type=1, x_start, y_start, x_end, y_end, factor on density, factor on celerity
  !!! If perturbation has the shape of a circle
  !!!!   Type=2, x_center, y_center, circle radius, None, factor on density, factor on celerity 
  double precision, dimension(NPERTURB_MODEL,7), parameter :: ADD_PERTURB_MODEL_INFO = &
                transpose(reshape( & ! (/ 1.d0, 12000.d0, 2000.d0, 15000.d0, 3000.d0,1.2d0,1.1d0, &
                   (/ 1.d0,  -10000.d0,  13000.d0,  61000.d0, 28000.0d0,1.d0,1.0d0/), (/7,NPERTURB_MODEL/))) 
                 !   (/ 2.d0,  7500.d0,  5000.d0,  1500.d0, 0.0d0,1.3d0,1.1d0/), (/7,NPERTURB_MODEL/))) 
                 
  
  logical, parameter :: add_windperturb_profile = .false.
  ! wind can be modeled by a gaussian
  ! expression if the model needs 3 parameters : the mean of the gaussian, the variance and the amplitude of the wind
  ! we add also two parameters to cancel the wind at the extremities
  integer, parameter :: ymin_wind = 8000 ! m
  integer, parameter :: ymax_wind = 28000 ! m
  integer, parameter :: jmin_wind = ymin_wind / DELTAY + 1
  integer, parameter :: jmax_wind = ymax_wind / DELTAY + 1
  double precision, parameter :: mean_gauss_wind = 18  ! km
  double precision, parameter :: sigma2_gauss_wind = 30 ! km^2
  double precision, parameter :: max_wind_factor = 10 ! meters/s 
   
  
 
! display information on the screen from time to time
  integer, parameter :: IT_DISPLAY = 1000 !NSTEP !  200! NSTEP !200
  ! save_normimage_overtime is a parameter to define the normalisation of the snapshots (only possible with method = 1)
  !! 1 : if want to save the image normalised over the maxmimum amplitude of the wave over all the propagation 
  !! 0 : if want to save the image normalised over the maximum amplitude of the current state 
  integer :: save_normimage_overtime = 0 
  !! 4 parameters to save the maximum amplitude of the waveform over all the propagation (if method = 1 and save_normimage_overtime = 1)
  double precision :: maxval_image_p = -1.0d0
  double precision :: maxval_image_rho = -1.0d0
  double precision :: maxval_image_vx = -1.0d0
  double precision :: maxval_image_vy = -1.0d0

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
  double precision, parameter :: ONE_OVER_DELTAX = 1.d0 / (DELTAX)
  double precision, parameter :: ONE_OVER_DELTAY = 1.d0 / (DELTAY)

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
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: & 
    pressure,rhop,vx,vy,                        &
    rhoa,pa,vax,vay,                            &
    kappa_unrelaxed_prior,rho0_prior,p0_prior,windx_prior,windy_prior,c0_prior, &
    kappa_unrelaxed_true,rho0_true, p0_true, windx_true, windy_true, &
    gamma_chimie,g                               


  ! for the source
  double precision :: t,t_demi,source_term, source_term_demi
  
  ! for receivers
  double precision xspacerec,yspacerec,distval,dist
  integer, dimension(NREC) :: ix_rec,iy_rec
  double precision, dimension(NREC) :: xrec,yrec
  integer :: myNREC
  logical :: save_sismos ! to save or not the seismograms
  
! for seismograms
  double precision, dimension(NSTEP,NREC) :: sisvx,sisvy,sispressure,sisrhop

  double precision, dimension(NSTEP,NREC) :: sispressure_true, sispressure_prior
  double precision :: normsq_pressure_true
  double precision, dimension(NREC) :: normsq_pressure_true_per_rec
  double precision :: regul_term_rho0_prior,regul_term_p0_prior,regul_term_windx_prior
  double precision :: normsq_rho0_prior, normsq_p0_prior, normsq_windx_prior, normsq_windy_prior

  ! Kernels shows the sensitivity to the model of :
  !! 0 : Full waveform
  !! 1 : Arrival time 
  integer, parameter :: observation = 1
  double precision, dimension(NSTEP) :: wr
  double precision, dimension(NREC,2), parameter :: REC_wr = transpose(reshape( &
      (/1315, 25, & !1280, 20, &
        1320, 25, & !1284, 20, &
        1325, 25, & !1288, 20, &
        1330, 25, & !1292, 20, &
        1335, 25, & !1296, 20, &
        1340, 25 /), (/2,NREC/)))  !1300, 20 /), (/2,NREC/)))
  
  !integer :: i_tmin = tmin / DELTAT
  !integer :: i_delta_tmin = delta_tmin / DELTAT
  double precision :: timeshift
  double precision, dimension(NSTEP,NREC) :: adjoint_source

  double precision, dimension(-1:NX_LOCAL+2, -1:NY_LOCAL+2) :: K_rho0, K_windx, K_windy, K_p0

  ! checkpointing
  integer, parameter :: NFRAMES = 100
  integer, parameter :: N_LOC_FRAMES = 8
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2,1:4,1:NFRAMES) :: FRAMES
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2,1:4,1:N_LOC_FRAMES) :: LOC_FRAMES

  ! PML parameters
  ! flags to add PML layers to the edges of the grid
  logical, parameter :: USE_PML_XMIN = .TRUE.
  logical, parameter :: USE_PML_XMAX = .TRUE.
  logical, parameter :: USE_PML_YMIN = .FALSE.
  logical, parameter :: USE_PML_YMAX = .TRUE.
  ! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 20
  ! power to compute d0 profile
  double precision, parameter :: NPOWER = 2.d0
  ! from Stephen Gedney's unpublished class notes for class EE699, lecture 8, slide 8-11
  double precision, parameter :: K_MAX_PML = 1.d0
  double precision, parameter :: ALPHA_MAX_PML = 2.d0*PI*(f0/2.d0) ! from Festa and Vilotte
  ! 1D arrays for the damping profiles
  double precision, dimension(1:NX) :: d_x,K_x,alpha_x,a_x,b_x,d_x_half,K_x_half,alpha_x_half,a_x_half,b_x_half,c_x,c_x_half, &
                                     one_over_K_x,one_over_K_x_half,one_over_Kdalpha_x,one_over_Kdalpha_x_half
  double precision, dimension(1:NY) :: d_y,K_y,alpha_y,a_y,b_y,d_y_half,K_y_half,alpha_y_half,a_y_half,b_y_half,c_y,c_y_half, &
                                     one_over_K_y,one_over_K_y_half,one_over_Kdalpha_y,one_over_Kdalpha_y_half

  double precision :: thickness_PML_x,thickness_PML_y,xoriginleft,xoriginright,yoriginbottom,yorigintop
  double precision :: Rcoef,d0_x,d0_y,xval,yval,abscissa_in_PML,abscissa_normalized
  
  
  ! PML memory variables
    double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) ::        &
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
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) ::                         &
      eq1_memory_drhoa_dx_adj,                                  &
      eq1_memory_dp0_dx_adj,                                    &
      eq1_memory_dgammap0pa_dx_adj, eq1_memory_dgammap0pa_dy_adj,         &
      eq1_memory_dwindx_dx_adj, eq1_memory_dwindy_dx_adj,       &
      eq1_memory_dvax_dx_adj, eq1_memory_dvax_dy_adj,           &
      eq1_memory_drhoa_dy_adj,                                  &
      eq1_memory_dp0_dy_adj,                                    &
      eq1_memory_dwindy_dy_adj, eq1_memory_dwindx_dy_adj,       &
      eq1_memory_dvay_dx_adj, eq1_memory_dvay_dy_adj,           &
      !
      eq2_memory_dpa_dx_adj,                               &
      eq2_memory_drhoa_dx_adj,                             &
      eq2_memory_dpa_dy_adj,                               &
      eq2_memory_drhoa_dy_adj,                             &
      eq2_memory_dwindx_dx_adj, eq2_memory_dwindy_dy_adj,       &
      eq2_memory_dvax_dx_adj, eq2_memory_dvay_dy_adj,           &
      eq2_memory_dwindy_dx_adj, eq2_memory_dwindx_dy_adj
  ! PML variables for checkpointing
      double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2,NFRAMES) ::                &
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

      double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2,N_LOC_FRAMES) ::                &
  ! for equation on rhop and pressure
      eq1_memory_dp0_dx_fw_loc_fr, eq1_memory_dp0_dy_fw_loc_fr,            &
      eq1_memory_drho0_dx_fw_loc_fr, eq1_memory_drho0_dy_fw_loc_fr,        &
      eq1_memory_dpressure_dx_fw_loc_fr, eq1_memory_dpressure_dy_fw_loc_fr,&
      eq1_memory_drhop_dx_fw_loc_fr, eq1_memory_drhop_dy_fw_loc_fr,        &
      eq1_memory_dvx_dx_fw_loc_fr, eq1_memory_dvy_dy_fw_loc_fr,            &
      eq1_memory_dwindx_dx_fw_loc_fr, eq1_memory_dwindy_dy_fw_loc_fr,      &
  ! for equation on vx
      eq2_memory_dpressure_dx_fw_loc_fr,                                   &
      eq2_memory_drho0_dx_fw_loc_fr, eq2_memory_drho0_dy_fw_loc_fr,        &
      eq2_memory_dvx_dx_fw_loc_fr, eq2_memory_dvx_dy_fw_loc_fr,            &
      eq2_memory_dwindx_dx_fw_loc_fr, eq2_memory_dwindx_dy_fw_loc_fr,      &
      eq2_memory_dwindy_dy_fw_loc_fr,                                      &
  ! for equation on vy
      eq3_memory_dpressure_dy_fw_loc_fr,                                   &
      eq3_memory_drho0_dy_fw_loc_fr, eq3_memory_drho0_dx_fw_loc_fr,        &
      eq3_memory_dvy_dy_fw_loc_fr, eq3_memory_dvy_dx_fw_loc_fr,            &
      eq3_memory_dwindy_dy_fw_loc_fr, eq3_memory_dwindy_dx_fw_loc_fr,      &
      eq3_memory_dwindx_dx_fw_loc_fr

 ! MPI variables
 ! array needed for MPI_RECV
  integer, dimension(MPI_STATUS_SIZE) :: message_status

 ! tag of the message to send
  integer, parameter :: message_tag = 0

 ! number of values to send or receive
  integer, parameter :: number_of_values_x = 2*(NY_LOCAL+4)
  integer, parameter :: number_of_values_y = 2*(NX_LOCAL+4)
  integer, parameter :: number_of_values_corner = 4

  integer :: row_Comm,ierr
  integer :: nb_procs,rank,i_rank,j_rank,code,rank_cut_plane,i_global,offset_i, &
                                               j_global, offset_j
  integer :: sender_right_shift,receiver_right_shift,sender_left_shift,receiver_left_shift,&
  sender_bottom_shift,receiver_bottom_shift,sender_top_shift,receiver_top_shift,&
             sender_right_top_shift, receiver_right_top_shift, sender_right_bottom_shift, receiver_right_bottom_shift,&
             sender_bottom_right_shift, receiver_bottom_right_shift, sender_bottom_left_shift, receiver_bottom_left_shift,&
             sender_top_left_shift, receiver_top_left_shift
  
  ! timer to count elapsed time
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  integer ihours,iminutes,iseconds,int_tCPU
  double precision :: time_start,time_end,tCPU


 integer, parameter :: Nflat = 3*NY_LOCAL 

 ! parameters of the inversion
 !! 1 density, celerity (wave speed), windx
 !! 2 density, pressure, windx
 !! 3 log density, log celerity, windx
 !! 35 log density, wave speed, windx
 !! 36 log density/density_prior, wave speed, windx
 !! 4 log density, log pressure, windx
 !! 5 log celerity, log pressure, windx
 integer, parameter :: parametrisation = 36

 ! contains the scaling to have x1, x2, and x3 of the inversion varying in the same way
 ! Depend on the chosen parameterization :
  ! if parametrisation == 1 then x1: density, x2: celerity, x3: windx (choose 1,100,1)
  ! if parametrisation == 2 then x1: density, x2: pressure, x3: windx (choose 1,1e5,1)
  ! if parametrisation == 3 then x1: log density, x2: log celerity, x3: windx (choose 1,1,1)
  ! if parametrisation == 4 then x1: log density, x2: log pressure, x3: windx (choose 1,1,1)
  ! if parametrisation == 5 then x1: log celerity, x2: log pressure, x3: windx (choose 1,1,1)
 ! Ref : Nocedal, (2006) Numerical Optimisation. 
 ! Scaling is defined in Scaling, page 26 (chapitre 2. Fundamentals of unconstrained optimization)
 double precision, dimension(3), parameter :: scale_model = (/1.0d0,100.0d0,100.0d0/)
 
 ! number of iterations to start with steepest descent direction
 integer :: steepest_nbiter_default = 5

 ! gradient
 !! 1 Fletcher Reeves
 !! 2 Polak Ribieres
 !! 3 Perry Shanno
 !! 4 Hager Zhang
 !! 5 Dai Kou
 !! 6 L-BFGS
 integer, parameter :: type_gradient = 6
 !! parameter for L-Bfgs
 integer,parameter :: mem_lbfgs = 5
 double precision, dimension(mem_lbfgs) :: RHO_list
 double precision, dimension(mem_lbfgs,1:Nflat) :: S_list, Y_list
 
 ! add a term of regularisation of the cost function of the inverse problem
 integer, parameter :: type_regul_term = 0 ! 0. None
                                           ! 1. Norm of current model - a priori model, 
                                           ! 2. (Not implemented - gradient),
                                           ! 3. (Not implemented - laplacian)  
 double precision, parameter :: regul_weight = 0.0000010d0
 double precision, dimension(1:Nflat) :: factor_regul_SRdist
 
 double precision :: alpha_start, alpha, alpha_prec, alpha_low, alpha_high
 
 integer :: count_grad, count_restart, count_f
 double precision :: reg_weight
  
 double precision, parameter :: c1 = 0.0001
 double precision, parameter :: c2 = 0.9
  
 ! define the factor of decreasing for the backtracking algorithm
 double precision, parameter :: rate = 0.8d0
 ! define the maximum number of iterations in the backtracking and linesearch/zoom algorithm
 integer, parameter :: maxiter_innerloop = 50

 ! define the maximum number of iterations in the main optimisation loop
 integer, parameter :: maxiter_outerloop = 100
 ! tolerance on the model x to stop the optimisation algorithm
 double precision :: tol_x = 1e-10

 ! alpha_max is maximum alpha that can be used in the line search procedure
 double precision :: alpha_max = 10
 
 double precision, dimension(1:Nflat) :: m0     ! a priori model
 double precision, dimension(1:Nflat) :: m1     ! a priori model
 
 double precision, dimension(1:Nflat) :: x,dfx
 double precision :: fx

 double precision, dimension(1:Nflat) :: x_old,dfx_old
 double precision :: fx_old
 double precision, dimension(1:Nflat) :: x_new,dfx_new
 double precision :: fx_new

 double precision, dimension(1:Nflat) :: x_low,dfx_low
 double precision :: fx_low
 double precision, dimension(1:Nflat) :: x_high,dfx_high
 double precision :: fx_high

 double precision, dimension(1:Nflat) :: r, r_old

 double precision :: dfx_r
 
 ! after updating model with the descent direction, can smooth the model
 !! 1 : mean filter (mean filter with a window of 5 elements)
 !! 2 : gaussian filter (gaussian filter with a window of 5 elements)
 !! 3 : median filter (median filter with a window of 9 elements)
 integer, parameter :: type_smoothing = 3

 double precision :: fx_data, fx_regul
end module parameters 
