module parameters

  use MPI
  implicit none

  !==========================================================================
  ! SECTION 1 — Mathematical and physical constants
  !==========================================================================
  double precision, parameter :: PI                  = 3.141592653589793238462643d0
  double precision, parameter :: ZERO                = 0.d0
  double precision, parameter :: HUGEVAL             = 1.d+30
  double precision, parameter :: TINYVAL             = 1e-16
  double precision, parameter :: STABILITY_THRESHOLD = 1.d+25
  double precision, parameter :: NPOWER              = 2.d0
  double precision, parameter :: K_MAX_PML           = 1.d0
  integer,          parameter :: NPOINTS_PML         = 20
  integer,          parameter :: message_tag         = 0
  integer,          parameter :: number_of_values_corner = 4

  !==========================================================================
  ! SECTION 2 — Parameters read from the parameter file
  !==========================================================================

  ! --- Method ---
  integer :: method
  logical :: validation
  integer :: which_kernel

  ! --- MPI decomposition ---
  integer :: NPROC_X
  integer :: NPROC_Y

  ! --- Grid ---
  integer          :: NX
  integer          :: NY
  double precision :: DELTAX
  double precision :: DELTAY

  ! --- Time stepping ---
  integer          :: NSTEP
  double precision :: DELTAT

  ! --- PML ---
  logical :: USE_PML_XMIN
  logical :: USE_PML_XMAX
  logical :: USE_PML_YMIN
  logical :: USE_PML_YMAX

  ! --- Source ---
  double precision :: f0
  integer          :: type_source
  integer          :: wavefront
  double precision :: xsource
  double precision :: ysource
  double precision :: SSF_Sigma

  ! --- Observations ---
  integer            :: observation_from_file
  character(len=200) :: path_obs_file
  integer            :: observation
  integer            :: window_waveform

  ! --- Receivers ---
  integer :: NREC_SET
  integer,          allocatable :: NREC_PER_SET(:)       ! (NREC_SET)
  double precision, allocatable :: REC_SET_INFO(:,:)     ! (NREC_SET, 4)
  double precision, allocatable :: REC_wr(:,:)           ! (NREC, 2)

  ! --- Atmospheric model ---
  logical            :: atmospheric_model_file
  character(len=200) :: atmospheric_file_name_true
  character(len=200) :: atmospheric_file_name_prior

  ! --- Physical parameters  ---
  double precision :: gamma_chemestry_value
  double precision :: cp_unrelaxed_true
  double precision :: density_true
  double precision :: windx_value_true
  double precision :: cp_unrelaxed_prior
  double precision :: density_prior
  double precision :: windx_value_prior

  ! --- Model perturbations ---
  integer :: NPERTURB_MODEL
  double precision, allocatable :: ADD_PERTURB_MODEL_INFO(:,:)  ! (NPERTURB_MODEL, 7)

  ! --- Wind perturbation ---
  logical          :: add_windperturb_profile
  integer          :: ymin_wind
  integer          :: ymax_wind
  double precision :: mean_gauss_wind
  double precision :: sigma2_gauss_wind
  double precision :: max_wind_factor

  ! --- Output and saving ---
  integer :: IT_DISPLAY
  integer :: save_normimage_overtime
  logical :: save_adjoint_source

  ! --- Checkpointing ---
  integer :: N_GLOB_FRAMES
  integer :: N_LOC_FRAMES

  ! --- Inversion ---
  integer          :: parametrisation
  double precision :: scale_model(3)
  integer          :: type_gradient
  integer          :: mem_lbfgs
  integer          :: type_regul_term
  double precision :: regul_weight
  integer          :: type_smoothing

  ! --- Optimization control ---
  integer          :: steepest_nbiter_default
  integer          :: maxiter_innerloop
  integer          :: maxiter_outerloop
  double precision :: tol_x
  double precision :: alpha_max
  double precision :: rate
  double precision :: c1
  double precision :: c2

  !==========================================================================
  ! SECTION 3 —  Derived parameters computed by init_parameters()
  !==========================================================================
  integer :: NPROC
  integer :: NX_LOCAL, NY_LOCAL
  integer :: Nflat
  integer :: number_of_values_x, number_of_values_y
  integer :: NREC
  integer :: ISOURCE, JSOURCE
  integer :: jmin_wind, jmax_wind

  double precision :: t0, a
  double precision :: factor = 1.0d0

  double precision :: ALPHA_MAX_PML

  double precision :: NINE_OVER_8_DELTAX,  NINE_OVER_8_DELTAY
  double precision :: ONE_OVER_24_DELTAX,  ONE_OVER_24_DELTAY
  double precision :: ONE_OVER_SIX_DELTAX, ONE_OVER_SIX_DELTAY
  double precision :: ONE_OVER_DELTAT,  ONE_OVER_DELTAX,  ONE_OVER_DELTAY

  !==========================================================================
  ! SECTION 4 —  Runtime variables allocated in init_parameters()
  !==========================================================================

  double precision, allocatable :: pressure(:,:), rhop(:,:), vx(:,:), vy(:,:)
  double precision, allocatable :: rhoa(:,:), pa(:,:), vax(:,:), vay(:,:)
  double precision, allocatable :: kappa_unrelaxed_prior(:,:), rho0_prior(:,:)
  double precision, allocatable :: p0_prior(:,:), windx_prior(:,:), windy_prior(:,:), c0_prior(:,:)
  double precision, allocatable :: kappa_unrelaxed_true(:,:), rho0_true(:,:)
  double precision, allocatable :: p0_true(:,:), windx_true(:,:), windy_true(:,:)
  double precision, allocatable :: gamma_chemestry(:,:), g(:,:)

  double precision :: t, t_demi, source_term, source_term_demi
  double precision :: xspacerec, yspacerec, distval, dist
  integer,          allocatable :: ix_rec(:), iy_rec(:)
  double precision, allocatable :: xrec(:), yrec(:)
  integer  :: myNREC
  logical  :: save_sismos

  double precision, allocatable :: sisvx(:,:), sisvy(:,:)
  double precision, allocatable :: sispressure(:,:), sisrhop(:,:)
  double precision, allocatable :: sispressure_true(:,:)
  double precision, allocatable :: sispressure_prior(:,:)
  double precision, allocatable :: sispressure_source(:)
  double precision, allocatable :: adjoint_source(:,:)
  double precision, allocatable :: wr(:)

  double precision :: normsq_pressure_true
  double precision, allocatable :: normsq_pressure_true_per_rec(:)
  double precision :: regul_term_rho0_prior, regul_term_p0_prior, regul_term_windx_prior
  double precision :: normsq_rho0_prior, normsq_p0_prior, normsq_windx_prior, normsq_windy_prior
  double precision :: timeshift

  double precision, allocatable :: K_rho0(:,:), K_windx(:,:), K_windy(:,:), K_p0(:,:), K(:,:)

  double precision, allocatable :: GLOB_FRAMES(:,:,:,:)
  double precision, allocatable :: LOC_FRAMES(:,:,:,:)

  double precision, allocatable :: d_x(:), K_x(:), alpha_x(:), a_x(:), b_x(:), c_x(:)
  double precision, allocatable :: d_x_half(:), K_x_half(:), alpha_x_half(:)
  double precision, allocatable :: a_x_half(:), b_x_half(:), c_x_half(:)
  double precision, allocatable :: one_over_K_x(:), one_over_K_x_half(:)
  double precision, allocatable :: one_over_Kdalpha_x(:), one_over_Kdalpha_x_half(:)
  double precision, allocatable :: d_y(:), K_y(:), alpha_y(:), a_y(:), b_y(:), c_y(:)
  double precision, allocatable :: d_y_half(:), K_y_half(:), alpha_y_half(:)
  double precision, allocatable :: a_y_half(:), b_y_half(:), c_y_half(:)
  double precision, allocatable :: one_over_K_y(:), one_over_K_y_half(:)
  double precision, allocatable :: one_over_Kdalpha_y(:), one_over_Kdalpha_y_half(:)

  double precision :: thickness_PML_x, thickness_PML_y
  double precision :: xoriginleft, xoriginright, yoriginbottom, yorigintop
  double precision :: Rcoef, d0_x, d0_y, xval, yval
  double precision :: abscissa_in_PML, abscissa_normalized

  double precision, allocatable :: &
    eq1_memory_dp0_dx_fw(:,:),         eq1_memory_dp0_dy_fw(:,:),          &
    eq1_memory_drho0_dx_fw(:,:),       eq1_memory_drho0_dy_fw(:,:),        &
    eq1_memory_dpressure_dx_fw(:,:),   eq1_memory_dpressure_dy_fw(:,:),    &
    eq1_memory_drhop_dx_fw(:,:),       eq1_memory_drhop_dy_fw(:,:),        &
    eq1_memory_dvx_dx_fw(:,:),         eq1_memory_dvy_dy_fw(:,:),          &
    eq1_memory_dwindx_dx_fw(:,:),      eq1_memory_dwindy_dy_fw(:,:),       &
    eq2_memory_dpressure_dx_fw(:,:),                                        &
    eq2_memory_drho0_dx_fw(:,:),       eq2_memory_drho0_dy_fw(:,:),        &
    eq2_memory_dvx_dx_fw(:,:),         eq2_memory_dvx_dy_fw(:,:),          &
    eq2_memory_dwindx_dx_fw(:,:),      eq2_memory_dwindx_dy_fw(:,:),       &
    eq2_memory_dwindy_dy_fw(:,:),                                           &
    eq3_memory_dpressure_dy_fw(:,:),                                        &
    eq3_memory_drho0_dy_fw(:,:),       eq3_memory_drho0_dx_fw(:,:),        &
    eq3_memory_dvy_dy_fw(:,:),         eq3_memory_dvy_dx_fw(:,:),          &
    eq3_memory_dwindy_dy_fw(:,:),      eq3_memory_dwindy_dx_fw(:,:),       &
    eq3_memory_dwindx_dx_fw(:,:)

  double precision, allocatable :: &
    eq1_memory_drhoa_dx_adj(:,:),                                           &
    eq1_memory_dp0_dx_adj(:,:),                                             &
    eq1_memory_dgammap0pa_dx_adj(:,:), eq1_memory_dgammap0pa_dy_adj(:,:),  &
    eq1_memory_dwindx_dx_adj(:,:),     eq1_memory_dwindy_dx_adj(:,:),      &
    eq1_memory_dvax_dx_adj(:,:),       eq1_memory_dvax_dy_adj(:,:),        &
    eq1_memory_drhoa_dy_adj(:,:),                                           &
    eq1_memory_dp0_dy_adj(:,:),                                             &
    eq1_memory_dwindy_dy_adj(:,:),     eq1_memory_dwindx_dy_adj(:,:),      &
    eq1_memory_dvay_dx_adj(:,:),       eq1_memory_dvay_dy_adj(:,:),        &
    eq2_memory_dpa_dx_adj(:,:),        eq2_memory_drhoa_dx_adj(:,:),       &
    eq2_memory_dpa_dy_adj(:,:),        eq2_memory_drhoa_dy_adj(:,:),       &
    eq2_memory_dwindx_dx_adj(:,:),     eq2_memory_dwindy_dy_adj(:,:),      &
    eq2_memory_dvax_dx_adj(:,:),       eq2_memory_dvay_dy_adj(:,:),        &
    eq2_memory_dwindy_dx_adj(:,:),     eq2_memory_dwindx_dy_adj(:,:)

  double precision, allocatable :: &
    eq1_memory_dp0_dx_fw_fr(:,:,:),       eq1_memory_dp0_dy_fw_fr(:,:,:),          &
    eq1_memory_drho0_dx_fw_fr(:,:,:),     eq1_memory_drho0_dy_fw_fr(:,:,:),        &
    eq1_memory_dpressure_dx_fw_fr(:,:,:), eq1_memory_dpressure_dy_fw_fr(:,:,:),    &
    eq1_memory_drhop_dx_fw_fr(:,:,:),     eq1_memory_drhop_dy_fw_fr(:,:,:),        &
    eq1_memory_dvx_dx_fw_fr(:,:,:),       eq1_memory_dvy_dy_fw_fr(:,:,:),          &
    eq1_memory_dwindx_dx_fw_fr(:,:,:),    eq1_memory_dwindy_dy_fw_fr(:,:,:),       &
    eq2_memory_dpressure_dx_fw_fr(:,:,:),                                           &
    eq2_memory_drho0_dx_fw_fr(:,:,:),     eq2_memory_drho0_dy_fw_fr(:,:,:),        &
    eq2_memory_dvx_dx_fw_fr(:,:,:),       eq2_memory_dvx_dy_fw_fr(:,:,:),          &
    eq2_memory_dwindx_dx_fw_fr(:,:,:),    eq2_memory_dwindx_dy_fw_fr(:,:,:),       &
    eq2_memory_dwindy_dy_fw_fr(:,:,:),                                              &
    eq3_memory_dpressure_dy_fw_fr(:,:,:),                                           &
    eq3_memory_drho0_dy_fw_fr(:,:,:),     eq3_memory_drho0_dx_fw_fr(:,:,:),        &
    eq3_memory_dvy_dy_fw_fr(:,:,:),       eq3_memory_dvy_dx_fw_fr(:,:,:),          &
    eq3_memory_dwindy_dy_fw_fr(:,:,:),    eq3_memory_dwindy_dx_fw_fr(:,:,:),       &
    eq3_memory_dwindx_dx_fw_fr(:,:,:)

  double precision, allocatable :: &
    eq1_memory_dp0_dx_fw_loc_fr(:,:,:),       eq1_memory_dp0_dy_fw_loc_fr(:,:,:),          &
    eq1_memory_drho0_dx_fw_loc_fr(:,:,:),     eq1_memory_drho0_dy_fw_loc_fr(:,:,:),        &
    eq1_memory_dpressure_dx_fw_loc_fr(:,:,:), eq1_memory_dpressure_dy_fw_loc_fr(:,:,:),    &
    eq1_memory_drhop_dx_fw_loc_fr(:,:,:),     eq1_memory_drhop_dy_fw_loc_fr(:,:,:),        &
    eq1_memory_dvx_dx_fw_loc_fr(:,:,:),       eq1_memory_dvy_dy_fw_loc_fr(:,:,:),          &
    eq1_memory_dwindx_dx_fw_loc_fr(:,:,:),    eq1_memory_dwindy_dy_fw_loc_fr(:,:,:),       &
    eq2_memory_dpressure_dx_fw_loc_fr(:,:,:),                                               &
    eq2_memory_drho0_dx_fw_loc_fr(:,:,:),     eq2_memory_drho0_dy_fw_loc_fr(:,:,:),        &
    eq2_memory_dvx_dx_fw_loc_fr(:,:,:),       eq2_memory_dvx_dy_fw_loc_fr(:,:,:),          &
    eq2_memory_dwindx_dx_fw_loc_fr(:,:,:),    eq2_memory_dwindx_dy_fw_loc_fr(:,:,:),       &
    eq2_memory_dwindy_dy_fw_loc_fr(:,:,:),                                                  &
    eq3_memory_dpressure_dy_fw_loc_fr(:,:,:),                                               &
    eq3_memory_drho0_dy_fw_loc_fr(:,:,:),     eq3_memory_drho0_dx_fw_loc_fr(:,:,:),        &
    eq3_memory_dvy_dy_fw_loc_fr(:,:,:),       eq3_memory_dvy_dx_fw_loc_fr(:,:,:),          &
    eq3_memory_dwindy_dy_fw_loc_fr(:,:,:),    eq3_memory_dwindy_dx_fw_loc_fr(:,:,:),       &
    eq3_memory_dwindx_dx_fw_loc_fr(:,:,:)

  ! MPI
  integer, dimension(MPI_STATUS_SIZE) :: message_status
  integer :: row_Comm, ierr
  integer :: nb_procs, rank, i_rank, j_rank, code
  integer :: rank_cut_plane, i_global, offset_i, j_global, offset_j
  integer :: sender_right_shift,  receiver_right_shift
  integer :: sender_left_shift,   receiver_left_shift
  integer :: sender_bottom_shift, receiver_bottom_shift
  integer :: sender_top_shift,    receiver_top_shift
  integer :: sender_right_top_shift,    receiver_right_top_shift
  integer :: sender_right_bottom_shift, receiver_right_bottom_shift
  integer :: sender_bottom_right_shift, receiver_bottom_right_shift
  integer :: sender_bottom_left_shift,  receiver_bottom_left_shift
  integer :: sender_top_left_shift,     receiver_top_left_shift

  ! Timer to count elapsed time
  character(len=8)  :: datein
  character(len=10) :: timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  integer :: ihours, iminutes, iseconds, int_tCPU
  double precision :: time_start, time_end, tCPU

  ! Inversion
  double precision, allocatable :: RHO_list(:)
  double precision, allocatable :: S_list(:,:), Y_list(:,:)
  double precision, allocatable :: factor_regul_SRdist(:)
  double precision, allocatable :: m0(:), m1(:)
  double precision, allocatable :: x(:), dfx(:)
  double precision, allocatable :: x_old(:), dfx_old(:)
  double precision, allocatable :: x_new(:), dfx_new(:)
  double precision, allocatable :: x_low(:), dfx_low(:)
  double precision, allocatable :: x_high(:), dfx_high(:)
  double precision, allocatable :: r(:), r_old(:)

  double precision :: fx, fx_old, fx_new, fx_low, fx_high
  double precision :: fx_data, fx_regul
  double precision :: dfx_r
  double precision :: alpha_start, alpha, alpha_prec, alpha_low, alpha_high
  double precision :: reg_weight
  integer :: count_grad, count_restart, count_f

  ! Miscellaneous
  double precision :: distance2, factor_ssf
  double precision :: maxval_image_p, maxval_image_rho, maxval_image_vx, maxval_image_vy
  integer, dimension(100) :: list_irec_alt
  integer :: dim_list_irec

contains

  !==========================================================================
  ! Read the parfile and populate all variables
  !==========================================================================
  subroutine read_parfile(parfile_path)
    character(len=*), intent(in) :: parfile_path

    integer            :: unit_par, ios, eq_pos, iset, ibuf
    character(len=512) :: line, key, val
    double precision   :: buf4(4), buf7(7)

    ! ---- First pass: read scalar parameters ----
    unit_par = 20
    open(unit=unit_par, file=trim(parfile_path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'ERROR: cannot open parfile: ', trim(parfile_path)
      stop
    end if
    do
      read(unit_par, '(A)', iostat=ios) line
      if (ios /= 0) exit                        ! end of file

      ! Remove inline comments
      line = strip_comment(line)
      line = adjustl(line)
      if (len_trim(line) == 0) cycle            ! blank line

      eq_pos = index(line, '=')
      if (eq_pos == 0) cycle
      key = to_lower(adjustl(trim(line(1:eq_pos-1))))
      val = adjustl(line(eq_pos+1:))
      call assign_scalar(key, val)
    end do
    close(unit_par)

    ! Validate required scalar parameters before allocation
    call check_parfile()

    ! Allocate arrays whose size depends on scalar parameters
    allocate(NREC_PER_SET(NREC_SET))
    allocate(REC_SET_INFO(NREC_SET, 4))
    REC_SET_INFO = 0.d0

    if (NPERTURB_MODEL > 0) then
      allocate(ADD_PERTURB_MODEL_INFO(NPERTURB_MODEL, 7))
      ADD_PERTURB_MODEL_INFO = 0.d0
    else
      allocate(ADD_PERTURB_MODEL_INFO(0, 7))
    end if

    ! ---- Second pass: read indexed arrays ----
    open(unit=unit_par, file=trim(parfile_path), status='old', action='read', iostat=ios)
    do
      read(unit_par, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = strip_comment(line)
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      eq_pos = index(line, '=')
      if (eq_pos == 0) cycle
      key = to_lower(adjustl(trim(line(1:eq_pos-1))))
      val = adjustl(line(eq_pos+1:))

      ! NREC_PER_SET = v1 [v2 ...]
      if (trim(key) == 'nrec_per_set') then
        read(val, *) NREC_PER_SET(1:NREC_SET)

      ! REC_SET_INFO_N = x1 y1 x2 y2
      else if (key(1:13) == 'rec_set_info_') then
        read(key(14:), *) iset
        if (iset >= 1 .and. iset <= NREC_SET) then
          read(val, *) buf4
          REC_SET_INFO(iset, :) = buf4
        end if

      ! PERTURB_N = type x1 y1 x2 y2 fac_rho fac_c
      else if (key(1:8) == 'perturb_') then
        read(key(9:), *) ibuf
        if (ibuf >= 1 .and. ibuf <= NPERTURB_MODEL) then
          read(val, *) buf7
          ADD_PERTURB_MODEL_INFO(ibuf, :) = buf7
        end if
      end if
    end do
    close(unit_par)

  end subroutine read_parfile


  !==========================================================================
  ! Assign one key=value pair to the correct variable
  !==========================================================================
  subroutine assign_scalar(key, val)
    character(len=*), intent(in) :: key, val

    select case (trim(key))
      case ('method')                      ; read(val,*) method
      case ('validation')                  ; validation   = parse_logical(val)
      case ('which_kernel')                ; read(val,*) which_kernel
      case ('nproc_x')                     ; read(val,*) NPROC_X
      case ('nproc_y')                     ; read(val,*) NPROC_Y
      case ('nx')                          ; read(val,*) NX
      case ('ny')                          ; read(val,*) NY
      case ('deltax')                      ; read(val,*) DELTAX
      case ('deltay')                      ; read(val,*) DELTAY
      case ('nstep')                       ; read(val,*) NSTEP
      case ('deltat')                      ; read(val,*) DELTAT
      case ('use_pml_xmin')                ; USE_PML_XMIN = parse_logical(val)
      case ('use_pml_xmax')                ; USE_PML_XMAX = parse_logical(val)
      case ('use_pml_ymin')                ; USE_PML_YMIN = parse_logical(val)
      case ('use_pml_ymax')                ; USE_PML_YMAX = parse_logical(val)
      case ('f0')                          ; read(val,*) f0
      case ('type_source')                 ; read(val,*) type_source
      case ('wavefront')                   ; read(val,*) wavefront
      case ('xsource')                     ; read(val,*) xsource
      case ('ysource')                     ; read(val,*) ysource
      case ('ssf_sigma')                   ; read(val,*) SSF_Sigma
      case ('observation_from_file')       ; read(val,*) observation_from_file
      case ('path_obs_file')               ; path_obs_file = trim(val)
      case ('observation')                 ; read(val,*) observation
      case ('window_waveform')             ; read(val,*) window_waveform
      case ('nrec_set')                    ; read(val,*) NREC_SET
      case ('atmospheric_model_file')      ; atmospheric_model_file    = parse_logical(val)
      case ('atmospheric_file_name_true')  ; atmospheric_file_name_true  = trim(val)
      case ('atmospheric_file_name_prior') ; atmospheric_file_name_prior = trim(val)
      case ('gamma_chemestry_value')       ; read(val,*) gamma_chemestry_value
      case ('cp_unrelaxed_true')           ; read(val,*) cp_unrelaxed_true
      case ('density_true')                ; read(val,*) density_true
      case ('windx_value_true')            ; read(val,*) windx_value_true
      case ('cp_unrelaxed_prior')          ; read(val,*) cp_unrelaxed_prior
      case ('density_prior')               ; read(val,*) density_prior
      case ('windx_value_prior')           ; read(val,*) windx_value_prior
      case ('nperturb_model')              ; read(val,*) NPERTURB_MODEL
      case ('add_windperturb_profile')     ; add_windperturb_profile = parse_logical(val)
      case ('ymin_wind')                   ; read(val,*) ymin_wind
      case ('ymax_wind')                   ; read(val,*) ymax_wind
      case ('mean_gauss_wind')             ; read(val,*) mean_gauss_wind
      case ('sigma2_gauss_wind')           ; read(val,*) sigma2_gauss_wind
      case ('max_wind_factor')             ; read(val,*) max_wind_factor
      case ('it_display')                  ; read(val,*) IT_DISPLAY
      case ('save_normimage_overtime')     ; read(val,*) save_normimage_overtime
      case ('save_adjoint_source')         ; save_adjoint_source = parse_logical(val)
      case ('n_glob_frames')                     ; read(val,*) N_GLOB_FRAMES
      case ('n_loc_frames')                ; read(val,*) N_LOC_FRAMES
      case ('parametrisation')             ; read(val,*) parametrisation
      case ('scale_model')                 ; read(val,*) scale_model
      case ('type_gradient')               ; read(val,*) type_gradient
      case ('mem_lbfgs')                   ; read(val,*) mem_lbfgs
      case ('type_regul_term')             ; read(val,*) type_regul_term
      case ('regul_weight')                ; read(val,*) regul_weight
      case ('type_smoothing')              ; read(val,*) type_smoothing
      case ('steepest_nbiter_default')     ; read(val,*) steepest_nbiter_default
      case ('maxiter_innerloop')           ; read(val,*) maxiter_innerloop
      case ('maxiter_outerloop')           ; read(val,*) maxiter_outerloop
      case ('tol_x')                       ; read(val,*) tol_x
      case ('alpha_max')                   ; read(val,*) alpha_max
      case ('rate')                        ; read(val,*) rate
      case ('c1')                          ; read(val,*) c1
      case ('c2')                          ; read(val,*) c2
      ! Composite keys are handled in read_parfile()
      case ('nrec_per_set')                ; return
      case default
        if (key(1:13) /= 'rec_set_info_' .and. &
            key(1:7)  /= 'rec_wr_'       .and. &
            key(1:8)  /= 'perturb_') then
          write(*,*)  'WARNING: unknown key "', trim(key), '" ignored.'
        end if
    end select
  end subroutine assign_scalar


  !==========================================================================
  ! Check if all require variables are given
  !==========================================================================
  subroutine check_parfile()
    logical :: ok
    ok = .true.

    ! Check that all required parameters have been properly initialized.
    ! Invalid or non-physical values indicate missing or incorrect entries
    ! in the parameter file.
    if (NX <= 0)      then ; write(*,*) 'ERROR parfile : missing or invalid NX'      ; ok = .false. ; end if
    if (NY <= 0)      then ; write(*,*) 'ERROR parfile : missing or invalid NY'      ; ok = .false. ; end if
    if (NPROC_X <= 0) then ; write(*,*) 'ERROR parfile : missing or invalid NPROC_X' ; ok = .false. ; end if
    if (NPROC_Y <= 0) then ; write(*,*) 'ERROR parfile : missing or invalid NPROC_Y' ; ok = .false. ; end if
    if (NSTEP <= 0)   then ; write(*,*) 'ERROR parfile : missing or invalid NSTEP'   ; ok = .false. ; end if
    if (DELTAX <= 0)  then ; write(*,*) 'ERROR parfile : missing or invalid DELTAX'  ; ok = .false. ; end if
    if (DELTAY <= 0)  then ; write(*,*) 'ERROR parfile : missing or invalid DELTAY'  ; ok = .false. ; end if
    if (DELTAT <= 0)  then ; write(*,*) 'ERROR parfile : missing or invalid DELTAT'  ; ok = .false. ; end if
    if (f0 <= 0)      then ; write(*,*) 'ERROR parfile : missing or invalid f0'      ; ok = .false. ; end if
    if (NREC_SET <= 0) then ; write(*,*) 'ERROR parfile : missing or invalid NREC_SET' ; ok = .false. ; end if
    if (method < 1 .or. method > 3) then
      write(*,*) 'ERROR parfile : method should be 1, 2 ou 3' ; ok = .false.
    end if

    if (.not. ok) stop 'Error: incomplete or invalid parameter file.'
  end subroutine check_parfile


  !==========================================================================
  ! Compute all derived quantities after parfile has been read
  !==========================================================================
  subroutine init_parameters(parfile_path)
    character(len=*), intent(in) :: parfile_path

    integer :: unit_par, ios, eq_pos, irec, i1
    character(len=512) :: line, key, val
    double precision   :: buf2(2)

    !
    NPROC    = NPROC_X * NPROC_Y
    NX_LOCAL = NX / NPROC_X
    NY_LOCAL = NY / NPROC_Y
    Nflat    = 3 * NY_LOCAL

    ! Number of receivers
    NREC     = sum(NREC_PER_SET)

    ! Source grid indices
    ISOURCE = int(xsource / DELTAX) + 1
    JSOURCE = int(ysource / DELTAY) + 1

    ! Wind-model grid indices
    jmin_wind = int(real(ymin_wind) / DELTAY) + 1
    jmax_wind = int(real(ymax_wind) / DELTAY) + 1

    ! Source time function constants
    t0 = 1.20d0 / f0
    a  = PI * PI * f0 * f0

    ! Finite-difference coefficients
    NINE_OVER_8_DELTAX  = 9.d0 / (8.0d0 * DELTAX)
    NINE_OVER_8_DELTAY  = 9.d0 / (8.0d0 * DELTAY)
    ONE_OVER_24_DELTAX  = 1.d0 / (24.d0 * DELTAX)
    ONE_OVER_24_DELTAY  = 1.d0 / (24.d0 * DELTAY)
    ONE_OVER_SIX_DELTAX = 1.d0 / (6.0d0 * DELTAX)
    ONE_OVER_SIX_DELTAY = 1.d0 / (6.0d0 * DELTAY)
    ONE_OVER_DELTAT     = 1.d0 / DELTAT
    ONE_OVER_DELTAX     = 1.d0 / DELTAX
    ONE_OVER_DELTAY     = 1.d0 / DELTAY

    ! MPI message sizes
    number_of_values_x = 2 * (NY_LOCAL + 4)
    number_of_values_y = 2 * (NX_LOCAL + 4)

    ! PML
    ALPHA_MAX_PML   = 2.d0 * PI * (f0 / 2.d0)

    ! Image normalization parameter
    maxval_image_p   = -1.0d0
    maxval_image_rho = -1.0d0
    maxval_image_vx  = -1.0d0
    maxval_image_vy  = -1.0d0

    ! Allocate arrays including halo cells
    i1 = -1 ! Lower array bound including two ghost cells

    ! --- Receivers ---
    allocate(REC_wr(NREC, 2))
    REC_wr = 0.d0

    ! Additional pass to read REC_wr_N entries
    unit_par = 21
    open(unit=unit_par, file=trim(parfile_path), status='old', action='read', iostat=ios)
    do
      read(unit_par, '(A)', iostat=ios) line
      if (ios /= 0) exit
      line = strip_comment(line)
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      eq_pos = index(line, '=')
      if (eq_pos == 0) cycle
      key = to_lower(adjustl(trim(line(1:eq_pos-1))))
      val = adjustl(line(eq_pos+1:))
      if (key(1:7) == 'rec_wr_') then
        read(key(8:), *) irec
        if (irec >= 1 .and. irec <= NREC) then
          read(val, *) buf2
          REC_wr(irec, :) = buf2
        end if
      end if
    end do
    close(unit_par)

    allocate(ix_rec(NREC), iy_rec(NREC))
    allocate(xrec(NREC), yrec(NREC))
    allocate(normsq_pressure_true_per_rec(NREC))

    ! --- Primary fields ---
    allocate(pressure(i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(rhop    (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(vx      (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(vy      (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(rhoa    (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(pa      (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(vax     (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(vay     (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(kappa_unrelaxed_prior(i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(rho0_prior           (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(p0_prior             (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(windx_prior          (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(windy_prior          (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(c0_prior             (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(kappa_unrelaxed_true (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(rho0_true            (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(p0_true              (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(windx_true           (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(windy_true           (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(gamma_chemestry      (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(g                    (i1:NX_LOCAL+2, i1:NY_LOCAL+2))

    ! --- Seismograms ---
    allocate(sisvx            (NSTEP, NREC))
    allocate(sisvy            (NSTEP, NREC))
    allocate(sispressure      (NSTEP, NREC))
    allocate(sisrhop          (NSTEP, NREC))
    allocate(sispressure_true (NSTEP, NREC))
    allocate(sispressure_prior(NSTEP, NREC))
    allocate(sispressure_source(NSTEP))
    allocate(adjoint_source   (NSTEP, NREC))
    allocate(wr(NSTEP))

    ! --- Sensitivity kernels ---
    allocate(K_rho0 (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(K_windx(i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(K_windy(i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(K_p0   (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(K      (i1:NX_LOCAL+2, i1:NY_LOCAL+2))

    ! --- Checkpointing ---
    allocate(GLOB_FRAMES    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, 4, N_GLOB_FRAMES))
    allocate(LOC_FRAMES(i1:NX_LOCAL+2, i1:NY_LOCAL+2, 4, N_LOC_FRAMES))

    ! --- One-dimensional PML profiles ---
    allocate(d_x(NX), K_x(NX), alpha_x(NX), a_x(NX), b_x(NX), c_x(NX))
    allocate(d_x_half(NX), K_x_half(NX), alpha_x_half(NX))
    allocate(a_x_half(NX), b_x_half(NX), c_x_half(NX))
    allocate(one_over_K_x(NX), one_over_K_x_half(NX))
    allocate(one_over_Kdalpha_x(NX), one_over_Kdalpha_x_half(NX))
    allocate(d_y(NY), K_y(NY), alpha_y(NY), a_y(NY), b_y(NY), c_y(NY))
    allocate(d_y_half(NY), K_y_half(NY), alpha_y_half(NY))
    allocate(a_y_half(NY), b_y_half(NY), c_y_half(NY))
    allocate(one_over_K_y(NY), one_over_K_y_half(NY))
    allocate(one_over_Kdalpha_y(NY), one_over_Kdalpha_y_half(NY))

    ! --- Forward PML checkpoint arrays ---
    allocate(eq1_memory_dp0_dx_fw       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dp0_dy_fw       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_drho0_dx_fw     (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_drho0_dy_fw     (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dpressure_dx_fw (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dpressure_dy_fw (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_drhop_dx_fw     (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_drhop_dy_fw     (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dvx_dx_fw       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dvy_dy_fw       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dwindx_dx_fw    (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dwindy_dy_fw    (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dpressure_dx_fw (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_drho0_dx_fw     (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_drho0_dy_fw     (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dvx_dx_fw       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dvx_dy_fw       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dwindx_dx_fw    (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dwindx_dy_fw    (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dwindy_dy_fw    (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq3_memory_dpressure_dy_fw (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq3_memory_drho0_dy_fw     (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq3_memory_drho0_dx_fw     (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq3_memory_dvy_dy_fw       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq3_memory_dvy_dx_fw       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq3_memory_dwindy_dy_fw    (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq3_memory_dwindy_dx_fw    (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq3_memory_dwindx_dx_fw    (i1:NX_LOCAL+2, i1:NY_LOCAL+2))

    ! --- Adjoint PML memory variables ---
    allocate(eq1_memory_drhoa_dx_adj        (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dp0_dx_adj          (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dgammap0pa_dx_adj   (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dgammap0pa_dy_adj   (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dwindx_dx_adj       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dwindy_dx_adj       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dvax_dx_adj         (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dvax_dy_adj         (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_drhoa_dy_adj        (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dp0_dy_adj          (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dwindy_dy_adj       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dwindx_dy_adj       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dvay_dx_adj         (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq1_memory_dvay_dy_adj         (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dpa_dx_adj          (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_drhoa_dx_adj        (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dpa_dy_adj          (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_drhoa_dy_adj        (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dwindx_dx_adj       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dwindy_dy_adj       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dvax_dx_adj         (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dvay_dy_adj         (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dwindy_dx_adj       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))
    allocate(eq2_memory_dwindx_dy_adj       (i1:NX_LOCAL+2, i1:NY_LOCAL+2))

    ! --- Forward PML checkpoint arrays ---
    allocate(eq1_memory_dp0_dx_fw_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq1_memory_dp0_dy_fw_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq1_memory_drho0_dx_fw_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq1_memory_drho0_dy_fw_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq1_memory_dpressure_dx_fw_fr (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq1_memory_dpressure_dy_fw_fr (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq1_memory_drhop_dx_fw_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq1_memory_drhop_dy_fw_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq1_memory_dvx_dx_fw_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq1_memory_dvy_dy_fw_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq1_memory_dwindx_dx_fw_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq1_memory_dwindy_dy_fw_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq2_memory_dpressure_dx_fw_fr (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq2_memory_drho0_dx_fw_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq2_memory_drho0_dy_fw_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq2_memory_dvx_dx_fw_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq2_memory_dvx_dy_fw_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq2_memory_dwindx_dx_fw_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq2_memory_dwindx_dy_fw_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq2_memory_dwindy_dy_fw_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq3_memory_dpressure_dy_fw_fr (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq3_memory_drho0_dy_fw_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq3_memory_drho0_dx_fw_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq3_memory_dvy_dy_fw_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq3_memory_dvy_dx_fw_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq3_memory_dwindy_dy_fw_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq3_memory_dwindy_dx_fw_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))
    allocate(eq3_memory_dwindx_dx_fw_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_GLOB_FRAMES))

    ! --- Forward PML local checkpoint arrays ---
    allocate(eq1_memory_dp0_dx_fw_loc_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq1_memory_dp0_dy_fw_loc_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq1_memory_drho0_dx_fw_loc_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq1_memory_drho0_dy_fw_loc_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq1_memory_dpressure_dx_fw_loc_fr (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq1_memory_dpressure_dy_fw_loc_fr (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq1_memory_drhop_dx_fw_loc_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq1_memory_drhop_dy_fw_loc_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq1_memory_dvx_dx_fw_loc_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq1_memory_dvy_dy_fw_loc_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq1_memory_dwindx_dx_fw_loc_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq1_memory_dwindy_dy_fw_loc_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq2_memory_dpressure_dx_fw_loc_fr (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq2_memory_drho0_dx_fw_loc_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq2_memory_drho0_dy_fw_loc_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq2_memory_dvx_dx_fw_loc_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq2_memory_dvx_dy_fw_loc_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq2_memory_dwindx_dx_fw_loc_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq2_memory_dwindx_dy_fw_loc_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq2_memory_dwindy_dy_fw_loc_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq3_memory_dpressure_dy_fw_loc_fr (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq3_memory_drho0_dy_fw_loc_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq3_memory_drho0_dx_fw_loc_fr     (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq3_memory_dvy_dy_fw_loc_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq3_memory_dvy_dx_fw_loc_fr       (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq3_memory_dwindy_dy_fw_loc_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq3_memory_dwindy_dx_fw_loc_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))
    allocate(eq3_memory_dwindx_dx_fw_loc_fr    (i1:NX_LOCAL+2, i1:NY_LOCAL+2, N_LOC_FRAMES))

    ! --- Inversion ---
    allocate(RHO_list(mem_lbfgs))
    allocate(S_list(mem_lbfgs, Nflat), Y_list(mem_lbfgs, Nflat))
    allocate(factor_regul_SRdist(Nflat))
    allocate(m0(Nflat), m1(Nflat))
    allocate(x(Nflat), dfx(Nflat))
    allocate(x_old(Nflat), dfx_old(Nflat))
    allocate(x_new(Nflat), dfx_new(Nflat))
    allocate(x_low(Nflat), dfx_low(Nflat))
    allocate(x_high(Nflat), dfx_high(Nflat))
    allocate(r(Nflat), r_old(Nflat))

  end subroutine init_parameters


  !==========================================================================
  ! Utility functions
  !==========================================================================
  function strip_comment(line) result(clean)
    character(len=*), intent(in) :: line
    character(len=len(line))     :: clean
    integer :: pos
    pos = index(line, '#')
    if (pos == 0) then ; clean = line
    else               ; clean = line(1:pos-1)
    end if
  end function strip_comment

  function to_lower(str) result(low)
    character(len=*), intent(in) :: str
    character(len=len(str))      :: low
    integer :: i, c
    low = str
    do i = 1, len(str)
      c = ichar(str(i:i))
      if (c >= 65 .and. c <= 90) low(i:i) = char(c + 32)
    end do
  end function to_lower

  function parse_logical(val) result(b)
    character(len=*), intent(in) :: val
    logical :: b
    character(len=len(val)) :: v
    v = to_lower(adjustl(val))
    select case (trim(v))
      case ('.true.', 'true', '1', 'yes') ; b = .true.
      case default                         ; b = .false.
    end select
  end function parse_logical

end module parameters
