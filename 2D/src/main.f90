!======================================================================
!> Main program.
!!
!! Initializes the simulation, configures the MPI domain decomposition,
!! builds the computational model, and executes one of the available
!! computation modes:
!!   - forward modeling,
!!   - sensitivity-kernel computation,
!!   - full-waveform inversion.
!======================================================================

program main

 use MPI
 use parameters
 implicit none
 
 integer :: irec
 double precision :: Courant_number

 character(len=200) :: parfile_path

! ------------------------------------------------------------------
! Read the input parameter file
! ------------------------------------------------------------------
  if (command_argument_count() >= 1) then
    call get_command_argument(1, parfile_path)
  else
    parfile_path = "parfile"   ! valeur par defaut
  end if

! Read simulation parameters, including those required to configure MPI decomposition
  call read_parfile(trim(parfile_path))
  call init_parameters(trim(parfile_path)) 


! ------------------------------------------------------------------
! Initialize MPI
! ------------------------------------------------------------------
  ! start MPI processes
  call MPI_INIT(code)
  
! Get the total number of MPI processes
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

! Get the rank of the current MPI process (0 to nb_procs-1)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)
  
! Determine neighboring MPI ranks used for data exchanges
  call get_neighbors()
  
  
! ------------------------------------------------------------------
! Initialize the simulation timer
! ------------------------------------------------------------------
  call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! Note: this simple computation does not account for month/year changes
  time_start = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
               60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0
  

! ------------------------------------------------------------------
! Print simulation information and create output directories
! ------------------------------------------------------------------
  if (rank == 0) then
  print *
  print *,'2D acoustic finite-differences code in density, velocity and pressure'
  print *

! display size of the model
  print *
  print *,'NX = ',NX
  print *,'NY = ',NY
  print *
  print *,'size of the model along X = ',(NX - 1) * DELTAX
  print *,'size of the model along Y = ',(NY - 1) * DELTAY
  print *
  print *,'Total number of grid points = ',NX * NY
  print *
 endif

! create output directories
  if (rank == 0 )then
    ! suppress old files (can be commented out if "call system" is missing in your compiler)
    call system('rm -rf OUTPUT/ MODELS')
    ! Create output directories
    call execute_command_line ('mkdir -p OUTPUT/')
    call execute_command_line ('mkdir -p MODELS/')
    if (method == 3) then
     call execute_command_line ('mkdir -p OUTPUT_INVERSION/')
    endif
  endif


! ------------------------------------------------------------------
! Initialize the Perfectly Matched Layers (PML)
! ------------------------------------------------------------------
!  Define the parameters for the Perfectly Matched Layers (PML)
  call computePML()


! ------------------------------------------------------------------
! Initialize the background atmospheric models
! ------------------------------------------------------------------
! Define the density, pressure, and wind fields for the true model and the prior model
 call init_backgrounds()

 
! ------------------------------------------------------------------
! Initialize the source and receiver geometry
! ------------------------------------------------------------------
! Define the position of the source and receivers
  call init_source_recvrs()


! ------------------------------------------------------------------
! Check the CFL stability condition
! ------------------------------------------------------------------
! Check the CFL stability condition for the explicit time integration scheme
! R. Courant, K. O. Friedrichs and H. Lewy (1928)
! For this O(2,4) scheme, when DELTAX == DELTAY the Courant number is given by Levander, Geophysics, vol. 53(11), p. 1427,
! equation (7) and is equal to 0.606 (it is thus smaller than that of the O(2,2) scheme, which is 1/sqrt(2) = 0.707,
! i.e. when switching to a fourth-order spatial scheme one needs a time step that is about 0.707 / 0.606 = 1.167 times smaller.
  if (DELTAX == DELTAY) then
    Courant_number = cp_unrelaxed_prior * DELTAT / DELTAX
    if (rank == 0) then
      print *,'Courant number is ',Courant_number
      print *,' (the maximum possible value is 0.606; in practice for accuracy reasons a value not larger than 0.30 is recommended)'
      print *
    endif
    if (Courant_number > 0.606) stop 'time step is too large, simulation will be unstable'
  endif
    
! ------------------------------------------------------------------
! Select the computation mode
! ------------------------------------------------------------------
! method = 1 : Forward modeling
! method = 2 : Sensitivity-kernel computation
! method = 3 : Full-waveform inversion
! ------------------------------------------------------------------

! ==================================================================
! Forward modeling
! ==================================================================
 if (method == 1) then
 
   ! Initialize all parameters required for the forward problem
   call reset_forward()

   ! Option to save normalized snapshots
   ! If save_normimage_overtime = 1: Save snapshots normalized by the maximum value of the waveform during all the propagation
   ! WARNING: This option needs two simulations:
   ! One to identify the maximum value, and another to save the normalized snapshots.
   if (save_normimage_overtime == 1) then
     save_sismos = .False.
     call forwardproblem(p0_true, rho0_true, windx_true, windy_true,  1, NSTEP, 2) 
     save_normimage_overtime = 0
   endif
    
   ! Initialize all the parameters to solve the forward problem (p',rho',vx',vy' and PML memory variables)
   call reset_forward()
   ! Parameter to save images and pressure records
   save_sismos = .True.
   ! Solve the forward problem
   call forwardproblem(p0_true, rho0_true, windx_true, windy_true,  1, NSTEP, 2) 
   
   
! ==================================================================
! Sensitivity-kernel computation
! ==================================================================
  else if (method == 2) then
 
   ! Generate or load observed data
   if (observation_from_file == 0) then
     ! Initialize all the parameters to solve the forward problem (p',rho',vx',vy' and PML memory variables)
     call reset_forward()
     ! Parameter to save images and pressure records
     save_sismos = .True.
     ! Solve the forward problem
     call forwardproblem(p0_true, rho0_true, windx_true, windy_true,  1, NSTEP, 2) 
   elseif (observation_from_file == 1) then
       save_sismos = .True. 
       call read_obs()
   elseif (observation_from_file == 2) then ! Observation is delay time from cross-correlation 
      if (observation == 0) then
        print *, "ERROR: if the given observation is a delay time, it should be sensitivity kernels of arrival time."
      endif
      if (NREC /= 1) then
        print *, "ERROR: if the observation is a delay time, only 1 receiver is possible (for now)."
      endif
      
      REC_wr(:,1) = 0 
      REC_wr(:,2) = NSTEP * DELTAT
   endif
    
   if (.not. (observation_from_file == 2)) then
       ! Get normalization information from observations
       call write_seismograms(sisvx,sisvy,sispressure,sisrhop,NSTEP,NREC,DELTAT,t0,2)
       call MPI_BARRIER(MPI_COMM_WORLD, code)
       sispressure_true(:,:) = sispressure(:,:) 
   
       do irec=1,NREC
        normsq_pressure_true_per_rec(irec) = DELTAT * DELTAX *DELTAY * sum(sispressure_true(:,irec)**2)
       enddo
   endif


   if (validation) then
   call compute_gradient()
   
   else 
   save_sismos = .True.
   ! Compute sensivity kernels
   call compute_kernel()
   call write_kernels()
endif	
 

! ==================================================================
! Full-waveform inversion
! ==================================================================
  elseif (method == 3 ) then
     
     ! TODO add an option to load seismograms directly
   
     ! Get an observation (for now, you can only create observations)
     ! Initialize all the parameters to solve the forward problem (p',rho',vx',vy' and PML memory variables)
     call reset_forward()
     ! Parameter to save images and pressure records
     save_sismos = .True.
     ! Solve the forward problem
     call forwardproblem(p0_true, rho0_true, windx_true, windy_true,  1, NSTEP, 2) 

     ! Compute normalization factors from the observations.
     sispressure_true(:,:) = sispressure(:,:) 
     do irec=1,NREC
       normsq_pressure_true_per_rec(irec) =  DELTAT * DELTAX * DELTAY * sum(sispressure_true(:,irec)**2)
     enddo

     ! Atmospheric model is saved in a vector form for inverse problem
     call priormodel2flatmodel(m0)
   
     call MPI_ALLREDUCE(maxval(m0(1:NY_LOCAL)), regul_term_rho0_prior, 1,MPI_DOUBLE_PRECISION,&
                  MPI_MAX,MPI_COMM_WORLD, code)
     call MPI_ALLREDUCE(maxval(m0(NY_LOCAL+1:2*NY_LOCAL)),regul_term_p0_prior,1,MPI_DOUBLE_PRECISION,&
                  MPI_MAX,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(maxval(m0(2*NY_LOCAL+1:Nflat)), regul_term_windx_prior,1, MPI_DOUBLE_PRECISION, &
                  MPI_MAX,MPI_COMM_WORLD, code)
     if (regul_term_windx_prior < TINYVAL) then
        regul_term_windx_prior = 1.0
     endif

    
     ! Solve the inverse problem
     call optimisation()
     
   else 
     print *, "Can only do 1: Forward, 2: Kernel or 3: Inversion computations"
     stop
 endif


! ------------------------------------------------------------------
! Finalize MPI
! ------------------------------------------------------------------  
  ! Close MPI program
  call MPI_FINALIZE(code)

 end program  main
