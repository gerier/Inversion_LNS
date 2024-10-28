

program main

 use MPI
 use parameters
 implicit none
 
 !integer i,j,
 integer :: irec
 double precision :: Courant_number

     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                      MPI Init                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! start MPI processes
  call MPI_INIT(code)
  
! get total number of MPI processes in variable nb_procs
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

! get the rank of our process from 0 (master) to nb_procs-1 (workers)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)
  
! get the neighboor of each process, to know for information exchanges
  call get_neighboors()
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       Init time to follow the evolution                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! this fails if we cross the end of the month
  time_start = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
               60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                    INITIALISATION                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


  if (rank == 0 )then
    ! suppress old files (can be commented out if "call system" is missing in your compiler)
    call system('if [ ! -f ../OtherStation/ ] ; then rm -r OUTPUT/ MODELS/; fi')
    ! create directory to store results
    call execute_command_line ('mkdir '//'OUTPUT/')
    call execute_command_line ('mkdir '//'MODELS/')
    if (method == 3) then
     call execute_command_line ('mkdir '//'OUTPUT_INVERSION/')
    endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                       PML INIT                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !  define the parameters for the Perfectly Match Layers (PML)
  call computePML()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                Background domain                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! define the matrix of density, pressure and wind for the true model and the a priori model
 call init_backgrounds()
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          source / receivers INFORMATION          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! define the position of the source and receivers
  call init_source_recvrs()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                   CFL Condition                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! check the Courant stability condition for the explicit time scheme
! R. Courant, K. O. Friedrichs and H. Lewy (1928)
! For this O(2,4) scheme, when DELTAX == DELTAY the Courant number is given by Levander, Geophysics, vol. 53(11), p. 1427,
! equation (7) and is equal to 0.606 (it is thus smaller than that of the O(2,2) scheme, which is 1/sqrt(2) = 0.707,
! i.e. when switching to a fourth-order spatial scheme one needs a time step that is about 0.707 / 0.606 = 1.167 times smaller.
  if (DELTAX == DELTAY) then
    Courant_number = cp_unrelaxed * DELTAT / DELTAX
    if (rank == 0) then
      print *,'Courant number is ',Courant_number
      print *,' (the maximum possible value is 0.606; in practice for accuracy reasons a value not larger than 0.30 is recommended)'
      print *
    endif
    if (Courant_number > 0.606) stop 'time step is too large, simulation will be unstable'
  endif
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!               OPTIONS OF APPLICATIONS                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  ! INITIALISATION

 !!! FORWARD PROPAGATION
 if (method == 1) then
 
   ! initialise all the parameters to solve the forward problem (p',rho',vx',vy' and PML memory variables)
   call reset_forward()

   ! option to save normalised image (on all the snapshots)
   ! if save_normimage_overtime = 1 : save snapshots normalised by the maximum alue of the waveform during all the propagation
   !! BE AWARE: this option needs two simulations: one to identify the maximum value, and a second one to save the snapshots 
   if (save_normimage_overtime == 1) then
     save_sismos = .False.
     call forwardproblem(p0_true, rho0_true, windx_true, windy_true,  1, NSTEP, 2) 
     save_normimage_overtime = 0
   endif
    
   ! initialise all the parameters to solve the forward problem (p',rho',vx',vy' and PML memory variables)
   call reset_forward()
   ! parameter to save images and pressure records
   save_sismos = .True.
   ! solve the forward problem
   call forwardproblem(p0_true, rho0_true, windx_true, windy_true,  1, NSTEP, 2) 
   
   
  !!! SENSITIVITY KERNELS 
  else if (method == 2) then
 
   ! TODO add an option to load seismograms directly
   
   ! get an observation (for now, you can only create observations)
   ! initialise all the parameters to solve the forward problem (p',rho',vx',vy' and PML memory variables)
   call reset_forward()
   ! parameter to save images and pressure records
   save_sismos = .True.
   ! solve the forward problem
   call forwardproblem(p0_true, rho0_true, windx_true, windy_true,  1, NSTEP, 2) 

   ! get normalisation informations from observations
   call write_seismograms(sisvx,sisvy,sispressure,sisrhop,NSTEP,NREC,DELTAT,t0,2)
   call MPI_BARRIER(MPI_COMM_WORLD, code)
   sispressure_true(:,:) = sispressure(:,:) 
   
   do irec=1,NREC
    normsq_pressure_true_per_rec(irec) = sum(sispressure_true(:,irec)**2)
   enddo
     
   ! compute kernel
   call compute_kernel()
   call write_kernels()
	
 

  !!! INVERSE PROBLEM
  elseif (method == 3 ) then
     
     ! TODO add an option to load seismograms directly
   
     ! get an observation (for now, you can only create observations)
     ! initialise all the parameters to solve the forward problem (p',rho',vx',vy' and PML memory variables)
     call reset_forward()
     ! parameter to save images and pressure records
     save_sismos = .True.
     ! solve the forward problem
     call forwardproblem(p0_true, rho0_true, windx_true, windy_true,  1, NSTEP, 2) 

     ! get normalisation informations from observations
     sispressure_true(:,:) = sispressure(:,:) 
     do irec=1,NREC
       normsq_pressure_true_per_rec(irec) = sum(sispressure_true(:,irec)**2)
     enddo

     ! atmospheric model is saved in a vector form for inverse problem
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

      
     ! init for regularisation term
     !if (type_regul_term > 0) then
     !  call init_factor_regul()
     !  call get_norm_apriori()
     !endif 
     
     ! resolution of the inverse problem
     call optimisation()
     
   else 
     print *, "Can only do 1: Forward, 2: Kernel or 3: Inversion computations"
     stop
 endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                          End                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  ! close MPI program
  call MPI_FINALIZE(code)

 end program  main
