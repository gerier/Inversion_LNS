

program main

 use MPI
 use parameters
 implicit none
 
 integer i,j,irec, ii, jj
 double precision :: Courant_number
 
  integer :: it, it_time, time_last_frame

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                      MPI Init                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! start MPI processes
  call MPI_INIT(code)
  
! get total number of MPI processes in variable nb_procs
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

! get the rank of our process from 0 (master) to nb_procs-1 (workers)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, code)
  
  
! check that code was compiled with the right number of slices
  if (nb_procs /= NPROC) then
    print *,'error in MPI number of slices: nb_procs,NPROC = ',nb_procs,NPROC,' but they should be equal'
    stop 'nb_procs must be equal to NPROC'
  endif

! we restrict ourselves to an even number of slices
! in order to have a cut plane in the middle of the mesh for visualization purposes
  if (mod(nb_procs,2) /= 0) stop 'nb_procs must be even'

! check that we can cut along Z in an exact number of slices
  if (mod(NX,nb_procs) /= 0) stop 'NX must be a multiple of nb_procs'

! check that a slice is at least as thick as a PML layer
  if (NX_LOCAL < NPOINTS_PML) stop 'NX_LOCAL must be greater than NPOINTS_PML'

  ! offset of this slice when we cut along Z
  offset_i = rank * NX_LOCAL
  
  ! slice number for the cut plane in the middle of the mesh
  rank_cut_plane = nb_procs/2 - 1
  
  
  
  ! we receive from the process on the left, and send to the process on the right
  sender_right_shift = rank - 1
  receiver_right_shift = rank + 1

! if we are the first process, there is no neighbor on the left
  if (USE_PML_XMIN) then
    if (rank == 0) sender_right_shift = MPI_PROC_NULL
  else
    if (rank == 0) sender_right_shift = nb_procs - 1
  endif
  
! if we are the last process, there is no neighbor on the right
  if (USE_PML_XMAX) then
    if (rank == nb_procs - 1) receiver_right_shift = MPI_PROC_NULL
  else
    if (rank == nb_procs - 1) receiver_right_shift = 0
  endif
  ! we receive from the process on the right, and send to the process on the left
  sender_left_shift = rank + 1
  receiver_left_shift = rank - 1

! if we are the first process, there is no neighbor on the left
  if (USE_PML_XMIN) then
    if (rank == 0) receiver_left_shift = MPI_PROC_NULL
  else
    if (rank == 0) receiver_left_shift = nb_procs - 1
  endif

! if we are the last process, there is no neighbor on the right
  if (USE_PML_XMAX) then
    if (rank == nb_procs - 1) sender_left_shift = MPI_PROC_NULL
  else
    if (rank == nb_procs - 1) sender_left_shift = 0
  endif
  
  i2begin = 1
  if (rank == 0) i2begin = 2 ! TODO modifier pour voir si necessaire la taille supp ?

  iminus1end = NX_LOCAL
  if (rank == nb_procs - 1) iminus1end = NX_LOCAL-1 ! TODO idem l.78
  
  
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
  print *,'2D acoustic finite-difference code in density, velocity and pressure'
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
    call execute_command_line ('mkdir '//'OUTPUT/')
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                       PML INIT                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call computePML()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                Background domain                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! compute the Lame parameter and density


!!! Compute the true model
  do j = 1,NY
    do i = 1,NX_LOCAL
    
      i_global = i + offset_i
      if ((i_global < IObs_start) .or. (i_global >= IObs_end) .or. (j < JObs_start) .or. (j >= JObs_end)) then 
        ! to have adiscontinuity in the background model
        ! attention, pour la comparaison
        !avec le code python, les bornes sont placés en Iobs-1 et non en Iobs
        rho0_true(i,j)            = density
        p0_true(i,j)              = density * cp_unrelaxed * cp_unrelaxed / gamma_chimie
        kappa_unrelaxed_true(i,j) = gamma_chimie * p0_true(i,j)
        windx_true(i,j)           = 0
        windy_true(i,j)           = 0
      else
        rho0_true(i,j)            = density  * obstacle_factor_rho
        p0_true(i,j)              = rho0_true(i,j) * (cp_unrelaxed)**2 * obstacle_factor_c2 / gamma_chimie
        kappa_unrelaxed_true(i,j) = gamma_chimie * p0_true(i,j)
        windx_true(i,j)           = 0
        windy_true(i,j)           = 0
       endif
       
      ! to have a wind profil  
      if (add_wind_profile .and. (j > 25 .and. j <= 175)) then 
        windy_true(i,j) = exp(- ((j-1)*0.1d0 - 10)**2 / 5 ) * 150
      endif
      
    enddo
  enddo


!!! Compute the a priori model
  do j = 1,NY
    do i = 1,NX_LOCAL
        ! to have a discontinuity in the background model
        ! attention, pour la comparaison
        !avec le code python, les bornes sont placés en Iobs-1 et non en Iobs
        rho0_prior(i,j)            = density
        kappa_unrelaxed_prior(i,j) = density * cp_unrelaxed * cp_unrelaxed
        p0_prior(i,j)              = density * cp_unrelaxed * cp_unrelaxed / gamma_chimie
        windx_prior(i,j)           = 0
        windy_prior(i,j)           = 0
    enddo
  enddo

        
  ! write the background model
  call write_background(rho0_true, kappa_unrelaxed_true, p0_true, windx_true, windy_true, gamma_chimie)
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          source / receivers INFORMATION          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (rank == 0) then
    ! print position of the source
    print *,'Position of the source:'
    print *
    print *,'x = ',xsource
    print *,'y = ',ysource
    print *

    ! define location of receivers
    print *,'There are ',nrec,' receivers'
    print *
  endif
  
  if (NREC > 1) then
! this is to avoid a warning with GNU gfortran at compile time about division by zero when NREC = 1
    myNREC = NREC
    xspacerec = (xfin-xdeb) / dble(myNREC-1)
    yspacerec = (yfin-ydeb) / dble(myNREC-1)
  else
    xspacerec = 0.d0
    yspacerec = 0.d0
  endif
  do irec=1,nrec
    xrec(irec) = xdeb + dble(irec-1)*xspacerec
    yrec(irec) = ydeb + dble(irec-1)*yspacerec
  enddo

! find closest grid point for each receiver
  do irec=1,nrec
    dist = HUGEVAL
    do j = 1,NY
    do i = 1,NX
      distval = sqrt((DELTAX*dble(i-1) - xrec(irec))**2 + (DELTAY*dble(j-1) - yrec(irec))**2)
      if (distval < dist) then
        dist = distval
        ix_rec(irec) = i
        iy_rec(irec) = j
      endif
    enddo
    enddo
    
    if (rank ==0) then
    print *,'receiver ',irec,' x_target,y_target = ',xrec(irec),yrec(irec)
    print *,'closest grid point found at distance ',dist,' in i,j = ',ix_rec(irec),iy_rec(irec)
    print *
    endif
  enddo
    
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
!!             Retrieve the observations                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  ! INITIALISATION
  call reset_forward()
  save_sismos = .True.
  
  !call forwardproblem(p0_true, rho0_true, windx_true, windy_true, kappa_unrelaxed_true,  1, NSTEP, 3) 
  call forwardproblem(p0_prior, rho0_prior,windx_prior,windy_prior,kappa_unrelaxed_prior,  1, NSTEP, 3) 
  sispressure_true(:,:) = sispressure(:,:)
  norm_pressure_true = maxval(abs(sispressure_true))**2
  
  do i=0,nb_procs-1
  if (rank == i) then
    OPEN(UNIT=12, FILE="OUTPUT/p_true.txt", position="append", ACTION="write")
    DO ii=1,NX_LOCAL
      WRITE(12,*) (pressure(ii,jj), jj=1,NY)
    END DO
    CLOSE(12)
  endif
  call mpi_barrier(mpi_comm_world,code)
  enddo


  call reset_kernel()
  save_sismos = .False.
  
  call save_frames()
   do it=0,NSTEP
   
   if (NSTEP-it == NSTEP-1 .or. modulo(NSTEP-it,NSTEP/NFRAMES) == NSTEP/NFRAMES-1 ) then
      print *, 1,it, NSTEP-it
      call save_local_frames(NSTEP-it)
   else if (NSTEP-it /= NSTEP) then
     
     call load_frame(NSTEP-it, time_last_frame)
     print *, 2,it, NSTEP-it, time_last_frame
     call forwardproblem(p0_prior,rho0_prior,windx_prior,windy_prior,kappa_unrelaxed_prior,time_last_frame, NSTEP-it,3)

   else
     print *, 3, it, NSTEP-it
   endif

    if (rank == 0) then
   OPEN(UNIT=1222, FILE="./OUTPUT/p_checkpoint.txt", position="append", ACTION="write")
   WRITE(1222,*) it, pressure(ix_rec(1) - offset_i,iy_rec(1))
   !print *, pressure(ix_rec(1),ix_rec(1)), ix_rec(1),ix_rec(1)
   CLOSE(1222)
   endif
    
   enddo
  
   
  ! close MPI program
  call MPI_FINALIZE(code)
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                  Compute gradient                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !call compute_kernel()
    
 end program  main
