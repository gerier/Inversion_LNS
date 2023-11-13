

program main

 use MPI
 use parameters
 implicit none
 
 integer i,j,irec, ii, jj
 double precision :: Courant_number

double precision, dimension(1:NX,1:NY) :: pressure_global
integer :: rk, ii_rk,jj_rk
character(len=100) :: file_name
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
  !if (mod(NPROC_X,2) /= 0) stop 'nb_procs_x must be even'
  !if (mod(NPROC_Y,2) /= 0) stop 'nb_procs_y must be even'
  !if (NPROC_X == 1) stop 'nb_procs_x must be not 1'
  !if (NPROC_Y == 1) stop 'nb_procs_y must be not 1'
  
! check that we can cut along Z in an exact number of slices
  if (mod(NX,NPROC_X) /= 0) stop 'NX must be a multiple of nb_procs_x'
  if (mod(NY,NPROC_Y) /= 0) stop 'NY must be a multiple of nb_procs_y'

! check that a slice is at least as thick as a PML layer
  if (NX_LOCAL < NPOINTS_PML) stop 'NX_LOCAL must be greater than NPOINTS_PML' ! TODO add condition of PML
  if (NY_LOCAL < NPOINTS_PML) stop 'NY_LOCAL must be greater than NPOINTS_PML' ! TODO idem
 
  ! rank on line or column 
  i_rank = modulo(rank,NPROC_X)
  j_rank = rank /NPROC_X
   
  ! offset of this slice when we cut along Z
  offset_i = i_rank * NX_LOCAL
  offset_j = j_rank * NY_LOCAL
   
  ! we receive from the process on the left, and send to the process on the right
  sender_right_shift = rank - 1
  receiver_right_shift = rank + 1
  sender_bottom_shift = rank + NPROC_X
  receiver_bottom_shift = rank - NPROC_X
  sender_left_shift = rank + 1
  receiver_left_shift = rank - 1
  sender_top_shift = rank - NPROC_X
  receiver_top_shift = rank + NPROC_X
  
  sender_right_top_shift = rank - 1 - NPROC_X
  receiver_right_top_shift = rank + 1 + NPROC_X
  sender_right_bottom_shift = rank - 1 + NPROC_X
  receiver_right_bottom_shift = rank + 1 - NPROC_X
  sender_bottom_right_shift = rank - 1 + NPROC_X
  receiver_bottom_right_shift = rank + 1 - NPROC_X
  sender_bottom_left_shift = rank + 1 + NPROC_X
  receiver_bottom_left_shift = rank - 1 - NPROC_X
  sender_top_left_shift = rank + 1 - NPROC_X
  receiver_top_left_shift = rank - 1 + NPROC_X
    
! if we are the first process, there is no neighbor on the left
  if (USE_PML_XMIN) then
    if (i_rank == 0) sender_right_shift = MPI_PROC_NULL
    if (i_rank == 0) receiver_left_shift = MPI_PROC_NULL
        
    if (i_rank == 0) sender_right_top_shift = MPI_PROC_NULL
    if (i_rank == 0) sender_right_bottom_shift = MPI_PROC_NULL
    if (i_rank == 0) sender_bottom_right_shift = MPI_PROC_NULL
    if (i_rank == 0) receiver_bottom_left_shift = MPI_PROC_NULL
    if (i_rank == 0) receiver_top_left_shift = MPI_PROC_NULL
  else
    if (i_rank == 0) sender_right_shift = NPROC_X * (j_rank+1) - 1
    if (i_rank == 0) receiver_left_shift = NPROC_X * (j_rank+1) - 1
    
    if (i_rank == 0) sender_right_top_shift = rank - 1
    if (i_rank == 0) sender_right_bottom_shift = rank + 2*NPROC_X - 1
    if (i_rank == 0) sender_bottom_right_shift = rank + 2*NPROC_X - 1
    if (i_rank == 0) receiver_bottom_left_shift = rank - 1
    if (i_rank == 0) receiver_top_left_shift = rank + 2*NPROC_X - 1
  endif
  
! if we are the last process, there is no neighbor on the right
  if (USE_PML_XMAX) then
    if (i_rank == NPROC_X - 1) receiver_right_shift = MPI_PROC_NULL
    if (i_rank == NPROC_X - 1) sender_left_shift = MPI_PROC_NULL   
    
    if (i_rank == NPROC_X - 1) sender_bottom_left_shift = MPI_PROC_NULL
    if (i_rank == NPROC_X - 1) receiver_right_top_shift = MPI_PROC_NULL
    if (i_rank == NPROC_X - 1) receiver_right_bottom_shift = MPI_PROC_NULL
    if (i_rank == NPROC_X - 1) receiver_bottom_right_shift = MPI_PROC_NULL
    if (i_rank == NPROC_X - 1) sender_top_left_shift = MPI_PROC_NULL
  else
    if (i_rank == NPROC_X - 1) receiver_right_shift = rank - i_rank
    if (i_rank == NPROC_X - 1) sender_left_shift = rank - i_rank
    
    if (i_rank == NPROC_X - 1) sender_bottom_left_shift = rank + 1
    if (i_rank == NPROC_X - 1) receiver_right_top_shift = rank + 1 
    if (i_rank == NPROC_X - 1) receiver_right_bottom_shift = rank - 2*NPROC_X + 1 
    if (i_rank == NPROC_X - 1) receiver_bottom_right_shift = rank - 2*NPROC_X + 1 
    if (i_rank == NPROC_X - 1) sender_top_left_shift = rank - 2*NPROC_X + 1
  endif

! if we are the first process, there is no neighbor on the left
  if (USE_PML_YMIN) then
    if (j_rank == 0) sender_top_shift = MPI_PROC_NULL
    if (j_rank == 0) receiver_bottom_shift = MPI_PROC_NULL
    
    if (j_rank == 0) sender_right_top_shift = MPI_PROC_NULL
    if (j_rank == 0) receiver_right_bottom_shift = MPI_PROC_NULL
    if (j_rank == 0) receiver_bottom_right_shift = MPI_PROC_NULL
    if (j_rank == 0) receiver_bottom_left_shift = MPI_PROC_NULL
    if (j_rank == 0) sender_top_left_shift = MPI_PROC_NULL
  else
    if (j_rank == 0) sender_top_shift = rank + (NPROC_X)
    if (j_rank == 0) receiver_bottom_shift = rank + (NPROC_X)
    
    if (j_rank == 0) sender_right_top_shift = i_rank - 1 + NPROC_X * (NPROC_Y - 1)
    if (j_rank == 0) receiver_right_bottom_shift = rank + 1 + NPROC_X * (NPROC_Y - 1) 
    if (j_rank == 0) receiver_bottom_right_shift = rank + 1 + NPROC_X * (NPROC_Y - 1)
    if (j_rank == 0) receiver_bottom_left_shift = i_rank + NPROC_X * (NPROC_Y - 1) -1
    if (j_rank == 0) sender_top_left_shift = i_rank + (NPROC_Y-1) * NPROC_X + 1 
  endif
  
! if we are the last process, there is no neighbor on the right
  if (USE_PML_YMAX) then
    if (j_rank == NPROC_Y - 1) receiver_top_shift = MPI_PROC_NULL
    if (j_rank == NPROC_Y - 1) sender_bottom_shift = MPI_PROC_NULL
    
    if (j_rank == NPROC_Y - 1) receiver_right_top_shift = MPI_PROC_NULL
    if (j_rank == NPROC_Y - 1) sender_bottom_left_shift = MPI_PROC_NULL
    if (j_rank == NPROC_Y - 1) sender_bottom_right_shift = MPI_PROC_NULL
    if (j_rank == NPROC_Y - 1) sender_right_bottom_shift = MPI_PROC_NULL
    if (j_rank == NPROC_Y - 1) receiver_top_left_shift = MPI_PROC_NULL
  else
    if (j_rank == NPROC_Y - 1) receiver_top_shift = rank - (NPROC_X)
    if (j_rank == NPROC_Y - 1) sender_bottom_shift = rank - (NPROC_X)
    
    if (j_rank == NPROC_Y - 1) receiver_right_top_shift = i_rank + 1
    if (j_rank == NPROC_Y - 1) sender_bottom_left_shift = rank - (NPROC_Y - 1)*NPROC_X + 1
    if (j_rank == NPROC_Y - 1) sender_bottom_right_shift = rank - (NPROC_Y - 1)*NPROC_X - 1
    if (j_rank == NPROC_Y - 1) sender_right_bottom_shift = rank - (NPROC_Y - 1)*NPROC_X - 1
    if (j_rank == NPROC_Y - 1) receiver_top_left_shift = i_rank - 1
  endif

  
  if (.not.(USE_PML_YMAX) .and. .not.(USE_PML_XMAX)) then
    if (i_rank == NPROC_X-1 .and. j_rank == NPROC_Y - 1) then
       sender_bottom_left_shift = 0
       receiver_right_top_shift = 0
    endif
  endif
  
  if (USE_PML_YMAX .or. USE_PML_XMAX) then
    if (i_rank == NPROC_X-1 .and. j_rank == NPROC_Y - 1) then
       sender_bottom_left_shift = MPI_PROC_NULL
       receiver_right_top_shift = MPI_PROC_NULL
    endif
  endif
  
 
  if (.not.(USE_PML_YMAX) .and. .not.(USE_PML_XMIN)) then
      if (i_rank == 0 .and. j_rank == NPROC_Y - 1) then
         sender_bottom_right_shift = NPROC_X-1   
         sender_right_bottom_shift = NPROC_X-1
         receiver_top_left_shift = NPROC_X-1
      endif   
  endif
  
    if (USE_PML_YMAX .or. USE_PML_XMIN) then
      if (i_rank == 0 .and. j_rank == NPROC_Y - 1) then
         sender_bottom_right_shift = MPI_PROC_NULL 
         sender_right_bottom_shift = MPI_PROC_NULL
         receiver_top_left_shift = MPI_PROC_NULL
      endif   
  endif
  
  if (.not.(USE_PML_YMIN) .and. .not.(USE_PML_XMIN)) then
      if (i_rank == 0 .and. j_rank == 0) then
       receiver_bottom_left_shift = NPROC-1
       sender_right_top_shift = NPROC-1 
    endif
  endif
  
    if (USE_PML_YMIN .or. USE_PML_XMIN) then
      if (i_rank == 0 .and. j_rank == 0) then
       receiver_bottom_left_shift = MPI_PROC_NULL
       sender_right_top_shift = MPI_PROC_NULL
    endif
  endif
  
  if (.not.(USE_PML_YMIN) .and. .not.(USE_PML_XMAX)) then
      if (i_rank == NPROC_X-1 .and. j_rank == 0) then
         receiver_bottom_right_shift = NPROC_X*(NPROC_Y-1)   
         receiver_right_bottom_shift = NPROC_X*(NPROC_Y-1) 
         sender_top_left_shift = NPROC_X*(NPROC_Y-1) 
      endif
  endif
  
    if (USE_PML_YMIN .or. USE_PML_XMAX) then
      if (i_rank == NPROC_X-1 .and. j_rank == 0) then
         receiver_bottom_right_shift = MPI_PROC_NULL 
         receiver_right_bottom_shift = MPI_PROC_NULL
         sender_top_left_shift = MPI_PROC_NULL
      endif
  endif
  
  
  
  
  do rk=0,NPROC-1
  if (rk == rank) then
    print *, rank, sender_right_shift, receiver_right_shift
    print *, rank, sender_left_shift, receiver_left_shift
    print *, rank, sender_bottom_shift, receiver_bottom_shift    
    print *, rank, sender_top_shift, receiver_top_shift
    print *, rank, sender_right_bottom_shift, receiver_right_bottom_shift
    print *, rank, sender_bottom_right_shift, receiver_bottom_right_shift
    print *, rank, sender_right_top_shift, receiver_right_top_shift
    print *, rank, sender_bottom_left_shift, receiver_bottom_left_shift
    print *, rank, sender_top_left_shift, receiver_top_left_shift
    endif
    enddo
  

  
  i2begin = 1
  if (j_rank == 0) i2begin = 2 ! TODO modifier pour voir si necessaire la taille supp ?

  iminus1end = NX_LOCAL
  if (j_rank == NPROC_Y - 1) iminus1end = NX_LOCAL-1 ! TODO idem l.78
  
  j2begin = 1
  if (i_rank == 0) j2begin = 2 ! TODO modifier pour voir si necessaire la taille supp ?

  jminus1end = NY_LOCAL
  if (i_rank == NPROC_X - 1) jminus1end = NY_LOCAL-1 ! TODO idem l.78
  
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                       PML INIT                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call computePML()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                Background domain                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! compute the Lame parameter and density


!!! Compute the true model
  do j = 1,NY_LOCAL
    do i = 1,NX_LOCAL
    
      i_global = i + offset_i
      j_global = j + offset_j 
      
      if ((i_global < IObs_start) .or. (i_global >= IObs_end) .or. (j_global < JObs_start) .or. (j_global >= JObs_end)) then 
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
      if (add_wind_profile .and. (j_global > ymin_wind .and. j_global <= ymax_wind)) then 
        windy_true(i,j) = exp(- ((j_global-1)*0.1d0 - mean_gauss_wind)**2 / sigma2_gauss_wind ) * max_wind_factor
      endif
      
    enddo
  enddo

  call send_receive_rightleft(windx_true)
  call send_receive_rightleft(windy_true)
  call send_receive_rightleft(rho0_true)
  call send_receive_rightleft(p0_true)
     
  call send_receive_topbottom(windx_true)
  call send_receive_topbottom(windy_true)
  call send_receive_topbottom(rho0_true)
  call send_receive_topbottom(p0_true)
     
  call send_receive_corners(windx_true) 
  call send_receive_corners(windy_true)
  call send_receive_corners(rho0_true)

  ! write the background model
  call write_background(rho0_true, kappa_unrelaxed_true, p0_true, windx_true, windy_true, gamma_chimie)


!!! Compute the a priori model
  do j = 1,NY_LOCAL
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

  call send_receive_rightleft(windx_prior)
  call send_receive_rightleft(windy_prior)
  call send_receive_rightleft(rho0_prior)
  call send_receive_rightleft(p0_prior)
     
  call send_receive_topbottom(windx_prior)
  call send_receive_topbottom(windy_prior)
  call send_receive_topbottom(rho0_prior)
  call send_receive_topbottom(p0_prior)
     
  call send_receive_corners(windx_prior) 
  call send_receive_corners(windy_prior)
  call send_receive_corners(rho0_prior)
       
    ! write the background model
  call write_background(rho0_prior, kappa_unrelaxed_prior, p0_prior, windx_prior, windy_prior, gamma_chimie)

 
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
     
  call forwardproblem(p0_true, rho0_true, windx_true, windy_true, kappa_unrelaxed_true,  1, NSTEP, 1) 
  !call forwardproblem(p0_prior, rho0_prior, windx_prior, windy_prior, kappa_unrelaxed_prior,  1, NSTEP, 3) 
  sispressure_true(:,:) = sispressure(:,:)
  
  norm_pressure_true = maxval(abs(sispressure_true))**2
      call MPI_ALLREDUCE(maxval(abs(sispressure_true))**2, &
        norm_pressure_true,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,code)

        
  print *, "je suis vraiment sorti, j'ai fini", rank
  
 
    write(file_name, "('./OUTPUT/p_true_',i6.6,'.txt')") rank
    OPEN(UNIT=12, FILE=file_name, ACTION="write")
    DO ii=1,NX_LOCAL
      WRITE(12,*) (pressure(ii,jj), jj=1,NY_LOCAL)
    END DO
    CLOSE(12)

  
  

    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                  Compute gradient                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  call compute_kernel()
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                          End                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      ! close MPI program
  call MPI_FINALIZE(code)
 end program  main
