

program main

 use MPI
 use parameters
 implicit none
 
 integer i,j,irec, ii, jj
 double precision :: Courant_number

double precision, dimension(1:NX,1:NY) :: pressure_global
integer :: rk, ii_rk,jj_rk
character(len=100) :: file_name

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
    call execute_command_line ('mkdir '//'MODELS/')
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                       PML INIT                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call computePML()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                Background domain                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call init_backgrounds()
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          source / receivers INFORMATION          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!!             Retrieve the observations                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  ! INITIALISATION
!  call reset_forward()
!  save_sismos = .True.
     
!  call forwardproblem(p0_true, rho0_true, windx_true, windy_true, kappa_unrelaxed_true,  1, NSTEP, 1) 
!  !call forwardproblem(p0_prior, rho0_prior, windx_prior, windy_prior, kappa_unrelaxed_prior,  1, NSTEP, 3) 
!  sispressure_true(:,:) = sispressure(:,:)
  
!  norm_pressure_true = maxval(abs(sispressure_true))**2
!      call MPI_ALLREDUCE(maxval(abs(sispressure_true))**2, &
!        norm_pressure_true,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,code)

        
!  print *, "je suis vraiment sorti, j'ai fini", rank
  
 
!    write(file_name, "('./OUTPUT/p_true_',i6.6,'.txt')") rank
!    OPEN(UNIT=12, FILE=file_name, ACTION="write")
!    DO ii=1,NX_LOCAL
!      WRITE(12,*) (pressure(ii,jj), jj=1,NY_LOCAL)
!    END DO
!    CLOSE(12)

    call reset_forward()
  save_sismos = .True.
  
  !call forwardproblem(p0_true, rho0_true, windx_true, windy_true, kappa_unrelaxed_true,  1, NSTEP, 2) 
  call forwardproblem(p0_prior, rho0_prior,windx_prior,windy_prior,kappa_unrelaxed_prior,  1, NSTEP, 2) 
  sispressure_true(:,:) = sispressure(:,:)
  norm_pressure_true = maxval(abs(sispressure_true))**2
  
  call reset_kernel()
  save_sismos = .False.
  
  call save_frames()
   do it=0,NSTEP
   
   if (NSTEP-it == NSTEP-1 .or. modulo(NSTEP-it,NSTEP/NFRAMES) == NSTEP/NFRAMES-1 ) then
      !print *, 1,it, NSTEP-it
      call save_local_frames(NSTEP-it)
   else if (NSTEP-it /= NSTEP) then
     
     call load_frame(NSTEP-it, time_last_frame)
     !print *, 2,it, NSTEP-it, time_last_frame
     call forwardproblem(p0_prior,rho0_prior,windx_prior,windy_prior,kappa_unrelaxed_prior,time_last_frame, NSTEP-it,2)

   else
     !print *, 3, it, NSTEP-it
   endif

    if (i_rank == ix_rec(1)/NX_LOCAL .and. j_rank == iy_rec(1)/NY_LOCAL) then
   OPEN(UNIT=1222, FILE="./OUTPUT/p_checkpoint.txt", position="append", ACTION="write")
   WRITE(1222,*) it, pressure(ix_rec(1)-offset_i,iy_rec(1)-offset_j)
   !print *, pressure(ix_rec(1),ix_rec(1)), ix_rec(1),ix_rec(1)
   CLOSE(1222)
   endif
    
   enddo
  

    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                  Compute gradient                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!  call compute_kernel()
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                          End                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      ! close MPI program
  call MPI_FINALIZE(code)
 end program  main
