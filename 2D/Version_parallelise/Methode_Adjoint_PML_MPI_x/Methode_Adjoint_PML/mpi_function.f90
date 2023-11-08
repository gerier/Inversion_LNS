
subroutine read_parameters()

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DEFINE MPI MAIN PROPERTIES
  NX_LOCAL = NX / NPROC
  number_of_values = (NY+4)*2
  ! offset of this slice when we cut along X
  offset_i = rank * NX_LOCAL
  ! slice number for the cut plane in the middle of the mesh
  rank_cut_plane = nb_procs/2 - 1
  if(nb_procs == 1) rank_cut_plane = 0
	
  
  
  ! check that code was compiled with the right number of slices
  if(nb_procs /= NPROC) then
    print *,'nb_procs,NPROC = ',nb_procs,NPROC
    stop 'nb_procs must be equal to NPROC'
  endif
  
  if(mod(NX,nb_procs) /= 0) then
         print *,'NX must be a multiple of nb_procs' 
        stop 'NX must be a multiple of nb_procs'
  endif
  
  ! check that a slice is at least as thick as a PML layer
  if(NX_LOCAL < NPOINTS_PML .AND. (USE_PML_XMIN .OR. USE_PML_XMAX) ) then
        print *,'NX_LOCAL must be greater than NPOINTS_PML'
        stop 'NX_LOCAL must be greater than NPOINTS_PML'
  endif
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! #ALLOCATIONS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TODO
  
  
end subroutine read_parameters





subroutine init_MPI_neighbor()
  
  use MPI
  
 
  implicit none
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DEFINE MPI NEIGHBORS FROM BOUNDARY COND.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! we receive from the process on the left, and send to the process on the right
  sender_right_shift = rank - 1
  receiver_right_shift = rank + 1

  ! if we are the first process, there is no neighbor on the left
  if(rank == 0) then
        sender_right_shift = MPI_PROC_NULL
        if(ADD_PERIODIC_BOUNDARY_X) sender_right_shift = nb_procs - 1
  endif

  ! if we are the last process, there is no neighbor on the right
  if(rank == nb_procs - 1) then
        receiver_right_shift = MPI_PROC_NULL
        if(ADD_PERIODIC_BOUNDARY_X) receiver_right_shift = 0
  endif

!---

! we receive from the process on the right, and send to the process on the left
  sender_left_shift = rank + 1
  receiver_left_shift = rank - 1

! if we are the first process, there is no neighbor on the left
  if(rank == 0) then
        receiver_left_shift = MPI_PROC_NULL
        if(ADD_PERIODIC_BOUNDARY_X) receiver_left_shift = nb_procs - 1
  endif
  

! if we are the last process, there is no neighbor on the right
  if(rank == nb_procs - 1)  then
        sender_left_shift = MPI_PROC_NULL
        if(ADD_PERIODIC_BOUNDARY_X) sender_left_shift = 0
  endif

  ! If periodic boundary conditions or not last/first process (right/left edges)
  ! => each computation along x is made from 1 to NX_LOCAL
  i2begin = 1
  if(rank == 0 .AND. .not. ADD_PERIODIC_BOUNDARY_X) i2begin = 2

  iminus1end = NX_LOCAL
  if(rank == nb_procs - 1 .AND. .not. ADD_PERIODIC_BOUNDARY_X) iminus1end = NX_LOCAL - 1

  end subroutine init_MPI_neighbor
  
  
  
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! GENERATE SEISMOGRAMS FILES
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_source()
  
  use AGW_par_2D,only: rank_source, rank, ISOURCE, NX_LOCAL, offset_i, &
        rad_source, DELTAX
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DETERMINE IF CURRENT PROCESS CONTAINS THE SOURCE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rank_source = -1
  if( (ISOURCE <= NX_LOCAL+offset_i &
        .AND. ISOURCE > offset_i) &
        .OR. (rad_source > 0 .AND. &
              abs(ISOURCE - offset_i) <= rad_source/DELTAX .OR. &
              abs(ISOURCE - (NX_LOCAL+offset_i)) <= rad_source/DELTAX &
              ) ) &
        rank_source = rank
        
  end subroutine init_source
  
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! SAVE WAVEFIELDS FOR SEISMOGRAMS 
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine store_sismograms(sisp, sisrho, sisAy, sisvy)
  
  use AGW_par_2D,only: it, nb_recs, nb_recs_local, recs_local_to_global, &
        i_interp, j_interp, rank,DELTAX, DELTAY, NX_LOCAL, NY,&
        xrec, yrec, SAVE_SISMO_AT_THE_END, INTERPOLATE_POS_STATION, &
        DELTAT, NSTEP, offset_i, &
        p, rho1, rho,ay, uy, dvy,velocity, ix_rec, iy_rec, &
        rank, rank_source, ISOURCE, JSOURCE
  
  implicit none
  
  ! Local variables
  integer :: irec_local, irec
  double precision :: temp_interp
  
  !double precision, dimension(-1:NX_LOCAL+2,-1:NY+2) :: &
  !      uy, vx, vy, ay, p, rho1
  double precision, dimension(1:NSTEP,1:nb_recs_local) :: &
        sisp, sisrho, sisAy, sisvy    

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! STORE SISMOGRAMS
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do irec_local = 1,nb_recs_local
    irec = recs_local_to_global(irec_local)
    
    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Pressure interpolation
    temp_interp = p(ix_rec(irec)-offset_i,iy_rec(irec))
    
    if(INTERPOLATE_POS_STATION) then
        temp_interp = 0d0
        call lagrange_interp_2d ( 1, 1, i_interp(:,irec_local), j_interp(:,irec_local),rank,DELTAX, DELTAY, NX_LOCAL, NY, p,&
                         1, xrec(irec), yrec(irec), temp_interp )
    endif
        
    sisp(it,irec_local)  = temp_interp
    
    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Density interpolation
    !temp_interp = rho1(ix_rec(irec)-offset_i,iy_rec(irec))
    ! Test value pressure
    temp_interp = rho(ix_rec(irec)-offset_i,iy_rec(irec))*velocity(ix_rec(irec)-offset_i,iy_rec(irec))*&
        dvy(1,ix_rec(irec)-offset_i,1)
    
    if(INTERPOLATE_POS_STATION) then
        temp_interp = 0d0
        call lagrange_interp_2d ( 1, 1, i_interp(:,irec_local), j_interp(:,irec_local),rank,DELTAX, DELTAY, NX_LOCAL, NY, rho1,&
                         1, xrec(irec), yrec(irec), temp_interp )
    endif
    
    sisrho(it,irec_local) = temp_interp
    
    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Pressure interpolation
    temp_interp = ay(ix_rec(irec)-offset_i,iy_rec(irec))
    
    if(INTERPOLATE_POS_STATION) then
        temp_interp = 0d0
        call lagrange_interp_2d ( 1, 1, i_interp(:,irec_local), j_interp(:,irec_local),rank,DELTAX, DELTAY, NX_LOCAL, NY, ay,&
                         1, xrec(irec), yrec(irec), temp_interp )
    endif
    
    sisay(it,irec_local) = temp_interp
    
    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Vertical velocity interpolation
    temp_interp = dvy(1,ix_rec(irec)-offset_i,iy_rec(irec))!uy(ix_rec(irec)-offset_i,iy_rec(irec))
    
    if(INTERPOLATE_POS_STATION) then
        temp_interp = 0d0
        call lagrange_interp_2d ( 1, 1, i_interp(:,irec_local), j_interp(:,irec_local),rank,DELTAX, DELTAY, NX_LOCAL, NY, uy,&
                         1, xrec(irec), yrec(irec), temp_interp )
    endif
        
    sisvy(it,irec_local) = temp_interp
    
  enddo
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! #SAVE SISMOGRAM IF STEP BY STEP SAVE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(nb_recs_local > 0 .AND. .not. SAVE_SISMO_AT_THE_END) &
    call write_seismograms(it, 1, sisrho(it,:), sisay(it,:), sisvy(it,:), sisp(it,:))
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! #SAVE SISMOGRAMS IN FILE IF AT THE END
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(nb_recs_local > 0 .AND. SAVE_SISMO_AT_THE_END .AND. it >= NSTEP) &
    call write_seismograms(it, NSTEP, sisrho(:,:), sisay(:,:), sisvy(:,:), sisp(:,:))
    
  end subroutine store_sismograms
  

