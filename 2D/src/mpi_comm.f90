!!!
! MPI communication routines for 2D domain decomposition.
! Each process exchanges halo cells with its neighboring processes.
!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                  MPI Communications
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine send_receive_rightleft(var)
  ! Two MPI communications:
  !   Send halo region to the left neighbor and receive from the right neighbor
  !   Send halo region to the right neighbor and receive from the left neighbor
  !  (bidirectional data exchange between left and right neighbors)
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_x, &
                     receiver_left_shift,sender_left_shift,receiver_right_shift,sender_right_shift,&
                     message_tag, message_status, code
  implicit none
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: var

  call MPI_SENDRECV(var(1:2,:),number_of_values_x,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,var(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values_x, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(var(NX_LOCAL-1:NX_LOCAL,:),number_of_values_x,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,var(-1:0,:),number_of_values_x, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)

endsubroutine send_receive_rightleft


subroutine send_receive_left(var)
  ! One MPI communication:
  ! Send halo region to the left neighbor and receive from the right neighbor
  ! (information transfer from right to left)
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_x, &
                     receiver_left_shift,sender_left_shift,&
                     message_tag, message_status, code
  implicit none

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: var

call MPI_SENDRECV(var(1:2,:),number_of_values_x,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,var(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values_x, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)

endsubroutine send_receive_left


subroutine send_receive_right(var)
  ! One MPI communication:
  ! Send halo region to the right neighbor and receive from the left neighbor
  ! (information transfer from left to right)
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_x, &
                     receiver_right_shift,sender_right_shift,&
                     message_tag, message_status, code
  implicit none

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: var

call MPI_SENDRECV(var(NX_LOCAL-1:NX_LOCAL,:),number_of_values_x,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,var(-1:0,:),number_of_values_x, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)

endsubroutine send_receive_right


subroutine send_receive_topbottom(var)
  ! Two MPI communications:
  !   Send halo region to the bottom neighbor and receive from the top neighbor
  !   Send halo region to the top neighbor and receive from the bottom neighbor
  ! (information transfer top <=> bottom)
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_y, &
                     receiver_top_shift,sender_top_shift,receiver_bottom_shift,sender_bottom_shift,&
                     message_tag, message_status, code
  implicit none

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: var

   call MPI_SENDRECV(var(:,NY_LOCAL-1:NY_LOCAL),number_of_values_y,MPI_DOUBLE_PRECISION, &
         receiver_top_shift,message_tag,var(:,-1:0),number_of_values_y, &
         MPI_DOUBLE_PRECISION,sender_top_shift,message_tag,MPI_COMM_WORLD,message_status,code)
   call MPI_SENDRECV(var(:,1:2),number_of_values_y,MPI_DOUBLE_PRECISION, &
         receiver_bottom_shift,message_tag,var(:,NY_LOCAL+1:NY_LOCAL+2),number_of_values_y, &
         MPI_DOUBLE_PRECISION,sender_bottom_shift,message_tag,MPI_COMM_WORLD,message_status,code)

endsubroutine send_receive_topbottom


subroutine send_receive_bottom(var)
  ! One MPI communication:
  ! Send halo region to the bottom neighbor and receive from the top neighbor
  ! (information transfer from top to bottom)
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_y, &
                     receiver_bottom_shift,sender_bottom_shift,&
                     message_tag, message_status, code
  implicit none

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: var

   call MPI_SENDRECV(var(:,1:2),number_of_values_y,MPI_DOUBLE_PRECISION, &
         receiver_bottom_shift,message_tag,var(:,NY_LOCAL+1:NY_LOCAL+2),number_of_values_y, &
         MPI_DOUBLE_PRECISION,sender_bottom_shift,message_tag,MPI_COMM_WORLD,message_status,code)

endsubroutine send_receive_bottom


subroutine send_receive_top(var)
  ! One MPI communication:
  ! Send halo region to the top neighbor and receive from the bottom neighbor
  ! (information transfer from bottom to top)
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_y, &
                     receiver_top_shift,sender_top_shift,&
                     message_tag, message_status, code
  implicit none

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: var

   call MPI_SENDRECV(var(:,NY_LOCAL-1:NY_LOCAL),number_of_values_y,MPI_DOUBLE_PRECISION, &
         receiver_top_shift,message_tag,var(:,-1:0),number_of_values_y, &
         MPI_DOUBLE_PRECISION,sender_top_shift,message_tag,MPI_COMM_WORLD,message_status,code)


endsubroutine send_receive_top


subroutine send_receive_righttop(var)
  ! One MPI communication:
  ! Send halo region to the right/top neighbor and receive from the left/bottom neighbor
  ! (diagonal information transfer from left/bottom to right/top)
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_corner, &
                     receiver_right_top_shift,sender_right_top_shift,&
                     message_tag, message_status, code
  implicit none

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: var

    call MPI_SENDRECV(var(NX_LOCAL-1:NX_LOCAL,NY_LOCAL-1:NY_LOCAL),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_right_top_shift,message_tag,var(-1:0,-1:0),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_right_top_shift,message_tag,MPI_COMM_WORLD,message_status,code)

endsubroutine send_receive_righttop


subroutine send_receive_leftbottom(var)
  ! One MPI communication:
  ! Send halo region to the left/bottom neighbor and receive from the right/top neighbor
  ! (diagonal information transfer from right/top to left/bottom)
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_corner, &
                     receiver_bottom_left_shift,sender_bottom_left_shift,&
                     message_tag, message_status, code
  implicit none

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: var

    call MPI_SENDRECV(var(1:2,1:2),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_bottom_left_shift,message_tag,var(NX_LOCAL+1:NX_LOCAL+2,NY_LOCAL+1:NY_LOCAL+2),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_bottom_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)

endsubroutine send_receive_leftbottom


subroutine send_receive_rightbottom(var)
  ! One MPI communication:
  ! Send halo region to the right/bottom neighbor and receive from the left/top neighbor
  ! (diagonal information transfer from left/top to right/bottom)
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_corner, &
                     receiver_right_bottom_shift,sender_right_bottom_shift,&
                     message_tag, message_status, code
  implicit none

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: var

    call MPI_SENDRECV(var(NX_LOCAL-1:NX_LOCAL,1:2),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_right_bottom_shift,message_tag,var(-1:0,NY_LOCAL+1:NY_LOCAL+2),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_right_bottom_shift,message_tag,MPI_COMM_WORLD,message_status,code)

endsubroutine send_receive_rightbottom


subroutine send_receive_lefttop(var)
  ! One MPI communication:
  ! Send halo region to the left/top neighbor and receive from the right/bottom neighbor
  ! (diagonal information transfer from right/bottom to left/top)
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_corner, &
                     receiver_top_left_shift,sender_top_left_shift,&
                     message_tag, message_status, code
  implicit none

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: var

    call MPI_SENDRECV(var(1:2,NY_LOCAL-1:NY_LOCAL),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_top_left_shift,message_tag,var(NX_LOCAL+1:NX_LOCAL+2,-1:0),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_top_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)

endsubroutine send_receive_lefttop


subroutine send_receive_corners(var)
  ! Four MPI communications:
  !   - diagonal information transfer from right/top to left/bottom
  !   - diagonal information transfer from left/top to right/bottom
  !   - diagonal information transfer from left/bottom to right/top 
  !   - diagonal information transfer from right/bottom to left/top
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_corner, &
                     receiver_right_top_shift,sender_right_top_shift,&
                     receiver_bottom_left_shift,sender_bottom_left_shift,&
                     receiver_top_left_shift,sender_top_left_shift,&
                     receiver_bottom_right_shift,sender_bottom_right_shift,&
                     message_tag, message_status,code
  implicit none

  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: var

    call MPI_SENDRECV(var(1:2,1:2),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_bottom_left_shift,message_tag,var(NX_LOCAL+1:NX_LOCAL+2,NY_LOCAL+1:NY_LOCAL+2),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_bottom_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)

    call MPI_SENDRECV(var(NX_LOCAL-1:NX_LOCAL,1:2),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_bottom_right_shift,message_tag,var(-1:0,NY_LOCAL+1:NY_LOCAL+2),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_bottom_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)

    call MPI_SENDRECV(var(NX_LOCAL-1:NX_LOCAL,NY_LOCAL-1:NY_LOCAL),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_right_top_shift,message_tag,var(-1:0,-1:0),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_right_top_shift,message_tag,MPI_COMM_WORLD,message_status,code)

    call MPI_SENDRECV(var(1:2,NY_LOCAL-1:NY_LOCAL),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_top_left_shift,message_tag,var(NX_LOCAL+1:NX_LOCAL+2,-1:0),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_top_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)

endsubroutine send_receive_corners


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                  Definition of MPI Neighbor Processes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_neighbors()
!==============================================================================
!> Initialize the MPI Cartesian neighborhood information.
!>
!> Checks the domain decomposition, computes the local subdomain offsets,
!> assigns the neighboring MPI ranks used for halo exchanges, and accounts
!> for the selected boundary conditions (periodic or PML).
!==============================================================================
use MPI
use parameters, only : nb_procs, NPROC, NPROC_X, NPROC_Y, NX, NY, NX_LOCAL, NY_LOCAL,       &
                       offset_i, offset_j, &
                       USE_PML_XMIN, USE_PML_XMAX, USE_PML_YMIN, USE_PML_YMAX, NPOINTS_PML, &
                       rank, i_rank, j_rank, code, ierr, row_comm, &
                       receiver_bottom_left_shift, receiver_bottom_right_shift, receiver_bottom_shift, &
                       receiver_left_shift, receiver_right_bottom_shift, receiver_right_shift, &
                       receiver_right_top_shift, receiver_top_left_shift, receiver_top_shift, &
                       sender_bottom_left_shift, sender_bottom_right_shift, sender_bottom_shift, &
                       sender_left_shift, sender_right_bottom_shift, sender_right_shift, &
                       sender_right_top_shift, sender_top_left_shift, sender_top_shift
implicit none

! Check that the code was compiled with the correct number of MPI processes
  if (nb_procs /= NPROC) then
    print *,'incorrect number of MPI processes: nb_procs,NPROC = ',nb_procs,NPROC,' but they should be equal'
    stop 'nb_procs must be equal to NPROC'
  endif

! Check that the domain can be decomposed into an exact number of slices
! along the X and Y directions
  if (mod(NX,NPROC_X) /= 0) stop 'NX must be a multiple of NPROC_X'
  if (mod(NY,NPROC_Y) /= 0) stop 'NY must be a multiple of NPROC_Y'

! Check that each local domain is at least as large as the PML thickness
  if (NX_LOCAL <= NPOINTS_PML .and. USE_PML_XMIN) stop 'NX_LOCAL must be larger than NPOINTS_PML'
  if (NX_LOCAL <= NPOINTS_PML .and. USE_PML_XMAX) stop 'NX_LOCAL must be larger than NPOINTS_PML'
  if (NY_LOCAL <= NPOINTS_PML .and. USE_PML_YMIN) stop 'NY_LOCAL must be larger than NPOINTS_PML'
  if (NY_LOCAL <= NPOINTS_PML .and. USE_PML_YMAX) stop 'NY_LOCAL must be larger than NPOINTS_PML'

! Rank coordinates in the process grid
!
! Example:
! Horizontal direction: x (i_rank)
! Vertical direction:   y (j_rank)
!
! ---------------
! |rank 2 |rank 3 |
! | (0,1) | (1,1) |
! ---------------
! |rank 0 |rank 1 |
! | (0,0) | (1,0) |
! ---------------
  i_rank = modulo(rank,NPROC_X)
  j_rank = rank /NPROC_X

  ! Split processes into row-based communicators
  call MPI_BARRIER(mpi_comm_world,code)
  call MPI_Comm_Split(MPI_COMM_WORLD, j_rank, rank, row_Comm,ierr)

  ! Offset of the current local domain in the global grid
  offset_i = i_rank * NX_LOCAL
  offset_j = j_rank * NY_LOCAL

  ! Define neighbors for communication along the horizontal and vertical directions
  sender_right_shift = rank - 1
  receiver_right_shift = rank + 1
  sender_bottom_shift = rank + NPROC_X
  receiver_bottom_shift = rank - NPROC_X
  sender_left_shift = rank + 1
  receiver_left_shift = rank - 1
  sender_top_shift = rank - NPROC_X
  receiver_top_shift = rank + NPROC_X

  ! Define neighbors for diagonal communications
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

! If this is the first column process, define neighbors according to boundary conditions
  if (USE_PML_XMIN .or. USE_PML_XMAX) then
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

! If this is the last column process, define neighbors according to boundary conditions
  if (USE_PML_XMAX .or. USE_PML_XMIN) then
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

! If this is the first row process, define neighbors according to boundary conditions
  if (USE_PML_YMIN .or. USE_PML_YMAX ) then
    if (j_rank == 0) sender_top_shift = MPI_PROC_NULL
    if (j_rank == 0) receiver_bottom_shift = MPI_PROC_NULL

    if (j_rank == 0) sender_right_top_shift = MPI_PROC_NULL
    if (j_rank == 0) receiver_right_bottom_shift = MPI_PROC_NULL
    if (j_rank == 0) receiver_bottom_right_shift = MPI_PROC_NULL
    if (j_rank == 0) receiver_bottom_left_shift = MPI_PROC_NULL
    if (j_rank == 0) sender_top_left_shift = MPI_PROC_NULL
  else
    if (j_rank == 0) sender_top_shift = i_rank + NPROC_X * (NPROC_Y - 1)
    if (j_rank == 0) receiver_bottom_shift = i_rank + NPROC_X * (NPROC_Y - 1)

    if (j_rank == 0) sender_right_top_shift = i_rank - 1 + NPROC_X * (NPROC_Y - 1)
    if (j_rank == 0) receiver_right_bottom_shift = rank + 1 + NPROC_X * (NPROC_Y - 1)
    if (j_rank == 0) receiver_bottom_right_shift = rank + 1 + NPROC_X * (NPROC_Y - 1)
    if (j_rank == 0) receiver_bottom_left_shift = i_rank + NPROC_X * (NPROC_Y - 1) -1
    if (j_rank == 0) sender_top_left_shift = i_rank + (NPROC_Y-1) * NPROC_X + 1
  endif

! If this is the last row process, define neighbors according to boundary conditions
  if (USE_PML_YMAX .or. USE_PML_YMIN) then
    if (j_rank == NPROC_Y - 1) receiver_top_shift = MPI_PROC_NULL
    if (j_rank == NPROC_Y - 1) sender_bottom_shift = MPI_PROC_NULL

    if (j_rank == NPROC_Y - 1) receiver_right_top_shift = MPI_PROC_NULL
    if (j_rank == NPROC_Y - 1) sender_bottom_left_shift = MPI_PROC_NULL
    if (j_rank == NPROC_Y - 1) sender_bottom_right_shift = MPI_PROC_NULL
    if (j_rank == NPROC_Y - 1) sender_right_bottom_shift = MPI_PROC_NULL
    if (j_rank == NPROC_Y - 1) receiver_top_left_shift = MPI_PROC_NULL
  else
    if (j_rank == NPROC_Y - 1) receiver_top_shift = i_rank
    if (j_rank == NPROC_Y - 1) sender_bottom_shift = i_rank

    if (j_rank == NPROC_Y - 1) receiver_right_top_shift = i_rank + 1
    if (j_rank == NPROC_Y - 1) sender_bottom_left_shift = rank - (NPROC_Y - 1)*NPROC_X + 1
    if (j_rank == NPROC_Y - 1) sender_bottom_right_shift = rank - (NPROC_Y - 1)*NPROC_X - 1
    if (j_rank == NPROC_Y - 1) sender_right_bottom_shift = rank - (NPROC_Y - 1)*NPROC_X - 1
    if (j_rank == NPROC_Y - 1) receiver_top_left_shift = i_rank - 1
  endif

! If a process is located at a domain corner,
! define diagonal neighbors according to boundary conditions

  ! if process on the right, top corner, different neighbors depending on boundary conditions
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

  ! if process on the left, top corner, different neighbors depending on boundary conditions
  if (i_rank == 0 .and. j_rank == NPROC_Y - 1) then
     if (.not.(USE_PML_YMAX) .and. .not.(USE_PML_XMIN) .and. .not.(USE_PML_YMIN) .and. .not.(USE_PML_XMAX)) then ! All periodics
         sender_bottom_right_shift = NPROC_X-1
         sender_right_bottom_shift = NPROC_X-1
         receiver_top_left_shift = NPROC_X-1
      else
         sender_bottom_right_shift = MPI_PROC_NULL
         sender_right_bottom_shift = MPI_PROC_NULL
         receiver_top_left_shift = MPI_PROC_NULL
      endif
  endif

  ! if process on the left, bottom corner, different neighbors depending on boundary conditions
  if (i_rank == 0 .and. j_rank == 0) then
      if (.not.(USE_PML_YMIN) .and. .not.(USE_PML_XMIN) .and. .not.(USE_PML_YMAX) .and. .not.(USE_PML_XMAX) ) then
       receiver_bottom_left_shift = NPROC-1
       sender_right_top_shift = NPROC-1
      else
       receiver_bottom_left_shift = MPI_PROC_NULL
       sender_right_top_shift = MPI_PROC_NULL
    endif
  endif

  ! if process on the right, top corner, different neighbors depending on boundary conditions
  if (i_rank == NPROC_X-1 .and. j_rank == 0) then
    if (.not.(USE_PML_YMIN) .and. .not.(USE_PML_XMAX)  .and. .not.(USE_PML_YMAX) .and. .not.(USE_PML_XMAX) ) then
         receiver_bottom_right_shift = NPROC_X*(NPROC_Y-1)
         receiver_right_bottom_shift = NPROC_X*(NPROC_Y-1)
         sender_top_left_shift = NPROC_X*(NPROC_Y-1)
     else
         receiver_bottom_right_shift = MPI_PROC_NULL
         receiver_right_bottom_shift = MPI_PROC_NULL
         sender_top_left_shift = MPI_PROC_NULL
      endif
  endif


endsubroutine get_neighbors
