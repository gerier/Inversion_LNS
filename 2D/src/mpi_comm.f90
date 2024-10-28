!!!
! Functions relative to software parallelisation through MPI 
!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                  Communications MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine send_receive_rightleft(var)
  ! 2 communications
  !   from right to left 
  !   from left to right 
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_x, &
                     receiver_left_shift,sender_left_shift,receiver_right_shift,sender_right_shift,&
                     message_tag, message_status, code
  implicit none
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
  call MPI_SENDRECV(var(1:2,:),number_of_values_x,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,var(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values_x, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
  call MPI_SENDRECV(var(NX_LOCAL-1:NX_LOCAL,:),number_of_values_x,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,var(-1:0,:),number_of_values_x, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
     
endsubroutine send_receive_rightleft


subroutine send_receive_left(var)
  ! 1 communication: 
  ! direction of information: to the left
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_x, &
                     receiver_left_shift,sender_left_shift,&
                     message_tag, message_status, code
  implicit none
  
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
call MPI_SENDRECV(var(1:2,:),number_of_values_x,MPI_DOUBLE_PRECISION, &
         receiver_left_shift,message_tag,var(NX_LOCAL+1:NX_LOCAL+2,:),number_of_values_x, &
         MPI_DOUBLE_PRECISION,sender_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
         
endsubroutine send_receive_left


subroutine send_receive_right(var)
  ! 1 communication: 
  ! direction of information: to the right
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_x, &
                     receiver_right_shift,sender_right_shift,&
                     message_tag, message_status, code
  implicit none
  
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
call MPI_SENDRECV(var(NX_LOCAL-1:NX_LOCAL,:),number_of_values_x,MPI_DOUBLE_PRECISION, &
         receiver_right_shift,message_tag,var(-1:0,:),number_of_values_x, &
         MPI_DOUBLE_PRECISION,sender_right_shift,message_tag,MPI_COMM_WORLD,message_status,code)
     
endsubroutine send_receive_right


subroutine send_receive_topbottom(var)
  ! 2 communications
  !   from top to bottom 
  !   from bottom to top 
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_y, &
                     receiver_top_shift,sender_top_shift,receiver_bottom_shift,sender_bottom_shift,&
                     message_tag, message_status, code
  implicit none
  
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
   call MPI_SENDRECV(var(:,NY_LOCAL-1:NY_LOCAL),number_of_values_y,MPI_DOUBLE_PRECISION, &
         receiver_top_shift,message_tag,var(:,-1:0),number_of_values_y, &
         MPI_DOUBLE_PRECISION,sender_top_shift,message_tag,MPI_COMM_WORLD,message_status,code)
   call MPI_SENDRECV(var(:,1:2),number_of_values_y,MPI_DOUBLE_PRECISION, &
         receiver_bottom_shift,message_tag,var(:,NY_LOCAL+1:NY_LOCAL+2),number_of_values_y, &
         MPI_DOUBLE_PRECISION,sender_bottom_shift,message_tag,MPI_COMM_WORLD,message_status,code)

endsubroutine send_receive_topbottom


subroutine send_receive_bottom(var)
  ! 1 communication: 
  ! direction of information: to the bottom
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_y, &
                     receiver_bottom_shift,sender_bottom_shift,&
                     message_tag, message_status, code
  implicit none
  
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
   call MPI_SENDRECV(var(:,1:2),number_of_values_y,MPI_DOUBLE_PRECISION, &
         receiver_bottom_shift,message_tag,var(:,NY_LOCAL+1:NY_LOCAL+2),number_of_values_y, &
         MPI_DOUBLE_PRECISION,sender_bottom_shift,message_tag,MPI_COMM_WORLD,message_status,code)

endsubroutine send_receive_bottom


subroutine send_receive_top(var)
  ! 1 communication: 
  ! direction of information: to the top
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_y, &
                     receiver_top_shift,sender_top_shift,&
                     message_tag, message_status, code
  implicit none
  
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
   call MPI_SENDRECV(var(:,NY_LOCAL-1:NY_LOCAL),number_of_values_y,MPI_DOUBLE_PRECISION, &
         receiver_top_shift,message_tag,var(:,-1:0),number_of_values_y, &
         MPI_DOUBLE_PRECISION,sender_top_shift,message_tag,MPI_COMM_WORLD,message_status,code)


endsubroutine send_receive_top


subroutine send_receive_righttop(var)
  ! 1 communication: 
  ! direction of information: diagonal, to the right/top
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_corner, &
                     receiver_right_top_shift,sender_right_top_shift,&
                     message_tag, message_status, code
  implicit none
  
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
    call MPI_SENDRECV(var(NX_LOCAL-1:NX_LOCAL,NY_LOCAL-1:NY_LOCAL),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_right_top_shift,message_tag,var(-1:0,-1:0),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_right_top_shift,message_tag,MPI_COMM_WORLD,message_status,code)
      
endsubroutine send_receive_righttop


subroutine send_receive_leftbottom(var)
  ! 1 communication: 
  ! direction of information: diagonal, to the left/bottom
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_corner, &
                     receiver_bottom_left_shift,sender_bottom_left_shift,&
                     message_tag, message_status, code
  implicit none
  
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
    call MPI_SENDRECV(var(1:2,1:2),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_bottom_left_shift,message_tag,var(NX_LOCAL+1:NX_LOCAL+2,NY_LOCAL+1:NY_LOCAL+2),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_bottom_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
      
endsubroutine send_receive_leftbottom


subroutine send_receive_rightbottom(var)
  ! 1 communication: 
  ! direction of information: diagonal, to the right/bottom
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_corner, &
                     receiver_right_bottom_shift,sender_right_bottom_shift,&
                     message_tag, message_status, code
  implicit none
  
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
    call MPI_SENDRECV(var(NX_LOCAL-1:NX_LOCAL,1:2),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_right_bottom_shift,message_tag,var(-1:0,NY_LOCAL+1:NY_LOCAL+2),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_right_bottom_shift,message_tag,MPI_COMM_WORLD,message_status,code)
      
endsubroutine send_receive_rightbottom


subroutine send_receive_lefttop(var)
  ! 1 communication: 
  ! direction of information: diagonal, to the left/top
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_corner, &
                     receiver_top_left_shift,sender_top_left_shift,&
                     message_tag, message_status, code
  implicit none
  
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
    call MPI_SENDRECV(var(1:2,NY_LOCAL-1:NY_LOCAL),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_top_left_shift,message_tag,var(NX_LOCAL+1:NX_LOCAL+2,-1:0),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_top_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
      
endsubroutine send_receive_lefttop


subroutine send_receive_corners(var)
  ! 4 communications: 
  !   - direction: diagonal, to the left/bottom
  !   - direction: diagonal, to the right/bottom
  !   - direction: diagonal, to the right/top
  !   - direction: diagonal, to the left/top
  
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_corner, &
                     receiver_right_top_shift,sender_right_top_shift,&
                     receiver_bottom_left_shift,sender_bottom_left_shift,&
                     receiver_top_left_shift,sender_top_left_shift,&
                     receiver_bottom_right_shift,sender_bottom_right_shift,&
                     message_tag, message_status,code
  implicit none
  
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
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
!!                  Definition of Neighbours Processus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_neighboors()
use MPI
use parameters
implicit none

! check that code was compiled with the right number of slices
  if (nb_procs /= NPROC) then
    print *,'error in MPI number of slices: nb_procs,NPROC = ',nb_procs,NPROC,' but they should be equal'
    stop 'nb_procs must be equal to NPROC'
  endif

  ! check that we can cut along Z in an exact number of slices
  if (mod(NX,NPROC_X) /= 0) stop 'NX must be a multiple of nb_procs_x'
  if (mod(NY,NPROC_Y) /= 0) stop 'NY must be a multiple of nb_procs_y'

! check that a slice is at least as thick as a PML layer
  if (NX_LOCAL < NPOINTS_PML .and. USE_PML_XMIN) stop 'NX_LOCAL must be greater than NPOINTS_PML' 
  if (NX_LOCAL < NPOINTS_PML .and. USE_PML_XMAX) stop 'NX_LOCAL must be greater than NPOINTS_PML' 
  if (NY_LOCAL < NPOINTS_PML .and. USE_PML_YMIN) stop 'NY_LOCAL must be greater than NPOINTS_PML' 
  if (NY_LOCAL < NPOINTS_PML .and. USE_PML_YMAX) stop 'NY_LOCAL must be greater than NPOINTS_PML' 
 
  ! rank on line or column 
  ! Example: 
  !  horizontal direction: x,i 
  !  vertical direction: y,j 
  !  ---------------      
  ! |rank 2 |rank 3 |
  ! | (0,1) | (1,1) |
  !  ---------------      
  ! |rank 0 |rank 1 |
  ! | (0,0) | (1,0) |
  !  ---------------      
  i_rank = modulo(rank,NPROC_X)
  j_rank = rank /NPROC_X
  
  ! split process in subgroup
  call MPI_BARRIER(mpi_comm_world,code)
  call MPI_Comm_Split(MPI_COMM_WORLD, j_rank, rank, row_Comm,ierr)
 
  ! offset of this slice when we cut along Z
  offset_i = i_rank * NX_LOCAL
  offset_j = j_rank * NY_LOCAL
   
  ! neighbours for communication to the right, to the left, to the bottom or to the top
  sender_right_shift = rank - 1
  receiver_right_shift = rank + 1
  sender_bottom_shift = rank + NPROC_X
  receiver_bottom_shift = rank - NPROC_X
  sender_left_shift = rank + 1
  receiver_left_shift = rank - 1
  sender_top_shift = rank - NPROC_X
  receiver_top_shift = rank + NPROC_X
  
  ! neighbours for communication in diagonal
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
    
! if we are the first column processus, different neighboor depending on boundary conditions
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
  
! if we are the last column processus, different neighboor depending on boundary conditions
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

! if we are the first row processus, different neighboor depending on boundary conditions
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
  
! if we are the last row processus, different neighboor depending on boundary conditions
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

  ! if processus on the right, top corner, other neighboor depending on boundary conditions
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
  
  ! if processus on the left, top corner, other neighboor depending on boundary conditions
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
  
  ! if processus on the left, bottom corner, other neighboor depending on boundary conditions
  if (i_rank == 0 .and. j_rank == 0) then
      if (.not.(USE_PML_YMIN) .and. .not.(USE_PML_XMIN) .and. .not.(USE_PML_YMAX) .and. .not.(USE_PML_XMAX) ) then
       receiver_bottom_left_shift = NPROC-1
       sender_right_top_shift = NPROC-1 
      else
       receiver_bottom_left_shift = MPI_PROC_NULL
       sender_right_top_shift = MPI_PROC_NULL
    endif
  endif
  
  ! if processus on the right, top corner, other neighboor depending on boundary conditions
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
  
  
endsubroutine get_neighboors


