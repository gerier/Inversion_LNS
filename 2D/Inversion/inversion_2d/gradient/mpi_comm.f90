subroutine send_receive_rightleft(var)
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
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_corner, &
                     receiver_right_top_shift,sender_right_top_shift,&
                     receiver_bottom_left_shift,sender_bottom_left_shift,&
                     message_tag, message_status, code
  implicit none
  
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
    call MPI_SENDRECV(var(1:2,1:2),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_bottom_left_shift,message_tag,var(NX_LOCAL+1:NX_LOCAL+2,NY_LOCAL+1:NY_LOCAL+2),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_bottom_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
      
endsubroutine send_receive_leftbottom

subroutine send_receive_lefttop(var)
  use mpi
  use parameters, only : NX_LOCAL,NY_LOCAL,number_of_values_corner, &
                     receiver_right_top_shift,sender_right_top_shift,&
                     receiver_top_left_shift,sender_top_left_shift,&
                     message_tag, message_status, code
  implicit none
  
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: var
  
    call MPI_SENDRECV(var(1:2,NY_LOCAL-1:NY_LOCAL),number_of_values_corner,MPI_DOUBLE_PRECISION, &
         receiver_top_left_shift,message_tag,var(NX_LOCAL+1:NX_LOCAL+2,-1:0),number_of_values_corner, &
         MPI_DOUBLE_PRECISION,sender_top_left_shift,message_tag,MPI_COMM_WORLD,message_status,code)
      
endsubroutine send_receive_lefttop

subroutine send_receive_corners(var)
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



subroutine get_neighboors()

use parameters

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
  
  
  !do rk=0,NPROC-1
  !if (rk == rank) then
  !  print *, rank, sender_right_shift, receiver_right_shift
  !  print *, rank, sender_left_shift, receiver_left_shift
  !  print *, rank, sender_bottom_shift, receiver_bottom_shift    
  !  print *, rank, sender_top_shift, receiver_top_shift
  !  print *, rank, sender_right_bottom_shift, receiver_right_bottom_shift
  !  print *, rank, sender_bottom_right_shift, receiver_bottom_right_shift
  !  print *, rank, sender_right_top_shift, receiver_right_top_shift
  !  print *, rank, sender_bottom_left_shift, receiver_bottom_left_shift
  !  print *, rank, sender_top_left_shift, receiver_top_left_shift
  !  endif
  !  enddo
  

  
  i2begin = 1
  if (j_rank == 0) i2begin = 2 ! TODO modifier pour voir si necessaire la taille supp ?

  iminus1end = NX_LOCAL
  if (j_rank == NPROC_Y - 1) iminus1end = NX_LOCAL-1 ! TODO idem l.78
  
  j2begin = 1
  if (i_rank == 0) j2begin = 2 ! TODO modifier pour voir si necessaire la taille supp ?

  jminus1end = NY_LOCAL
  if (i_rank == NPROC_X - 1) jminus1end = NY_LOCAL-1 ! TODO idem l.78
  
endsubroutine get_neighboors


