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


