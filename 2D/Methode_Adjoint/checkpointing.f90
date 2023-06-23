
subroutine load_last_frame(it_time, it_frame, pressure, rhop, vx, vy)
  use parameters, only : NX, NY, NSTEP, NFRAMES, FRAMES, ZERO
  implicit none
  
  double precision, dimension(0:NX+1,0:NY+1) :: pressure, rhop, vx, vy
  integer :: it_frame, it_time, step_between_frames, frame
  logical :: condition 

  step_between_frames = NSTEP / NFRAMES
  condition = .False.
  frame = 1
  
  do while ( (.not. condition) .and. frame < NFRAMES) 
    if ( frame * step_between_frames > it_time ) then
      condition = .True.
      it_frame = 1
      pressure(:,:) = ZERO
      rhop(:,:)     = ZERO
      vx(:,:)       = ZERO
      vy(:,:)       = ZERO
    else if ((frame * step_between_frames <= it_time) .and. &
      ((frame + 1) * step_between_frames > it_time)) then 
      condition = .True.
      it_frame = frame * step_between_frames + 1 
      pressure(:,:) = FRAMES(:,:,1,frame)
      rhop(:,:)     = FRAMES(:,:,2,frame)
      vx(:,:)       = FRAMES(:,:,3,frame)
      vy(:,:)       = FRAMES(:,:,4,frame)
    else 
      frame = frame + 1
    endif
  enddo

endsubroutine load_last_frame



subroutine save_frames()
  
 use parameters
 implicit none
 
 integer :: it_frame
 integer :: it_start, it_end
 integer :: step_between_frames
 
 step_between_frames = NSTEP / NFRAMES - 1
 
 it_end = 0
 
 do it_frame=1,NFRAMES
   it_start = it_end + 1
   it_end   = it_start + step_between_frames
   call forwardproblem(pressure, rhop, vx, vy, p0, rho, v0x, v0y, kappa_unrelaxed, it_start, it_end, 3) 
   FRAMES(:,:,1,it_frame) = pressure(:,:)
   FRAMES(:,:,2,it_frame) = rhop(:,:)
   FRAMES(:,:,3,it_frame) = vx(:,:)
   FRAMES(:,:,4,it_frame) = vy(:,:)

 enddo
 
endsubroutine save_frames

