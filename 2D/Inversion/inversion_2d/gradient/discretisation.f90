

subroutine compute_decentered_forward_dU(u_m, u, u_p, u_pp, derivative, ONE_OVER_6_DSPATIAL)
   implicit none
   double precision :: u_m, u, u_p, u_pp 
   double precision :: derivative             ! output 
   double precision :: ONE_OVER_6_DSPATIAL

   derivative = - (u_pp - 6.d0* u_p + u *3.d0 + u_m*2.d0) * ONE_OVER_6_DSPATIAL

end subroutine compute_decentered_forward_dU

subroutine compute_decentered_backward_dU(u_mm, u_m, u, u_p, derivative, ONE_OVER_6_DSPATIAL)
   implicit none
   double precision :: u_mm, u_m, u, u_p 
   double precision :: derivative             ! output 
   double precision :: ONE_OVER_6_DSPATIAL

   derivative = (u_mm - 6.d0* u_m + u *3.d0 + u_p*2.d0) * ONE_OVER_6_DSPATIAL

end subroutine compute_decentered_backward_dU



subroutine compute_centered_dU(u_mm,u_m,u_p,u_pp, derivative, NINE_OVER_8_DSPATIAL, ONE_OVER_24_DSPATIAL)
   implicit none

   double precision :: u_mm, u_m, u_p, u_pp 
   double precision :: derivative             ! output 
   double precision :: NINE_OVER_8_DSPATIAL, ONE_OVER_24_DSPATIAL

  ! be aware, the mesh is as follow:
  ! + --x-- +
  ! |       |
  ! |       o
  ! |       |
  ! + --x-- +
  ! with + node where pressure is discretised, x where vx is discretised, o where vy is discretised
  ! Be aware when you do a centered derivative. 
  ! If you want derivative of vx in pressure point: u_mm=u(i-1), u_m=u(i), u_p=u(i+1), u_pp(i+2)
  ! If you want derivative of vy in pressure point: u_mm=u(i-2), u_m=u(i-1), u_p=u(i), u_pp(i+1)
  ! If you want derivative of pressure in vx point: u_mm=u(i-2), u_m=u(i-1), u_p=u(i), u_pp(i+1)
  ! If you want derivative of pressure in vy point: u_mm=u(i-1), u_m=u(i), u_p=u(i+1), u_pp(i+2) 
  
  derivative = 0d0  
  derivative = (u_p - u_m) * NINE_OVER_8_DSPATIAL + (u_mm - u_pp) * ONE_OVER_24_DSPATIAL


end subroutine compute_centered_dU



subroutine get_index_boundarycondition(i,j, Im2, Im1,Ip1, Ip2, Jm2, Jm1,Jp1, Jp2)

  use parameters, only : USE_PML_XMIN, USE_PML_YMIN, NX, NY, &
                         rank,i_rank,j_rank,NPROC_X,NPROC_Y,NX_LOCAL,NY_LOCAL
  implicit none
  
  integer :: &
    i,j, &! input
    Im2, Im1,Ip1, Ip2, Jm2, Jm1,Jp1, Jp2
  
  
    if (USE_PML_YMIN) then
      if (j_rank == 0) then
         if (j == 1) then
	        Jm1 = 0
	        Jm2 = -1
         elseif (j == 2) then
	        Jm2 = 0
	 endif
      endif
      if (j_rank == NPROC_Y-1) then
         if (j == NY_LOCAL) then
	        Jp1 = NY_LOCAL+1
	        Jp2 = NY_LOCAL+2
         else if (j == NY_LOCAL-1) then
               Jp2 = NY_LOCAL+1
         endif
      endif
    !else 
     
      !if (j == 1) then
	!      Jm1 = NY-1
	!      Jm2 = NY-2
       !elseif (j == 2) then
	!      Jm2 = NY-1
       !else if (j == NY) then
	!      Jp1 = 2
	!      Jp2 = 3
       !else if (j == NY-1) then
       !      Jp2 = 2
       !endif
    
     endif
     
     
     if (USE_PML_XMIN) then
       if (i_rank == 0) then
         if (i == 1) then
	   Im1 = 0
	   Im2 = -1
         elseif (i == 2) then
	   Im2 = 0
         endif
       endif
       if (i_rank == NPROC_X-1) then
         if (i == NX_LOCAL) then
	   Ip1 = NX_LOCAL+1
	   Ip2 = NX_LOCAL+2
         else if (i == NX_LOCAL-1) then
	   Ip2 = NX_LOCAL+1
         endif
       endif  
          
 !    else
   
       !if (rank==0) then
         !if (i == 1) then
	 !  Im1 = NX_LOCAL-1
	 !  Im2 = NX_LOCAL-2
         !elseif (i == 2) then
	  ! Im2 = NX_LOCAL-1
	  !endif
       !else if (rank==nb_procs-1) then
       !  if (i == NX_LOCAL) then
!	   Ip1 = 2
!	   Ip2 = 3
!        else if (i == NX_LOCAL-1) then
	!   Ip2 = 2
       ! endif  
       !endif
       
    endif
  
endsubroutine get_index_boundarycondition
