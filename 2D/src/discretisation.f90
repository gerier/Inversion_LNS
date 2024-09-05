!!!
!  Functions for Finite Differences Method
!!!


subroutine compute_decentered_forward_dU(u_m, u, u_p, u_pp, derivative, ONE_OVER_6_DSPATIAL)
   ! DF decentered forward, order 4
   implicit none
   double precision :: u_m, u, u_p, u_pp 
   double precision :: derivative             ! output 
   double precision :: ONE_OVER_6_DSPATIAL

   derivative = - (u_pp - 6.d0* u_p + u *3.d0 + u_m*2.d0) * ONE_OVER_6_DSPATIAL

end subroutine compute_decentered_forward_dU


subroutine compute_decentered_backward_dU(u_mm, u_m, u, u_p, derivative, ONE_OVER_6_DSPATIAL)
   ! DF decentered backward, order 4
   implicit none
   double precision :: u_mm, u_m, u, u_p 
   double precision :: derivative             ! output 
   double precision :: ONE_OVER_6_DSPATIAL

   derivative = (u_mm - 6.d0* u_m + u *3.d0 + u_p*2.d0) * ONE_OVER_6_DSPATIAL

end subroutine compute_decentered_backward_dU


subroutine compute_decentered_forward_dU_o2(u, u_p, u_pp, derivative, ONE_OVER_2_DSPATIAL)
   ! DF decentered forward, order 2
   implicit none
   double precision :: u, u_p, u_pp 
   double precision :: derivative             ! output 
   double precision :: ONE_OVER_2_DSPATIAL

   derivative = - (u_pp - 4.d0* u_p + u *3.d0) * ONE_OVER_2_DSPATIAL

end subroutine compute_decentered_forward_dU_o2


subroutine compute_decentered_backward_dU_o2(u_mm, u_m, u, derivative, ONE_OVER_2_DSPATIAL)
   ! DF decentered backward, order 2
   implicit none
   double precision :: u_mm, u_m, u 
   double precision :: derivative             ! output 
   double precision :: ONE_OVER_2_DSPATIAL

   derivative = (u_mm - 4.d0* u_m + u *3.d0) * ONE_OVER_2_DSPATIAL

end subroutine compute_decentered_backward_dU_o2


subroutine compute_centered_dU(u_mm,u_m,u_p,u_pp, derivative, NINE_OVER_8_DSPATIAL, ONE_OVER_24_DSPATIAL)
   ! DF backcentered, order 4
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
