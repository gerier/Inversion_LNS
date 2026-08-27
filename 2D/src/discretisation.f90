!==============================================================================
!  Finite-difference derivative operators
!
!  This file contains spatial derivative operators used in the finite-
!  difference solver. Both centered and decentered finite-difference schemes
!  are implemented for second- and fourth-order accuracy.
!
!  The operators are designed for staggered-grid discretizations commonly
!  used in wave propagation solvers.
!
!  Available routines:
!    - compute_centered_dU
!    - compute_decentered_forward_dU
!    - compute_decentered_backward_dU
!    - compute_decentered_forward_dU_o2
!    - compute_decentered_backward_dU_o2
!
!==============================================================================


subroutine compute_decentered_forward_dU(u_m, u, u_p, u_pp, derivative, ONE_OVER_6_DSPATIAL)
   !-----------------------------------------------------------------------
! Compute a fourth-order forward decentered finite-difference derivative.
!
! The derivative is evaluated using four consecutive points located
! in the positive spatial direction.
!
! Inputs:
!   u_m, u, u_p, u_pp : values of the field at neighboring points
!   ONE_OVER_6_DSPATIAL : inverse of 6 times the spatial step
!
! Output:
!   derivative : approximated spatial derivative
!
!-----------------------------------------------------------------------
   implicit none
   double precision, intent(in)  :: u_m, u, u_p, u_pp 
   double precision, intent(out) :: derivative 
   double precision, intent(in)  :: ONE_OVER_6_DSPATIAL

   derivative = - (u_pp - 6.0d0*u_p + 3.0d0*u + 2.0d0*u_m) * ONE_OVER_6_DSPATIAL

end subroutine compute_decentered_forward_dU


subroutine compute_decentered_backward_dU(u_mm, u_m, u, u_p, derivative, ONE_OVER_6_DSPATIAL)
!-----------------------------------------------------------------------
! Compute a fourth-order backward decentered finite-difference derivative.
!
! The derivative is evaluated using four consecutive points located
! in the negative spatial direction.
!
! Inputs:
!   u_mm, u_m, u, u_p : values of the field at neighboring points
!   ONE_OVER_6_DSPATIAL : inverse of 6 times the spatial step
!
! Output:
!   derivative : approximated spatial derivative
!
!-----------------------------------------------------------------------
   implicit none
   double precision, intent(in)  :: u_mm, u_m, u, u_p      
   double precision, intent(out) :: derivative           
   double precision, intent(in)  :: ONE_OVER_6_DSPATIAL

   derivative = (u_mm - 6.0d0*u_m + 3.0d0*u + 2.0d0*u_p) * ONE_OVER_6_DSPATIAL

end subroutine compute_decentered_backward_dU


subroutine compute_decentered_forward_dU_o2(u, u_p, u_pp, derivative, ONE_OVER_2_DSPATIAL)
!-----------------------------------------------------------------------
! Compute a second-order forward decentered finite-difference derivative.
!
! This lower-order stencil is used near boundaries where a fourth-order
! approximation cannot be applied.
!
! Inputs:
!   u, u_p, u_pp : values of the field at neighboring points
!   ONE_OVER_2_DSPATIAL : inverse of 2 times the spatial step
!
! Output:
!   derivative : approximated spatial derivative
!
!-----------------------------------------------------------------------
   implicit none
   double precision, intent(in)  :: u, u_p, u_pp 
   double precision, intent(out) :: derivative   
   double precision, intent(in)  :: ONE_OVER_2_DSPATIAL

   derivative = - (u_pp - 4.0d0*u_p + 3.0d0*u) * ONE_OVER_2_DSPATIAL

end subroutine compute_decentered_forward_dU_o2


subroutine compute_decentered_backward_dU_o2(u_mm, u_m, u, derivative, ONE_OVER_2_DSPATIAL)
!-----------------------------------------------------------------------
! Compute a second-order backward decentered finite-difference derivative.
!
! This stencil is used near boundaries where a fourth-order approximation
! is not available.
!
! Inputs:
!   u_mm, u_m, u : values of the field at neighboring points
!   ONE_OVER_2_DSPATIAL : inverse of 2 times the spatial step
!
! Output:
!   derivative : approximated spatial derivative
!
!-----------------------------------------------------------------------
   implicit none
   double precision, intent(in)  :: u_mm, u_m, u 
   double precision, intent(out) :: derivative
   double precision, intent(in)  :: ONE_OVER_2_DSPATIAL

   derivative = (u_mm - 4.0d0*u_m + 3.0d0*u) * ONE_OVER_2_DSPATIAL

end subroutine compute_decentered_backward_dU_o2


subroutine compute_centered_dU(u_mm,u_m,u_p,u_pp, derivative, NINE_OVER_8_DSPATIAL, ONE_OVER_24_DSPATIAL)
!-----------------------------------------------------------------------
! Compute a fourth-order centered finite-difference derivative.
!
! The derivative is computed using a symmetric four-point stencil.
! This routine is designed for the staggered-grid discretization used
! in the finite-difference solver.
!
! Inputs:
!   u_mm, u_m, u_p, u_pp :
!       Field values surrounding the derivative location.
!
!   NINE_OVER_8_DSPATIAL :
!       9/(8*dx) coefficient of the centered fourth-order scheme.
!
!   ONE_OVER_24_DSPATIAL :
!       1/(24*dx) coefficient of the correction term.
!
! Output:
!   derivative : Approximated spatial derivative.
!
! Note:
!   The stencil indexing depends on the staggered-grid location:
!
!        +---x---+
!        |       |
!        |   o   |
!        |       |
!        +---x---+
!
!   + : pressure nodes
!   x : vx nodes
!   o : vy nodes
!
! Examples:
! - For the derivative of vx in pressure point: u_mm=u(i-1), u_m=u(i),   u_p=u(i+1), u_pp(i+2)
! - For the derivative of vy in pressure point: u_mm=u(i-2), u_m=u(i-1), u_p=u(i),   u_pp(i+1)
! - For the derivative of pressure in vx point: u_mm=u(i-2), u_m=u(i-1), u_p=u(i),   u_pp(i+1)
! - For the derivative of pressure in vy point: u_mm=u(i-1), u_m=u(i),   u_p=u(i+1), u_pp(i+2) 
!----------------------------------------------------------------------- 
   implicit none

   double precision, intent(in)  :: u_mm, u_m, u_p, u_pp 
   double precision, intent(out) :: derivative      
   double precision, intent(in)  :: NINE_OVER_8_DSPATIAL, ONE_OVER_24_DSPATIAL
   
  derivative = (u_p - u_m) * NINE_OVER_8_DSPATIAL + (u_mm - u_pp) * ONE_OVER_24_DSPATIAL
end subroutine compute_centered_dU
