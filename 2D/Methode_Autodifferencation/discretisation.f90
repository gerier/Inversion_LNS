!subroutine compute_centered_derivative(u_mm, u_m, u_p, u_pp, derivative, ONE_OVER_24_DSPATIAL, NINE_OVER_8_DSPATIAL)
!   double precision :: u_mm, u_m, u_p, u_pp    ! input u(i-1), u(i-1/2), u(i+1/2), u(i+1)
!   double precision :: derivative             ! output 
!   double precision :: ONE_OVER_24_DSPATIAL,NINE_OVER_8_DSPATIAL
!   derivative = (u_p - u_m) * NINE_OVER_8_DSPATIAL + (u_mm - u_pp) * ONE_OVER_24_DSPATIAL
!end subroutine compute_centered_derivative


!subroutine compute_decentered_derivative(u_mm, u_m, u, u_p, derivative, ONE_OVER_6_DSPATIAL)
!   double precision :: u_mm, u_m, u_p, u_pp    ! input u(i-2), u(i-1), u(i), u(i+1)
!   double precision :: derivative             ! output 
!   double precision :: ONE_OVER_6_DSPATIAL
!   derivative = (u_mm  - 6.d0 * u_m  + u * 3.d0  + u_p  *2.d0)  * ONE_OVER_SIX_DELTAX
!end subroutine compute_centered_derivative

!subroutine compute_decentered_dU_dx_in_i(U, i, j, derivative, ONE_OVER_6_DSPATIAL, NX, NY)
!   implicit none
!   integer :: NX, NY
!   integer :: i,j
!   double precision, dimension (0:NX+1,0:NY+1) :: U               ! input 
!   double precision :: derivative             ! output 
!   double precision :: ONE_OVER_6_DSPATIAL
!
!   if (i == 2) then
!     derivative = (U(NX-1,j)  - 6.d0 * U(1,j)  + U(2,j) * 3.d0  + U(3,j)  *2.d0)  * ONE_OVER_6_DSPATIAL
!   elseif (i==1) then
!     derivative = (U(NX-2,j)  - 6.d0 * U(NX-1,j)  + U(1,j) * 3.d0  + U(2,j)  *2.d0)  * ONE_OVER_6_DSPATIAL
!   elseif (i==NX) then
!     derivative = (U(NX-2,j)  - 6.d0 * U(NX-1,j)  + U(NX,j) * 3.d0  + U(2,j)  *2.d0)  * ONE_OVER_6_DSPATIAL
!   else
!     derivative = (U(i-2,j)  - 6.d0 * U(i-1,j)  + U(i,j) * 3.d0  + U(i+1,j)  *2.d0)  * ONE_OVER_6_DSPATIAL
!   endif
!   
!end subroutine compute_decentered_dU_dx_in_i
!
!
subroutine compute_decentered_dU_dz_in_i(U, i, j, derivative, ONE_OVER_6_DSPATIAL,NX,NY,wind)
   implicit none
   integer :: NX, NY
   integer :: i,j
   double precision, dimension (0:NX+1,0:NY+1) :: U               ! input 
   double precision :: derivative             ! output 
   double precision :: ONE_OVER_6_DSPATIAL
   double precision :: wind
   derivative = 0

   if (j == 2) then
     derivative = (U(i,NY-1)  - 6.d0 * U(i,1)     + 3.d0 * U(i,j)  + 2.d0 * U(i,3))    * ONE_OVER_6_DSPATIAL
   elseif (j==1) then
     derivative = (U(i,NY-2)  - 6.d0 * U(i,NY-1)  + 3.d0 * U(i,j)  + 2.d0 * U(i,2))    * ONE_OVER_6_DSPATIAL
   elseif (j==NY) then
     derivative = (U(i,NY-2)  - 6.d0 * U(i,NY-1)  + 3.d0 *  U(i,j) + 2.d0 * U(i,2))    * ONE_OVER_6_DSPATIAL
   else
     derivative = (U(i,j-2)   - 6.d0 * U(i,j-1)   + 3.d0 * U(i,j)  + 2.d0 * U(i,j+1))  * ONE_OVER_6_DSPATIAL
   endif
   
end subroutine compute_decentered_dU_dz_in_i

!subroutine compute_decentered_dU_dz_in_i(U, i, j, derivative, ONE_OVER_6_DSPATIAL,NX,NY,wind)
!   implicit none
!   integer :: NX, NY
!   integer :: i,j
!   double precision, dimension (0:NX+1,0:NY+1) :: U               ! input 
!   double precision :: derivative             ! output 
!   double precision :: ONE_OVER_6_DSPATIAL
!   double precision :: wind, coef
!   integer :: coef_int
!
!   coef = sign(1d0,wind)
!   coef_int = INT(coef)
!
!   if (coef >= 0 .and. j == 2) then
!     derivative = (U(i,NY-1)  - 6.d0 * U(i,1)  + U(i,2) * 3.d0  + U(i,3)  *2.d0) * ONE_OVER_6_DSPATIAL
!   elseif (coef >= 0 .and. j==1) then
!     derivative = (U(i,NY-2)  - 6.d0 * U(i,NY-1)  + U(i,1) * 3.d0  + U(i,2) *2.d0)  * ONE_OVER_6_DSPATIAL
!   elseif (coef >= 0 .and. j==NY) then
!     derivative = (U(i,NY-2)  - 6.d0 * U(i,NY-1)  + U(i,NY) * 3.d0  + U(i,2) *2.d0)  * ONE_OVER_6_DSPATIAL
!   elseif (coef < 0 .and. j == 1) then
!     derivative = - (U(i,3)  - 6.d0 * U(i,2)  + U(i,1) * 3.d0  + U(i,NY)  *2.d0) * ONE_OVER_6_DSPATIAL
!   elseif (coef < 0 .and. j==NY) then
!     derivative = - (U(i,3)  - 6.d0 * U(i,2)  + U(i,NY) * 3.d0  + U(i,NY-1) *2.d0)  * ONE_OVER_6_DSPATIAL
!   elseif (coef < 0 .and. j==NY-1) then
!     derivative = - (U(i,2)  - 6.d0 * U(i,NY)  + U(i,NY-1) * 3.d0  + U(i,NY-2) *2.d0)  * ONE_OVER_6_DSPATIAL
!
!   else
!     derivative = coef* (U(i,j-2*coef_int) - 6.d0* U(i,j-coef_int) + U(i,j)* 3.d0 + U(i,j+coef_int)*2.d0) * ONE_OVER_6_DSPATIAL
!   endif
!
!end subroutine compute_decentered_dU_dz_in_i


subroutine compute_decentered_dU_dx_in_i(U, i, j, derivative, ONE_OVER_6_DSPATIAL, NX, NY, wind)
   implicit none
   integer :: NX, NY
   integer :: i,j
   double precision, dimension (0:NX+1,0:NY+1) :: U               ! input 
   double precision :: derivative             ! output 
   double precision :: ONE_OVER_6_DSPATIAL
   double precision :: wind, coef
   integer :: coef_int

   coef = sign(1d0,wind)
   coef_int = INT(coef)

   if (coef >= 0 .and. i == 2) then
     derivative = (U(NX-1,j)  - 6.d0 * U(1,j)  + U(2,j) * 3.d0  + U(3,j)  *2.d0) * ONE_OVER_6_DSPATIAL
   elseif (coef >= 0 .and. i==1) then
     derivative = (U(NX-2,j)  - 6.d0 * U(NX-1,j)  + U(1,j) * 3.d0  + U(2,j) *2.d0)  * ONE_OVER_6_DSPATIAL
   elseif (coef >=0 .and. i==NX) then
     derivative = (U(NX-2,j)  - 6.d0 * U(NX-1,j)  + U(NX,j) * 3.d0  + U(2,j) *2.d0)  * ONE_OVER_6_DSPATIAL
   elseif (coef < 0 .and. i == 1) then
     derivative = (U(3,j)  - 6.d0 * U(2,j)  + U(1,j) * 3.d0  + U(NX-1,j)  *2.d0) * ONE_OVER_6_DSPATIAL
   elseif (coef < 0 .and. i==NX) then
     derivative = (U(3,j)  - 6.d0 * U(2,j)  + U(NX,j) * 3.d0  + U(NX-1,j) *2.d0) * ONE_OVER_6_DSPATIAL
   elseif (coef < 0 .and. i==NX-1) then
     derivative = (U(2,j)  - 6.d0 * U(1,j)  + U(NX-1,j) * 3.d0  + U(NX-2,j) *2.d0)  * ONE_OVER_6_DSPATIAL
   else
     derivative = coef* (U(i-2*coef_int,j) - 6.d0* U(i-coef_int,j) + U(i,j) *3.d0 + U(i+coef_int,j)*2.d0) * ONE_OVER_6_DSPATIAL
   endif

end subroutine compute_decentered_dU_dx_in_i


subroutine compute_centered_dU_dx_in_i(U, i, j, derivative, NINE_OVER_8_DSPATIAL, ONE_OVER_24_DSPATIAL,NX,NY)
   implicit none
   integer :: NX,NY
   integer :: i,j
   double precision, dimension (0:NX+1,0:NY+1) :: U   ! input 
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
  
  if (i == 1) then
    u_mm = U(NX-1,j)
    u_m  = U(1,j)
    u_p  = U(2,j)
    u_pp = U(3,j)  
  elseif (i ==2 ) then
    u_mm = U(1,j)
    u_m  = U(2,j)
    u_p  = U(3,j)
    u_pp = U(4,j)  
  elseif (i== NX) then
    u_mm = U(NX-1,j)
    u_m  = U(NX,j)
    u_p  = U(2,j)
    u_pp = U(3,j)  
  elseif (i == NX-1) then 
    u_mm = U(NX-2,j)
    u_m  = U(NX-1,j)
    u_p  = U(NX,j)
    u_pp = U(2,j)      
  else 
    u_mm = U(i-1,j)
    u_m  = U(i,j)
    u_p  = U(i+1,j)
    u_pp = U(i+2,j)  
  endif
  
  derivative = (u_p - u_m) * NINE_OVER_8_DSPATIAL + (u_mm - u_pp) * ONE_OVER_24_DSPATIAL

end subroutine compute_centered_dU_dx_in_i


subroutine compute_centered_dU_dz_in_i(U, i, j, derivative, NINE_OVER_8_DSPATIAL, ONE_OVER_24_DSPATIAL,NX,NY)
   implicit none
   integer :: NX,NY
   integer :: i,j
   double precision, dimension (0:NX+1,0:NY+1) :: U   ! input 
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
  
  if (j == 1) then
    u_mm = U(i,NY-1)
    u_m  = U(i,1)
    u_p  = U(i,2)
    u_pp = U(i,3)  
  elseif (j == 2) then
    u_mm = U(i,1)
    u_m  = U(i,2)
    u_p  = U(i,3)
    u_pp = U(i,4) 
  elseif (j == NY) then
    u_mm = U(i,NY-1)
    u_m  = U(i,NY)
    u_p  = U(i,2)
    u_pp = U(i,3)  
  elseif (j == NY-1) then 
    u_mm = U(i,NY-2)
    u_m  = U(i,NY-1)
    u_p  = U(i,NY)
    u_pp = U(i,2)        
  else 
    u_mm = U(i,j-1)
    u_m  = U(i,j)
    u_p  = U(i,j+1)
    u_pp = U(i,j+2)  
  endif
  
  derivative = (u_p - u_m) * NINE_OVER_8_DSPATIAL + (u_mm - u_pp) * ONE_OVER_24_DSPATIAL

end subroutine compute_centered_dU_dz_in_i
