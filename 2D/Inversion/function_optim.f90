
module parameters 

 integer,parameter :: Nflat = 2
 
 double precision, dimension(Nflat) :: x0
 
 double precision, dimension(Nflat) :: x,dfx
 double precision :: fx

 double precision, dimension(Nflat) :: x_old,dfx_old
 double precision :: fx_old

 double precision, dimension(Nflat) :: x_new,dfx_new
 double precision :: fx_new

 double precision, dimension(Nflat) :: x_low,dfx_low
 double precision :: fx_low

 double precision, dimension(Nflat) :: x_high,dfx_high
 double precision :: fx_high

 double precision, dimension(Nflat) :: r, r_old

 double precision :: dfx_r
 
 double precision :: alpha_start, alpha, alpha_prec, alpha_low, alpha_high
 
 integer :: count_grad, count_restart, count_alpha
 double precision :: reg_weight
  
 double precision, parameter :: c1 = 0.0001
 double precision, parameter :: c2 = 0.9
  
 double precision, parameter :: rate = 0.8d0
 integer, parameter :: maxiter_backtracking = 20

 integer, parameter :: maxiter = 100
 double precision :: tol_x = 1e-10
 double precision :: alpha_max = 10
endmodule parameters


subroutine f(x,fx)

  use parameters, only : Nflat
  implicit none
  double precision, dimension(Nflat) :: x
  double precision :: fx
  
  fx = x(1) - x(2) + 2*x(1)*x(2) + 2*x(1)**2 + x(2)**2 + 1.5d0

endsubroutine f


subroutine df(x,dfx)

  use parameters, only : Nflat
  implicit none
  double precision, dimension(Nflat) :: x
  double precision, dimension(Nflat) :: dfx
  
  dfx(1) = 1 + 2*x(2) + 4*x(1)
  dfx(2) = - 1 + 2*x(1) + 2*x(2)

endsubroutine df


