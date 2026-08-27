
!! =====================================================================
!!  Optimisation routines
!! =====================================================================

subroutine update(x,alpha,p,x_new)
!-----------------------------------------------------------------------
! Update the optimization model.
!
! Compute the new model according to:
!   x_new = x + alpha * p
! where p is the current search direction.
!-----------------------------------------------------------------------
 use parameters, only : Nflat,MPI_COMM_WORLD, code
 implicit none
 double precision, dimension(1:Nflat) :: x, p, x_new
 double precision :: alpha

 if (alpha == 0) then
   x_new = x
 else
   x_new = x + alpha * p
   call MPI_BARRIER(MPI_COMM_WORLD, code)

   ! Optional smoothing can be applied here to reduce artifacts
   ! caused by nonlinear effects near the source and adjoint source locations.
   ! call flatmodel2priormodel(x)
   ! TODO
 endif
endsubroutine update


subroutine optimisation() ! TODO
!-----------------------------------------------------------------------
! Perform the nonlinear optimization procedure.
!
! This routine minimizes the objective function using either:
!   - Steepest descent,
!   - Nonlinear conjugate-gradient methods,
!   - Limited-memory BFGS (L-BFGS).
!
! A line-search strategy is used to determine the step length and
! optimization history is stored after each iteration.
!-----------------------------------------------------------------------
 use mpi
 use parameters, only : m0, x, fx, dfx, r, dfx_r, Nflat, &
                        x_old, fx_old, dfx_old, r_old, x_new, fx_new, dfx_new, &
                        alpha, alpha_start, &
                        rank, code, &
                        maxiter_outerloop, type_gradient, steepest_nbiter_default, mem_lbfgs, &
                        Y_list, S_list, RHO_list, &
                        count_f, count_grad, count_restart
 implicit none
 integer :: iter
 logical notFound
 double precision :: dfxdfx_old_r_old, norm2_dfxdfx_old,y_s

 count_grad = 0
 count_f = 0
 count_restart = 0


 ! Initialization
 x(:) = m0(:)
 call f(x,fx)
 call df(x,dfx)
 call save_info_inversion(0)

 x_old(:) = 1e6 * x(:)
 fx_old = 1e10
 dfx_old(:) = 0

 iter = 1
 notFound = .False.

 if (type_gradient == 6) then
    steepest_nbiter_default = mem_lbfgs
 endif

 do while (iter < maxiter_outerloop .and. (abs(fx-fx_old)>1e-14) )

  if (rank == 0) then
    print *, "[Iteration ", iter, "]"
  endif

  if (iter <= steepest_nbiter_default) then
  ! Compute the steepest descent direction
  ! and determine the step length using backtracking

    r(:) = - dfx(:)
    call scalar_product(dfx,r,dfx_r,Nflat)

    ! Store displacement (S) and gradient difference (Y)
    ! vectors required for the L-BFGS update.
    if (type_gradient == 6 .and. iter /= 1) then
      Y_list(iter-1,:) = dfx - dfx_old
      S_list(iter-1,:) = x - x_old
      call scalar_product(Y_list(iter-1,:),S_list(iter-1,:),y_s,Nflat)
      RHO_list(iter-1) = 1 / y_s
    endif

   ! Initialize the trial step length
    if (iter == 1) then
      alpha_start = 1.0d0
    else
      alpha_start = min(1.d0, 1.01d0 * 2.d0 * (fx - fx_old) / dfx_r)
    endif
    if (alpha_start < 0) then
    	alpha_start = 1.0d0
    endif

    ! Find a step length satisfying the Armijo condition using backtracking.
    call backtracking()
    if (rank == 0) then
      print *, "Backtracking : ", alpha
      print *, 'New fx :', fx_new, '(', rank, ')'
    endif

  else ! Compute the nonlinear conjugate-gradient or L-BFGS search direction.

    ! Check whether the previous search direction remains valid.
    ! Restart with steepest descent if the conjugate direction becomes unreliable.
    call scalar_product(dfx-dfx_old,r_old,dfxdfx_old_r_old,Nflat)
    call scalar_product(dfx-dfx_old,dfx-dfx_old,norm2_dfxdfx_old,Nflat)
    if( ((abs(dfxdfx_old_r_old))<1e-8 .and. (type_gradient>=3 .and. type_gradient<=5)) &
       .or. sqrt(norm2_dfxdfx_old)<1e-15) then
       ! TODO The restart criterion should be revisited.
       ! Investigate whether other terms contribute to this instability.
       ! Possible issue with sign conventions.
       r(:) = - dfx(:)
       count_restart = count_restart+1
    else
       call get_descent_direction(type_gradient)
    endif

    ! Initialize the first trial step length.
    ! Nocedal and Wright (2006), eq. 3.60
    call scalar_product(dfx,r,dfx_r,Nflat)
    alpha_start = min(1.0d0, 1.01d0 * 2.d0 * (fx - fx_old) / dfx_r)

    if (rank == 0) then
     print *, "Search direction dot Gradient =", dfx_r
    endif
    if (fx - fx_old >= 0) then
        r = -dfx
        call scalar_product(dfx,r,dfx_r,Nflat)
        alpha_start = 1.d0
        count_restart = count_restart+1
    endif

    if (alpha_start <= 0) then
      print *,  "Alpha start < 0: the previous iteration accepted an increase in the cost function."
      alpha_start = 1.d0
    endif

    !  Find a step length satisfying the strong Wolfe conditions.
    call line_search()
    if (rank == 0) then
    print *, "Previous cost function:", fx
    print *, "New cost function: ", fx_new
    endif
  endif

  ! Update optimization variables and history vectors.
  x_old(:) = x(:)
  fx_old = fx
  dfx_old = dfx(:)
  r_old = r(:)

  x(:) = x_new(:)
  fx = fx_new
  dfx(:) = dfx_new(:)

  call f(x,fx_new)
  call save_info_inversion(iter)
  call MPI_BARRIER(MPI_COMM_WORLD, code)


  iter = iter + 1

  if (rank ==0) then
    print *, "==> Improvement in cost function:", fx -fx_old
  endif

 enddo

endsubroutine optimisation



!! =====================================================================
!!  Line search routines
!! =====================================================================

subroutine get_descent_direction(conj_number)
!-----------------------------------------------------------------------
! Compute a new optimization search direction.
!
! Depending on the selected method, this routine implements:
!   1 - Fletcher-Reeves nonlinear conjugate gradient
!   2 - Polak-Ribière nonlinear conjugate gradient
!   3 - Perry-Shanno nonlinear conjugate gradient
!   4 - Hager-Zhang nonlinear conjugate gradient
!   5 - Dai-Kou nonlinear conjugate gradient
!   6 - Limited-memory BFGS (L-BFGS)
! The algorithms are described in Dener A. et al. (2019)
!-----------------------------------------------------------------------

  use parameters, only : dfx, dfx_old, r, r_old,x,x_old,Nflat, mem_lbfgs, RHO_list, S_list, Y_list
  implicit none
  double precision :: beta,eta,gama,tau,rho
  integer :: conj_number,i
  double precision, dimension(1:Nflat) :: y,s, h0_diag,q
  double precision :: dfx_dfx, dfx_old_dfx_old, dfx_y, y_r_old, dfx_r_old,y_y,y_s, r_old_r_old, dfx_s,s_q,y_r
  double precision, dimension(mem_lbfgs) :: a

  if (conj_number == 1) then ! Fletcher-Reeves
    call scalar_product(dfx,dfx,dfx_dfx,Nflat) !
    call scalar_product(dfx_old,dfx_old,dfx_old_dfx_old,Nflat) !
    beta = dfx_dfx / dfx_old_dfx_old
    r = - dfx + beta * r_old
  else if (conj_number == 2) then ! Polak-Ribieres
    y = dfx - dfx_old
    s = x - x_old
    call scalar_product(dfx,y,dfx_y,Nflat) !
    call scalar_product(dfx_old,dfx_old,dfx_old_dfx_old,Nflat) !
    beta = dfx_y / dfx_old_dfx_old
    r = - dfx + beta * r_old
  else if (conj_number == 3) then ! Perry-Shanno
    y = dfx - dfx_old
    s = x - x_old
    call scalar_product(dfx,y,dfx_y,Nflat) !
    call scalar_product(y,y,y_y,Nflat) !
    call scalar_product(y,r_old,y_r_old,Nflat) !
    call scalar_product(dfx,r_old,dfx_r_old,Nflat)
    call scalar_product(y,s,y_s,Nflat) !

    beta = (dfx_y - 2*y_y / y_r_old * dfx_r_old) / y_r_old
    gama = dfx_r_old / y_r_old
    tau = y_s / y_y
    r = tau * (-dfx + gama * y + beta * r_old)
  else if (conj_number == 4) then ! Hager-Zhang
    y = dfx - dfx_old
    s = x - x_old

    call scalar_product(dfx,y,dfx_y,Nflat) !
    call scalar_product(y,y,y_y,Nflat) !
    call scalar_product(y,r_old,y_r_old,Nflat) !
    call scalar_product(dfx_old,dfx_old,dfx_old_dfx_old,Nflat) !
    call scalar_product(r_old,r_old,r_old_r_old,Nflat) !
    call scalar_product(dfx,r_old,dfx_r_old,Nflat) !

    beta = (dfx_y - 2*y_y / y_r_old * dfx_r_old) / dfx_old_dfx_old
    eta = - 1 / sqrt(r_old_r_old) / min(0.01, sqrt(dfx_old_dfx_old))
    beta = max(beta,eta)
    r = - dfx + beta * r_old
  else if (conj_number == 5) then  ! Dai-Kou
    y = dfx - dfx_old
    s = x - x_old

    call scalar_product(dfx,y,dfx_y,Nflat) !
    call scalar_product(dfx,s,dfx_s,Nflat) !
    call scalar_product(y,y,y_y,Nflat) !
    call scalar_product(y,r_old,y_r_old,Nflat) !
    call scalar_product(y,s,y_s,Nflat) !
    call scalar_product(r_old,r_old,r_old_r_old,Nflat) !
    call scalar_product(dfx,r_old,dfx_r_old,Nflat) !

    beta = dfx_y/y_r_old - y_y/y_s * dfx_s / y_r_old
    beta = max(beta,0.5*dfx_r_old/r_old_r_old)
    r = - dfx + beta * r_old


   else if (conj_number == 6) then  ! L-BFGS
    q = dfx
    y = dfx - dfx_old
    s = x - x_old
    call scalar_product(y,s,y_s,Nflat)
    rho = 1 / y_s

    Y_list = cshift(Y_list, 1, 1)
    S_list = cshift(S_list, 1, 1)
    RHO_list = cshift(RHO_list, 1)

    Y_list(mem_lbfgs,:) = y
    S_list(mem_lbfgs,:) = s
    RHO_list(mem_lbfgs) = rho

    a(:) = 0
    do i=1,mem_lbfgs
      call scalar_product(S_list(mem_lbfgs-i+1,:),q,s_q,Nflat) !
      a(mem_lbfgs-i+1) = RHO_list(mem_lbfgs-i+1) * s_q
      q = q - a(mem_lbfgs-i+1) * Y_list(mem_lbfgs-i+1,:)
    enddo

    call scalar_product(S_list(mem_lbfgs,:),Y_list(mem_lbfgs,:),y_s,Nflat)
    call scalar_product(Y_list(mem_lbfgs,:),Y_list(mem_lbfgs,:),y_y,Nflat)
    h0_diag = y_s / y_y
    r = h0_diag * q
    do i=1,mem_lbfgs
      call scalar_product(Y_list(i,:),r,y_r,Nflat) !
      beta = RHO_list(i) * y_r
      r = r + S_list(i,:) * (a(i)-beta)
    enddo
    r = -r
  endif

endsubroutine get_descent_direction




subroutine line_search()
!-----------------------------------------------------------------------
! Perform a line search satisfying the strong Wolfe conditions.
!
! The trial step length is increased until either the Wolfe conditions
! are satisfied or a bracketing interval is found. In the latter case,
! the zoom procedure is called to refine the solution.
!
! Algorithm 3.5 from Nocedal and Wright (2006), p60.
!-----------------------------------------------------------------------
 use parameters, only : alpha_max, x,r, fx, dfx_r, alpha, alpha_low, alpha_high, alpha_start,&
                        x_new, fx_new, dfx_new, Nflat,HUGEVAL
 implicit none

 double precision :: alpha_prec, fx_prec
 logical stop_cond
 logical :: wolfe1, wolfe2
 double precision :: dfx_new_r

 ! Initialization
 alpha_prec = 0
 fx_prec = HUGEVAL
 alpha = alpha_start

 stop_cond = .False.

 do while (.not.stop_cond)

   call update(x,alpha,r,x_new)
   call f(x_new, fx_new)


   if (.not. wolfe1(alpha,fx,dfx_r,fx_new) .or. fx_new >= fx_prec) then ! not Armijo condition
     ! In this case, alpha_prec is assumed to satisfy the Armijo condition.
     ! Search for a better step length (in the interval [alpha_prec, alpha])
     alpha_high = alpha
     alpha_low = alpha_prec
     call zoom()
     stop_cond = .True.

   else !  Armijo condition satisfied
     print *, "Line search :",alpha
     call df(x_new, dfx_new)
     call scalar_product(dfx_new, r, dfx_new_r, Nflat)
     if (wolfe2()) then ! Check the strong Wolfe curvature condition.
       stop_cond = .True. ! Strong Wolfe conditions are satisfied.
     else
       if (dfx_new_r>0) then
          ! In this case, phi(alpha) = f(x+alpha*r) is increasing
          ! A larger step length cannot further decrease phi.
          ! Since both alpha_prec and alpha satisfy the Armijo condition,
          ! search for a better step length within this interval.
          alpha_high = alpha_prec
          alpha_low = alpha
          call zoom()
          stop_cond = .True.
       else
          ! Search for a larger step length than alpha_prec
          ! that still satisfies the Armijo condition.
          alpha_prec = alpha
          fx_prec = fx_new
          alpha = (alpha_max+alpha)/2

          if (abs(alpha - alpha_max) < 1e-6) then
            ! if alpha is very close to alpha_max, then accept alpha
            stop_cond = .True.
          endif ! alpha == alpha_max

       endif ! < dfx_new , r >   > 0

     endif ! wolfe2

   endif ! wolfe1

 enddo

endsubroutine line_search




subroutine bestAlpha(alpha_hist,size_alpha_hist)
!-----------------------------------------------------------------------
! Select the best tested step length.
!
! Among all candidate step lengths stored during the line search,
! choose the one producing the smallest objective function value f(x+alpha*r)
! and update the trial model accordingly.
!-----------------------------------------------------------------------
use parameters, only : x, r, alpha, x_new,fx_new,maxiter_innerloop,rank,code, MPI_COMM_WORLD
implicit none
double precision :: f_best,best
integer :: i
integer :: size_alpha_hist
double precision, dimension(maxiter_innerloop,2) :: alpha_hist

 call MPI_BARRIER(MPI_COMM_WORLD, code)
 best = alpha_hist(1,1)
 f_best = alpha_hist(1,2)
 do i=2,size_alpha_hist
   if (alpha_hist(i,2) < f_best) then
     f_best = alpha_hist(i,2)
     best = alpha_hist(i,1)
   endif
 enddo

 ! update
 alpha = best
 call update(x, alpha, r, x_new)
 fx_new = f_best

print *, "fx new best ", rank, fx_new
endsubroutine bestAlpha

subroutine backtracking()
!-----------------------------------------------------------------------
! Perform an Armijo backtracking line search.
!
! Starting from the initial trial step length, the step is repeatedly
! reduced until the sufficient decrease (Armijo) condition is satisfied.
! If no acceptable step is found, the best tested value is selected.
!
! Algorithm 3.1 from Nocedal and Wright (2006), p37.
!-----------------------------------------------------------------------
use parameters, only : alpha,fx,fx_new,rate,maxiter_innerloop,x,r,x_new,fx_new,dfx_new,alpha_start,dfx_r
implicit none
 logical :: sufficientdecrease = .False.
 integer :: iter
 double precision, dimension(maxiter_innerloop,2) :: alpha_hist
 integer :: size_alpha_hist
 logical :: wolfe1

    ! Initialization
    alpha = alpha_start

    alpha_hist(:,:) = 0
    size_alpha_hist = 0

    sufficientdecrease = .False.
    iter = 1

    do while ((.not. sufficientdecrease) .and. (iter < maxiter_innerloop))

        ! Update
        call update(x, alpha, r, x_new)
        call f(x_new, fx_new)

        alpha_hist(iter,1) = alpha
        alpha_hist(iter,2) = fx_new
        size_alpha_hist = size_alpha_hist + 1

        if (wolfe1(alpha,fx,dfx_r,fx_new)) then ! Armijo condition
            sufficientdecrease = .True.
        else
            ! Decrease alpha
            alpha = rate * alpha
        endif
        iter = iter + 1

    enddo

    if (iter == maxiter_innerloop) then
        ! Backtracking failed to find a step length satisfying the Armijo condition.
        print *, "ERROR : No acceptable step length found during backtracking."
        print *, "iter = ", iter
        ! Select the step length that produced the smallest objective value.
        call bestAlpha(alpha_hist,size_alpha_hist)
        print *, "Choose alpha = ", alpha
    endif

    call df(x_new,dfx_new)

endsubroutine backtracking


subroutine zoom()
!-----------------------------------------------------------------------
! Refine the step length within a bracketing interval.
!
! The interval is successively reduced using quadratic and cubic
! interpolation until the strong Wolfe conditions are satisfied or the
! interval becomes sufficiently small.
!
! Algorithm 3.6 from Nocedal and Wright (2006), p61.
!-----------------------------------------------------------------------
 use parameters, only : Nflat, alpha, alpha_low, alpha_high,maxiter_innerloop,&
                       x,fx,dfx_r,r,x_new,fx_new,dfx_new,rank,TINYVAL,HUGEVAL
 implicit none

 double precision, dimension(1:Nflat) :: x_low, x_high, dfx_low, dfx_high
 double precision :: fx_low, fx_high

 double precision, dimension(maxiter_innerloop,2) :: alpha_hist
 integer :: size_alpha_hist

 logical :: stop_cond
 integer :: iter
 logical :: wolfe1, wolfe2

 double precision :: dfx_low_r, dfx_high_r, dfx_new_r

 ! Initialisation
 stop_cond = .False.
 iter = 1
 alpha_hist(:,:) = HUGEVAL
 size_alpha_hist=0

 call update(x,alpha_low,r, x_low)
 call f(x_low,fx_low)
 call df(x_low,dfx_low)

 call update(x,alpha_high,r,x_high)
 call f(x_high,fx_high)
 call df(x_high,dfx_high)

 !if (.not.wolfe1(alpha_low,fx,dfx_r,fx_low)) then ! not Armijo condition
 !  print *, "WARNING: alpha_low do not respect sufficient decrease condition"
 !  call backtracking()
 !  stop_cond = .True.
 !endif

 do while (.not. stop_cond)

   ! Interpolate the step length by minimizing a quadratic or cubic model
   ! built from the current interval bounds (alpha_low and alpha_high).
   if (iter == 1) then
     ! On the first iteration, estimate alpha by minimizing the quadratic
     ! interpolant defined by alpha_low and alpha_high.
     call scalar_product(dfx_low,r,dfx_low_r, Nflat)
     call quadratic(alpha_low,fx_low,dfx_low_r,alpha_high, fx_high)
   else
     if (modulo(iter,10) < TINYVAL) then
        if (alpha_low /= 0) then
        alpha = exp( (log(alpha_low) + log(alpha_high))/2)
        else
          alpha = (alpha_low + alpha_high)/2
        endif
     else
      ! On subsequent iterations, estimate alpha by minimizing the cubic
      ! interpolant defined by alpha_low and alpha_high.
      call scalar_product(dfx_low,r,dfx_low_r, Nflat)
      call scalar_product(dfx_high,r,dfx_high_r, Nflat)
      call cubic(alpha_low,fx_low,dfx_low_r,alpha_high, fx_high,dfx_high_r)
     endif
   endif

   ! If the interpolated step length lies outside the interval,
   ! replace it by a safer value inside the interval (mid-point).
   if (alpha <= alpha_low .or. alpha >= alpha_high .or. alpha <= 0 .or. isnan(alpha)) then
     if (alpha_low /= 0 .and. alpha_high /= 0) then
     ! To accelerate convergence, when alpha_low and alpha_high are > 0,
     ! use their geometric mean.
        alpha = exp( (log(alpha_low) + log(alpha_high))/2)
     else
        ! choose the middle point between alpha_low and alpha_high
        alpha = (alpha_low + alpha_high) / 2
      endif
   endif

   if (rank == 0) then
     print*, "Zoom: ", alpha, ' (', alpha_low,",",alpha_high,")", fx_low,dfx_low_r,fx_high, dfx_high_r
   endif

   ! update
   call update(x,alpha,r,x_new)
   call f(x_new, fx_new)
   call df(x_new,dfx_new)
   call scalar_product(dfx_new,r,dfx_new_r, Nflat)

   ! Store tested step lengths and corresponding objective values.
   alpha_hist(iter,1) = alpha
   alpha_hist(iter,2) = fx_new
   size_alpha_hist = size_alpha_hist + 1


   if (wolfe1(alpha,fx,dfx_r,fx_new)) then ! Armijo condition


     if (wolfe2()) then ! Curvature condition
       stop_cond = .True.
     else

       if (dfx_new_r * (alpha_high - alpha_low) >=0) then
         ! in that case, phi(alpha) = f(x+alpha*r) is increasing
         ! we want to look for alpha_st in [alpha_low, alpha] to decrease more f(x+alpha_st*r)
         alpha_high = alpha_low
         x_high = x_low
         fx_high = fx_low
         dfx_high = dfx_low

       endif

       ! If the previous condition was satisfied:
       ! since  alpha respect the decrease condition, we can reduce our interval to "[alpha_low, alpha_high] = [alpha, alpha_low]"
       ! then in line 375, we permut alpha_low and alpha_high to respect the inequality alpha_low < alpha_high
       ! lines are coded in that manner because the code was first inspired by Nocedal and Wright (2006)
       ! if the previous condition was not satisfied:
       ! Since alpha satisfies the sufficient decrease condition,
       ! we can reduce our interval and increase the lower bound
       ! and look for alpha_st in [alpha, alpha_high]
       alpha_low = alpha
       x_low = x_new
       fx_low = fx_new
       dfx_low = dfx_new


       if (abs(alpha_low - alpha_high)<1e-15 .or. maxiter_innerloop <= size_alpha_hist) then
          ! If the interval becomes too small,
          ! select the best tested step length.
          call bestAlpha(alpha_hist,size_alpha_hist)
          call df(x_new,dfx_new)
          stop_cond = .True.
       endif


       !if (alpha_low < 1e-16) then
       ! If alpha_low is too small, stop the algorithm and select an alpha respecting the Armijo condition
       ! call bestAlpha(alpha_hist,size_alpha_hist)
       ! stop_cond = .True.
       !endif

     endif ! Wolfe2

   else

       ! in that case, we can reduce the interval by decreasing the upper bound
       alpha_high = alpha
       x_high = x_new
       fx_high = fx_new
       dfx_high = dfx_new


       if (abs(alpha_low - alpha_high)<1e-15 .or. maxiter_innerloop <= size_alpha_hist) then
         ! if alpha_high too close of alpha_low, stop the algorithm and select the best alpha from the tested ones
         call bestAlpha(alpha_hist,size_alpha_hist)
         call df(x_new,dfx_new)
         stop_cond = .True.
       endif


   endif ! wolfe1

   iter = iter + 1
 enddo

endsubroutine zoom



subroutine quadratic(xa,fxa,dfxa,xb,fxb)
!-----------------------------------------------------------------------
! Estimate the trial step length by quadratic interpolation.
!
! A quadratic model of the objective function is built from two points
! and the derivative at the first point.
!
! Nocedal and Wright (2006), Eq. (3.58), p. 79.
!-----------------------------------------------------------------------
 ! xa   : point a
 ! fxa  : f evaluated in a
 ! dfxa : gradient of f evaluated in a
 ! xb   : point b
 ! fxb  : f evaluated in b
 use parameters, only : alpha
 implicit none
 double precision :: fxa,fxb
 double precision :: xa,xb,dfxa
 double precision :: denom

 denom = (fxb - fxa - (xb - xa) * dfxa ) / ( (xb - xa)**2)
 alpha = xa-dfxa/(2*denom)
endsubroutine quadratic


subroutine cubic(xa,fa,fpa,xb,fb,fpb)
!-----------------------------------------------------------------------
! Estimate the trial step length by cubic interpolation.
!
! A cubic model of the objective function is built from the objective
! values and directional derivatives at two points.
!
! Nocedal and Wright (2006), Eq. (3.59), p. 79.
!-----------------------------------------------------------------------
 ! xa  : point a
 ! fxa : f evaluated in a
 ! fpa : gradient of f evaluated in a
 ! xb  : point b
 ! fxb : f evaluated in b
 ! fpb : gradient of f evaluated in b
 use parameters, only : alpha
 implicit none
 double precision :: fa,fb
 double precision :: denom, num, d1, d2
 double precision :: xa,xb,fpa,fpb

 d1 = fpb + fpa - 3*(fb-fa)/(xb-xa)
 d2 = abs(xa-xb)/(xa-xb) * sqrt(d1**2 - fpb*fpa)
 num = fpa + d2 - d1
 denom = fpa - fpb + 2*d2
 alpha = xa - (xa-xb) * num/denom
endsubroutine cubic


logical function wolfe1(alpha,fx,dfx_r,fx_new)
!-----------------------------------------------------------------------
! Evaluate the Armijo (sufficient decrease) condition.
!
! Returns .true. if the current trial step satisfies the first Wolfe
! condition.
!-----------------------------------------------------------------------
 use parameters, only: c1
 implicit none
 double precision ::  alpha, fx, fx_new, dfx_r
 wolfe1 = fx_new <= fx + alpha * c1 * dfx_r
 return
end


function wolfe2()
!-----------------------------------------------------------------------
! Evaluate the strong Wolfe curvature condition.
!
! Returns .true. if the directional derivative at the trial point
! satisfies the second Wolfe condition.
!-----------------------------------------------------------------------
 use parameters, only :  Nflat,dfx_r,c2,r ,dfx_new
 implicit none
 logical :: wolfe2
 double precision :: dfx_new_r
 call scalar_product(dfx_new, r, dfx_new_r, Nflat)
 wolfe2 = abs(dfx_new_r) <= c2 * abs(dfx_r)
endfunction wolfe2

!! =====================================================================
!!  Routines to smooth a vertical atmospheric profile 
!! =====================================================================

subroutine gaussian_filter(u,mask_dim)
!-----------------------------------------------------------------------
! Apply a one-dimensional Gaussian smoothing filter.
!
! The filter is currently applied along a single model dimension using
! a five-point Gaussian stencil.
!-----------------------------------------------------------------------
  use parameters, only : NX_LOCAL, NY_LOCAL
  implicit none
  integer :: mask_dim
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: u
  double precision, dimension(-1:NY_LOCAL+2) :: u_old
  double precision:: u_aux
  integer :: j
  double precision, dimension(5) :: mask_1d = (/1., 4., 6., 4., 1./) /16

  if (mask_dim == 1) then
    u_old(:) = u(1,:)
    do j=3,NY_LOCAL-2
      u_aux = sum(u_old(j-2:j+2) * mask_1d(:))
      u(:,j) = u_aux
    enddo
    u(:,1) = u(1,3)
    u(:,2) = u(1,3)
    u(:,NY_LOCAL-1) = u(1,NY_LOCAL-2)
    u(:,NY_LOCAL) = u(1,NY_LOCAL-2)
  endif

endsubroutine gaussian_filter


subroutine mean_filter(u,mask_dim)
!-----------------------------------------------------------------------
! Apply a one-dimensional moving-average filter.
!
! A five-point averaging stencil is applied along the selected model
! dimension.
!-----------------------------------------------------------------------
  use parameters, only : NX_LOCAL, NY_LOCAL
  implicit none
  integer :: mask_dim
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: u
  double precision, dimension(-1:NY_LOCAL+2) :: u_old
  double precision:: u_aux
  integer :: j

  if (mask_dim == 1) then
    u_old(:) = u(1,:)
    do j=3,NY_LOCAL-2
      u_aux = sum(u_old(j-2:j+2)) /5
      u(:,j) = u_aux
    enddo
    u(:,1) = u(1,3)
    u(:,2) = u(1,3)
    u(:,NY_LOCAL-1) = u(1,NY_LOCAL-2)
    u(:,NY_LOCAL) = u(1,NY_LOCAL-2)
  endif

endsubroutine mean_filter

subroutine median_filter(u,mask_dim)
!-----------------------------------------------------------------------
! Apply a one-dimensional median filter.
!
! A five-point median filter is applied along the selected model
! dimension in order to reduce isolated numerical artifacts.
!-----------------------------------------------------------------------
  use parameters, only : NX_LOCAL, NY_LOCAL
  implicit none
  integer :: mask_dim
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: u
  double precision, dimension(-1:NY_LOCAL+2) :: u_old
  integer :: j
  double precision, dimension(5) :: imageA

  if (mask_dim == 1) then
    u_old(:) = u(1,:)
    do j=3,NY_LOCAL-2
      imageA = u_old(j-2:j+2)
      call sort(imageA,5)
      u(:,j) = imageA(3)
    enddo
  endif

endsubroutine median_filter


subroutine smoothing(x) ! TODO
!-----------------------------------------------------------------------
! Apply smoothing to the inversion model.
!
! The optimization variables are converted to physical parameters,
! filtered according to the selected smoothing method, and converted
! back to the inversion parameterization.
!-----------------------------------------------------------------------
  use parameters, only : Nflat,c0_prior,rho0_prior,p0_prior,windx_prior,gamma_chemestry, &
                             type_smoothing
  implicit none
  double precision, dimension(1:Nflat), intent(inout) :: x
  call flatmodel2priormodel(x)

  if (type_smoothing == 1) then
  ! Apply a moving-average filter.
    call mean_filter(rho0_prior,1)
    call mean_filter(p0_prior,1)
    call mean_filter(windx_prior,1)
  elseif (type_smoothing == 2) then
  ! Apply a Gaussian filter.
    call gaussian_filter(rho0_prior,1)
    call gaussian_filter(p0_prior,1)
    call gaussian_filter(windx_prior,1)
  elseif (type_smoothing == 3) then
  ! Apply a median filter.
    call median_filter(rho0_prior,1)
    call median_filter(p0_prior,1)
    call median_filter(windx_prior,1)

  else
    if (type_smoothing /= 0) then
     print *, "ERROR: Type smoothing unknown"
    endif
  endif

  c0_prior(:,:)   = sqrt(gamma_chemestry(:,:) * p0_prior(:,:) / rho0_prior(:,:))

  call priormodel2flatmodel(x)

endsubroutine smoothing

!! =====================================================================
!!  Other routines
!! =====================================================================

subroutine scalar_product(u,v,res,dims)
!-----------------------------------------------------------------------
! Compute the global scalar product of two distributed vectors.
!
! Each MPI process computes the local contribution and an MPI reduction
! is performed to obtain the global dot product.
!-----------------------------------------------------------------------
 use MPI
 use parameters, only : code
 implicit none 
 integer :: dims
 double precision, dimension(dims) :: u,v
 double precision :: res, res_local

 res_local = dot_product(u,v)

 call MPI_BARRIER(MPI_COMM_WORLD, code)
 call MPI_ALLREDUCE(res_local, res, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)

endsubroutine scalar_product


subroutine get_alpha_low_high(alpha_1, alpha_2, alpha_low, alpha_high)
!-----------------------------------------------------------------------
! Order two trial step lengths.
!
! Returns alpha_low = min(alpha_1, alpha_2) and
! alpha_high = max(alpha_1, alpha_2).
!-----------------------------------------------------------------------
 implicit none
 double precision :: alpha_1, alpha_2, alpha_low, alpha_high
 if (alpha_1 < alpha_2) then
   alpha_low = alpha_1
   alpha_high = alpha_2
 else
   alpha_low = alpha_2
   alpha_high = alpha_1
 endif
end subroutine


subroutine sort(array, array_size)
!-----------------------------------------------------------------------
! Sort an array in ascending order.
!
! This routine is used by the median filter to determine the median
! value of a local stencil.
!-----------------------------------------------------------------------
implicit none
 integer :: array_size, i,j
 double precision, dimension(array_size), intent(inout) :: array
 double precision :: box

 do i=1,int(array_size/2)+1
  do j=i+1,array_size
   if (array(j) <=array(i)) then
    box = array(i)
    array(i)=array(j)
    array(j)=box
   endif
  enddo
 enddo

endsubroutine sort
