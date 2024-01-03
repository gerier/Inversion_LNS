subroutine scalar_product(u,v,res,dims)
 use MPI
 use parameters, only : code
 integer :: dims
 double precision, dimension(dims) :: u,v
 double precision :: res, res_local
 
 res_local = dot_product(u,v)

 call MPI_BARRIER(MPI_COMM_WORLD, code)
 call MPI_ALLREDUCE(res_local, res, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)

endsubroutine scalar_product


subroutine optimisation() ! TODO

 use parameters
 implicit none
 integer :: iter
 logical notFound 
 double precision :: dfxdfx_old_r_old
 
 count_grad = 0
 count_alpha = 0
 count_restart = 0
 
 ! initialisation 
 x(:) = x0(:) 
 call f(x,fx)
 call df(x,dfx) 
 
 x_old(:) = 1e6 * x(:)
 fx_old = 1e10
 dfx_old(:) = 0
 
 iter = 1
 notFound = .False.
 
 do while (iter < maxiter .and. (abs(fx-fx_old)>1e-14) )
  print *, "[Iteration ", iter, "]"
 
  if (iter <= 2) then  ! compute the steepest descent gradient

    r(:) = - dfx(:)
    call scalar_product(dfx,r,dfx_r,Nflat) ! dfx_r = dot_product(dfx, r)

   ! initialise a first step
    if (iter == 1.0d0) then
      alpha_start = 1.0d0
    else 
      alpha_start = min(1.d0, 1.01d0 * 2.d0 * (fx - fx_old) / dfx_r)
    endif
    
    ! find a step that respects only the Armijo condition
    call backtracking()

    
  else ! compute the conjugate gradient
  
    ! use the conjugate gradient
    call scalar_product(dfx-dfx_old,r_old,dfxdfx_old_r_old,Nflat) 
    if ((abs(dfxdfx_old_r_old)) < 1e-8 ) then ! TODO pas logique <dfx,r> ! .or. dot_product(dfx,r) > 0 )
       r(:) = - dfx(:)
       count_restart = count_restart+1
    else 
       call conjugateGradient(5)
    endif
    call scalar_product(dfx,r,dfx_r,Nflat) ! dfx_r = dot_product(dfx, r)
    
    ! initialise a first step
    ! Nocedal, 2006, eq. 3.60
    alpha_start = min(1.0d0, 1.01d0 * 2.d0 * (fx - fx_old) / dfx_r)
    if (alpha_start < 0) then
      print *, "Alpha start < 0 :  we have accepted to increase the mistfit in the previous iteration"
      alpha_start = 1.d0
    endif
    
    ! find a step that respects the Strong Wolfe conditions
    call line_search() 
    print *, fx_new
  endif
 
  ! update iterates
  iter = iter + 1 
  
  x_old(:) = x(:)
  fx_old = fx
  dfx_old = dfx(:)
  r_old = r(:)
  
  x(:) = x_new(:)
  fx = fx_new 
  dfx(:) = dfx_new(:)
 
 enddo
 
endsubroutine optimisation



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine conjugateGradient(conj_number)
  ! Dener A. et al., 2019
  use parameters, only : dfx, dfx_old, r, r_old,x,x_old,Nflat
  implicit none
  double precision :: beta,eta,gama,tau
  integer :: conj_number
  double precision, dimension(Nflat) :: y,s
  double precision :: dfx_dfx, dfx_old_dfx_old, dfx_y, y_r_old, dfx_r_old,y_y,y_s, r_old_r_old, dfx_s
  
  if (conj_number == 1) then ! Fletcher-Reeves
    !beta =  dot_product(dfx,dfx) / dot_product(dfx_old,dfx_old)
    call scalar_product(dfx,dfx,dfx_dfx,Nflat) !
    call scalar_product(dfx_old,dfx_old,dfx_old_dfx_old,Nflat) !
    beta = dfx_dfx / dfx_old_dfx_old
    r = - dfx + beta * r_old
  else if (conj_number == 2) then ! Polak-Ribieres
    y = dfx - dfx_old
    s = x - x_old
    !beta = dot_product(dfx,y) / dot_product(dfx_old,dfx_old)
    call scalar_product(dfx,y,dfx_y,Nflat) !
    call scalar_product(dfx_old,dfx_old,dfx_old_dfx_old,Nflat) !
    beta = dfx_dfx / dfx_old_dfx_old
    r = - dfx + beta * r_old
  else if (conj_number == 3) then ! Perry-Shanno
    y = dfx - dfx_old
    s = x - x_old
    call scalar_product(dfx,y,dfx_y,Nflat) !
    call scalar_product(y,y,y_y,Nflat) !
    call scalar_product(y,r_old,y_r_old,Nflat) !
    call scalar_product(dfx,r_old,dfx_r_old,Nflat) 
    call scalar_product(y,s,y_s,Nflat) !
    
    ! beta = (dot_product(dfx,y) - 2*dot_product(y,y) / dot_product(y,r_old) * dot_product(dfx,r_old)) / dot_product(y,r_old)
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
    
    beta = dfx_y/y_r_old - y_y/y_s * dfx_s / y_r_old
    beta = max(beta,0.5*dfx_r_old/r_old_r_old)
    r = - dfx + beta * r_old
  endif
  
  
 
endsubroutine conjugateGradient




subroutine line_search() 
 ! algorithm 3.5 p60, Nocedal, 2006
 use parameters, only : alpha_max, x,r, fx, dfx_r,c1,c2, alpha, alpha_low, alpha_high, alpha_start,&
                        x_new, fx_new, dfx_new, Nflat
 implicit none
 
 double precision :: alpha_prec
 logical stop_cond
 logical :: wolfe1, wolfe2
 double precision :: dfx_new_r
 
 ! initialisation
 alpha_prec = 1.0d-10
 alpha = alpha_start
 
 stop_cond = .False.
 
 do while (.not.stop_cond)
 
   call update(x,alpha,r,x_new)
   call f(x_new, fx_new)
 
    
   if (.not. wolfe1(alpha,fx,dfx_r,fx_new)) then ! not Armijo condition 
     ! in that case, we suppose that alpha_prec respect the Armijo condition
     ! we look for a better alpha in the intervalle [alpha_prec, alpha] 
     call get_alpha_low_high(alpha_prec, alpha, alpha_low, alpha_high)
     call zoom()  ! TODO
     stop_cond = .True.
     
   else !  Armijo condition satisfied
     call df(x_new, dfx_new) 
     call scalar_product(dfx_new, r, dfx_new_r, Nflat)
     if (wolfe2()) then ! check curvature condition
       stop_cond = .True. ! Strong Wolfe condition are satisifed 
     else
       if (dfx_new_r>0) then
          ! in that case, phi(alpha) = f(x+apha*r) is increasing
          ! it is not possible to find a larger alpha to minimise phi
          ! alpha_prec and alpha satisfy the Armijo condition, we look for the better in that intervalle
          call get_alpha_low_high(alpha_prec, alpha, alpha_low, alpha_high)
          call zoom()  ! TODO
          stop_cond = .True.
       else 
          ! in that case, look for a larger alpha (than alpha_prec) respecting at least the Armijo condition 
          alpha_prec = alpha
          alpha = (alpha_max+alpha)/2
          if (alpha == alpha_max) then
            stop_cond = .True.
          endif ! alpha == alpha_max
          
       endif ! < dfx_new , r >   > 0
       
     endif ! wolfe2
   
   endif ! wolfe1
 
 enddo 
 
 
endsubroutine line_search


subroutine get_alpha_low_high(alpha_1, alpha_2, alpha_low, alpha_high)
 double precision :: alpha_1, alpha_2, alpha_low, alpha_high
 ! subroutine that define the lowest/largest between alpha_1 and alpha_2
 if (alpha_1 < alpha_2) then
   alpha_low = alpha_1
   alpha_high = alpha_2
 else 
   alpha_low = alpha_2
   alpha_high = alpha_1
 endif 
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bestAlpha(alpha_hist,size_alpha_hist) ! Comment faire le array ? ne connait pas la taille 
use parameters, only : x, r, alpha, x_new, fx_new,maxiter_backtracking
implicit none
double precision :: f_best,best
integer :: i 
integer :: size_alpha_hist
double precision, dimension(maxiter_backtracking,2) :: alpha_hist

 ! subroutine to find the alpha that minimise f(x+alpha*r)
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
 
endsubroutine bestAlpha

subroutine backtracking()
! algorithm 3.1 p37 from Nocedal, 2006
use parameters, only : alpha, Nflat,fx,fx_new,rate,maxiter_backtracking,x,r,x_new,fx_new,dfx_new,alpha_start,dfx_r
implicit none
 logical :: sufficientdecrease = .False.
 integer :: iter
 double precision, dimension(maxiter_backtracking,2) :: alpha_hist
 integer :: size_alpha_hist
 logical :: wolfe1
 
    ! initialisation
    alpha = alpha_start 
    
    alpha_hist(:,:) = 0
    size_alpha_hist = 0
    
    sufficientdecrease = .False.
    iter = 1
 
    do while ((.not. sufficientdecrease) .and. (iter < maxiter_backtracking)) 
        print *, "dans backtracking, iter : ", iter, " alpha = ", alpha, ", fx_new = ", fx_new
        ! update
        call update(x, alpha, r, x_new)
        call f(x_new, fx_new)

        alpha_hist(iter,1) = alpha
        alpha_hist(iter,2) = fx_new
        size_alpha_hist = size_alpha_hist + 1 

        if (wolfe1(alpha,fx,dfx_r,fx_new)) then ! Armijo condition
            sufficientdecrease = .True.
        else 
            ! decrease alpha 
            alpha = rate * alpha
        endif
        iter = iter + 1

    enddo
     
    if (iter == maxiter_backtracking) then
        print *, "ERROR : No acceptable step find in backtracking"
        print *, "iter = ", iter
        call bestAlpha(alpha_hist,size_alpha_hist)
        print *, "Choose alpha = ", alpha
    endif

    call df(x_new,dfx_new) 

endsubroutine backtracking






subroutine zoom()
 ! algorithm 3.6 p61 from Nocedal, 2006
 use parameters, only : c1, c2, Nflat, alpha, alpha_low, alpha_high,maxiter_backtracking,&
                       x, fx, dfx_r, r,x_new,fx_new,dfx_new
 implicit none
 
 double precision, dimension(Nflat) :: x_low, x_high, dfx_low, dfx_high, x_old 
 double precision :: fx_low, fx_high, fx_old

 double precision, dimension(maxiter_backtracking,2) :: alpha_hist
 integer :: size_alpha_hist
 
 logical :: stop_cond 
 integer :: iter 
 logical :: wolfe1, wolfe2
 
 double precision :: dfx_low_r, dfx_high_r, dfx_new_r
 
 ! initailisation
 stop_cond = .False.
 iter = 1
 size_alpha_hist=0
 
 call update(x,alpha_low,r, x_low)
 call f(x_low,fx_low)
 call df(x_low,dfx_low)
 
 call update(x,alpha_high,r,x_high)
 call f(x_high,fx_high)
 call df(x_high,dfx_high)
 
 if (.not.wolfe1(alpha_low,fx,dfx_r,fx_low)) then ! not Armijo condition
   print *, "WARNING: alpha_low do not respect sufficient decrease condition"
   call backtracking()
   stop_cond = .True.
 endif

 do while (.not. stop_cond)
   !print *, alpha, alpha_low, alpha_high
   
   ! interpolate alpha as the minimum of a cubic or quadratic function passing though alpha_low and alpha_high
   if (iter == 0) then
     call scalar_product(dfx_low,r,dfx_low_r, Nflat)
     !call quadratic(alpha_low,fx_low,dot_product(dfx_low,r),alpha_high, fx_high)
     call quadratic(alpha_low,fx_low,dfx_low_r,alpha_high, fx_high)
   else
     call scalar_product(dfx_low,r,dfx_low_r, Nflat)
     call scalar_product(dfx_high,r,dfx_high_r, Nflat)
     ! call cubic(alpha_low,fx_low,dot_product(dfx_low,r),alpha_high, fx_high,alpha_old,fx_old)
     !call cubic(alpha_low,fx_low,dot_product(dfx_low,r),alpha_high, fx_high,dot_product(dfx_high,r))
     call cubic(alpha_low,fx_low,dfx_low_r,alpha_high, fx_high,dfx_high_r)
   endif
 
   ! if alpha is given outside of the intervalle, we choose the middle of the intervalle
   if (alpha < alpha_low .or. alpha >= alpha_high) then
     alpha = (alpha_low + alpha_high) / 2 
   endif
 
   ! update
   call update(x,alpha,r,x_new)
   call f(x_new, fx_new)
   call df(x_new,dfx_new)
   call scalar_product(dfx_new,r,dfx_new_r, Nflat)
   
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
     
       ! if the previous condition was satisfied:
       ! since  alpha respect the decrease condition, we can reduce our intervalle to "[alpha_low, alpha_high] = [alpha, alpha_low]"
       ! then in line 375, we permut alpha_low and alpha_high to respect the inequality alpha_low < alpha_high
       ! lines are coded in that manner because the code was first inspired by Nocedal 2006
       ! if the previous condition was not satisfied:
       ! since alpha respect the decrease condition, we can reduce our intervalle and increase the lower bound
       ! and look for alpha_st in [alpha, alpha_high] 
       alpha_low = alpha
       x_low = x_new
       fx_low = fx_new
       dfx_low = dfx_new
       
       if (abs(alpha_low - alpha_high)<1e-15) then
          ! if alpha_high too close of alpha_low, stop the algorithm and select an apha respecting the Armijo condition 
          call bestAlpha(alpha_hist,size_alpha_hist)
          stop_cond = .True.
       endif
   
       
       if (alpha_high < alpha_low) then
         ! reorder alpha_high and alph_low
         call permut(alpha_high, alpha_low, x_high, x_low, fx_high, fx_low, dfx_high, dfx_low)
         stop_cond = .True.
       end if 
     
       if (alpha_low < 1e-16) then
        ! if alpha_low is too small, stop the algorithm and select an alpha respecting the Armijo condition
        call bestAlpha(alpha_hist,size_alpha_hist)
        stop_cond = .True.
       endif
       
     endif ! Wolfe2
     
   else
     
       ! in that case, we can reduce the intervalle by decreasing the upper bound
       alpha_high = alpha
       x_high = x_new
       fx_high = fx_new
       dfx_high = dfx_new
       
       if (abs(alpha_low - alpha_high)<1e-15) then
         ! if alpha_high too close of alpha_low, stop the algorithm and select the best apha from the tested ones
         call bestAlpha(alpha_hist,size_alpha_hist)
       endif
     
   
   endif ! wolfe1
   
   iter = iter + 1
 enddo

endsubroutine zoom

subroutine permut(alpha_high, alpha_low, x_high, x_low, fx_high, fx_low, dfx_high, dfx_low)

 use parameters, only : Nflat
 implicit none
 double precision :: alpha_high, alpha_low, fx_high, fx_low
 double precision, dimension(Nflat) :: x_high, x_low, dfx_high, dfx_low

 double precision :: aux
 double precision, dimension(Nflat) :: tab_aux
 
 aux = alpha_high
 alpha_high = alpha_low
 alpha_low = aux
 
 aux = fx_high
 fx_high = fx_low
 fx_low = aux

 tab_aux(:) = x_high(:)
 x_high(:) = x_low(:)
 x_low(:) = tab_aux(:)
 
 tab_aux(:) = dfx_high(:)
 dfx_high(:) = dfx_low(:)
 dfx_low(:) = tab_aux(:)

endsubroutine permut


!function lissage() ! TODO

!endfunction lissage()


subroutine update(x,alpha,p,x_new) 
 use parameters, only : Nflat
 implicit none
 double precision, dimension(Nflat) :: x, p, x_new
 double precision :: alpha 
 x_new = x + alpha * p
endsubroutine update



subroutine quadratic(xa,fxa,dfxa,xb,fxb)
 use parameters, only : alpha, Nflat    ! TODO define in parameters
 implicit none
 ! Nocedal, 2006, p.79 eq 3.58
 double precision :: fxa,fxb
 double precision :: a,b
 double precision :: xa,xb,dfxa

 a = fxb-fxa-dfxa*xb 
 b = dfxa*xb**2
 alpha = -b/(2*a)
endsubroutine quadratic


subroutine cubic_2(xa,fa,fpa,xb,fb,xc,fc)
 use parameters, only : alpha, Nflat    ! TODO define in parameters
 implicit none
 ! https://github.com/scipy/scipy/blob/v1.11.3/scipy/optimize/_linesearch.py#L183-L322, function _cubicmin
 double precision :: fa,fb,fc
 double precision :: denom, radical
 double precision :: xa,xb,xc,fpa
 double precision, dimension(2,2) :: d1
 double precision, dimension(2) :: d2, d3

 denom = (xb-xa)**2 * (xc-xa)**2 * (xb-xc)

 d1(1,1) = (xc-xa)**2
 d1(1,2) = (xa-xb)**2
 d1(2,1) = (xa-xc)**3
 d1(2,2) = (xb-xa)**3 

 d2(1) = fb - fa - fpa * (xb-xa)
 d2(2) = fc - fa - fpa * (xc-xa)

 d3 = matmul(d1,d2)
 d3(:) = d3(:) / denom 
 radical = d3(2)*d3(2)-3*d3(1)*fpa
 alpha = xa + (-d3(2) + sqrt(radical)) / (3*d3(1))
endsubroutine cubic_2

subroutine cubic(xa,fa,fpa,xb,fb,fpb)
 use parameters, only : alpha, Nflat    ! TODO define in parameters
 implicit none
 ! Nocedal, 2006, p.79 eq 3.59
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
 use parameters, only: Nflat,c1
 double precision ::  alpha, fx, fx_new, dfx_r
 ! Armijo condition: decreasing condition     
 wolfe1 = fx_new <= fx + alpha * c1 * dfx_r 
 return 
end


function wolfe2()
 use parameters, only :  Nflat,dfx_r,c2,r ,dfx_new   
 implicit none
 logical :: wolfe2
 double precision :: dfx_new_r
 ! 2nd String Wolfe condition: Curvature condition
 call scalar_product(dfx_new, r, dfx_new_r, Nflat)
 !wolfe2 = abs(dot_product(dfx_new,r)) <= c2 * abs(dfx_r)
 wolfe2 = abs(dfx_new_r) <= c2 * abs(dfx_r)
endfunction wolfe2

