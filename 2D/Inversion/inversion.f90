

subroutine optimisation() ! TODO

 use parameters
 implicit none
 integer :: iter
 logical notFound 
 
 count_grad = 0
 count_alpha = 0
 count_restart = 0
 
 x(:) = x0(:) 
 call f(x,fx)
 call df(x,dfx) 
 
 x_old(:) = 1e6 * x(:)
 fx_old = 0
 dfx_old(:) = 0
 
 iter = 1
 notFound = .False.
 
 do while (iter < maxiter)
  if (iter <= 2) then
  
    ! compute the steepest descent gradient
    r(:) = - dfx(:)
    dfx_r = dot_product(dfx, r)
   
   ! line search
    if (iter == 1) then
      alpha_start = 1 
    else 
      alpha_start = min(1., 1.01 * 2 * (fx - fx_old) / dfx_r)
    endif
    call backtracking()

  else 
 
    ! compute the conjugate gradient
    if (dot_product((dfx - dfx_old),r_old) < 1e-8) then
       r(:) = - dfx(:)
       count_restart = count_restart+1
    else 
       call conjugateGradient(5)
    endif
    dfx_r = dot_product(dfx, r)
    
    ! line search
    alpha_start = min(1., 1.01 * 2 * (fx - fx_old) / dfx_r)
    if (alpha_start < 0) then
      print *, "Alpha start < 0 :  we have accepted to increase the mistfit in the previous iteration"
    endif
    
    call backtracking() ! TODO : modify to have line search
    
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
  use parameters, only : dfx, dfx_old, r, r_old,x,x_old,Nflat
  implicit none
  double precision :: beta,eta,gama,tau
  integer :: conj_number
  double precision, dimension(Nflat) :: y,s
  
  if (conj_number == 1) then 
    beta = dot_product(dfx,dfx) / dot_product(dfx_old,dfx_old)
    r = - dfx + beta * r_old
  else if (conj_number == 2) then
    y = dfx - dfx_old
    s = x - x_old
    beta = dot_product(dfx,y) / dot_product(dfx_old,dfx_old)
    r = - dfx + beta * r_old
  else if (conj_number == 3) then
    y = dfx - dfx_old
    s = x - x_old
    beta = (dot_product(dfx,y) - 2*dot_product(y,y) / dot_product(y,r_old) * dot_product(dfx,r_old)) / dot_product(y,r_old)
    gama = dot_product(dfx,r_old) / dot_product(y,r_old)
    tau = dot_product(y,s) / dot_product(y,y)
    r = tau * (-dfx + gama * y + beta * r_old)
  else if (conj_number == 4) then
    y = dfx - dfx_old
    s = x - x_old
    beta = (dot_product(dfx,y) - 2*dot_product(y,y) / dot_product(y,r_old) * dot_product(dfx,r_old)) / dot_product(dfx_old,dfx_old)
    eta = - 1 / sqrt( dot_product(r_old,r_old)) / min(0.01, sqrt(dot_product(dfx_old,dfx_old)))
    beta = max(beta,eta)
    r = - dfx + beta * r_old
  else if (conj_number == 5) then
    y = dfx - dfx_old
    s = x - x_old
    beta = dot_product(dfx,y)/dot_product(y,r_old) - dot_product(y,y)/dot_product(s,y) * dot_product(dfx,s) / dot_product(y,r_old)
    beta = max(beta,0.5*dot_product(dfx,r_old)/dot_product(r_old,r_old))
    r = - dfx + beta * r_old
  endif
  
  
 
endsubroutine conjugateGradient

subroutine line_search_wolfe ! TODO

 implicit none
 
endsubroutine line_search_wolfe


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bestAlpha(alpha_hist,size_alpha_hist) ! Comment faire le array ? ne connait pas la taille 
use parameters, only : x, r, alpha, x_new, fx_new,maxiter_backtracking
implicit none
double precision :: f_best
integer :: best, i 
integer :: size_alpha_hist
double precision, dimension(maxiter_backtracking,2) :: alpha_hist

 best = 0
 f_best = alpha_hist(2,1)

 do i=2,size_alpha_hist 
   if (alpha_hist(i,2) < f_best) then ! TODO
     f_best = alpha_hist(i,2)
     best = alpha_hist(i,1)
   endif
 enddo
 alpha = best
 call update(x, alpha, r, x_new)
 fx_new = f_best
 
endsubroutine bestAlpha

subroutine backtracking() 
use parameters, only : alpha, Nflat,fx,fx_new,rate,maxiter_backtracking,x,r,x_new,fx_new,dfx_new,alpha_start
implicit none
 logical :: sufficientdecrease = .False.
 integer :: iter
 double precision, dimension(maxiter_backtracking,2) :: alpha_hist
 integer :: size_alpha_hist
 logical :: wolfe1
 
    alpha = alpha_start 
    
    alpha_hist(:,:) = 0
    size_alpha_hist = 0
    
    sufficientdecrease = .False.
    iter = 1
 
    do while ((.not. sufficientdecrease) .and. (iter < maxiter_backtracking)) 
        call update(x, alpha, r, x_new)
        call f(x_new, fx_new)

        print *, fx_new
        alpha_hist(iter,1) = alpha
        alpha_hist(iter,2) = fx_new
        size_alpha_hist = size_alpha_hist + 1 

        if (wolfe1()) then
            sufficientdecrease = .True.
        else 
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


logical function wolfe1()
 use parameters, only : dfx_r,c1,alpha,fx,fx_new     
 wolfe1 = fx_new <= fx + alpha * c1 * dfx_r 
 return 
end


function wolfe2()
 use parameters, only :  Nflat,dfx_r,c2,r ,dfx_new   
 implicit none
 logical :: wolfe2

 wolfe2 = abs(dot_product(dfx_new,r)) <= c2 * abs(dfx_r)
endfunction wolfe2



!subroutine zoom(alpha_low, alpha_high,)
!
! use parameters, only : c1, c2, f, gradf, x, fx, dfx_r, r, Nflat
! implicit none
! double precision :: alpha_low, alpha_high
 
! double precision, dimension(Nflat) :: x_low, x_high, dfx_low, dfx_high 
! double precision :: fx_low, fx_high

 
! x_low = update(x,alpha_low,r)
! fx_low = f(x_low)
! dfx_low = gradf(x_low)
 
 
! x_high = update(x,alpha_high,r)
! fx_high = f(x_high)
! dfx_high = gradf(x_high)
 
 
! array ! TODO
 
! if (.not.wolfe1()) then
 
! endif
 

!endsubroutine zoom



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
 double precision :: fxa,fxb
 double precision :: a,b
 double precision :: xa,xb,dfxa

 a = fxb-fxa-dfxa*xb 
 b = dfxa*xb**2
 alpha = -b/(2*a)
endsubroutine quadratic


subroutine cubic(xa,fa,fpa,xb,fb,xc,fc)
 use parameters, only : alpha, Nflat    ! TODO define in parameters
 implicit none
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
endsubroutine cubic


