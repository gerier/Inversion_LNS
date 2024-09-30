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
 double precision :: dfxdfx_old_r_old, norm2_dfxdfx_old,y_s
 
 count_grad = 0
 count_f = 0
 count_restart = 0
 

 ! initialisation 
 x(:) = m0(:) 
 call save_info_inversion(0)
 call f(x,fx)
 call df(x,dfx) 

 x_old(:) = 1e6 * x(:)
 fx_old = 1e10
 dfx_old(:) = 0
 
 iter = 1
 notFound = .False.

 if (type_gradient == 6) then
    steepest_nbiter_default = mem_lbfgs
 endif
 
 do while (iter < maxiter .and. (abs(fx-fx_old)>1e-14) )
 
  if (rank == 0) then
    print *, "[Iteration ", iter, "]"
  endif
  
  if (iter <= steepest_nbiter_default) then  ! compute the steepest descent gradient and backtracking step

    r(:) = - dfx(:)
    call scalar_product(dfx,r,dfx_r,Nflat) 
   
    ! if L-BFGS is chosen, save informations for next iterations
    if (type_gradient == 6 .and. iter /= 1) then
      Y_list(iter-1,:) = dfx - dfx_old
      S_list(iter-1,:) = x - x_old
      call scalar_product(Y_list(iter-1,:),S_list(iter-1,:),y_s,Nflat)
      RHO_list(iter-1) = 1 / y_s
    endif

   ! initialise a first step
    if (iter == 1.0d0) then
      alpha_start = 1.0d0
    else 
      alpha_start = min(1.d0, 1.01d0 * 2.d0 * (fx - fx_old) / dfx_r)
    endif
    if (alpha_start < 0) then
    	alpha_start = 1.0d0
    endif
    
    ! find a step that respects only the Armijo condition (using backtracking)
    call backtracking()
    if (rank == 0) then
      print *, "Backtracking : ", alpha
      print *, 'New fx :', fx_new, '(', rank, ')'
    endif
 
  else ! compute the conjugate gradient
  
    ! use the conjugate gradient
    ! first, check if descent direction can be applied (test if descent direction will not be equal to infinity)
    call scalar_product(dfx-dfx_old,r_old,dfxdfx_old_r_old,Nflat)
    call scalar_product(dfx-dfx_old,dfx-dfx_old,norm2_dfxdfx_old,Nflat) 
    if (((abs(dfxdfx_old_r_old))<1e-8 .and. (type_gradient>=3 .and. type_gradient<=5)) .or. sqrt(norm2_dfxdfx_old)<1e-15) then ! TODO pas logique <dfx,r> ! .or. dot_product(dfx,r) > 0 )
       r(:) = - dfx(:)
       count_restart = count_restart+1
    else 
       call get_descent_direction(type_gradient)
    endif
    
    ! initialise a first step
    ! Nocedal, 2006, eq. 3.60
    call scalar_product(dfx,r,dfx_r,Nflat)
    alpha_start = min(1.0d0, 1.01d0 * 2.d0 * (fx - fx_old) / dfx_r)
    
    if (rank == 0) then   
     print *, "Direction descent x Gradient", dfx_r  
    endif
    if (fx - fx_old >= 0) then
        r = -dfx
        call scalar_product(dfx,r,dfx_r,Nflat)
        alpha_start = 1.d0
        count_restart = count_restart+1
    endif

    if (alpha_start <= 0) then
      print *, "Alpha start < 0 :  we have accepted to increase the mistfit in the previous iteration"
      alpha_start = 1.d0
    endif
    
    ! find a step that respects the Strong Wolfe conditions
    call line_search() 
    if (rank == 0) then
    print *, "Previous cost function :", fx
    print *, "New cost function : ", fx_new
    endif
  endif
 
  ! update iterates
  x_old(:) = x(:)
  fx_old = fx
  dfx_old = dfx(:)
  r_old = r(:)
  
  x(:) = x_new(:)
  fx = fx_new    
  dfx(:) = dfx_new(:)
 
  call save_info_inversion(iter)
  call MPI_BARRIER(MPI_COMM_WORLD, code)


  iter = iter + 1  
  
  if (rank ==0) then
    print *, "==> Improvement in misfit:", fx -fx_old
  endif
  
 enddo
 
endsubroutine optimisation



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_descent_direction(conj_number)
  ! Dener A. et al., 2019
  use parameters, only : dfx, dfx_old, r, r_old,x,x_old,Nflat, mem_lbfgs, RHO_list, S_list, Y_list
  implicit none
  double precision :: beta,eta,gama,tau,rho
  integer :: conj_number,i
  double precision, dimension(Nflat) :: y,s, h0_diag,q
  double precision :: dfx_dfx, dfx_old_dfx_old, dfx_y, y_r_old, dfx_r_old,y_y,y_s, r_old_r_old, dfx_s,s_q,y_r
  double precision, dimension(mem_lbfgs) :: a
  
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
 ! algorithm 3.5 p60, Nocedal, 2006
 use parameters, only : alpha_max, x,r, fx, dfx_r,c1,c2, alpha, alpha_low, alpha_high, alpha_start,&
                        x_new, fx_new, dfx_new, Nflat,HUGEVAL
 implicit none
 
 double precision :: alpha_prec, fx_prec
 logical stop_cond
 logical :: wolfe1, wolfe2
 double precision :: dfx_new_r
 
 ! initialisation
 alpha_prec = 0
 fx_prec = HUGEVAL
 alpha = alpha_start
 
 stop_cond = .False.
 
 do while (.not.stop_cond)
 
   call update(x,alpha,r,x_new)
   call f(x_new, fx_new)
 
    
   if (.not. wolfe1(alpha,fx,dfx_r,fx_new) .or. fx_new >= fx_prec) then ! not Armijo condition 
     ! in that case, we suppose that alpha_prec respect the Armijo condition
     ! we look for a better alpha in the intervalle [alpha_prec, alpha] 
     !call get_alpha_low_high(alpha_prec, alpha, alpha_low, alpha_high)
     alpha_high = alpha
     alpha_low = alpha_prec
     call zoom()  
     stop_cond = .True.
     
   else !  Armijo condition satisfied
     print *, "LIne search :",alpha 
     call df(x_new, dfx_new) 
     call scalar_product(dfx_new, r, dfx_new_r, Nflat)
     if (wolfe2()) then ! check curvature condition
       stop_cond = .True. ! Strong Wolfe condition are satisifed 
     else
       if (dfx_new_r>0) then
          ! in that case, phi(alpha) = f(x+apha*r) is increasing
          ! it is not possible to find a larger alpha to minimise phi
          ! alpha_prec and alpha satisfy the Armijo condition, we look for the better in that intervalle
          !call get_alpha_low_high(alpha_prec, alpha, alpha_low, alpha_high)
          alpha_high = alpha_prec
          alpha_low = alpha
          call zoom()  
          stop_cond = .True.
       else 
          ! in that case, look for a larger alpha (than alpha_prec) respecting at least the Armijo condition 
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
use parameters, only : x, r, alpha, x_new,fx_new,maxiter_backtracking,rank
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
 
print *, "fx new best ", rank, fx_new
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
        ! Backtracking fails to find alpha satisfying Armijo condition 
        print *, "ERROR : No acceptable step find in backtracking"
        print *, "iter = ", iter
        ! Select the alpha given the best results (the smallest cost function)
        call bestAlpha(alpha_hist,size_alpha_hist)
        print *, "Choose alpha = ", alpha
    endif

    call df(x_new,dfx_new) 

endsubroutine backtracking






subroutine zoom()
 ! algorithm 3.6 p61 from Nocedal, 2006
 use parameters, only : c1, c2, Nflat, alpha, alpha_low, alpha_high,maxiter_backtracking,&
                       x,fx,dfx_r,r,x_new,fx_new,dfx_new,rank,TINYVAL,HUGEVAL
 implicit none
 
 double precision, dimension(Nflat) :: x_low, x_high, dfx_low, dfx_high
 double precision :: fx_low, fx_high

 double precision, dimension(maxiter_backtracking,2) :: alpha_hist
 integer :: size_alpha_hist
 
 logical :: stop_cond 
 integer :: iter 
 logical :: wolfe1, wolfe2
 
 double precision :: dfx_low_r, dfx_high_r, dfx_new_r
 
 ! initailisation
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
   !print *, alpha, alpha_low, alpha_high
   
   ! interpolate alpha as the minimum of a cubic or quadratic function passing though alpha_low and alpha_high
   if (iter == 0) then
     ! to start, choose alpha as the minimum of the quadratic function interpolated the curve, between alpha low and alpha_high 
     call scalar_product(dfx_low,r,dfx_low_r, Nflat)
     call quadratic(alpha_low,fx_low,dfx_low_r,alpha_high, fx_high)
   else
     if (modulo(iter,10) < TINYVAL) then
        alpha = (alpha_low + alpha_high)/2
     else
      ! then, choose alpha as the minimum of the cubic function interpolated the curve, between alpha low and alpha_high 
      call scalar_product(dfx_low,r,dfx_low_r, Nflat)
      call scalar_product(dfx_high,r,dfx_high_r, Nflat)
      call cubic(alpha_low,fx_low,dfx_low_r,alpha_high, fx_high,dfx_high_r)
     endif
   endif
 
   ! if alpha is given outside of the intervalle, we choose the middle of the intervalle
   if (alpha <= alpha_low .or. alpha >= alpha_high .or. alpha <= 0 .or. isnan(alpha)) then
     if (alpha_low /= 0 .and. alpha_high /= 0) then
        ! to accelerate, if alpha_low and alpha_high > 0, then use the middle point between the logarithm of alpha_low and alpha_high 
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
   
   ! init array to records the tested step alpha
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
          call df(x_new,dfx_new) 
          stop_cond = .True.
       endif
   

       !if (alpha_low < 1e-16) then
       ! ! if alpha_low is too small, stop the algorithm and select an alpha respecting the Armijo condition
       ! call bestAlpha(alpha_hist,size_alpha_hist)
       ! stop_cond = .True.
       !endif

     endif ! Wolfe2
     
   else

       ! in that case, we can reduce the intervalle by decreasing the upper bound
       alpha_high = alpha
       x_high = x_new
       fx_high = fx_new
       dfx_high = dfx_new
      
       
       if (abs(alpha_low - alpha_high)<1e-15 .or. maxiter_backtracking <= size_alpha_hist) then
         ! if alpha_high too close of alpha_low, stop the algorithm and select the best apha from the tested ones
         call bestAlpha(alpha_hist,size_alpha_hist)
         call df(x_new,dfx_new) 
         stop_cond = .True.
       endif
     
   
   endif ! wolfe1
   
   iter = iter + 1
 enddo

endsubroutine zoom


subroutine gaussian_filter(u,mask_dim)
  use parameters, only : NX_LOCAL, NY_LOCAL
  integer :: mask_dim
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: u
  double precision, dimension(-1:NY_LOCAL+2) :: u_old
  double precision:: u_aux
  integer :: j
  
 ! gaussian window (defined as a gaussian kernel image)                                       
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
  use parameters, only : NX_LOCAL, NY_LOCAL
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
  use parameters, only : NX_LOCAL, NY_LOCAL
  integer :: mask_dim
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2), intent(inout) :: u
  double precision, dimension(-1:NY_LOCAL+2) :: u_old
  integer :: j
  double precision, dimension(5) :: imageA
                                                                       
  if (mask_dim == 1) then
    u_old(:) = u(1,:)
    do j=1,NY_LOCAL
      imageA = u_old(j-2:j+2)
      call sort(imageA,5)
      u(:,j) = imageA(3)
   
    enddo
  endif
 
endsubroutine median_filter



subroutine sort(array, array_size)

 integer :: array_size
 double precision, dimension(array_size), intent(inout) :: array

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

subroutine smoothing(x) ! TODO

  use parameters, only : Nflat,c0_prior,rho0_prior,p0_prior,windx_prior,gamma_chimie, &
                             NX_LOCAL,NY_LOCAL,type_smoothing
  implicit none
  double precision, dimension(Nflat), intent(inout) :: x
  call flatmodel2priormodel(x)
  ! apply a mean, gaussian or median filter to smooth the solution x 
  
  if (type_smoothing == 1) then
  ! mean filter is applied
    call mean_filter(rho0_prior,1)
    call mean_filter(p0_prior,1)
    call mean_filter(windx_prior,1)
  elseif (type_smoothing == 2) then
  ! gaussian filter is applied
    call gaussian_filter(rho0_prior,1)
    call gaussian_filter(p0_prior,1)
    call gaussian_filter(windx_prior,1)
  elseif (type_smoothing == 3) then
  ! gaussian filter is applied
    call median_filter(rho0_prior,1)
    call median_filter(p0_prior,1)
    call median_filter(windx_prior,1)  
    
  else 
    if (type_smoothing /= 0) then
     print *, "ERROR: Type smoothing unknown"
    endif
  endif

  c0_prior(:,:)   = sqrt(gamma_chimie(:,:) * p0_prior(:,:) / rho0_prior(:,:))
  call priormodel2flatmodel(x)
 
endsubroutine smoothing


subroutine update(x,alpha,p,x_new) 
 use parameters, only : Nflat,MPI_COMM_WORLD, code
 implicit none
 double precision, dimension(Nflat) :: x, p, x_new
 double precision :: alpha 
 ! update the current solution witht the descent direction and the line search step 
 if (alpha == 0) then
   x_new = x
 else 
   x_new = x + alpha * p
   call MPI_BARRIER(MPI_COMM_WORLD, code)
 
   ! smooth the solution because of artefacts due to on-linearity at the source and at the adjoint sources
   call smoothing(x_new)
   call MPI_BARRIER(MPI_COMM_WORLD, code)
 endif
endsubroutine update


subroutine quadratic(xa,fxa,dfxa,xb,fxb)
 ! xa   : point a
 ! fxa  : f evaluated in a
 ! dfxa : gradient of f evaluated in a
 ! xb   : point b
 ! fxb  : f evaluated in b
 use parameters, only : alpha
 implicit none
 ! The curve between hthe two points is a quadratic interpretation and return the minimum value "alpha"
 ! Nocedal, 2006, p.79 eq 3.58
 double precision :: fxa,fxb
 double precision :: xa,xb,dfxa
 double precision :: denom
 
 denom = (fxb - fxa - (xb - xa) * dfxa ) / ( (xb - xa)**2)
 alpha = xa-dfxa/(2*denom)
endsubroutine quadratic


subroutine cubic(xa,fa,fpa,xb,fb,fpb)
 ! xa  : point a
 ! fxa : f evaluated in a
 ! fpa : gradient of f evaluated in a
 ! xb  : point b
 ! fxb : f evaluated in b
 ! fpb : gradient of f evaluated in b
 use parameters, only : alpha   
 implicit none
 ! The curve between hthe two points is a cubic interpretation and return the minimum value "alpha"
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

