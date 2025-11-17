subroutine get_timeshift(s1, s2, i_tmin,i_delta_tmin, shift)

 use parameters, only : DELTAT, NSTEP
 implicit none
 double precision, dimension(NSTEP) :: s1,s2
 double precision, dimension(NSTEP) :: cp_s1, cp_s2
 double precision :: shift
 integer :: i_tmin, i_delta_tmin
 double precision, allocatable, dimension(:) :: corr
 integer :: index_max, id_timeshift,i


    index_max = i_tmin + i_delta_tmin
    
    cp_s1(:) = 0.0d0
    cp_s1(i_tmin:index_max) = s1(i_tmin:index_max)
    cp_s2(:) = 0.0d0
    cp_s2(i_tmin:index_max) = s2(i_tmin:index_max)

    allocate(corr(1:2*i_delta_tmin+1))
    corr(:) = 0.0d0
    
    cp_s2 = cshift(cp_s2, i_delta_tmin+1)    
    do i=1,2*i_delta_tmin+1
	cp_s2 = cshift(cp_s2, -1) 
        corr(i) = sum(s1(i_tmin:index_max) * cp_s2(i_tmin:index_max))
    enddo

    id_timeshift = maxloc(corr,1)

    shift = (id_timeshift-1 - i_delta_tmin) * DELTAT
 
    deallocate(corr) 

endsubroutine get_timeshift




!subroutine compute_centered_dt(u_m,u_p, derivative, ONE_OVER_2_DTEMP)
!   ! DF centered, order 2
!   implicit none!
!
!   double precision :: u_m, u_p 
!   double precision :: derivative             ! output 
!   double precision :: ONE_OVER_2_DTEMP

!  derivative = 0.0d0  
!  derivative = (u_p - u_m) * ONE_OVER_2_DTEMP
!end subroutine compute_centered_dt


subroutine integral_time_prod(s1, s2, res)
 
 use parameters, only : wr, DELTAT, NSTEP
 implicit none
 double precision, dimension(1:NSTEP) :: s1,s2
 double precision :: res
 
 res = DELTAT * sum( s1 * s2 * wr)
 
endsubroutine integral_time_prod 
 
 
 
 
subroutine get_Tr_adjoint_source()

use parameters, only : NSTEP, DELTAT, ONE_OVER_DELTAT, sispressure_prior, sispressure_true, NREC, &
                       adjoint_source, wr, timeshift,t0, PI, REC_wr
implicit none

double precision, dimension(NREC) :: Nr
double precision, dimension(NSTEP,NREC) :: dtsispressure
integer :: irec, it_t0, it, i_tmin, i_delta_tmin

dtsispressure(:,:) = 0.d0
dtsispressure(2:NSTEP-1,:) = (sispressure_prior(3:NSTEP,:) - sispressure_prior(1:NSTEP-2,:)) * 0.5 * ONE_OVER_DELTAT

Nr(:) = 0.0d0
do irec=1,NREC

   wr(:) = 0.d0
   i_tmin = REC_wr(irec,1) / DELTAT
   i_delta_tmin = REC_wr(irec,2) / DELTAT
   wr(i_tmin: i_tmin+i_delta_tmin) = 1.0d0
   it_t0 = t0/DELTAT
   do it=1,it_t0
     wr(i_tmin - it) = 0.5 + 0.5 * cos(PI*it*DELTAT/t0)
     wr(it + i_tmin + i_delta_tmin) = 0.5 + 0.5 * cos(PI*it*DELTAT/t0)
   enddo

   call integral_time_prod(dtsispressure(:,irec), dtsispressure(:,irec), Nr(irec))
   Nr(irec) = -1 * Nr(irec)
   call get_timeshift(sispressure_prior(:,irec), sispressure_true(:,irec), i_tmin, i_delta_tmin, timeshift)  ! TODO check the sign, if bad sign, invere prior and true
   adjoint_source(:,irec) = timeshift * dtsispressure(:,irec) * wr(:)/Nr(irec)

enddo

endsubroutine get_Tr_adjoint_source
