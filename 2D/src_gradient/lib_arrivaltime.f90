
subroutine get_timeshift(s1, s2, i_tmin, i_delta_tmin, shift)

  use parameters, only : DELTAT, NSTEP
  implicit none

  double precision, dimension(NSTEP), intent(in) :: s1, s2
  double precision :: shift

  integer, intent(in) :: i_tmin, i_delta_tmin

  double precision, dimension(NSTEP) :: cp_s1, cp_s2
  double precision, allocatable :: corr(:)

  integer :: index_max, id_timeshift
  integer :: lag
  double precision, dimension(NSTEP) :: tmp

  index_max = i_tmin + i_delta_tmin

  ! Fenêtrage
  cp_s1 = 0.0d0
  cp_s2 = 0.0d0

  cp_s1(i_tmin:index_max) = s1(i_tmin:index_max)
  cp_s2(i_tmin:index_max) = s2(i_tmin:index_max)

  ! Allocation
  allocate(corr(-i_delta_tmin:i_delta_tmin))
  corr = 0.0d0

  ! Corrélation croisée (on décale s1 au lieu de s2)
  do lag = -i_delta_tmin, i_delta_tmin
     tmp = cshift(cp_s1, lag)
     corr(lag) = sum(tmp(i_tmin:index_max) * cp_s2(i_tmin:index_max))
  end do

  ! Trouver max
  id_timeshift = maxloc(corr, dim=1)

  ! Conversion index -> lag réel
  lag = id_timeshift - (i_delta_tmin + 1)

  shift = lag * DELTAT

  deallocate(corr)

end subroutine get_timeshift



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
 use MPI
use parameters, only : NSTEP, DELTAT, ONE_OVER_DELTAT, sispressure_prior, sispressure_true, NREC, &
                       adjoint_source, wr, timeshift,t0, PI, REC_wr, NY_LOCAL, NX_LOCAL, j_rank, i_rank, &
                       iy_rec,ix_rec, observation_from_file, sispressure_source, isource, jsource, rank, &
                       ierr, MPI_COMM_WORLD
implicit none

double precision, dimension(NREC) :: Nr
double precision, dimension(NSTEP,NREC) :: dtsispressure
double precision, dimension(NSTEP,NREC) :: d2tsispressure
double precision :: aux
integer :: irec, it_t0, it, i_tmin, i_delta_tmin
integer :: i
integer :: rank_source

dtsispressure(:,:) = 0.d0
dtsispressure(2:NSTEP-1,:) = (sispressure_prior(3:NSTEP,:) - sispressure_prior(1:NSTEP-2,:)) * 0.5 * ONE_OVER_DELTAT
d2tsispressure(2:NSTEP-1,:) = -(sispressure_prior(3:NSTEP,:) - sispressure_prior(1:NSTEP-2,:)) * 0.5 * ONE_OVER_DELTAT


rank_source = (isource-1)/NX_LOCAL + ((jsource-1)/NY_LOCAL)*NY_LOCAL
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Bcast(sispressure_source, NSTEP, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)


Nr(:) = 0.0d0
do irec=1,NREC

   ! create a flag for each timestamp : 1=in the cross-correlation window, 0=outside
   wr(:) = 0.d0
   i_tmin = REC_wr(irec,1) / DELTAT
   i_delta_tmin = REC_wr(irec,2) / DELTAT
   wr(i_tmin: i_tmin+i_delta_tmin) = 1.0d0
   ! the cross-correlation window start with the first part of a hann window, and end with the second half of the hann window 
   ! this allows to have a signal that have 0 energy at the start and end of the window  
   it_t0 = t0/DELTAT
   do it=1,it_t0
     wr(i_tmin - it) = 0.5 + 0.5 * cos(PI*it*DELTAT/t0)
     wr(it + i_tmin + i_delta_tmin) = 0.5 + 0.5 * cos(PI*it*DELTAT/t0)
   enddo

   if (i_rank == (ix_rec(irec)-1)/NX_LOCAL .and. j_rank == (iy_rec(irec)-1)/NY_LOCAL) then
     call integral_time_prod(dtsispressure(:,irec), dtsispressure(:,irec), Nr(irec))
     Nr(irec) = -1 * Nr(irec)
     ! if observation is a delay time, compute the delay time from source to receiver 
     if (observation_from_file == 2) then
       ! TODO: if arrival time has undergo a refraction, hilbert transform to do on source
       call get_timeshift(sispressure_prior(:,irec), sispressure_source, i_tmin, i_delta_tmin, timeshift) 
     ! if no observation, use an arbitrary delay time 
     elseif (observation_from_file == 3) then
        if (NREC /= 1) then
          print *, "ERROR: should have only 1 receiver if no observation"
          stop
        endif
        timeshift = 1.0
     else
        call get_timeshift(sispressure_prior(:,irec), sispressure_true(:,irec), i_tmin, i_delta_tmin, timeshift)  ! TODO check the sign, if bad sign, invere prior and true
     endif
   print *, "Timeshift between observation and synthetic: ", timeshift
   adjoint_source(:,irec) = timeshift * d2tsispressure(:,irec) * wr(:)/Nr(irec) 
   endif

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
enddo

endsubroutine get_Tr_adjoint_source
