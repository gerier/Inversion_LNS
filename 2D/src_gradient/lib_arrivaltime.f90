!==============================================================================
!> @file traveltime_adjoint.f90
!!
!! @brief Utilities for traveltime-misfit measurements and adjoint-source
!!        construction.
!!
!! This file gathers the routines required to estimate traveltime shifts
!! between waveforms and to build the corresponding adjoint sources.
!!
!!
!! Main routines:
!!   - interp1                : One-dimensional linear interpolation.
!!   - get_timeshift          : Traveltime estimation by discrete
!!                              cross-correlation.
!!   - get_timeshift_interp   : Sub-sample traveltime estimation using
!!                              interpolation.
!!   - integral_time_prod     : Weighted time integral of two signals.
!!   - get_Tr_adjoint_source  : Construction of the traveltime adjoint
!!                              source.
!!
!! References:
!!   Tromp, J., Tape, C., & Liu, Q. (2005),
!!   "Seismic tomography, adjoint methods, time reversal and banana-doughnut
!!   kernels", Geophysical Journal International, 160(1), 195-216.
!!
!==============================================================================


!==============================================================================
!                              Interpolation
!==============================================================================

subroutine interp1( xData, yData, xVal, yVal, lenData, lenVal )
!======================================================================
!> Perform one-dimensional piecewise-linear interpolation.
!!
!! This routine interpolates a tabulated function defined by
!! (xData, yData) onto a new set of abscissae xVal using piecewise
!! linear interpolation.
!! The xData and xVal arrays are assumed to be sorted in
!! ascending order.
!!
!! Arguments:
!!   xData   : Original x-coordinates.
!!   yData   : Original function values.
!!   xVal    : Target interpolation points.
!!   yVal    : Interpolated values.
!!   lenData : Number of input samples.
!!   lenVal  : Number of interpolation points.
!======================================================================
  implicit none

  integer, intent(in) :: lenData, lenVal 
  double precision, dimension(lenData), intent(in) :: xData, yData
  double precision, dimension(lenVal), intent(in) ::  xVal
  double precision, dimension(lenVal), intent(out) :: yVal
  
  integer :: inputIndex, dataIndex
 
  !inputIndex = 0
  !do dataIndex = 1, lenData
  !    do while ((inputIndex <= lenVal) .and. (Xval(inputIndex) < xData(dataIndex+1)))
  !      inputIndex = inputIndex + 1 
  !      yVal(inputIndex) = yData(dataIndex) + (Xval(inputIndex) - xData(dataIndex)) &
  !              * (YData(dataIndex+1) - YData(dataIndex)) / (XData(dataIndex+1) - xData(dataIndex))
  !    enddo
  !end do
  dataIndex = 1
  do inputIndex = 1, lenVal

     ! Advance dataIndex until the current interpolation point
     ! lies within the interval [xData(dataIndex), xData(dataIndex+1)].
     ! The xVal array is assumed to be sorted in ascending order.
     do while ( (dataIndex < lenData - 1) .and. (xVal(inputIndex) > xData(dataIndex+1)) )
        dataIndex = dataIndex + 1
     end do

     yVal(inputIndex) = yData(dataIndex) + (xVal(inputIndex) - xData(dataIndex)) &
             * (yData(dataIndex+1) - yData(dataIndex)) / (xData(dataIndex+1) - xData(dataIndex))

  end do
end subroutine


!==============================================================================
!                         Traveltime estimation
!==============================================================================

subroutine get_timeshift(s1, s2, i_tmin, i_delta_tmin, shift)
!======================================================================
!> Estimate the delay between two signals.
!!
!! The delay is obtained by maximizing the discrete cross-correlation
!! between the two signals over a prescribed time window.
!!
!! Arguments:
!!   s1            : First signal.
!!   s2            : Second signal.
!!   i_tmin        : Start index of the correlation window.
!!   i_delta_tmin  : Number of samples in the correlation window.
!!   shift         : Estimated time shift.
!======================================================================
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

  ! Window the input signals. 
  cp_s1 = 0.0d0
  cp_s2 = 0.0d0

  cp_s1(i_tmin:index_max) = s1(i_tmin:index_max)
  cp_s2(i_tmin:index_max) = s2(i_tmin:index_max)

  ! Allocate temporary arrays.
  allocate(corr(-i_delta_tmin:i_delta_tmin))
  corr = 0.0d0

  ! Compute the cross-correlation.
  do lag = -i_delta_tmin, i_delta_tmin
     tmp = cshift(cp_s1, lag)
     corr(lag) = sum(tmp(i_tmin:index_max) * cp_s2(i_tmin:index_max))
  end do

  ! Find the maximum correlation.
  id_timeshift = maxloc(corr, dim=1)

  ! Convert the array index to the corresponding lag.
  lag = id_timeshift - (i_delta_tmin + 1)

  shift = lag * DELTAT

  deallocate(corr)

end subroutine get_timeshift


subroutine get_timeshift_interp(s1, s2, i_tmin, i_delta_tmin, shift)
!======================================================================
!> Estimate the delay between two signals using time-interpolated signals.
!!
!! The signals are first interpolated onto a refined time grid using
!! linear interpolation. The time shift is then estimated by maximizing
!! the cross-correlation on the refined grid, providing sub-sample
!! accuracy.
!!
!! Arguments:
!!   s1            : First signal.
!!   s2            : Second signal.
!!   i_tmin        : Start index of the correlation window.
!!   i_delta_tmin  : Number of samples in the correlation window.
!!   shift         : Estimated time shift.
!======================================================================
  use parameters, only : DELTAT, NSTEP
  implicit none
  double precision, dimension(NSTEP), intent(in) :: s1, s2
  double precision, intent(out) :: shift
  integer, intent(in) :: i_tmin, i_delta_tmin

  double precision, parameter :: REFINE = 100.0d0
  double precision :: DELTAT_FINE
  integer :: NSTEP_FINE

  ! Define the coarse time window.
  integer :: i_start, i_end, n_coarse_win

  ! Arrays defined on the coarse time window.
  double precision, allocatable :: t_coarse_win(:), s1_coarse_win(:), s2_coarse_win(:)

  ! Arrays defined on the refined time window.
  double precision, allocatable :: t_fine_win(:), s1_fine_win(:), s2_fine_win(:)

  ! Temporary arrays used for the cross-correlation.
  double precision, allocatable :: cp_s1(:), cp_s2(:), tmp(:), corr(:)

  ! Indices for the refined grid 
  integer :: i_tmin_fine, i_delta_tmin_fine, index_max_fine
  integer :: i, lag, id_timeshift

  ! ! Define the refined time step
  DELTAT_FINE = DELTAT / REFINE

  ! Define the coarse time window, including a margin on both sides.
  i_start      = max(1, i_tmin - i_delta_tmin)
  i_end        = min(NSTEP, i_tmin + 2 * i_delta_tmin)
  n_coarse_win = i_end - i_start + 1
  NSTEP_FINE   = int(n_coarse_win * REFINE)

  allocate(t_coarse_win(n_coarse_win))
  allocate(s1_coarse_win(n_coarse_win))
  allocate(s2_coarse_win(n_coarse_win))
  allocate(t_fine_win(NSTEP_FINE))
  allocate(s1_fine_win(NSTEP_FINE))
  allocate(s2_fine_win(NSTEP_FINE))

  ! Extract the coarse-grid time window.
  do i = 1, n_coarse_win
    t_coarse_win(i)  = (i_start + i - 2) * DELTAT   ! Absolute time.
    s1_coarse_win(i) = s1(i_start + i - 1)
    s2_coarse_win(i) = s2(i_start + i - 1)
  end do

  ! Build the refined time vector.
  do i = 1, NSTEP_FINE
    t_fine_win(i) = t_coarse_win(1) + (i - 1) * DELTAT_FINE
  end do

  ! Interpolate both signals onto the refined grid.
  call interp1(t_coarse_win, s1_coarse_win, t_fine_win, s1_fine_win, n_coarse_win, NSTEP_FINE)
  call interp1(t_coarse_win, s2_coarse_win, t_fine_win, s2_fine_win, n_coarse_win, NSTEP_FINE)

  ! Convert the global window index into the corresponding
  ! local index on the refined time grid.
  ! i_tmin is absolute, i_start is the starting index of the window
  ! → local offset = (i_tmin - i_start) * REFINE
  i_tmin_fine       = int((i_tmin - i_start) * REFINE) + 1
  i_delta_tmin_fine = int(i_delta_tmin * REFINE)
  index_max_fine    = i_tmin_fine + i_delta_tmin_fine

  ! Compute the cross-correlation on the refined grid.
  allocate(cp_s1(NSTEP_FINE))
  allocate(cp_s2(NSTEP_FINE))
  allocate(tmp(NSTEP_FINE))
  allocate(corr(-i_delta_tmin_fine:i_delta_tmin_fine))

  cp_s1 = 0.0d0
  cp_s2 = 0.0d0
  cp_s1(i_tmin_fine:index_max_fine) = s1_fine_win(i_tmin_fine:index_max_fine)
  cp_s2(i_tmin_fine:index_max_fine) = s2_fine_win(i_tmin_fine:index_max_fine)

  corr = 0.0d0
  do lag = -i_delta_tmin_fine, i_delta_tmin_fine
    tmp = cshift(cp_s1, lag)
    corr(lag) = sum(tmp(i_tmin_fine:index_max_fine) * cp_s2(i_tmin_fine:index_max_fine))
  end do

  ! Convert the array index to the corresponding lag.
  id_timeshift = maxloc(corr, dim=1)
  lag = id_timeshift - (i_delta_tmin_fine + 1)
  shift = lag * DELTAT_FINE

  deallocate(t_coarse_win, s1_coarse_win, s2_coarse_win)
  deallocate(t_fine_win, s1_fine_win, s2_fine_win)
  deallocate(cp_s1, cp_s2, tmp, corr)

end subroutine get_timeshift_interp


subroutine integral_time_prod(s1, s2, res)
 !======================================================================
!> Compute the weighted time integral of two signals.
!!
!! This routine evaluates
!!
!!    ∫ s1(t) s2(t) w(t) dt
!!
!! using the discrete time step DELTAT.
!!
!! Arguments:
!!   s1  : First signal.
!!   s2  : Second signal.
!!   res : Weighted time integral.
!======================================================================
 use parameters, only : wr, DELTAT, NSTEP
 implicit none
 double precision, dimension(1:NSTEP) :: s1,s2
 double precision :: res
 
 res = DELTAT * sum( s1 * s2 * wr)
 
endsubroutine integral_time_prod 
 
 
 
!==============================================================================
!                       Adjoint-source construction
!==============================================================================
 
subroutine get_Tr_adjoint_source()
!======================================================================
!> Build the adjoint source associated with traveltime measurements.
!!
!! This routine computes the adjoint source corresponding to the
!! traveltime misfit. The time shift between observed and synthetic
!! signals is estimated through cross-correlation, and the adjoint
!! source is constructed following the standard traveltime-adjoint
!! formulation.
!======================================================================
use MPI
use parameters, only : NSTEP, DELTAT, ONE_OVER_DELTAT, sispressure_prior, sispressure_true, NREC, &
                       adjoint_source, wr, timeshift,t0, PI, REC_wr, NY_LOCAL, NX_LOCAL, j_rank, i_rank, &
                       iy_rec,ix_rec, observation_from_file, sispressure_source, isource, jsource, &
                       ierr, MPI_COMM_WORLD
implicit none

double precision, dimension(NREC) :: Nr
double precision, dimension(NSTEP,NREC) :: dtsispressure
double precision, dimension(NSTEP,NREC) :: dtsispressure_opposite

integer :: irec, it_t0, it, i_tmin, i_delta_tmin
integer :: rank_source

dtsispressure(:,:) = 0.d0
dtsispressure(2:NSTEP-1,:) = (sispressure_prior(3:NSTEP,:) - sispressure_prior(1:NSTEP-2,:)) * 0.5 * ONE_OVER_DELTAT
dtsispressure_opposite(2:NSTEP-1,:) = -(sispressure_prior(3:NSTEP,:) - sispressure_prior(1:NSTEP-2,:)) * 0.5 * ONE_OVER_DELTAT


rank_source = (isource-1)/NX_LOCAL + ((jsource-1)/NY_LOCAL)*NY_LOCAL
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Bcast(sispressure_source, NSTEP, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
call MPI_Barrier(MPI_COMM_WORLD, ierr)


Nr(:) = 0.0d0
do irec=1,NREC

   ! Create a weighting function equal to one inside the correlation window.
   wr(:) = 0.d0
   i_tmin = REC_wr(irec,1) / DELTAT
   i_delta_tmin = REC_wr(irec,2) / DELTAT
   wr(i_tmin: i_tmin+i_delta_tmin) = 1.0d0
   ! Apply half-Hann tapers at both ends of the correlation window.
   ! This ensures that the window smoothly reaches zero at both ends,
   ! thereby reducing edge effects.
   it_t0 = t0/DELTAT
   do it=1,it_t0
     wr(i_tmin - it) = 0.5 + 0.5 * cos(PI*it*DELTAT/t0)
     wr(it + i_tmin + i_delta_tmin) = 0.5 + 0.5 * cos(PI*it*DELTAT/t0)
   enddo

   if (i_rank == (ix_rec(irec)-1)/NX_LOCAL .and. j_rank == (iy_rec(irec)-1)/NY_LOCAL) then
     call integral_time_prod(dtsispressure(:,irec), dtsispressure(:,irec), Nr(irec))
     Nr(irec) = -1 * Nr(irec)
     ! If observation is a delay time, compute the delay time from source to receiver 
     if (observation_from_file == 2) then
       ! TODO: if the arrival has undergone refraction,
       ! compute the Hilbert transform of the source signal.
       call get_timeshift(sispressure_prior(:,irec) * wr(:), sispressure_source* wr(:), &
                          i_tmin - it_t0, i_delta_tmin + 2*it_t0, timeshift) 
     ! If no observations are available, prescribe an arbitrary time shift.
     elseif (observation_from_file == 3) then
        if (NREC /= 1) then
          print *, "ERROR: should have only 1 receiver if no observation"
          stop
        endif
        timeshift = 1.0
     else
        call get_timeshift(sispressure_prior(:,irec) * wr(:), sispressure_true(:,irec) * wr(:), &
                           i_tmin - it_t0, i_delta_tmin + 2*it_t0, timeshift) 
     endif
   print *, "Estimated time shift between observed and synthetic signals: ", timeshift
   adjoint_source(:,irec) = timeshift * dtsispressure_opposite(:,irec) * wr(:)/Nr(irec) 
   endif

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
enddo

endsubroutine get_Tr_adjoint_source
