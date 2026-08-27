subroutine compute_gradient()
!==============================================================================
! Compute a finite-difference approximation of the sensitivity kernel (small-perturbation method).
!
! This routine estimates the sensitivity kernel associated with the selected
! model parameter by perturbing the prior model locally and evaluating the
! resulting variation of the objective function.
!
! For each perturbation:
!   - the forward simulation is recomputed,
!   - the objective function is evaluated,
!   - a finite-difference approximation of the Fréchet derivative is stored.
!
! The resulting kernel is written to disk for post-processing.
!
!==============================================================================
  use mpi
  use parameters, only: p0_prior, rho0_prior, windx_prior, windy_prior, K, which_kernel,  &
                        NSTEP, NX_LOCAL, NY_LOCAL, NX, NSTEP, DELTAX, DELTAY, observation, &
                        save_sismos, sispressure, sispressure_true, &
                        TINYVAL, ZERO,&
                        NPOINTS_PML, &
                        code, rank, i_rank, j_rank, offset_i, offset_j
  implicit none

  integer :: i,j, ii,jj, iglob, jglob,iii,jjj
  double precision :: delta, obj_function_apriori, obj_function_perturb
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: rho0_perturb, &  
                    p0_perturb, windx_perturb, windy_perturb
  double precision:: sigma_perturb
  integer :: n_points_sigma

  character(len=100) :: file_name

  double precision :: norm_frechet_deriv

  ! For arrival-time based misfits, the perturbation is spatially smoothed
  ! using a Gaussian kernel in order to increase the timeshift due to the perturbation in model
  ! sigma_perturb controls the spatial extension of the perturbation.
  sigma_perturb = 1.5
  n_points_sigma = int(3 * sigma_perturb)

  ! Define the amplitude of the finite perturbation used for the
  ! numerical approximation of the Frechet derivative.
  ! The perturbation amplitude is chosen as a fraction of the reference
  ! model value for pressure and density, while a fixed amplitude is used
  ! for wind components.
  if (which_kernel == 1) then ! for pressure sensitivity kernel 
    delta =  0.05 * p0_prior(1,1)
  elseif (which_kernel == 2) then ! for density sensitivity kernel 
   delta =  0.05 * rho0_prior(1,1)
  elseif (which_kernel == 3 .or. which_kernel == 4) then ! for wind sensitivity kernel 
   delta = 5 
  endif
   
  ! Initialize the simulation variables.
  call reset_forward()
  
  K(:,:)  = ZERO
  
  p0_perturb(:,:) = p0_prior(:,:)
  rho0_perturb(:,:) = rho0_prior(:,:)
  windx_perturb(:,:) = windx_prior(:,:)
  windy_perturb(:,:) = windy_prior(:,:)
  
  
  ! Compute the reference wavefield using the prior (background) model.
  ! This simulation provides the reference objective function against which
  ! perturbed simulations will be compared.

  ! Compute the wavefield associated with the prior model.
  call forwardproblem(p0_prior,rho0_prior,windx_prior,windy_prior,1,NSTEP, 1) 

  ! Evaluate the objective function for the prior model.
  if (observation == 0) then
    call compute_objective_function_fullwaveform(obj_function_apriori, sispressure, sispressure_true)
  else
    call compute_objective_function_arrivaltime(obj_function_apriori, sispressure, sispressure_true)
  endif

  ! Loop over all model parameters for which the sensitivity kernel is computed.
  ! Each grid point is perturbed independently and the corresponding variation
  ! of the objective function provides an estimate of the local sensitivity.
  do jglob=6,150
   do iglob= NPOINTS_PML,NX-NPOINTS_PML
    
    ! Convert global indices to local MPI indices.
    i = iglob - offset_i
    j = jglob - offset_j
    
    ! Apply a local perturbation to the selected model parameter:
    ! pressure, density, or wind component.
    if (which_kernel == 1) then
      p0_perturb(:,:) = p0_prior(:,:)
      
      if (observation == 0) then
        p0_perturb(i,j) = p0_prior(i,j) + delta
        norm_frechet_deriv = 1.0d0
      else 
        norm_frechet_deriv = 1.0d0   
        do ii=-n_points_sigma,n_points_sigma
          do jj =-n_points_sigma,n_points_sigma
          if ( i_rank == (iglob+ii-1)/NX_LOCAL .and. j_rank == (jglob+jj-1)/NY_LOCAL) then 
            iii = iglob + ii - i_rank * NX_LOCAL
            jjj = jglob + jj - j_rank * NY_LOCAL
            p0_perturb(iii,jjj) = p0_prior(iii,jjj) + delta * exp( -1.0 * (ii**2 + jj**2)*0.5/sigma_perturb**2)
            endif
          enddo
         enddo
         call MPI_BARRIER(MPI_COMM_WORLD, code)
       endif

    elseif (which_kernel == 2) then
      rho0_perturb(:,:) = rho0_prior(:,:) 
  
      if (observation == 0) then
        rho0_perturb(i,j) = rho0_prior(i,j) + delta
        norm_frechet_deriv = 1.0d0
      else 
        norm_frechet_deriv = 1.0d0  
        do ii=-n_points_sigma,n_points_sigma
          do jj =-n_points_sigma,n_points_sigma
          if ( i_rank == (iglob+ii-1)/NX_LOCAL .and. j_rank == (jglob+jj-1)/NY_LOCAL) then  
            iii = iglob + ii - i_rank * NX_LOCAL
            jjj = jglob + jj - j_rank * NY_LOCAL
            rho0_perturb(iii,jjj) = rho0_prior(iii,jjj) + delta * exp( -1.0 * (ii**2 + jj**2)*0.5/sigma_perturb**2)
          endif
          enddo
        enddo
       call MPI_BARRIER(MPI_COMM_WORLD, code)
      endif
    elseif (which_kernel == 3) then
      windx_perturb(:,:) = windx_prior(:,:)

      if (observation == 0) then
         windx_perturb(i,j) = windx_prior(i,j) + delta
         norm_frechet_deriv = 1.0d0
      else 
        norm_frechet_deriv = 1.0d0
        do ii=-n_points_sigma,n_points_sigma
          do jj =-n_points_sigma,n_points_sigma
          if ( i_rank == (iglob+ii-1)/NX_LOCAL .and. j_rank == (jglob+jj-1)/NY_LOCAL) then  
            iii = iglob + ii - i_rank * NX_LOCAL
            jjj = jglob + jj - j_rank * NY_LOCAL
            windx_perturb(iii,jjj) = windx_prior(iii,jjj) + delta * exp( -1.0 * (ii**2 + jj**2)*0.5/sigma_perturb**2)
          endif
          enddo
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD, code)
      endif
    endif

    call MPI_BARRIER(MPI_COMM_WORLD, code)
   
    ! Exchange ghost-cell values after modifying the model parameters.
    ! This ensures that all MPI subdomains have consistent values before
    ! computing spatial derivatives during the forward simulation.
    call send_receive_rightleft(windx_perturb)
    call send_receive_rightleft(windy_perturb)
    call send_receive_rightleft(rho0_perturb)
    call send_receive_rightleft(p0_perturb)
     
    call send_receive_topbottom(windx_perturb)
    call send_receive_topbottom(windy_perturb)
    call send_receive_topbottom(rho0_perturb)
    call send_receive_topbottom(p0_perturb)
     
    call send_receive_corners(windx_perturb) 
    call send_receive_corners(windy_perturb)
    call send_receive_corners(rho0_perturb)
 
 

    ! Enable receiver recording for the perturbed simulation.
    ! The resulting synthetic data are used to evaluate the objective function.
    save_sismos= .TRUE.
    
    ! Reset all dynamic variables before running a new forward simulation.
    ! The previous wavefield must not contaminate the perturbed simulation.
    call reset_forward()
    
    ! Run the forward simulation for the perturbed model.
    call forwardproblem(p0_perturb,rho0_perturb,windx_perturb,windy_perturb, 1,NSTEP,1)
 
 
    ! Evaluate the objective function for the perturbed model.
    if (observation == 0) then
      call compute_objective_function_fullwaveform(obj_function_perturb, sispressure, sispressure_true)
    else
      call compute_objective_function_arrivaltime(obj_function_perturb, sispressure, sispressure_true)
    endif
    
    
    if (i_rank == (iglob-1)/NX_LOCAL .and. j_rank == (jglob-1)/NY_LOCAL) then 
      print *, obj_function_perturb-obj_function_apriori,obj_function_perturb, iglob,jglob
      ! Compute the finite-difference approximation of the sensitivity kernel:
      !     df/dm ≈ (f(m + dm) - f(m)) / dm
      !
      ! The kernel quantifies the sensitivity of the objective function
      ! with respect to the selected model parameter.
      K(i,j) = K(i,j) + (obj_function_perturb-obj_function_apriori)/(delta/norm_frechet_deriv)
    endif
    
    ! Periodically save the current kernel for monitoring.
    if ((mod(j,10) == 0)) then
     write(file_name, "('./OUTPUT/Kernel_true_',i6.6,'.txt')") rank
     open(unit=12, file=file_name, action="write")
     do ii=1,NX_LOCAL
       write(12,*) (K(ii,jj), jj=1,NY_LOCAL)
     enddo
     close(12)  
    endif
    
  enddo
 enddo
    
 ! Multiply K by ΔxΔy to account for the control volume (cell area in 2D) in the discrete formulation.
 K(:,:) = DELTAX * DELTAY * K(:,:)
    
 ! Save the final finite-difference kernel.
 write(file_name, "('./OUTPUT/Kernel_true_',i6.6,'.txt')") rank
 open(unit=12, file=file_name, action="write")
 do ii=1,NX_LOCAL
  write(12,*) (K(ii,jj), jj=1,NY_LOCAL)
 enddo
 close(12)

endsubroutine compute_gradient



subroutine compute_objective_function_fullwaveform(fm, signal1, signal2)
!==============================================================================
! Evaluate the waveform misfit objective function.
!
! The objective function corresponds to a normalized L2 distance between
! simulated and observed pressure waveforms:
!
! J = 1/2 * sum_r ||p_sim - p_obs||², with r a receiver 
!
! The normalization by the observed signal energy makes the contribution
! of receivers with different amplitudes comparable.
!
! Inputs:
!   signal1
!       Simulated pressure waveforms.
!   signal2
!       Observed pressure waveforms.
!
! Output:
!   fm
!       Value of the waveform objective function.
!
!==============================================================================
 use mpi
 use parameters, only : NX_LOCAL, NY_LOCAL, NSTEP, &
                        NREC, ix_rec, iy_rec, &
                        normsq_pressure_true_per_rec, &
                        code, i_rank, j_rank, &
                        DELTAT, DELTAX, DELTAY, TINYVAL
  implicit none

  integer :: irec
  double precision, intent(out) :: fm
  double precision, dimension(NSTEP,NREC), intent(in) :: signal1, signal2
  double precision, dimension(NREC) :: fm_local_per_rec

  fm = 0.d0
  fm_local_per_rec(:) = 0.0d0
   do irec=1,nrec
     if (i_rank == (ix_rec(irec))/NX_LOCAL .and. j_rank == (iy_rec(irec))/NY_LOCAL) then
       ! Avoid numerical instability for receivers with negligible signal energy.
       if ( normsq_pressure_true_per_rec(irec) > TINYVAL) then
   	fm_local_per_rec(irec) = sum((signal1(:,irec) - signal2(:,irec))**2) &
                    / normsq_pressure_true_per_rec(irec) 
       endif
     endif
  enddo

 call MPI_BARRIER(MPI_COMM_WORLD, code)
 ! Sum the local receiver contributions over all MPI processes.
 call MPI_ALLREDUCE(sum(fm_local_per_rec), fm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)

     fm = fm * 0.5d0 * DELTAT * DELTAX * DELTAY 
     
endsubroutine compute_objective_function_fullwaveform



subroutine compute_objective_function_arrivaltime(fm, signal1, signal2)
!==============================================================================
! Evaluate the arrival-time objective function.
!
! The arrival-time misfit is based on the time shift between simulated
! and observed signals inside the selected receiver window.
! The objective function is defined as:
!
! J = 1/2 * sum_r(T_r^sim - T_r^obs) with r a receiver
!
! Inputs:
!   signal1
!       Simulated pressure waveforms.
!   signal2
!       Observed pressure waveforms.
!
! Output:
!   fm
!       Value of the arrival-time objective function.
!
!==============================================================================
 use mpi 
 use parameters, only : NX_LOCAL, NY_LOCAL, DELTAT, PI, t0, &
                        REC_wr, wr, NREC, NSTEP, ix_rec, iy_rec, &
                        code, i_rank, j_rank
  implicit none

  integer :: irec, i_tmin, i_delta_tmin, it, it_t0
  double precision, intent(out) :: fm
  double precision :: shift
  double precision, dimension(NSTEP,NREC), intent(in) :: signal1, signal2
  double precision, dimension(NREC) :: fm_local_per_rec

  fm = 0.d0
  fm_local_per_rec(:) = 0.0d0

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
   
   if (i_rank == (ix_rec(irec))/NX_LOCAL .and. j_rank == (iy_rec(irec))/NY_LOCAL) then
       ! Window for the cross-correlation
       ! TODO add the cos taper  
       i_tmin = REC_wr(irec,1) / DELTAT
       i_delta_tmin = REC_wr(irec,2) / DELTAT
       ! Estimate the sub-sample travel-time shift using interpolation
       ! within the selected correlation window.
       call get_timeshift_interp(signal1(:,irec)*wr(:), signal2(:,irec)*wr(:), i_tmin-it_t0, i_delta_tmin+2*it_t0, shift)
       fm_local_per_rec(irec) = shift**2
   endif
 enddo

 call MPI_BARRIER(MPI_COMM_WORLD, code)
 call MPI_ALLREDUCE(sum(fm_local_per_rec), fm, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)

     fm = fm * 0.5d0 
   
endsubroutine compute_objective_function_arrivaltime
