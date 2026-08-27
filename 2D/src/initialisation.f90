
subroutine init_backgrounds()
!======================================================================
!> Initialize the background atmospheric models.
!!
!! This routine initializes the true and prior atmospheric models either
!! from external atmospheric profile files or from the parameters
!! specified in the input parameter file. It also:
!!   - applies user-defined perturbations to the true model,
!!   - applies boundary conditions,
!!   - exchanges halo values between MPI processes,
!!   - computes the prior wave speed,
!!   - writes the background model to disk.
!======================================================================
use mpi
use parameters, only : rho0_true, p0_true, kappa_unrelaxed_true, g, windx_true, windy_true, & 
                       cp_unrelaxed_true, gamma_chemestry, &
                       rho0_prior, p0_prior, kappa_unrelaxed_prior, windx_prior, windy_prior, &
                       c0_prior, cp_unrelaxed_prior, &
                       density_true, windx_value_true, gamma_chemestry_value, &
                       density_prior, windx_value_prior, &
                       nperturb_model, add_perturb_model_info, add_windperturb_profile, &
                       sigma2_gauss_wind,  jmin_wind, jmax_wind, max_wind_factor, mean_gauss_wind, &
                       atmospheric_model_file, atmospheric_file_name_prior, atmospheric_file_name_true, &
                       method, NX_LOCAL, NY_LOCAL, DELTAX, DELTAY, NY, &
                       NPROC_Y, code, j_rank, i_global, j_global, offset_i, offset_j, &
                       USE_PML_YMIN
implicit none
integer :: i,j, ind, i_perturb_model

double precision :: dist2circ, x_circ, y_circ
integer :: i_min,i_max,j_min,j_max

double precision :: z, rho,  c, p, grav, w_P, gamma_ratio, &
                       dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, &
                       dummy7, dummy8, dummy9, dummy10, dummy11

 
 if (atmospheric_model_file) then ! Define the model from an external file
  
  ! ------------------------------------------------------------------
  ! True model
  ! ------------------------------------------------------------------
  ! Load the atmospheric model from file
  open(2,file=atmospheric_file_name_true, status='old', action='read')
    read (2,*) ! to skip first line, with metadata
    read (2,*) ! to skip first line, with metadata
    read (2,*) ! to skip first line, with metadata
   
    ! Read the model line by line (one altitude level per line)
    do j = 2, NY
      read(2,*)  z, rho, dummy1, c, p, dummy2, grav, dummy4, dummy5, dummy6, dummy7, &
                 dummy8, dummy9, w_P, dummy10, dummy11, gamma_ratio
                
      if (j_rank == (j-1)/NY_LOCAL) then
        ind= j - offset_j 
        rho0_true(:,ind)            = rho
        gamma_chemestry(:,ind)      =  gamma_ratio
        p0_true(:,ind)              = p
        kappa_unrelaxed_true(:,ind) = p * gamma_ratio
        g(:,ind)                    = grav
        windx_true(:,ind)           = w_P
        windy_true(:,ind)           = 0.d0
      endif
    enddo
  close(2)

  if (method > 1) then
  ! ------------------------------------------------------------------
  ! Prior model
  ! ------------------------------------------------------------------
    ! Load atmospheric model from file
    open(2,file=atmospheric_file_name_prior, status='old', action='read')
     read (2,*) ! to skip first line, with metadata
     read (2,*) ! to skip first line, with metadata
     read (2,*) ! to skip first line, with metadata
   
     ! Read the model line by line (one altitude level per line)
     do j = 2, NY
        read(2,*)  z, rho, dummy1, c, p, dummy2, dummy3, dummy4, dummy5, dummy6, dummy7, &
                  dummy8, dummy9, w_P, dummy10, dummy11, gamma_ratio
                
        if (j_rank == (j-1)/NY_LOCAL) then
          ind=j - offset_j 
          rho0_prior(:,ind)            = rho
          p0_prior(:,ind)              = p
          kappa_unrelaxed_prior(:,ind) = p * gamma_ratio
          windx_prior(:,ind)           = w_P
          windy_prior(:,ind)           = 0.d0
        endif
     enddo
    close(2)
  endif 

 else ! Define the model from the parameters specified in the parameter file
 
  ! ------------------------------------------------------------------
  ! True model
  ! ------------------------------------------------------------------
  gamma_chemestry(:,:)         = gamma_chemestry_value
  rho0_true(:,:)            = density_true
  p0_true(:,:)              = rho0_true(:,:) * (cp_unrelaxed_true*cp_unrelaxed_true) / gamma_chemestry(:,:)
  windx_true(:,:)           = windx_value_true
  windy_true(:,:)           = 0.d0
  g(:,:)                    = 9.81d0
  
  ! ------------------------------------------------------------------
  ! Prior model
  ! ------------------------------------------------------------------
  rho0_prior(:,:)            = density_prior
  p0_prior(:,:)              = density_prior * (cp_unrelaxed_prior * cp_unrelaxed_prior) / gamma_chemestry(:,:)
  windx_prior(:,:)           = windx_value_prior
  windy_prior(:,:)           = 0.0d0
  
endif  


! Apply user-defined perturbations to the true model
 do j = 1,NY_LOCAL
    do i = 1,NX_LOCAL
    
      i_global = i + offset_i
      j_global = j + offset_j 
      
      if (NPERTURB_MODEL > 0) then 
        do i_perturb_model=1,NPERTURB_MODEL
        
          ! TODO: smooth the perturbation boundaries to avoid sharp discontinuities.
          
          ! Rectangular perturbation  
          if (ADD_PERTURB_MODEL_INFO(i_perturb_model,1) == 1)then
            i_min = ADD_PERTURB_MODEL_INFO(i_perturb_model,2) / DELTAX + 1
            j_min = ADD_PERTURB_MODEL_INFO(i_perturb_model,3) / DELTAY + 1
            i_max = ADD_PERTURB_MODEL_INFO(i_perturb_model,4) / DELTAX + 1
            j_max = ADD_PERTURB_MODEL_INFO(i_perturb_model,5) / DELTAY + 1      
            if (i_global >= i_min .and. i_global < i_max .and. j_global >= j_min .and. j_global < j_max ) then
                rho0_true(i,j)            = rho0_true(i,j) * ADD_PERTURB_MODEL_INFO(i_perturb_model,6)
              p0_true(i,j)              = rho0_true(i,j) * cp_unrelaxed_true*ADD_PERTURB_MODEL_INFO(i_perturb_model,7)&
                    *cp_unrelaxed_true*ADD_PERTURB_MODEL_INFO(i_perturb_model,7) / gamma_chemestry(i,j)
            endif
          
          ! Circular perturbation
          elseif (ADD_PERTURB_MODEL_INFO(i_perturb_model,1) == 2) then
            ! Center and radius of the circular perturbation
            x_circ = ADD_PERTURB_MODEL_INFO(i_perturb_model,2)
            y_circ = ADD_PERTURB_MODEL_INFO(i_perturb_model,3)
            dist2circ = sqrt( (((i_global+0.5)*DELTAX) - x_circ)**2 + (((j_global+0.5)*DELTAY) - y_circ)**2)
            ! Apply the perturbation inside the circle
            if ( dist2circ < ADD_PERTURB_MODEL_INFO(i_perturb_model,4) ) then 
              rho0_true(i,j)  = rho0_true(i,j) * ADD_PERTURB_MODEL_INFO(i_perturb_model,6)
              p0_true(i,j)    = rho0_true(i,j) * cp_unrelaxed_prior*ADD_PERTURB_MODEL_INFO(i_perturb_model,7) &
                           *cp_unrelaxed_prior*ADD_PERTURB_MODEL_INFO(i_perturb_model,7) / gamma_chemestry(i,j)
            endif
          endif
        enddo
      endif
      
      ! Add a prescribed wind profile
       
      if (add_windperturb_profile .and. (j_global > jmin_wind .and. j_global <= jmax_wind)) then 
        windx_true(i,j) = exp(- ((j_global-1)*DELTAY/1e3 - mean_gauss_wind)**2 / sigma2_gauss_wind ) * max_wind_factor
      endif
          
    enddo
  enddo

  ! TODO: apply spatial smoothing to perturbations to avoid sharp discontinuities
  
  
  ! BOUNDARY CONDITIONS
  
  ! Boundary conditions for the true model at the bottom
  if (j_rank == 0 .and. (.not. USE_PML_YMIN)) then
     do j=-1,1 
       rho0_true(:,j) = rho0_true(:,3-j)
       p0_true(:,j) = p0_true(:,3-j)
       gamma_chemestry(:,j) = gamma_chemestry(:,3-j)
       g(:,j) = -g(:,3-j)
       windx_true(:,j) = - windx_true(:,3-j)
     enddo
  endif

  ! Boundary conditions for the true model at the top
  if (j_rank == NPROC_Y-1) then
    do j=1,2
    windx_true(:,NY_LOCAL+j) = windx_true(:,NY_LOCAL)
    rho0_true(:,NY_LOCAL+j) = rho0_true(:,NY_LOCAL)
    p0_true(:,NY_LOCAL+j) = p0_true(:,NY_LOCAL)
    gamma_chemestry(:,NY_LOCAL+j) = gamma_chemestry(:,NY_LOCAL)
    g(:,NY_LOCAL+j) = g(:,NY_LOCAL)
    enddo
  endif
  
  ! Boundary conditions for the prior model at the bottom
  if (j_rank == 0) then
     do j=-1,1 
        rho0_prior(:,j) = rho0_prior(:,3-j)
        p0_prior(:,j) = p0_prior(:,3-j)
        windx_prior(:,j) = - windx_prior(:,3-j)
     enddo
  endif  
  
  ! Boundary conditions for the prior model at the top
  ! TODO: optional:  replace the constant extension by a linear extrapolation to reduce artificial gradients.
  if (j_rank == NPROC_Y-1) then
    do j=1,2
    windx_prior(:,NY_LOCAL+j) = windx_prior(:,NY_LOCAL)
    rho0_prior(:,NY_LOCAL+j) = rho0_prior(:,NY_LOCAL)
    p0_prior(:,NY_LOCAL+j) = p0_prior(:,NY_LOCAL)
    enddo
  endif

  ! Ensure all MPI processes have completed model initialization.
  call MPI_BARRIER(MPI_COMM_WORLD, code)

  ! Exchange halo values of the true model with neighboring MPI processes
  call send_receive_rightleft(windx_true)
  call send_receive_rightleft(windy_true)
  call send_receive_rightleft(rho0_true)
  call send_receive_rightleft(p0_true)
  call send_receive_rightleft(gamma_chemestry)
     
  call send_receive_topbottom(windx_true)
  call send_receive_topbottom(windy_true)
  call send_receive_topbottom(rho0_true)
  call send_receive_topbottom(p0_true)
  call send_receive_topbottom(gamma_chemestry)
     
  call send_receive_corners(windx_true) 
  call send_receive_corners(windy_true)
  call send_receive_corners(rho0_true)
  call send_receive_corners(p0_true)
  
 
  ! Exchange halo values of the prior model with neighboring MPI processes
  call send_receive_rightleft(windx_prior)
  call send_receive_rightleft(windy_prior)
  call send_receive_rightleft(rho0_prior)
  call send_receive_rightleft(p0_prior)
     
  call send_receive_topbottom(windx_prior)
  call send_receive_topbottom(windy_prior)
  call send_receive_topbottom(rho0_prior)
  call send_receive_topbottom(p0_prior)
     
  call send_receive_corners(windx_prior) 
  call send_receive_corners(windy_prior)
  call send_receive_corners(rho0_prior)
  call send_receive_corners(p0_prior)   
  
  ! Ensure all MPI processes have completed model initialization.
  call MPI_BARRIER(MPI_COMM_WORLD, code)
    
  ! Compute the prior-model wave speed
  c0_prior(:,:) = sqrt(gamma_chemestry(:,:)*p0_prior(:,:)/rho0_prior(:,:))

  ! Save the initialized background model for visualization and post-processing.
  call write_background(rho0_true, kappa_unrelaxed_true, p0_true, windx_true, windy_true, gamma_chemestry,g)

  
endsubroutine init_backgrounds
 

subroutine init_source_recvrs() 
!======================================================================
!> Initialize the source and receiver locations.
!!
!! This routine computes the receiver coordinates from the receiver
!! configuration, identifies the closest grid point associated with each
!! receiver, and prints the source/receiver information.
!======================================================================
use parameters, only : xsource, ysource, &
                       NREC, NREC_SET, NREC_PER_SET, REC_SET_INFO, myNREC, xrec, yrec, ix_rec, iy_rec, &
                       xspacerec, yspacerec, dist, distval, &
                       NX, NY, DELTAX, DELTAY, HUGEVAL, rank
implicit none
integer :: i,j, i_recset, irec, irec_loc


  if (rank == 0) then
    ! Print the position of the source
    print *,'Position of the source:'
    print *
    print *,'x = ',xsource
    print *,'y = ',ysource
    print *

    ! Print receiver information
    print *,'There are ',nrec,' receivers'
    print *
  endif
  
  irec = 0 
  do i_recset=1,NREC_SET! irec=1,nrec
    if (NREC_PER_SET(i_recset) > 1) then
      ! This avoids a GNU gfortran warning about division by zero when NREC = 1.
      myNREC = NREC_PER_SET(i_recset)
      xspacerec = (REC_SET_INFO(i_recset,3)-REC_SET_INFO(i_recset,1)) / dble(myNREC-1)
      yspacerec = (REC_SET_INFO(i_recset,4)-REC_SET_INFO(i_recset,2)) / dble(myNREC-1)
    else
      xspacerec = 0.d0
      yspacerec = 0.d0
    endif
   
    do irec_loc=1,NREC_PER_SET(i_recset)
     irec = irec + 1
     xrec(irec) = REC_SET_INFO(i_recset,1) + dble(irec_loc-1)*xspacerec
     yrec(irec) = REC_SET_INFO(i_recset,2) + dble(irec_loc-1)*yspacerec
    enddo
  enddo

  ! Find the closest grid point for each receiver
  do irec=1,nrec
    dist = HUGEVAL
    do j = 1,NY
    do i = 1,NX
      distval = sqrt((DELTAX*dble(i-1) - xrec(irec))**2 + (DELTAY*dble(j-1) - yrec(irec))**2)
      if (distval < dist) then
        dist = distval
        ix_rec(irec) = i
        iy_rec(irec) = j
      endif
    enddo
    enddo
    
    if (rank ==0) then
    print *,'receiver ',irec,' x_target,y_target = ',xrec(irec),yrec(irec)
    print *,'closest grid point found at distance ',dist,' in i,j = ',ix_rec(irec),iy_rec(irec)
    print *
    endif
  enddo

  
endsubroutine init_source_recvrs


subroutine init_factor_regul()
!======================================================================
!> Initialize the regularization weighting factors.
!!
!! Defines spatial weighting coefficients used during inversion to
!! penalize model updates in the vicinity of the source region.
!======================================================================
 use parameters , only : factor_regul_SRdist, distance2, NY_LOCAL, NX_LOCAL
 implicit none

 integer :: i,j
 
 do j=1,NY_LOCAL
  do i=1,NX_LOCAL
    factor_regul_SRdist(j) = 1 + 50 * exp(- distance2/1e5)
    factor_regul_SRdist(NY_LOCAL + j) = 1 + 50 * exp(- distance2/1e5)
    factor_regul_SRdist(2*NY_LOCAL + j) = 1 + 50 * exp(- distance2/1e5)
  enddo
 enddo

    
endsubroutine init_factor_regul



subroutine get_norm_apriori()
!======================================================================
!> Compute the global L2 norms of the prior model.
!!
!! Computes the distributed L2 norms of the prior pressure, density and
!! wind fields using MPI reductions. These norms are later used to
!! normalize regularization terms during inversion.
!======================================================================
  use mpi
  use parameters, only : normsq_p0_prior, normsq_rho0_prior, normsq_windx_prior, code, &
                         p0_prior, rho0_prior, windx_prior
  implicit none
 
  call MPI_ALLREDUCE(sum( p0_prior(1,:)**2), normsq_p0_prior, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)
  call MPI_ALLREDUCE(sum( rho0_prior(1,:)**2), normsq_rho0_prior, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)
  call MPI_ALLREDUCE(sum( windx_prior(1,:)**2), normsq_windx_prior, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)
  !call MPI_ALLREDUCE(sum( windy_prior(1,:)**2), normsq_windy_prior, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)
     
  if (normsq_windx_prior == 0.0d0) then  
     normsq_windx_prior = 1.0d0
  endif
  !if (normsq_windy_prior == 0.0d0) then  
  !   normsq_windy_prior = 0 
  !endif
   
endsubroutine get_norm_apriori


subroutine read_obs()
!======================================================================
!> Read observed pressure seismograms.
!!
!! Reads one pressure time series per receiver from disk and stores the
!! observations in the global array `sispressure_true`. Both whitespace-
!! and comma-separated input files are supported.
!======================================================================
use parameters, only : NREC, NSTEP, sispressure_true, path_obs_file
implicit none

double precision :: t, press
integer :: it, irec, ios, k
character(len=100) :: name_file_rec, path_obs_file_rec
character(len=200) :: line

do irec = 1, NREC

   write(name_file_rec,"(i6.6,'.dat')") irec
   path_obs_file_rec = trim(path_obs_file) // trim(name_file_rec)

   open(2, file=path_obs_file_rec, status='old', action='read')

   do it = 1, NSTEP

      read(2,'(A)', iostat=ios) line
      if (ios /= 0) exit

      ! Replace commas with spaces to allow free-format reading.
      do k = 1, len_trim(line)
         if (line(k:k) == ',') line(k:k) = ' '
      end do

      read(line, *, iostat=ios) t, press

      sispressure_true(it, irec) = press

   end do

   close(2)

end do

end subroutine read_obs
