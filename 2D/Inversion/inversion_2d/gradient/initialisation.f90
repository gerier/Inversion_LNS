
subroutine init_backgrounds()

use parameters

integer :: i,j, i_recset, irec_loc

double precision :: dist2circ, x_circ, y_circ
integer :: i_min,i_max,j_min,j_max

!!! Compute the true model
          rho0_true(:,:)            = density 
          p0_true(:,:)              = rho0_true(i,j) * (cp_unrelaxed)**2 / gamma_chimie
          kappa_unrelaxed_true(:,:) = gamma_chimie * p0_true(:,:)
          windx_true(:,:)           = 0
          windy_true(:,:)           = 0
          
  do j = 1,NY_LOCAL
    do i = 1,NX_LOCAL
    
      i_global = i + offset_i
      j_global = j + offset_j 
      
      !if ((i_global < IObs_start) .or. (i_global >= IObs_end) .or. (j_global < JObs_start) .or. (j_global >= JObs_end)) then 
      !  ! to have adiscontinuity in the background model
      !  ! attention, pour la comparaison
      !  !avec le code python, les bornes sont placés en Iobs-1 et non en Iobs
      !  rho0_true(i,j)            = density
      !  p0_true(i,j)              = density * cp_unrelaxed * cp_unrelaxed / gamma_chimie
      !  kappa_unrelaxed_true(i,j) = gamma_chimie * p0_true(i,j)
      !  windx_true(i,j)           = 0
      !  windy_true(i,j)           = 0
      
      do i_perturb_model=1,NPERTURB_MODEL
        
        if (ADD_PERTURB_MODEL_INFO(i_perturb_model,1) == 1)then
          i_min = ADD_PERTURB_MODEL_INFO(i_perturb_model,2) / DELTAX + 1
          j_min = ADD_PERTURB_MODEL_INFO(i_perturb_model,3) / DELTAY + 1
          i_max = ADD_PERTURB_MODEL_INFO(i_perturb_model,4) / DELTAX + 1
          j_max = ADD_PERTURB_MODEL_INFO(i_perturb_model,5) / DELTAY + 1      
          if (i_global >= i_min .and. i_global < i_max .and. j_global >= j_min .and. j_global < j_max ) then
               rho0_true(i,j)            = density * ADD_PERTURB_MODEL_INFO(i_perturb_model,6)
    		p0_true(i,j)              = density * (cp_unrelaxed*ADD_PERTURB_MODEL_INFO(i_perturb_model,7))**2  / gamma_chimie
    		kappa_unrelaxed_true(i,j) = gamma_chimie * p0_true(i,j)
          endif
        
        elseif (ADD_PERTURB_MODEL_INFO(i_perturb_model,1) == 2) then
          ! is the point in a circle 
          x_circ = ADD_PERTURB_MODEL_INFO(i_perturb_model,2)
          y_circ = ADD_PERTURB_MODEL_INFO(i_perturb_model,3)
          dist2circ = sqrt( (((i_global+0.5)*DELTAX) - x_circ)**2 + (((j_global+0.5)*DELTAY) - y_circ)**2)
          if ( dist2circ < ADD_PERTURB_MODEL_INFO(i_perturb_model,4) ) then
    		rho0_true(i,j)            = density * ADD_PERTURB_MODEL_INFO(i_perturb_model,6)
    		p0_true(i,j)              = density * (cp_unrelaxed*ADD_PERTURB_MODEL_INFO(i_perturb_model,7))**2  / gamma_chimie
    		kappa_unrelaxed_true(i,j) = gamma_chimie * p0_true(i,j)
          endif
        endif
      enddo
      
      ! to have a wind profil  
      if (add_wind_profile .and. (j_global > ymin_wind .and. j_global <= ymax_wind)) then 
        windy_true(i,j) = exp(- ((j_global-1)*0.1d0 - mean_gauss_wind)**2 / sigma2_gauss_wind ) * max_wind_factor
      endif
      
    enddo
  enddo

  call MPI_BARRIER(MPI_COMM_WORLD, code)

  call send_receive_rightleft(windx_true)
  call send_receive_rightleft(windy_true)
  call send_receive_rightleft(rho0_true)
  call send_receive_rightleft(p0_true)
     
  call send_receive_topbottom(windx_true)
  call send_receive_topbottom(windy_true)
  call send_receive_topbottom(rho0_true)
  call send_receive_topbottom(p0_true)
     
  call send_receive_corners(windx_true) 
  call send_receive_corners(windy_true)
  call send_receive_corners(rho0_true)

  ! write the background model
  call write_background(rho0_true, kappa_unrelaxed_true, p0_true, windx_true, windy_true, gamma_chimie)


!!! Compute the a priori model
  do j = 1,NY_LOCAL
    do i = 1,NX_LOCAL
        ! to have a discontinuity in the background model
        ! attention, pour la comparaison
        !avec le code python, les bornes sont placés en Iobs-1 et non en Iobs
        rho0_prior(i,j)            = density
        kappa_unrelaxed_prior(i,j) = density * cp_unrelaxed * cp_unrelaxed
        p0_prior(i,j)              = density * cp_unrelaxed * cp_unrelaxed / gamma_chimie
        windx_prior(i,j)           = 0
        windy_prior(i,j)           = 0
    enddo
  enddo

  call MPI_BARRIER(MPI_COMM_WORLD, code)

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
       
    ! write the background model
  !call write_background(rho0_prior, kappa_unrelaxed_prior, p0_prior, windx_prior, windy_prior, gamma_chimie)

  c0_prior(:,:) = sqrt(gamma_chimie * p0_prior(:,:) / rho0_prior(:,:)) 

endsubroutine init_backgrounds
 


subroutine init_source_recvrs() 

use parameters

  if (rank == 0) then
    ! print position of the source
    print *,'Position of the source:'
    print *
    print *,'x = ',xsource
    print *,'y = ',ysource
    print *

    ! define location of receivers
    print *,'There are ',nrec,' receivers'
    print *
  endif
  
  irec = 0 
  do i_recset=1,NREC_SET! irec=1,nrec
  
    if (NREC_PER_SET(i_recset) > 1) then
! this is to avoid a warning with GNU gfortran at compile time about division by zero when NREC = 1
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

! find closest grid point for each receiver
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

 use parameters !, only : factor_regul_SRdist, NX_LOCAL, NY_LOCAL, DELTAX, DELTAY,  wavefront,&
                !        type_source, f0, ISOURCE, JSOURCE,offset_i, offset_j, distance2
 implicit none
 character(len=100) :: file_name
 integer :: i,j, ii,jj
 double precision :: distance2_wind

 do j=1,NY_LOCAL
  do i=1,NX_LOCAL
    if (type_source == 1 .and. wavefront == 1) then
      distance2 = ((i + offset_i - Isource) * DELTAX)**2  
    elseif (type_source == 1 .and. wavefront == 2) then
     distance2 = ((j + offset_j - Jsource) * DELTAY)**2 
    elseif (type_source ==2 .or. type_source == 3) then
       distance2 = ((i + offset_i - Isource) * DELTAX)**2 + ((j + offset_j - Jsource) * DELTAY)**2
    endif
     factor_regul_SRdist(i+(j-1)*NX_LOCAL) = 1 + 50 * exp(- distance2/1e5)
     factor_regul_SRdist(prod_NXNY_LOCAL + i+(j-1)*NX_LOCAL) = 1 + 50 * exp(- distance2/1e5)
  enddo
  if (type_source == 1 .and. wavefront == 1) then
    distance2_wind = 0
  else 
    distance2_wind =  ((j + offset_j - Jsource) * DELTAY)**2  
  endif
  factor_regul_SRdist(2*prod_NXNY_LOCAL + j) = 1 + 50 * exp(- distance2_wind*f0/1e5)
 enddo

    
endsubroutine init_factor_regul
