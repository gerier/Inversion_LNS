
subroutine init_backgrounds()

use parameters
implicit none
integer :: i,j, ind, i_perturb_model

double precision :: dist2circ, x_circ, y_circ
integer :: i_min,i_max,j_min,j_max

double precision :: z, rho,  c, p, grav, w_P, gamma_ratio, &
                       dummy1, dummy2, dummy3, dummy4, dummy5, dummy6, &
                       dummy7, dummy8, dummy9, dummy10, dummy11

 
 if (atmospheric_model_file) then ! Define model from external file
  
  !!! 
  ! TRUE MODEL
  !!!
  
  ! Load atmospheric model from file
  open(2,file=atmospheric_file_name_true, status='old', action='read')
    read (2,*) ! to skip first line, with metadata
    read (2,*) ! to skip first line, with metadata
    read (2,*) ! to skip first line, with metadata
   
    ! load model line by line (or altitude by altitude)
    do j = 2, NY
      read(2,*)  z, rho, dummy1, c, p, dummy2, grav, dummy4, dummy5, dummy6, dummy7, &
                 dummy8, dummy9, w_P, dummy10, dummy11, gamma_ratio
                
      if (j_rank == (j-1)/NY_LOCAL) then
        ind= j - offset_j 
        rho0_true(:,ind)            = rho
        gamma_chimie(:,ind)         =  gamma_ratio
        p0_true(:,ind)              = p
        kappa_unrelaxed_true(:,ind) = p * gamma_ratio
        g(:,ind)                    = grav
        windx_true(:,ind)           = w_P
        windy_true(:,ind)           = 0.d0
      endif
    enddo
  close(2)

  if (method > 1) then
    !!! 
    ! PRIOR MODEL
    !!!

    ! Load atmospheric model from file
    open(2,file=atmospheric_file_name_prior, status='old', action='read')
     read (2,*) ! to skip first line, with metadata
     read (2,*) ! to skip first line, with metadata
     read (2,*) ! to skip first line, with metadata
   
     ! load model line by line (or altitude by altitude)
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

 else ! Define model from informations of parameters.f90
 
  !!! 
  ! True model
  !!!

  gamma_chimie(:,:)         = gamma_chimie_value
  rho0_true(:,:)            = density_true
  p0_true(:,:)              = rho0_true(:,:) * (cp_unrelaxed_true*cp_unrelaxed_true) / gamma_chimie(:,:)
  windx_true(:,:)           = windx_value_true
  windy_true(:,:)           = 0.d0
  g(:,:)                    = 0.d0
  
  !!! 
  ! Priori model
  !!!
  
  rho0_prior(:,:)            = density_prior
  p0_prior(:,:)              = density_prior * (cp_unrelaxed_prior * cp_unrelaxed_prior) / gamma_chimie(:,:)
  windx_prior(:,:)           = windx_value_prior
  windy_prior(:,:)           = 0.0d0
  
endif  


 ! PERTURBATION TO TRUE MODEL
 do j = 1,NY_LOCAL
    do i = 1,NX_LOCAL
    
      i_global = i + offset_i
      j_global = j + offset_j 
      
      if NPERTURB_MODEL > 0 : 
        do i_perturb_model=1,NPERTURB_MODEL
        
          ! TODO: add blur on the shape of the perturbation to have a smoother model
          
          ! perturbation of rectangular shape
          if (ADD_PERTURB_MODEL_INFO(i_perturb_model,1) == 1)then
            i_min = ADD_PERTURB_MODEL_INFO(i_perturb_model,2) / DELTAX + 1
            j_min = ADD_PERTURB_MODEL_INFO(i_perturb_model,3) / DELTAY + 1
            i_max = ADD_PERTURB_MODEL_INFO(i_perturb_model,4) / DELTAX + 1
            j_max = ADD_PERTURB_MODEL_INFO(i_perturb_model,5) / DELTAY + 1      
            if (i_global >= i_min .and. i_global < i_max .and. j_global >= j_min .and. j_global < j_max ) then
                rho0_true(i,j)            = rho0_true(i,j) * ADD_PERTURB_MODEL_INFO(i_perturb_model,6)
              p0_true(i,j)              = rho0_true(i,j) * cp_unrelaxed_true*ADD_PERTURB_MODEL_INFO(i_perturb_model,7)&
                    *cp_unrelaxed_true*ADD_PERTURB_MODEL_INFO(i_perturb_model,7) / gamma_chimie(i,j)
            endif
          
          ! perturbation of circular shape
          elseif (ADD_PERTURB_MODEL_INFO(i_perturb_model,1) == 2) then
            ! is the point in a circle 
            x_circ = ADD_PERTURB_MODEL_INFO(i_perturb_model,2)
            y_circ = ADD_PERTURB_MODEL_INFO(i_perturb_model,3)
            dist2circ = sqrt( (((i_global+0.5)*DELTAX) - x_circ)**2 + (((j_global+0.5)*DELTAY) - y_circ)**2)
            if ( dist2circ < ADD_PERTURB_MODEL_INFO(i_perturb_model,4) ) then
          rho0_true(i,j)            = rho0_true(i,j) * ADD_PERTURB_MODEL_INFO(i_perturb_model,6)
          p0_true(i,j)              = rho0_true(i,j) * cp_unrelaxed_prior*ADD_PERTURB_MODEL_INFO(i_perturb_model,7) &
          *cp_unrelaxed_prior*ADD_PERTURB_MODEL_INFO(i_perturb_model,7) / gamma_chimie(i,j)
            endif
          endif
        enddo
      endif
      
      ! to have a wind profil 
       
      if (add_windperturb_profile .and. (j_global > jmin_wind .and. j_global <= jmax_wind)) then 
        windx_true(i,j) = exp(- ((j_global-1)*DELTAX/1e3 - mean_gauss_wind)**2 / sigma2_gauss_wind ) * max_wind_factor
      endif
          
    enddo
  enddo


  ! BOUNDARY CONDITIONS
  
  ! boundary conditions at the bottom
  if (j_rank == 0 .and. (.not. USE_PML_YMIN)) then
     do j=-1,1 
       rho0_true(:,j) = rho0_true(:,3-j)
       p0_true(:,j) = p0_true(:,3-j)
       gamma_chimie(:,j) = gamma_chimie(:,3-j)
       g(:,j) = -g(:,3-j)
       windx_true(:,j) = - windx_true(:,3-j)
     enddo
  endif

  ! boundary conditions at the top
  ! TODO: extension of the model should be linear and no constant, to avoid strong variation of derivatives
  if (j_rank == NPROC_Y-1) then
    do j=1,2
    windx_true(:,NY_LOCAL+j) = windx_true(:,NY_LOCAL)
    rho0_true(:,NY_LOCAL+j) = rho0_true(:,NY_LOCAL)
    p0_true(:,NY_LOCAL+j) = p0_true(:,NY_LOCAL)
    gamma_chimie(:,NY_LOCAL+j) = gamma_chimie(:,NY_LOCAL)
    g(:,NY_LOCAL+j) = g(:,NY_LOCAL)
    enddo
  endif
  
  
    ! boundary conditions at the bottom
  if (j_rank == 0) then
     do j=-1,1 
        rho0_prior(:,j) = rho0_prior(:,3-j)
        p0_prior(:,j) = p0_prior(:,3-j)
        windx_prior(:,j) = - windx_prior(:,3-j)
     enddo
  endif  
  
  ! boundary conditions at the top
  ! TODO: extension of the model should be linear and no constant, to avoid strong variation of derivatives
  if (j_rank == NPROC_Y-1) then
    do j=1,2
    windx_prior(:,NY_LOCAL+j) = windx_prior(:,NY_LOCAL)
    rho0_prior(:,NY_LOCAL+j) = rho0_prior(:,NY_LOCAL)
    p0_prior(:,NY_LOCAL+j) = p0_prior(:,NY_LOCAL)
    enddo
  endif


  call MPI_BARRIER(MPI_COMM_WORLD, code)

  ! send information to neighboors
  call send_receive_rightleft(windx_true)
  call send_receive_rightleft(windy_true)
  call send_receive_rightleft(rho0_true)
  call send_receive_rightleft(p0_true)
  call send_receive_rightleft(gamma_chimie)
     
  call send_receive_topbottom(windx_true)
  call send_receive_topbottom(windy_true)
  call send_receive_topbottom(rho0_true)
  call send_receive_topbottom(p0_true)
  call send_receive_topbottom(gamma_chimie)
     
  call send_receive_corners(windx_true) 
  call send_receive_corners(windy_true)
  call send_receive_corners(rho0_true)
  call send_receive_corners(p0_true)
  
 
  ! send information to neighboors
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
  
  call MPI_BARRIER(MPI_COMM_WORLD, code)
    
  ! compute prior wave speed 
  c0_prior(:,:) = sqrt(gamma_chimie(:,:)*p0_prior(:,:)/rho0_prior(:,:))


  ! write the background model
  call write_background(rho0_true, kappa_unrelaxed_true, p0_true, windx_true, windy_true, gamma_chimie,g)

  
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
 ! function to avoid the source region in the inverse problem (update with a coefficient the model in this area) 
 use parameters 
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
 use parameters
 implicit none
 
  call MPI_ALLREDUCE(sum(  p0_prior(:,:)**2), normsq_p0_prior, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)
  call MPI_ALLREDUCE(sum(  rho0_prior(:,:)**2), normsq_rho0_prior, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)
  call MPI_ALLREDUCE(sum(  windx_prior(:,:)**2), normsq_windx_prior, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)
  call MPI_ALLREDUCE(sum(  windy_prior(1,:)**2), normsq_windy_prior, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)
     
   if (normsq_windx_prior == 0.0d0) then  
     normsq_windx_prior = 1.0d0
   endif
   if (normsq_windy_prior == 0.0d0) then  
     normsq_windy_prior = 0 
   endif
   
endsubroutine get_norm_apriori
