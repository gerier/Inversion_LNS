!==============================================================================
!  CPML initialization for finite-difference wave simulations
!
!  This file contains routines used to initialize the convolutional
!  perfectly matched layer (CPML) absorbing boundary conditions.
!
!  The routine computes the spatial profiles and auxiliary coefficients:
!
!    - damping profiles              d
!    - coordinate stretching factors K
!    - frequency shift parameters    alpha
!    - recursive convolution terms   a_x, a_y, b_x, b_y, c_x, c_y
!
!  Profiles are computed for both regular and staggered grid locations.
!
!  The implementation follows:
!
!    - INRIA Research Report RR-3471, Section 6.1
!    - S. Gedney, EE699 Lecture Notes, Lecture 8
!
!==============================================================================


subroutine computePML()
!-----------------------------------------------------------------------
! Initialize the Convolutional Perfectly Matched Layer (CPML) absorbing boundary.
!
! This routine computes the spatial CPML profiles along the X and Y
! directions and initializes the coefficients used in the time-stepping
! equations.
!
! The CPML parameters are computed using polynomial damping profiles:
!
!   d(x)     : damping profile
!   K(x)     : complex coordinate stretching coefficient
!   alpha(x) : frequency shift parameter
!
! Both regular grid points and staggered half-grid points are considered.
!
! References:
!   - INRIA report RR-3471, Section 6.1
!   - S. Gedney, EE699 lecture notes, Lecture 8
!
!-----------------------------------------------------------------------
 use parameters
 implicit none
 integer :: i,j

! CPML layer thickness in meters.
  thickness_PML_x = NPOINTS_PML * DELTAX
  thickness_PML_y = NPOINTS_PML * DELTAY

! Target reflection coefficient at the outer CPML boundary (INRIA report section 6.1)
! http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  Rcoef = 0.001d0

! Check that the polynomial damping order is valid.
  if (NPOWER < 1) stop 'NPOWER must be greater than or equal to 1'

! Compute the maximum damping coefficient d0 from INRIA report section 6.1 
! http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x = - (NPOWER + 1) * cp_unrelaxed_prior * log(Rcoef) / (2.d0 * thickness_PML_x)
  d0_y = - (NPOWER + 1) * cp_unrelaxed_prior * log(Rcoef) / (2.d0 * thickness_PML_y)

  if (rank == 0) then
  print *,'d0_x = ',d0_x
  print *,'d0_y = ',d0_y
  print *
  endif

  ! Initialization
  d_x(:) = ZERO
  d_x_half(:) = ZERO
  K_x(:) = 1.d0
  K_x_half(:) = 1.d0
  alpha_x(:) = ZERO
  alpha_x_half(:) = ZERO
  a_x(:) = ZERO
  a_x_half(:) = ZERO

  d_y(:) = ZERO
  d_y_half(:) = ZERO
  K_y(:) = 1.d0
  K_y_half(:) = 1.d0
  alpha_y(:) = ZERO
  alpha_y_half(:) = ZERO
  a_y(:) = ZERO
  a_y_half(:) = ZERO

! -----------------------------------------------------------------------------
! Precompute CPML coefficients in the X direction
! -----------------------------------------------------------------------------

! Position of the CPML interface along the x-axis.
! The CPML on the right starts at the domain boundary minus the CPML thickness.
  xoriginleft = thickness_PML_x
  xoriginright = (NX-1)*DELTAX - thickness_PML_x

  do i = 1,NX

! abscissa of current grid point along the damping profile
    xval = DELTAX * dble(i-1)

!---------- left edge----------
    if (USE_PML_XMIN) then

! Compute CPML coefficients at regular grid locations.
      abscissa_in_PML = xoriginleft - xval
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(i) = d0_x * abscissa_normalized**NPOWER
        ! K and alpha profiles follow Gedney's CPML formulation:
        ! S. Gedney, EE699 Lecture Notes, Lecture 8, slide 8-2.
        K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! Compute CPML coefficients at staggered half-grid locations.
      abscissa_in_PML = xoriginleft - (xval + DELTAX/2.d0)
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(i) = d0_x * abscissa_normalized**NPOWER
        ! K and alpha profiles follow Gedney's CPML formulation:
        ! S. Gedney, EE699 Lecture Notes, Lecture 8, slide 8-2.
        K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

!---------- right edge ----------
    if (USE_PML_XMAX) then

! Compute CPML coefficients at regular grid locations.
      abscissa_in_PML = xval - xoriginright
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x(i) = d0_x * abscissa_normalized**NPOWER
        ! K and alpha profiles follow Gedney's CPML formulation:
        ! S. Gedney, EE699 Lecture Notes, Lecture 8, slide 8-2.
        K_x(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! Compute CPML coefficients at staggered half-grid locations.
      abscissa_in_PML = xval + DELTAX/2.d0 - xoriginright
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_x
        d_x_half(i) = d0_x * abscissa_normalized**NPOWER
        ! K and alpha profiles follow Gedney's CPML formulation:
        ! S. Gedney, EE699 Lecture Notes, Lecture 8, slide 8-2.
        K_x_half(i) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_x_half(i) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

! Numerical correction: alpha must remain positive.
    if (alpha_x(i) < ZERO) alpha_x(i) = ZERO
    if (alpha_x_half(i) < ZERO) alpha_x_half(i) = ZERO

! Compute CPML coefficients
! Avoid division by zero outside the CPML region.
    b_x(i) = exp(- (d_x(i) / K_x(i) + alpha_x(i)) * DELTAT)
    b_x_half(i) = exp(- (d_x_half(i) / K_x_half(i) + alpha_x_half(i)) * DELTAT)

    if (abs(d_x(i)) > 1.d-6) a_x(i) = d_x(i) * (b_x(i) - 1.d0) / (K_x(i) * (d_x(i) + K_x(i) * alpha_x(i)))
    if (abs(d_x_half(i)) > 1.d-6) a_x_half(i) = d_x_half(i) * &
      (b_x_half(i) - 1.d0) / (K_x_half(i) * (d_x_half(i) + K_x_half(i) * alpha_x_half(i)))

    if (abs(d_x(i)) > 1.d-6) c_x(i) = - d_x(i) / (K_x(i) * (d_x(i) + K_x(i) * alpha_x(i)))
    if (abs(d_x_half(i)) > 1.d-6) c_x_half(i) = - d_x_half(i) / &
       (K_x_half(i) * (d_x_half(i) + K_x_half(i) * alpha_x_half(i)))

  enddo
 
! -----------------------------------------------------------------------------
! Precompute CPML coefficients in the Y direction
! -----------------------------------------------------------------------------

! Position of the CPML interface along the y-axis.
! The CPML on the top starts at the domain boundary minus the CPML thickness.
  yoriginbottom = thickness_PML_y
  yorigintop = (NY-1)*DELTAY - thickness_PML_y

  do j = 1,NY

! abscissa of current grid point along the damping profile
    yval = DELTAY * dble(j-1)

!---------- bottom edge--------------
    if (USE_PML_YMIN) then

! Compute CPML coefficients at regular grid locations.
      abscissa_in_PML = yoriginbottom - yval
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y(j) = d0_y * abscissa_normalized**NPOWER
        ! K and alpha profiles follow Gedney's CPML formulation:
        ! S. Gedney, EE699 Lecture Notes, Lecture 8, slide 8-2.
        K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! Compute CPML coefficients at staggered half-grid locations.
      abscissa_in_PML = yoriginbottom - (yval + DELTAY/2.d0)
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_half(j) = d0_y * abscissa_normalized**NPOWER
        ! K and alpha profiles follow Gedney's CPML formulation:
        ! S. Gedney, EE699 Lecture Notes, Lecture 8, slide 8-2.
        K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

!---------- top edge --------------
    if (USE_PML_YMAX) then

! Compute CPML coefficients at regular grid locations.
      abscissa_in_PML = yval - yorigintop
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y(j) = d0_y * abscissa_normalized**NPOWER
        ! K and alpha profiles follow Gedney's CPML formulation:
        ! S. Gedney, EE699 Lecture Notes, Lecture 8, slide 8-2.
        K_y(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

! Compute CPML coefficients at staggered half-grid locations.
      abscissa_in_PML = yval + DELTAY/2.d0 - yorigintop
      if (abscissa_in_PML >= ZERO) then
        abscissa_normalized = abscissa_in_PML / thickness_PML_y
        d_y_half(j) = d0_y * abscissa_normalized**NPOWER
        ! K and alpha profiles follow Gedney's CPML formulation:
        ! S. Gedney, EE699 Lecture Notes, Lecture 8, slide 8-2.
        K_y_half(j) = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
        alpha_y_half(j) = ALPHA_MAX_PML * (1.d0 - abscissa_normalized)
      endif

    endif

! Compute CPML coefficients
! Avoid division by zero outside the CPML region.
    b_y(j) = exp(- (d_y(j) / K_y(j) + alpha_y(j)) * DELTAT)
    b_y_half(j) = exp(- (d_y_half(j) / K_y_half(j) + alpha_y_half(j)) * DELTAT)

    if (abs(d_y(j)) > 1.d-6) a_y(j) = d_y(j) * (b_y(j) - 1.d0) / (K_y(j) * (d_y(j) + K_y(j) * alpha_y(j)))
    if (abs(d_y_half(j)) > 1.d-6) a_y_half(j) = d_y_half(j) * &
      (b_y_half(j) - 1.d0) / (K_y_half(j) * (d_y_half(j) + K_y_half(j) * alpha_y_half(j)))

    if (abs(d_y(j)) > 1.d-6) c_y(j) = - d_y(j) / (K_y(j) * (d_y(j) + K_y(j) * alpha_y(j)))
    if (abs(d_y_half(j)) > 1.d-6) c_y_half(j) = - d_y_half(j) / &
        (K_y_half(j) * (d_y_half(j) + K_y_half(j) * alpha_y_half(j)))
      
  enddo

! Precompute inverse coefficients to avoid divisions in the time loop.
! Multiplications are significantly cheaper than divisions on most architectures.
  one_over_K_x(:) = 1.d0 / K_x(:)
  one_over_K_x_half(:) = 1.d0 / K_x_half(:)
  one_over_K_y(:) = 1.d0 / K_y(:)
  one_over_K_y_half(:) = 1.d0 / K_y_half(:)


! -----------------------------------------------------------------------------
! Precompute CPML coefficients in X and Y directions
! -----------------------------------------------------------------------------
  ! Initialization
  one_over_Kdalpha_x(:) = 1.0d0
  one_over_Kdalpha_y(:) = 1.0d0
  one_over_Kdalpha_x_half(:) = 1.0d0
  one_over_Kdalpha_y_half(:) = 1.0d0
  

! Compute CPML coefficients
! Avoid division by zero outside the CPML region.
  do i=1,NX
    if (abs(d_x(i)) > 1e-6) one_over_Kdalpha_x(i) = 1.d0 / (K_x(i) + d_x(i) / alpha_x(i))
    if (abs(d_x_half(i)) > 1e-6)  one_over_Kdalpha_x_half(i) = 1.d0 / (K_x_half(i) + d_x_half(i) / alpha_x_half(i))
  enddo
  do j=1,NY
    if (abs(d_y(j)) > 1e-6)  one_over_Kdalpha_y(j) = 1.d0 / (K_y(j) + d_y(j) / alpha_y(j))
    if (abs(d_y_half(j)) > 1e-6) one_over_Kdalpha_y_half(j) = 1.d0 / (K_y_half(j) + d_y_half(j) / alpha_y_half(j))
  enddo
 

endsubroutine computePML
