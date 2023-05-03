module parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Paramaters 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! total number of grid points in each direction of the grid
  integer, parameter :: NY = 201
  integer, parameter :: NX = 201

! size of a grid cell
  double precision, parameter :: DELTAX = 100.d0
  double precision, parameter :: DELTAY = DELTAX

! P-velocity and density
! the unrelaxed value is the value at frequency = 0 (the relaxed value would be the value at frequency = +infinity)
  double precision, parameter :: cp_unrelaxed = 347.763977d0 
  double precision, parameter :: density = 1.13837624d0 
  double precision, parameter :: gamma_chimie = 1.4d0

! Time step in seconds.
! The CFL stability number for the O(2,2) algorithm is 1 / sqrt(2) = 0.707
! i.e. one must choose  cp * deltat / deltax < 0.707.
! For the O(2,4) algorithm used here it is a bit more restrictive,
! it is cp * deltat / deltax < 0.606  (see Levander 1988 eq (7)).
! However this only ensures that the scheme is stable. To have a scheme that is both stable and accurate,
! for O(2,4) some numerical tests show that one needs to take about half of that,
! i.e. choose deltat so that cp * deltat / deltax is equal to about 0.30 or so. (or any value below; but not above).
! Since the time scheme is only second order, this also depends on how many time steps are performed in total
! (i.e. what the value of NSTEP below is); for large values of NSTEP, of course numerical errors will start to accumulate.
  double precision, parameter :: DELTAT = 0.05d0

! total number of time steps
  integer, parameter :: NSTEP = 600

! parameters for the source
  double precision, parameter :: f0 = 0.1d0
  double precision, parameter :: t0 = 1.20d0 / f0
  double precision, parameter :: factor = 1.d0

! source (in pressure, thus at a gridpoint rather than half a grid cell away)
  double precision, parameter :: xsource = 10000.d0
  double precision, parameter :: ysource = 10000.d0
  integer, parameter :: ISOURCE = xsource / DELTAX + 1
  integer, parameter :: JSOURCE = ysource / DELTAY + 1
  double precision, parameter :: SSF_Sigma = 500.d0
  ! spread the source spatial function
  double precision :: distance2, factor_ssf
 
  double precision, parameter :: obs_xstart = 0.d0
  double precision, parameter :: obs_xend   = 21000.d0
  double precision, parameter :: obs_ystart = 11000.d0
  double precision, parameter :: obs_yend   = 21000.d0
  integer, parameter :: IObs_start = obs_xstart / DELTAX + 1
  integer, parameter :: IObs_end   = obs_xend   / DELTAX + 1
  integer, parameter :: JObs_start = obs_ystart / DELTAY + 1
  integer, parameter :: JObs_end   = obs_yend   / DELTAY + 1
  double precision, parameter :: obsfactor_rho = 1.0d0
  double precision, parameter :: obsfactor_c2 = 1.2d0
  
! receivers
  integer, parameter :: NREC = 1 !201
!! DK DK I use 2301 here instead of 2300 in order to fall exactly on a grid point
  double precision, parameter :: xdeb = 12500.d0   ! first receiver x in meters
  double precision, parameter :: ydeb = 10000.d0   ! first receiver y in meters
  double precision, parameter :: xfin = 12500.d0   ! last receiver x in meters
  double precision, parameter :: yfin = 10000.d0   ! last receiver y in meters

! display information on the screen from time to time
  integer, parameter :: IT_DISPLAY = 200

! compute some constants once and for all for the fourth-order spatial scheme
! These coefficients are given for instance by Levander, Geophysics, vol. 53(11), p. 1436, equation (A-2)
  double precision, parameter :: NINE_OVER_8_DELTAX = 9.d0 / (8.d0*DELTAX)
  double precision, parameter :: NINE_OVER_8_DELTAY = 9.d0 / (8.d0*DELTAY)
  double precision, parameter :: ONE_OVER_24_DELTAX = 1.d0 / (24.d0*DELTAX)
  double precision, parameter :: ONE_OVER_24_DELTAY = 1.d0 / (24.d0*DELTAY)
  ! compute some constants once and for all for the second-order spatial scheme
  double precision, parameter :: ONE_OVER_SIX_DELTAX = 1.d0 / (6.d0*DELTAX)
  double precision, parameter :: ONE_OVER_SIX_DELTAY = 1.d0 / (6.d0*DELTAY)
  ! Some constants for the derivation in time
  double precision, parameter :: ONE_OVER_DELTAT = 1.d0 / (DELTAT)

! value of PI
  double precision, parameter :: PI = 3.141592653589793238462643d0

! zero
  double precision, parameter :: ZERO = 0.d0

! large value for maximum
  double precision, parameter :: HUGEVAL = 1.d+30

! threshold above which we consider that the code became unstable
  double precision, parameter :: STABILITY_THRESHOLD = 1.d+25

! main arrays
! in order to be able to use a fourth-order spatial operator on the edges of the model
! here we define the arrays with size (0:NX+1,0:NY+1) instead of size (NX,NY) as in the second-order case
  double precision, dimension(0:NX+1,0:NY+1) :: & 
    pressure,rhop,vx,vy,vx_old, vy_old,         &
    kappa_unrelaxed,rho,p0,v0x,v0y,             &
    rhoa,pa,vax,vay,vax_old, vay_old,           &
    kappa_unrelaxed_obs,rho_obs, p0_obs, v0x_obs, v0y_obs                           


  ! for the source
  double precision :: a,t,source_term
  
  ! for receivers
  double precision xspacerec,yspacerec,distval,dist
  integer, dimension(NREC) :: ix_rec,iy_rec
  double precision, dimension(NREC) :: xrec,yrec
  integer :: myNREC

! for seismograms
  double precision, dimension(NSTEP,NREC) :: sisvx,sisvy,sispressure,sisrhop

  double precision, dimension(NSTEP,NREC) :: obspressure, calcpressure
  double precision norm_obs

  double precision, dimension(0:NX+1, 0:NY+1) :: Krho, Kvx, Kvy, Kp  

  ! checkpointing
  integer, parameter :: NFRAMES = 50
  double precision, dimension(0:NX+1,0:NY+1,1:4,1:NFRAMES) :: FRAMES
  logical :: use_checkpoint ! to save or not the seismograms
  
  ! precision machine
  double precision :: TINYVAL = 1e-16


end module parameters 
