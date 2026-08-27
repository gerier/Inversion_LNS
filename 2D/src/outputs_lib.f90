!==========================================================
! Routines for writing models, seismograms, kernels and images
!==========================================================

subroutine write_background(rho0, kappa_unrelaxed, p0, windx, windy, gamma_chemestry,g)
!==========================================================
! Write the background model to disk.
!
! Save the reference atmospheric state (density, pressure,
! wind and thermodynamic parameters) to disk for
! post-processing and visualization.
!==========================================================
    use mpi
    use parameters, only : NX_LOCAL,NY_LOCAL,DELTAY,NPROC_Y,rank,code
    implicit none
    integer j,rk
    double precision, dimension(-1:NX_LOCAL+2, -1:NY_LOCAL+2) :: rho0, kappa_unrelaxed, p0, windx, windy,g
    double precision :: gamma_chemestry
    
    do rk=0,NPROC_Y-1
    if (rank == rk) then
      open(unit=201,file='./OUTPUT/true_background_model.dat',status='unknown', position='append',action='write')
      do j = 1,NY_LOCAL
        write(201,*) (j-1)*DELTAY*1.0d-3, rho0(1,j), kappa_unrelaxed(1,j), p0(1,j), windx(1,j), windy(1,j), gamma_chemestry,g(1,j) 
      enddo
      close(201)
    endif
    call mpi_barrier(mpi_comm_world,code)
    enddo
    
endsubroutine write_background
  
  

subroutine write_one_seismogram(basename, irec, data, nt, DELTAT, t0, time_shift)
!==========================================================
! Write a single seismogram (one receiver, one field) to disk.
!==========================================================
  implicit none
  character(len=*), intent(in) :: basename
  integer, intent(in) :: irec, nt
  double precision, intent(in) :: data(nt)
  double precision, intent(in) :: DELTAT, t0, time_shift

  integer :: it
  character(len=100) :: file_name

  write(file_name,"(a,'_',i3.3,'.dat')") trim(basename), irec
  open(unit=11,file=file_name,status='unknown')
  do it=1,nt
    write(11,*) sngl(dble(it-1)*DELTAT - t0 + time_shift), ' ', sngl(data(it))
  enddo
  close(11)

end subroutine write_one_seismogram

subroutine write_seismograms(sisvx,sisvy,sispressure,sisrhop,nt,nrec,DELTAT,t0,type_number)
!==========================================================
! Write synthetic or observed seismograms.
!
! Save pressure, density and velocity time series recorded
! at receiver locations.
!
! type_number = 1: synthetic data
! type_number = 2: observed data
!==========================================================
  use MPI
  use parameters, only : ix_rec,iy_rec,NX_LOCAL,NY_LOCAL,i_rank,j_rank, code
  implicit none

  integer nt,nrec
  double precision DELTAT,t0
  double precision sisvx(nt,nrec)
  double precision sisvy(nt,nrec)
  double precision sispressure(nt,nrec)
  double precision sisrhop(nt,nrec)
  integer irec
  integer :: type_number
  character(len=10) :: suffix

  call mpi_barrier(mpi_comm_world, code)

  if (type_number == 1) then
    suffix = ''
  else if (type_number == 2) then
    suffix = '_obs'
  endif

  do irec=1,nrec
    if (i_rank == (ix_rec(irec)-1)/NX_LOCAL .and. j_rank == (iy_rec(irec)-1)/NY_LOCAL) then
    
      ! Pressure and density are staggered in time by DELTAT/2 with respect to velocity
      ! (cf. eq. 13, Robertsson, Blanch and Symes, Geophysics, vol. 59(9), 1994)
      call write_one_seismogram('./OUTPUT/pressure_file'//trim(suffix), irec, &
                                 sispressure(:,irec), nt, DELTAT, t0, DELTAT/2.d0)

      call write_one_seismogram('./OUTPUT/density_file'//trim(suffix), irec, &
                                 sisrhop(:,irec), nt, DELTAT, t0, DELTAT/2.d0)

      call write_one_seismogram('./OUTPUT/Vx_file'//trim(suffix), irec, &
                                 sisvx(:,irec), nt, DELTAT, t0, 0.d0)

      call write_one_seismogram('./OUTPUT/Vy_file'//trim(suffix), irec, &
                                 sisvy(:,irec), nt, DELTAT, t0, 0.d0)
    
    endif
  enddo
    
end subroutine write_seismograms



subroutine write_seismograms_inversion(sispressure,nt,nrec,DELTAT,t0, it_inv)
!==========================================================
! Write inversion seismograms
!
! Save pressure seismograms associated with the current
! inversion iteration.
!==========================================================
  use MPI
  use parameters, only : ix_rec,iy_rec,NX_LOCAL,NY_LOCAL,i_rank,j_rank, code
  implicit none

  integer nt,nrec,it_inv
  double precision DELTAT,t0


  double precision sispressure(nt,nrec)


  integer irec,it

  character(len=100) file_name

  call mpi_barrier(mpi_comm_world, code)

  do irec=1,nrec
    if (i_rank == (ix_rec(irec)-1)/NX_LOCAL .and. j_rank == (iy_rec(irec)-1)/NY_LOCAL) then
    
      ! pressure
      write(file_name,"('./OUTPUT_INVERSION/pressure_file_',i3.3,'_',i6.6,'.dat')") irec, it_inv
      open(unit=11,file=file_name,status='unknown')
      do it=1,nt
! In the scheme of eq (13) of Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)
! Pressure is defined at time t + DELTAT/2, i.e. staggered in time with respect to velocity.
! Account for this DELTAT/2 shift when writing the seismograms.
        write(11,*) sngl(dble(it-1)*DELTAT - t0 + DELTAT/2.d0),' ',sngl(sispressure(it,irec))
      enddo
      close(11)
    endif
  enddo

end subroutine write_seismograms_inversion



subroutine create_color_image(image_data_2D,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,field_number,type_number)
!==========================================================
! Create a color image of a scalar field.
!
! Generate a PNM image of a 2-D field with source,
! receivers and PML boundaries overlaid. 
! Conversion to GIF using ImageMagick is available but 
! currently disabled.
!==========================================================
  use parameters, only :rank, maxval_image_p,maxval_image_rho,maxval_image_vx, maxval_image_vy
  implicit none

! non linear display to enhance small amplitudes for graphics
  double precision, parameter :: POWER_DISPLAY = 0.30d0

! amplitude threshold above which we draw the color point
  double precision, parameter :: cutvect = 0.01d0

! Background color for amplitudes below the threshold.
  logical, parameter :: WHITE_BACKGROUND = .true.

! Size of the source cross and receiver markers (pixels).
  integer, parameter :: width_cross = 5, thickness_cross = 1, size_square = 3

  integer NX,NY,it,field_number,type_number,ISOURCE,JSOURCE,NPOINTS_PML,nrec
  logical USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX

! in order to be able to use a fourth-order spatial operator on the edges of the model
! here we define the array with size (0:NX+1,0:NY+1) instead of size (NX,NY) as in the second-order case
  double precision, dimension(0:NX+1,0:NY+1) :: image_data_2D
  double precision :: maxval_image

  integer, dimension(nrec) :: ix_rec,iy_rec

  integer :: ix,iy,irec

  character(len=100) :: file_name,system_command

  integer :: R, G, B

  double precision :: normalized_value,max_amplitude

! table of base names / kernel names indexed by field_number (1=Vx,2=Vy,3=pressure,4=density)
  character(len=10), dimension(4), parameter :: field_base   = &
       (/ 'Vx        ', 'Vy        ', 'pressure  ', 'density   ' /)
  character(len=10), dimension(4), parameter :: field_kernel = &
       (/ 'Kvx       ', 'Kvy       ', 'Kp        ', 'Krho      ' /)

  character(len=20) :: base_name
  
  
! open image file and create system command to convert image to more convenient format
! use the "convert" command from ImageMagick http://www.imagemagick.org

! Get the file name in function of the case.
! special case: field_number == 6 only exists for type_number == 1 (per-rank pressure dump)
  if (type_number == 1 .and. field_number == 6) then
    write(file_name,"('./OUTPUT/pressure_image_r',i6.6,'_',i6.6,'.pnm')") rank, it
    write(system_command, &
      "('convert pressure_image',i6.6,'.pnm pressure_image',i6.6,'.gif ; rm pressure_image',i6.6,'.pnm')") it,it,it
    maxval_image = maxval_image_p
    return
  endif
  
   if (field_number < 1 .or. field_number > 4) return

  select case (type_number)
    case (1)
      base_name = trim(field_base(field_number))
    case (2)
      base_name = trim(field_base(field_number))//'_obs'
    case (3)
      base_name = trim(field_kernel(field_number))
    case (4)
      base_name = trim(field_base(field_number))//'_adj'
    case default
      return
  end select

! maxval_image n'est utilisé que pour les champs "physiques" (types 1 et 2)
  if (type_number == 1 .or. type_number == 2) then
    select case (field_number)
      case (1); maxval_image = maxval_image_vx
      case (2); maxval_image = maxval_image_vy
      case (3); maxval_image = maxval_image_p
      case (4); maxval_image = maxval_image_rho
    end select
  endif

  write(file_name,"('./OUTPUT/',a,'_image',i6.6,'.pnm')") trim(base_name), it
  write(system_command, &
    "('convert ',a,'_image',i6.6,'.pnm ',a,'_image',i6.6,'.gif ; rm ',a,'_image',i6.6,'.pnm')") &
    trim(base_name),it,trim(base_name),it,trim(base_name),it
    
! Create the image.    
  open(unit=27, file=file_name, status='unknown')
  write(27,"('P3')") ! write image in PNM P3 format

  write(27,*) NX,NY ! write image size
  write(27,*) '255' ! maximum value of each pixel color

! Determine the amplitude used for normalization.
  if (maxval_image < 0d0 .or. type_number >= 3 ) then
    max_amplitude = maxval(abs(image_data_2D))
  else
    max_amplitude = maxval_image
  endif

! image starts in upper-left corner in PNM format
  do iy=NY,1,-1
    do ix=1,NX

! Normalize the field to the interval [-1,1].
! keeping in mind that amplitude can be negative
    normalized_value = image_data_2D(ix,iy) / max_amplitude

! Clip normalized values to avoid numerical edge effects.
    if (normalized_value < -1.d0) normalized_value = -1.d0
    if (normalized_value > 1.d0) normalized_value = 1.d0

! Draw the seismic source as an orange cross.
    if ((ix >= ISOURCE - width_cross .and. ix <= ISOURCE + width_cross .and. &
        iy >= JSOURCE - thickness_cross .and. iy <= JSOURCE + thickness_cross) .or. &
       (ix >= ISOURCE - thickness_cross .and. ix <= ISOURCE + thickness_cross .and. &
        iy >= JSOURCE - width_cross .and. iy <= JSOURCE + width_cross)) then
      R = 255
      G = 157
      B = 0

! Draw a black frame around the image.
  else if (ix <= 2 .or. ix >= NX-1 .or. iy <= 2 .or. iy >= NY-1) then
      R = 0
      G = 0
      B = 0

! Highlight the boundaries of the PML layers.
  else if ((USE_PML_XMIN .and. ix == NPOINTS_PML) .or. &
          (USE_PML_XMAX .and. ix == NX - NPOINTS_PML) .or. &
          (USE_PML_YMIN .and. iy == NPOINTS_PML) .or. &
          (USE_PML_YMAX .and. iy == NY - NPOINTS_PML)) then
      R = 255
      G = 150
      B = 0

! Display values below the threshold using the background color.
    else if (abs(image_data_2D(ix,iy)) <= max_amplitude * cutvect) then

! Use a black or white background for points that are below the threshold
      if (WHITE_BACKGROUND) then
        R = 255
        G = 255
        B = 255
      else
        R = 0
        G = 0
        B = 0
      endif

! Positive amplitudes are shown in red and negative amplitudes in blue.
    else if (normalized_value >= 0.d0) then
      R = nint(255.d0*normalized_value**POWER_DISPLAY)
      G = 0
      B = 0
    else
      R = 0
      G = 0
      B = nint(255.d0*abs(normalized_value)**POWER_DISPLAY)
    endif

! Draw receiver locations as dark green squares.
  do irec = 1,nrec
    if ((ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
        iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square) .or. &
       (ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
        iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square)) then
! Use dark green color
      R = 30
      G = 180
      B = 60
    endif
  enddo

! Write the RGB value of the current pixel.
    write(27,"(i3,' ',i3,' ',i3)") R,G,B

    enddo
  enddo

! Close file
  close(27)

! call the system to convert image to Gif (can be commented out if "call system" is missing in your compiler)
! call system(system_command)

  end subroutine create_color_image
  
 
 
subroutine write_one_kernel(basename, K, NX_LOCAL, NY_LOCAL, rank)
!==========================================================
! Write a single 2-D kernel array to disk, one file per MPI rank.
!==========================================================
  implicit none
  character(len=*), intent(in) :: basename
  integer, intent(in) :: NX_LOCAL, NY_LOCAL, rank
  double precision, intent(in) :: K(-1:NX_LOCAL+2,-1:NY_LOCAL+2)

  character(len=100) :: file_name

  real(kind=4) :: kernel_real_save(1:NX_LOCAL,1:NY_LOCAL)
  
  kernel_real_save= real(K(1:NX_LOCAL, 1:NY_LOCAL), kind=4)
  
  write(file_name, "('./OUTPUT/',a,'_',i6.6,'.bin')") trim(basename), rank
  open(unit=12, file=file_name, form="unformatted", access="stream", action="write")
  write(12) kernel_real_save
  close(12)
  
end subroutine write_one_kernel

subroutine write_kernels()
!==========================================================
! Write the sensitivity kernels to disk.
!
! Export the computed sensitivity kernels for all inverted
! parameters to text files.
!==========================================================
   use parameters, only : rank, K_p0, K_rho0, K_windx, K_windy, NX_LOCAL, NY_LOCAL
   implicit none

  call write_one_kernel('Kp0',    K_p0,    NX_LOCAL, NY_LOCAL, rank)
  call write_one_kernel('Krho0',  K_rho0,  NX_LOCAL, NY_LOCAL, rank)
  call write_one_kernel('Kwindx', K_windx, NX_LOCAL, NY_LOCAL, rank)
  call write_one_kernel('Kwindy', K_windy, NX_LOCAL, NY_LOCAL, rank)

endsubroutine write_kernels  
  

subroutine gather_and_generate_image(vx_loc,vy_loc,pressure_loc,rhop_loc,it,type_number)
!==========================================================
! Gather distributed fields and generate images
!
! Collect MPI subdomains on the root process and generate
! visualization images of the complete computational domain.
!==========================================================
  use MPI
  use parameters, only: NX, NY, NX_LOCAL, NY_LOCAL, &
                        NPROC, NPROC_X, code, rank,  &
                        ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
                        NPOINTS_PML, USE_PML_XMIN, USE_PML_XMAX, USE_PML_YMIN, USE_PML_YMAX
  implicit none
  
  ! Local variables
  integer :: it
  integer ii,jj,offset_ii,offset_jj,offset,rk
  integer :: num_pixel, image_split_size, i, j, type_number
  double precision, dimension(NX_LOCAL*NY_LOCAL) :: data_vx_image_loc, data_vy_image_loc, data_p_image_loc, &
        data_rho_image_loc
  double precision, dimension(NX*NY) :: data_vx_image_global, data_vy_image_global, data_p_image_global, &
        data_rho_image_global
  double precision, dimension(0:NX+1,0:NY+1) :: vx_all, vy_all, p_all, rhop_all
  double precision, dimension(-1:NX_LOCAL+2,-1:NY_LOCAL+2) :: vx_loc,vy_loc,pressure_loc,rhop_loc
   
    num_pixel = 1
    do i = 1,NX_LOCAL
        do j = 1,NY_LOCAL    
                data_vx_image_loc(num_pixel) = vx_loc(i, j)
                data_vy_image_loc(num_pixel) = vy_loc(i, j)
                data_p_image_loc(num_pixel)  = pressure_loc(i, j)
                data_rho_image_loc(num_pixel)  = rhop_loc(i, j)
                num_pixel = num_pixel + 1  
        enddo
    enddo
    
    
    image_split_size = NY_LOCAL * NX_LOCAL
   
    ! Gather the horizontal velocity field.
    call MPI_GATHER(data_vx_image_loc(:), image_split_size, MPI_DOUBLE_PRECISION, &
                 data_vx_image_global(:), image_split_size,MPI_DOUBLE_PRECISION,&
                 0, MPI_COMM_WORLD, code)
                 
    ! Gather the vertical velocity field.
    call MPI_GATHER(data_vy_image_loc(:), image_split_size, MPI_DOUBLE_PRECISION, &
                 data_vy_image_global(:), image_split_size,MPI_DOUBLE_PRECISION,&
                 0, MPI_COMM_WORLD, code)
    
    ! Gather the pressure field.
    call MPI_GATHER(data_p_image_loc(:), image_split_size, MPI_DOUBLE_PRECISION, &
                 data_p_image_global(:), image_split_size,MPI_DOUBLE_PRECISION,&
                 0, MPI_COMM_WORLD, code)
    
    ! Gather the density field.
    call MPI_GATHER(data_rho_image_loc(:), image_split_size, MPI_DOUBLE_PRECISION, &
                 data_rho_image_global(:), image_split_size,MPI_DOUBLE_PRECISION,&
                 0, MPI_COMM_WORLD, code)
       
              
    if(rank == 0) then
    
       do rk=0,NPROC-1
         offset_ii = modulo(rk,NPROC_X) * NX_LOCAL 
         offset_jj = rk / NPROC_X * NY_LOCAL
         offset = NX_LOCAL*NY_LOCAL*rk
         
         do ii = 1,NX_LOCAL
            do jj = 1,NY_LOCAL
                 vx_all(ii+offset_ii,jj+offset_jj) = data_vx_image_global(offset +(ii-1)*NY_LOCAL + jj)
                 vy_all(ii+offset_ii,jj+offset_jj) = data_vy_image_global(offset +(ii-1)*NY_LOCAL + jj)
                 p_all(ii+offset_ii,jj+offset_jj)  = data_p_image_global(offset +(ii-1)*NY_LOCAL + jj)
                 rhop_all(ii+offset_ii,jj+offset_jj)  = data_rho_image_global(offset + (ii-1)*NY_LOCAL + jj)
            enddo
        enddo
      enddo
  
        call create_color_image(vx_all,NX,NY,&
            it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,1,type_number)

        call create_color_image(vy_all,NX,NY,&
              it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,2,type_number)

        call create_color_image(p_all,NX,NY,&
              it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,3,type_number)
  
        call create_color_image(rhop_all,NX,NY,&
              it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,4,type_number)
              
  endif
end subroutine gather_and_generate_image
  
 
 subroutine write_inversion_array(basename, A, nx, ny, it, rank)
!==========================================================
! Write a 2-D array (or its first nx rows) to disk, tagged
! with iteration number and MPI rank.
!==========================================================
  implicit none
  character(len=*), intent(in) :: basename
  integer, intent(in) :: nx, ny, it, rank
  double precision, intent(in) :: A(-1:nx+2,-1:ny+2)

  character(len=100) :: file_name
 
  real(kind=4) :: a_real_save(1:nx,1:ny)
 
  a_real_save = real(A(1:nx, 1:ny), kind=4)
  
  write(file_name, "('./OUTPUT_INVERSION/',a,'_',i6.6,'_',i6.6,'.bin')") trim(basename), it, rank
  open(unit=12, file=file_name, form="unformatted", access="stream", action="write")
  write(12) a_real_save 
  close(12)

end subroutine write_inversion_array

subroutine save_info_inversion(it)
!==========================================================
! Save inversion results and diagnostics
!
! Save inversion statistics, current model parameters,
! sensitivity kernels and synthetic seismograms for the
! current optimization iteration.
!==========================================================
 use mpi
 use parameters, only : K_windx, K_p0, K_rho0, sispressure, &
                        windx_prior, p0_prior, rho0_prior, &
                        NX_LOCAL, NY_LOCAL, DELTAT, nstep, nrec, t0, &
                        fx_data, fx_regul, dfx_r, fx, fx_old, dfx, alpha, count_f, count_grad, count_restart, &
                        rank, i_rank, code
 implicit none
 integer :: it
 double precision :: dfx_dfx, dfx_dfx_local

 ! Initialize the inversion log file.
 if (it == 1 .and. rank == 0) then
   OPEN(UNIT=1222, FILE="./OUTPUT_INVERSION/iterations.txt", position="append", ACTION="write")
   WRITE(1222,*) "It, Misfit(data), Misfit(regul), Misfit(tot), Norm(D_Misfit), DF.r,",&
      "Improvement, alpha, count f, count grad, count restart"
   CLOSE(1222)
 endif
 
  ! Compute the gradient norm and append inversion statistics.
  dfx_dfx_local = dot_product(dfx,dfx)
  call MPI_BARRIER(MPI_COMM_WORLD, code)
  call MPI_ALLREDUCE(dfx_dfx_local, dfx_dfx, 1, MPI_DOUBLE_PRECISION, MPI_SUM,  MPI_COMM_WORLD, code)

 if (rank == 0) then 
  OPEN(UNIT=1222, FILE="./OUTPUT_INVERSION/iterations.txt", position="append", ACTION="write")
  WRITE(1222,*) it, fx_data, fx_regul,fx, sqrt(dfx_dfx), dfx_r, fx-fx_old, alpha, count_f, count_grad, count_restart
  CLOSE(1222)
 endif
 
 ! Write solution and gradient to disk (format: image)
 !call gather_and_generate_image(K_windx,K_windy,K_p0,K_rho0,it,3)
 !call gather_and_generate_image(windx_prior,windy_prior,p0_prior,rho0_prior,it,1)
 
 ! Write the current sensitivity kernels to disk.
  call write_inversion_array('Kwindx', K_windx, NX_LOCAL, NY_LOCAL, it, rank)
  call write_inversion_array('Kp0',    K_p0,    NX_LOCAL, NY_LOCAL, it, rank)
  call write_inversion_array('Krho0',  K_rho0,  NX_LOCAL, NY_LOCAL, it, rank)
    
 ! Write the current model parameters to disk.
  if (i_rank == 0) then
    call write_inversion_array('windx', windx_prior, 1, NY_LOCAL, it, rank)
    call write_inversion_array('p0',    p0_prior,    1, NY_LOCAL, it, rank)
    call write_inversion_array('rho0',  rho0_prior,  1, NY_LOCAL, it, rank)
  endif
 
 call write_seismograms_inversion(sispressure,nstep,nrec,DELTAT,t0, it)

end subroutine save_info_inversion

