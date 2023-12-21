!----
!----  save the seismograms in ASCII text format
!----


  subroutine write_background(rho0, kappa_unrelaxed, p0, windx, windy, gamma_chimie)
  
    use mpi
    use parameters, only : NX_LOCAL,NY_LOCAL,rank,NPROC_Y,code
    implicit none
    integer j,rk,ii,jj
    double precision, dimension(-1:NX_LOCAL+2, -1:NY_LOCAL+2) :: rho0, kappa_unrelaxed, p0, windx, windy
    double precision :: gamma_chimie
    character(len=100) :: file_name
    
    do rk=0,NPROC_Y-1
    if (rank == rk) then
      open(unit=201,file='./OUTPUT/true_background_model.dat',status='unknown', position='append',action='write')
      do j = 1,NY_LOCAL
        write(201,*) (j-1)*0.1d0, rho0(1,j), kappa_unrelaxed(1,j), p0(1,j), windx(1,j), windy(1,j), gamma_chimie 
      enddo
      close(201)
    endif
    call mpi_barrier(mpi_comm_world,code)
    enddo
    
    write(file_name, "('./MODELS/rho0_true_',i6.6,'.txt')") rank
    OPEN(UNIT=12, FILE=file_name, ACTION="write")
    DO ii=1,NX_LOCAL
      WRITE(12,*) (rho0(ii,jj), jj=1,NY_LOCAL)
    END DO
    CLOSE(12)
    
        write(file_name, "('./MODELS/p0_true_',i6.6,'.txt')") rank
    OPEN(UNIT=12, FILE=file_name, ACTION="write")
    DO ii=1,NX_LOCAL
      WRITE(12,*) (p0(ii,jj), jj=1,NY_LOCAL)
    END DO
    CLOSE(12)
    
        write(file_name, "('./MODELS/windx_true_',i6.6,'.txt')") rank
    OPEN(UNIT=12, FILE=file_name, ACTION="write")
    DO ii=1,NX_LOCAL
      WRITE(12,*) (windx(ii,jj), jj=1,NY_LOCAL)
    END DO
    CLOSE(12)
    
            write(file_name, "('./MODELS/windy_true_',i6.6,'.txt')") rank
    OPEN(UNIT=12, FILE=file_name, ACTION="write")
    DO ii=1,NX_LOCAL
      WRITE(12,*) (windy(ii,jj), jj=1,NY_LOCAL)
    END DO
    CLOSE(12)
    
    
  endsubroutine write_background
  
  
  
  subroutine write_seismograms(sisvx,sisvy,sispressure,sisrhop,nt,nrec,DELTAT,t0)

  use MPI
  use parameters, only : ix_rec,iy_rec,NX_LOCAL,NY_LOCAL,rank,i_rank,j_rank
  implicit none

  integer nt,nrec
  double precision DELTAT,t0

  double precision sisvx(nt,nrec)
  double precision sisvy(nt,nrec)
  double precision sispressure(nt,nrec)
  double precision sisrhop(nt,nrec)

  integer irec,it

  character(len=100) file_name

 

  do irec=1,nrec
    if (i_rank == ix_rec(irec)/NX_LOCAL .and. j_rank == iy_rec(irec)/NY_LOCAL) then
    
      ! pressure
      write(file_name,"('./OUTPUT/pressure_file_',i3.3,'.dat')") irec
      open(unit=11,file=file_name,status='unknown')
      do it=1,nt
! in the scheme of eq (13) of Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)
! pressure is defined at time t + DELTAT/2, i.e. staggered in time with respect to velocity.
! Here we must thus take this shift of DELTAT/2 into account to save the seismograms at the right time
        write(11,*) sngl(dble(it-1)*DELTAT - t0 + DELTAT/2.d0),' ',sngl(sispressure(it,irec))
      enddo
      close(11)
      
      ! density
      write(file_name,"('./OUTPUT/density_file_',i3.3,'.dat')") irec
      open(unit=11,file=file_name,status='unknown')
      do it=1,nt
        write(11,*) sngl(dble(it-1)*DELTAT - t0),' ',sngl(sisrhop(it,irec))
      enddo
      close(11)
      
      ! X component of velocity
      write(file_name,"('./OUTPUT/Vx_file_',i3.3,'.dat')") irec
      open(unit=11,file=file_name,status='unknown')
      do it=1,nt
        write(11,*) sngl(dble(it-1)*DELTAT - t0),' ',sngl(sisvx(it,irec))
      enddo
      close(11)
      
      ! Y component of velocity
      write(file_name,"('./OUTPUT/Vy_file_',i3.3,'.dat')") irec
      open(unit=11,file=file_name,status='unknown')
      do it=1,nt
        write(11,*) sngl(dble(it-1)*DELTAT - t0),' ',sngl(sisvy(it,irec))
      enddo
      close(11)
      
    endif
  enddo



  end subroutine write_seismograms

!----
!----  routine to create a color image of a given vector component
!----  the image is created in PNM format and then converted to GIF
!----

  subroutine create_color_image(image_data_2D,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,field_number,type_number)

use parameters, only :rank
  implicit none

! non linear display to enhance small amplitudes for graphics
  double precision, parameter :: POWER_DISPLAY = 0.30d0

! amplitude threshold above which we draw the color point
  double precision, parameter :: cutvect = 0.01d0

! use black or white background for points that are below the threshold
  logical, parameter :: WHITE_BACKGROUND = .true.

! size of cross and square in pixels drawn to represent the source and the receivers
  integer, parameter :: width_cross = 5, thickness_cross = 1, size_square = 3

  integer NX,NY,it,field_number,type_number,ISOURCE,JSOURCE,NPOINTS_PML,nrec
  logical USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX

! in order to be able to use a fourth-order spatial operator on the edges of the model
! here we define the array with size (0:NX+1,0:NY+1) instead of size (NX,NY) as in the second-order case
  double precision, dimension(0:NX+1,0:NY+1) :: image_data_2D

  integer, dimension(nrec) :: ix_rec,iy_rec

  integer :: ix,iy,irec

  character(len=100) :: file_name,system_command

  integer :: R, G, B

  double precision :: normalized_value,max_amplitude

! open image file and create system command to convert image to more convenient format
! use the "convert" command from ImageMagick http://www.imagemagick.org

  if (type_number == 1) then
    if (field_number == 1) then
      write(file_name,"('./OUTPUT/Vx_image',i6.6,'.pnm')") it
      write(system_command,"('convert Vx_image',i6.6,'.pnm Vx_image',i6.6,'.gif ; rm Vx_image',i6.6,'.pnm')") it,it,it
    else if (field_number == 2) then
      write(file_name,"('./OUTPUT/Vy_image',i6.6,'.pnm')") it
      write(system_command,"('convert Vy_image',i6.6,'.pnm Vy_image',i6.6,'.gif ; rm Vy_image',i6.6,'.pnm')") it,it,it
    else if (field_number == 3) then
      write(file_name,"('./OUTPUT/pressure_image',i6.6,'.pnm')") it
      write(system_command,"('convert pressure_image',i6.6,'.pnm pressure_image',i6.6,'.gif ; rm pressure_image',i6.6,'.pnm')") &
                               it,it,it
   else if (field_number == 4) then
      write(file_name,"('./OUTPUT/density_image',i6.6,'.pnm')") it
      write(system_command,"('convert density_image',i6.6,'.pnm density_image',i6.6,'.gif ; rm density_image',i6.6,'.pnm')") &
                               it,it,it
                               
   else if (field_number == 6) then
    write(file_name,"('./OUTPUT/pressure_image_r',i6.6,'_',i6.6,'.pnm')") rank, it
    write(system_command,"('convert pressure_image',i6.6,'.pnm pressure_image',i6.6,'.gif ; rm pressure_image',i6.6,'.pnm')") &
                               it,it,it                            
   endif 
                              
 else if (type_number == 2) then                                             
   if (field_number == 1) then
      write(file_name,"('./OUTPUT/Vx_obs_image',i6.6,'.pnm')") it
      write(system_command,"('convert Vx_obs_image',i6.6,'.pnm Vx_obs_image',i6.6,'.gif ; rm Vx_obs_image',i6.6,'.pnm')") it,it,it
    else if (field_number == 2) then
      write(file_name,"('./OUTPUT/Vy_obs_image',i6.6,'.pnm')") it
      write(system_command,"('convert Vy_obs_image',i6.6,'.pnm Vy_obs_image',i6.6,'.gif ; rm Vy_obs_image',i6.6,'.pnm')") it,it,it
    else if (field_number == 3) then
      write(file_name,"('./OUTPUT/pressure_obs_image',i6.6,'.pnm')") it
      write(system_command,&
      "('convert image',i6.6,'_pressure_obs.pnm image',i6.6,'_pressure_obs.gif ; rm pressure_obs_image',i6.6,'.pnm')") it,it,it
    endif  
                          
  else if (type_number == 3) then              
    if (field_number == 4) then
      write(file_name,"('./OUTPUT/Krho_image',i6.6,'.pnm')") it
      write(system_command,"('convert Krho_image',i6.6,'.pnm Krho_image',i6.6,'.gif ; rm Krho_image',i6.6,'.pnm')") it,it,it 
                                                           
    else if (field_number == 3) then
      write(file_name,"('./OUTPUT/Kp_image',i6.6,'.pnm')") it
      write(system_command,"('convert Kp_image',i6.6,'.pnm Kp_image',i6.6,'.gif ; rm Kp_image',i6.6,'.pnm')") it,it,it
    else if (field_number == 1) then
      write(file_name,"('./OUTPUT/Kvx_image',i6.6,'.pnm')") it
      write(system_command,"('convert Kvx_image',i6.6,'.pnm Kvx_image',i6.6,'.gif ; rm Kvx_image',i6.6,'.pnm')") it,it,it 
    else if (field_number == 2) then
      write(file_name,"('./OUTPUT/Kvy_image',i6.6,'.pnm')") it
      write(system_command,"('convert Kvy_image',i6.6,'.pnm Kvy_image',i6.6,'.gif ; rm Kvy_image',i6.6,'.pnm')") it,it,it
    endif
    
  else if (type_number == 4) then    
                           
    if (field_number == 1) then
      write(file_name,"('./OUTPUT/Vx_adj_image',i6.6,'.pnm')") it
      write(system_command,"('convert Vx_adj_image',i6.6,'.pnm Vx_adj_image',i6.6,'.gif ; rm Vx_adj_image',i6.6,'.pnm')") it,it,it
    else if (field_number == 2) then
      write(file_name,"('./OUTPUT/Vy_adj_image',i6.6,'.pnm')") it
      write(system_command,"('convert Vy_adj_image',i6.6,'.pnm Vy_adj_image',i6.6,'.gif ; rm Vy_adj_image',i6.6,'.pnm')") it,it,it
    else if (field_number == 3) then
      write(file_name,"('./OUTPUT/pressure_adj_image',i6.6,'.pnm')") it
      write(system_command,&
         "('convert pressure_adj_image',i6.6,'.pnm pressure_adj_image',i6.6,'.gif ; rm pressure_adj_image',i6.6,'.pnm')")it,it,it 
    else if (field_number == 4) then
      write(file_name,"('./OUTPUT/density_adj_image',i6.6,'.pnm')") it
      write(system_command,&
         "('convert density_adj_image',i6.6,'.pnm density_adj_image',i6.6,'.gif ; rm density_adj_image',i6.6,'.pnm')")it,it,it 
    endif
  endif

  open(unit=27, file=file_name, status='unknown')
  write(27,"('P3')") ! write image in PNM P3 format

  write(27,*) NX,NY ! write image size
  write(27,*) '255' ! maximum value of each pixel color

! compute maximum amplitude
  max_amplitude = maxval(abs(image_data_2D))

! image starts in upper-left corner in PNM format
  do iy=NY,1,-1
    do ix=1,NX

! define data as vector component normalized to [-1:1] and rounded to nearest integer
! keeping in mind that amplitude can be negative
    normalized_value = image_data_2D(ix,iy) / max_amplitude

! suppress values that are outside [-1:+1] to avoid small edge effects
    if (normalized_value < -1.d0) normalized_value = -1.d0
    if (normalized_value > 1.d0) normalized_value = 1.d0

! draw an orange cross to represent the source
    if ((ix >= ISOURCE - width_cross .and. ix <= ISOURCE + width_cross .and. &
        iy >= JSOURCE - thickness_cross .and. iy <= JSOURCE + thickness_cross) .or. &
       (ix >= ISOURCE - thickness_cross .and. ix <= ISOURCE + thickness_cross .and. &
        iy >= JSOURCE - width_cross .and. iy <= JSOURCE + width_cross)) then
      R = 255
      G = 157
      B = 0

! display two-pixel-thick black frame around the image
  else if (ix <= 2 .or. ix >= NX-1 .or. iy <= 2 .or. iy >= NY-1) then
      R = 0
      G = 0
      B = 0

! display edges of the PML layers
  else if ((USE_PML_XMIN .and. ix == NPOINTS_PML) .or. &
          (USE_PML_XMAX .and. ix == NX - NPOINTS_PML) .or. &
          (USE_PML_YMIN .and. iy == NPOINTS_PML) .or. &
          (USE_PML_YMAX .and. iy == NY - NPOINTS_PML)) then
      R = 255
      G = 150
      B = 0

! suppress all the values that are below the threshold
    else if (abs(image_data_2D(ix,iy)) <= max_amplitude * cutvect) then

! use a black or white background for points that are below the threshold
      if (WHITE_BACKGROUND) then
        R = 255
        G = 255
        B = 255
      else
        R = 0
        G = 0
        B = 0
      endif

! represent regular image points using red if value is positive, blue if negative
    else if (normalized_value >= 0.d0) then
      R = nint(255.d0*normalized_value**POWER_DISPLAY)
      G = 0
      B = 0
    else
      R = 0
      G = 0
      B = nint(255.d0*abs(normalized_value)**POWER_DISPLAY)
    endif

! draw a green square to represent the receivers
  do irec = 1,nrec
    if ((ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
        iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square) .or. &
       (ix >= ix_rec(irec) - size_square .and. ix <= ix_rec(irec) + size_square .and. &
        iy >= iy_rec(irec) - size_square .and. iy <= iy_rec(irec) + size_square)) then
! use dark green color
      R = 30
      G = 180
      B = 60
    endif
  enddo

! write color pixel
    write(27,"(i3,' ',i3,' ',i3)") R,G,B

    enddo
  enddo

! close file
  close(27)

! call the system to convert image to Gif (can be commented out if "call system" is missing in your compiler)
! call system(system_command)

  end subroutine create_color_image
  
  
  
  !!!!!!!!!!!!!!!!!!!!!
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! GATHER MPI DATA TO GENERATE IMAGE
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gather_and_generate_image(vx_loc,vy_loc,pressure_loc,rhop_loc,it,type_number)
  
  use MPI
  use parameters !,only: &
        !pressure, rhop, vx, vy,&
        !NX, NX_LOCAL, NY, rank, &
        !DELTAT
  
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
   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! GATHER AMPLITUDE DATA ON UX,UY,P BEFORE CREATING PICTURE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
   
    ! GATHER UX
    call MPI_GATHER(data_vx_image_loc(:), image_split_size, MPI_DOUBLE_PRECISION, &
                 data_vx_image_global(:), image_split_size,MPI_DOUBLE_PRECISION,&
                 0, MPI_COMM_WORLD, code)
                 
    ! GATHER UY
    call MPI_GATHER(data_vy_image_loc(:), image_split_size, MPI_DOUBLE_PRECISION, &
                 data_vy_image_global(:), image_split_size,MPI_DOUBLE_PRECISION,&
                 0, MPI_COMM_WORLD, code)
    
    ! GATHER P
    call MPI_GATHER(data_p_image_loc(:), image_split_size, MPI_DOUBLE_PRECISION, &
                 data_p_image_global(:), image_split_size,MPI_DOUBLE_PRECISION,&
                 0, MPI_COMM_WORLD, code)
    
    ! GATHER RHO
    call MPI_GATHER(data_rho_image_loc(:), image_split_size, MPI_DOUBLE_PRECISION, &
                 data_rho_image_global(:), image_split_size,MPI_DOUBLE_PRECISION,&
                 0, MPI_COMM_WORLD, code)
    
    
     call create_color_image(pressure(0:NX_LOCAL+1,0:NY_LOCAL+1),NX_LOCAL,NY_LOCAL,&
      it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
      NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,6,type_number)
              
              
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
  
  

