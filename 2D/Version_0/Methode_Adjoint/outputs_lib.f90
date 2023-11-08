!----
!----  save the seismograms in ASCII text format
!----

  subroutine write_seismograms(sisvx,sisvy,sispressure,sisrhop,nt,nrec,DELTAT,t0)

  implicit none

  integer nt,nrec
  double precision DELTAT,t0

  double precision sisvx(nt,nrec)
  double precision sisvy(nt,nrec)
  double precision sispressure(nt,nrec)
  double precision sisrhop(nt,nrec)

  integer irec,it

  character(len=100) file_name

 
! pressure
  do irec=1,nrec
    write(file_name,"('./OUTPUT/pressure_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
! in the scheme of eq (13) of Robertsson, Blanch and Symes, Geophysics, vol. 59(9), pp 1444-1456 (1994)
! pressure is defined at time t + DELTAT/2, i.e. staggered in time with respect to velocity.
! Here we must thus take this shift of DELTAT/2 into account to save the seismograms at the right time
      write(11,*) sngl(dble(it-1)*DELTAT - t0 + DELTAT/2.d0),' ',sngl(sispressure(it,irec))
    enddo
    close(11)
  enddo

! density
  do irec=1,nrec
    write(file_name,"('./OUTPUT/density_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT - t0),' ',sngl(sisrhop(it,irec))
    enddo
    close(11)
  enddo

! X component of velocity
  do irec=1,nrec
    write(file_name,"('./OUTPUT/Vx_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT - t0),' ',sngl(sisvx(it,irec))
    enddo
    close(11)
  enddo

! Y component of velocity
  do irec=1,nrec
    write(file_name,"('./OUTPUT/Vy_file_',i3.3,'.dat')") irec
    open(unit=11,file=file_name,status='unknown')
    do it=1,nt
      write(11,*) sngl(dble(it-1)*DELTAT - t0),' ',sngl(sisvy(it,irec))
    enddo
    close(11)
  enddo

  end subroutine write_seismograms

!----
!----  routine to create a color image of a given vector component
!----  the image is created in PNM format and then converted to GIF
!----

  subroutine create_color_image(image_data_2D,NX,NY,it,ISOURCE,JSOURCE,ix_rec,iy_rec,nrec, &
              NPOINTS_PML,USE_PML_XMIN,USE_PML_XMAX,USE_PML_YMIN,USE_PML_YMAX,field_number)

  implicit none

! non linear display to enhance small amplitudes for graphics
  double precision, parameter :: POWER_DISPLAY = 0.30d0

! amplitude threshold above which we draw the color point
  double precision, parameter :: cutvect = 0.01d0

! use black or white background for points that are below the threshold
  logical, parameter :: WHITE_BACKGROUND = .true.

! size of cross and square in pixels drawn to represent the source and the receivers
  integer, parameter :: width_cross = 5, thickness_cross = 1, size_square = 3

  integer NX,NY,it,field_number,ISOURCE,JSOURCE,NPOINTS_PML,nrec
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
  if (field_number == 1) then
    write(file_name,"('./OUTPUT/image',i6.6,'_Vx.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vx.pnm image',i6.6,'_Vx.gif ; rm image',i6.6,'_Vx.pnm')") it,it,it
  else if (field_number == 2) then
    write(file_name,"('./OUTPUT/image',i6.6,'_Vy.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vy.pnm image',i6.6,'_Vy.gif ; rm image',i6.6,'_Vy.pnm')") it,it,it
  else if (field_number == 3) then
    write(file_name,"('./OUTPUT/image',i6.6,'_pressure.pnm')") it
    write(system_command,"('convert image',i6.6,'_pressure.pnm image',i6.6,'_pressure.gif ; rm image',i6.6,'_pressure.pnm')") &
                               it,it,it
 else if (field_number == 11) then
    write(file_name,"('./OUTPUT/image',i6.6,'_Vx_obs.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vx_obs.pnm image',i6.6,'_Vx_obs.gif ; rm image',i6.6,'_Vx_obs.pnm')") it,it,it
  else if (field_number == 22) then
    write(file_name,"('./OUTPUT/image',i6.6,'_Vy_obs.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vy_obs.pnm image',i6.6,'_Vy_obs.gif ; rm image',i6.6,'_Vy_obs.pnm')") it,it,it
  else if (field_number == 33) then
    write(file_name,"('./OUTPUT/image',i6.6,'_pressure_obs.pnm')") it
    write(system_command,&
           "('convert image',i6.6,'_pressure_obs.pnm image',i6.6,'_pressure_obs.gif ; rm image',i6.6,'_pressure_obs.pnm')") &
                               it,it,it
 else if (field_number == 4) then
    write(file_name,"('./OUTPUT/image',i6.6,'_Krho.pnm')") it
    write(system_command,"('convert image',i6.6,'_Krho.pnm image',i6.6,'_Krho.gif ; rm image',i6.6,'_Krho.pnm')") &
                               it,it,it                             
  else if (field_number == 5) then
    write(file_name,"('./OUTPUT/image',i6.6,'_Kp.pnm')") it
    write(system_command,"('convert image',i6.6,'_Kp.pnm image',i6.6,'_Kp.gif ; rm image',i6.6,'_Kp.pnm')") &
                               it,it,it
  else if (field_number == 6) then
    write(file_name,"('./OUTPUT/image',i6.6,'_Kvx.pnm')") it
    write(system_command,"('convert image',i6.6,'_Kvx.pnm image',i6.6,'_Kvx.gif ; rm image',i6.6,'_Kvx.pnm')") &
                               it,it,it 
  else if (field_number == 7) then
    write(file_name,"('./OUTPUT/image',i6.6,'_Kvy.pnm')") it
    write(system_command,"('convert image',i6.6,'_Kvy.pnm image',i6.6,'_Kvy.gif ; rm image',i6.6,'_Kvy.pnm')") &
                               it,it,it
    else if (field_number == 8) then
    write(file_name,"('./OUTPUT/image',i6.6,'_Vx_adj.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vx_adj.pnm image',i6.6,'_Vx_adj.gif ; rm image',i6.6,'_Vx_adj.pnm')") it,it,it
  else if (field_number == 9) then
    write(file_name,"('./OUTPUT/image',i6.6,'_Vy_adj.pnm')") it
    write(system_command,"('convert image',i6.6,'_Vy_adj.pnm image',i6.6,'_Vy_adj.gif ; rm image',i6.6,'_Vy_adj.pnm')") it,it,it
  else if (field_number == 10) then
    write(file_name,"('./OUTPUT/image',i6.6,'_pressure_adj.pnm')") it
    write(system_command,&
           "('convert image',i6.6,'_pressure_adj.pnm image',i6.6,'_pressure_adj.gif ; rm image',i6.6,'_pressure_adj.pnm')") &
                               it,it,it 
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

