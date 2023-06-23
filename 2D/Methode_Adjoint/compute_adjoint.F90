subroutine compute_adjoint(pa, rhoa, vax, vay, it_step)

  use parameters, only : p0, rho, v0x, v0y,&
                       pressure, rhop, vx, vy, &
                       sispressure, sisrhop, sisvx, sisvy, NREC,IT_DISPLAY, ix_rec, iy_rec, &
                       NX, NY, NINE_OVER_8_DELTAX,ONE_OVER_24_DELTAX, &
                       NINE_OVER_8_DELTAY,ONE_OVER_24_DELTAY,ONE_OVER_SIX_DELTAX,ONE_OVER_SIX_DELTAY, &
                       DELTAT, NSTEP, t,   & 
                       ZERO, gamma_chimie
  implicit none

  integer it_step

  double precision, dimension(0:NX+1,0:NY+1) :: pa, rhoa, vax, vay

  double precision, dimension(0:NX+1,0:NY+1) :: adjoint_source_term
  
  double precision :: &
    value_drhoarho_dx,& 
    value_drho_dx,    &
    value_v0xdvax_dx,  &
    value_v0ydvax_dy,  &
    value_dv0x_dx,    &
    value_dv0y_dx,    &
    value_dvax_dx,     &
    value_dvay_dx,     &
    value_dvax_dy,     &
    value_dvay_dy,     &
    value_drhoarho_dy,&
    value_drho_dy,    &
    value_v0xdvay_dx,  &
    value_v0ydvay_dy,  &
    value_dv0x_dy,    &
    value_dv0y_dy,    &
    value_dp0pa_dx,   &
    value_dp0pa_dy,   &
    value_dp0_dx,     &
    value_dp0_dy
    
  double precision :: &
    value_drhoav0x_dx,&
    value_drhoav0y_dy,&
    value_vax_v0x_dv0x_dx,&
    value_vax_v0y_dv0x_dy,&
    value_vay_v0x_dv0y_dx,&
    value_vay_v0y_dv0y_dy,&
    value_dpav0x_dx,  &
    value_dpav0y_dy

  double precision :: &
  value_dv0x_dy_next, value_dv0x_dy_prec, &
  value_dv0y_dx_next, value_dv0y_dx_prec, &
  value_dvax_dy_next, value_dvax_dy_prec, &
  value_dvay_dx_next, value_dvay_dx_prec

  double precision ::           &
    rho_half_x, rho_half_y,     &
    rhoa_half_x, rhoa_half_y,   &
    pa_half_x, pa_half_y,       &
    vax_half_x_half_y, vay_half_x_half_y, &
    vax_half_x,vay_half_y,                & 
    v0x_half_x_half_y, v0y_half_x_half_y, &
    v0x_half_x,v0y_half_y
    
   double precision, dimension(0:NX+1,0:NY+1) :: aux_rhoa, aux_pa, aux_vax, aux_vay
 
   double precision, dimension(0:NX+1,0:NY+1) :: p0pa, rhoarho, &
               pav0x, pav0y, &
               rhoav0x, rhoav0y, &
               pap
               
          
 
  integer :: i,j
  integer :: Im1, Ip1, Jp1, Jm1


  ! initialisation
 

  call compute_adjoint_source(adjoint_source_term, it_step)

  
  do j=1,NY
   do i=1,NX
    
    !! Boundary condition
    Im1 = i-1
    Ip1 = i+1
    Jm1 = j-1
    Jp1 = j+1
    if (i == NX) then
      Ip1 = 2
    elseif (i == 1) then
      Im1 = NX-1
    endif
    if (j == NY) then
      Jp1 = 2
    elseif (j == 1) then
      Jm1 = NY-1
    endif    
    
    
    aux_vax(i,j) = vax(i,j)
    aux_vay(i,j) = vay(i,j)
    
    rhoarho(i,j) = rho(i,j) * rhoa(i,j)
    p0pa(i,j) = p0(i,j) * pa(i,j)
    
  
    enddo 
   enddo
  
  

  
  do j=1,NY
   do i=1,NX
    
    !! Boundary condition
    Im1 = i-1
    Ip1 = i+1
    Jm1 = j-1
    Jp1 = j+1
    if (i == NX) then
      Ip1 = 2
    elseif (i == 1) then
      Im1 = NX-1
    endif
    if (j == NY) then
      Jp1 = 2
    elseif (j == 1) then
      Jm1 = NY-1
    endif    
    
    ! compute va
    rhoa_half_x = 0.5d0 * (rhoa(i,j) + rhoa(Im1,j))
    rho_half_x  = 0.5d0 * (rho(i,j)  + rho(Im1,j))
    pa_half_x   = 0.5d0 * (pa(i,j)   + pa(Im1,j))
    vay_half_x_half_y = 0.25d0 * (aux_vay(i,j) + aux_vay(Im1,j) + aux_vay(i,Jm1) + aux_vay(Im1,Jm1))  
    v0y_half_x_half_y = 0.25d0 * (v0y(i,j)     + v0y(Im1,j)     + v0y(i,Jm1)     + v0y(Im1,Jm1)) 
    
    rhoa_half_y = 0.5d0 * (rhoa(i,j) + rhoa(i,Jp1))
    rho_half_y  = 0.5d0 * (rho(i,j)  + rho(i,Jp1))
    pa_half_y   = 0.5d0 * (pa(i,j)   + pa(i,Jp1))
    vax_half_x_half_y = 0.25d0 * (aux_vax(i,j) + aux_vax(Ip1,j) + aux_vax(i,Jp1) + aux_vax(Ip1,Jp1))  
    v0x_half_x_half_y = 0.25d0 * (v0x(i,j)     + v0x(Ip1,j)     + v0x(i,Jp1)     + v0x(Ip1,Jp1))  
    
    
    call compute_centered_dU_dx_in_i(rhoarho, Im1, j, value_drhoarho_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dx_in_i(rho, Im1, j, value_drho_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dx_in_i(p0pa, Im1, j, value_dp0pa_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dx_in_i(p0, Im1, j, value_dp0_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    
    call compute_decentered_dU_dx_in_i(aux_vax, i, j, value_dvax_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x(i,j))
    call compute_decentered_dU_dz_in_i(aux_vax, i, j, value_dvax_dy, ONE_OVER_SIX_DELTAY,NX,NY,v0y(i,j))
    
    call compute_decentered_dU_dx_in_i(v0x, i, j, value_dv0x_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0x(i,j))
    call compute_centered_dU_dx_in_i(v0y, Im1, j, value_dv0x_dy_prec, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dx_in_i(v0y, Im1, Jm1, value_dv0x_dy_next, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    value_dv0y_dx = 0.5d0 * ( value_dv0y_dx_prec + value_dv0y_dx_next) 
    
    value_v0xdvax_dx = v0x(i,j) * value_dvax_dx
    value_v0ydvax_dy = v0y_half_x_half_y * value_dvax_dy
  
   
    call compute_centered_dU_dz_in_i(rhoarho, i, j, value_drhoarho_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    call compute_centered_dU_dz_in_i(rho, i, j, value_drho_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    call compute_centered_dU_dz_in_i(p0pa, i, j, value_dp0pa_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    call compute_centered_dU_dz_in_i(pressure, i, j, value_dp0_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    
    call compute_decentered_dU_dz_in_i(aux_vay, i, j, value_dvay_dy, ONE_OVER_SIX_DELTAY,NX,NY,v0y(i,j))
    call compute_decentered_dU_dx_in_i(aux_vay, i, j, value_dvay_dx, ONE_OVER_SIX_DELTAX,NX,NY,v0y(i,j))
    
    call compute_decentered_dU_dz_in_i(v0y, i, j, value_dv0y_dy, ONE_OVER_SIX_DELTAX,NX,NY,v0y(i,j))
    call compute_centered_dU_dz_in_i(v0x, i, j, value_dv0x_dy_prec, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    call compute_centered_dU_dz_in_i(v0x, Ip1, j, value_dv0x_dy_next, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    value_dv0x_dy = 0.5d0 * ( value_dv0x_dy_prec + value_dv0x_dy_next) 
    
    value_v0xdvay_dx = v0x_half_x_half_y * value_dvay_dx
    value_v0ydvay_dy = v0y(i,j) * value_dvay_dy

    ! x component    
    vax(i,j) = vax(i,j) + (value_drhoarho_dx - rhoa_half_x * value_drho_dx) * DELTAT / rho_half_x
    !vax(i,j) = vax(i,j) + rho_half_x * (value_v0xdvax_dx + value_v0ydvax_dy) * DELTAT / rho_half_x
    !vax(i,j) = vax(i,j) - rho_half_x * (aux_vax(i,j) * value_dv0x_dx + vay_half_x_half_y * value_dv0y_dx) * DELTAT / rho_half_x
    vax(i,j) = vax(i,j) + (gamma_chimie * value_dp0pa_dx - pa_half_x * value_dp0_dx) * DELTAT / rho_half_x
    
    ! y component
    vay(i,j) = vay(i,j) + (value_drhoarho_dy - rhoa_half_y * value_drho_dy) * DELTAT / rho_half_y
    !vay(i,j) = vay(i,j) + rho_half_y * (value_v0xdvay_dx + value_v0ydvay_dy) * DELTAT / rho_half_y
    !vay(i,j) = vay(i,j) - rho_half_y * (vax_half_x_half_y * value_dv0x_dy + aux_vay(i,j) * value_dv0y_dy) * DELTAT / rho_half_y
    vay(i,j) = vay(i,j) + (gamma_chimie * value_dp0pa_dy - pa_half_y * value_dp0_dy) * DELTAT / rho_half_y
    
  enddo
  enddo
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do j=1,NY
   do i=1,NX
    
    !! Boundary condition
    Im1 = i-1
    Ip1 = i+1
    Jm1 = j-1
    Jp1 = j+1
    if (i == NX) then
      Ip1 = 2
    elseif (i == 1) then
      Im1 = NX-1
    endif
    if (j == NY) then
      Jp1 = 2
    elseif (j == 1) then
      Jm1 = NY-1
    endif 
   
   aux_rhoa(i,j) = rhoa(i,j)
   aux_pa(i,j) = pa(i,j)
   
   v0x_half_x = 0.5d0 * (v0x(i,j) + v0x(Ip1,j))
   v0y_half_y = 0.5d0 * (v0y(i,j) + v0y(Im1,j))
    
   rhoav0x(i,j) = rhoa(i,j) * v0x_half_x
   rhoav0y(i,j) = rhoa(i,j) * v0y_half_y
   pav0x(i,j) = pa(i,j) * v0x_half_x
   pav0y(i,j) = pa(i,j) * v0y_half_y

   enddo
  enddo
  
  do j=1,NY
   do i=1,NX
    
    !! Boundary condition
    Im1 = i-1
    Ip1 = i+1
    Jm1 = j-1
    Jp1 = j+1
    if (i == NX) then
      Ip1 = 2
    elseif (i == 1) then
      Im1 = NX-1
    endif
    if (j == NY) then
      Jp1 = 2
    elseif (j == 1) then
      Jm1 = NY-1
    endif 
    
    vax_half_x = 0.5d0 * (vax(i,j) + vax(Ip1,j))
    vay_half_y = 0.5d0 * (vay(i,j) + vay(Im1,j))
    v0x_half_x = 0.5d0 * (v0x(i,j) + v0x(Ip1,j))
    v0y_half_y = 0.5d0 * (v0y(i,j) + v0y(Im1,j))
    
    call compute_centered_dU_dx_in_i(vax, i, j, value_dvax_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dz_in_i(vay, i, Jm1, value_dvay_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
    call compute_centered_dU_dx_in_i(v0x, i, j, value_dv0x_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dz_in_i(v0y, i, Jm1, value_dv0y_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
        
    call compute_decentered_dU_dx_in_i(rhoav0x, i, j, value_drhoav0x_dx,ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
    call compute_decentered_dU_dz_in_i(rhoav0y, i, j, value_drhoav0y_dy,ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y)
    call compute_decentered_dU_dx_in_i(pav0x, i, j, value_dpav0x_dx,ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
    call compute_decentered_dU_dz_in_i(pav0y, i, j, value_dpav0y_dy,ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y)
    
    call compute_centered_dU_dx_in_i(v0x, i, j, value_dv0x_dx, NINE_OVER_8_DELTAX, ONE_OVER_24_DELTAX,NX,NY)
    call compute_centered_dU_dz_in_i(v0y, i, Jm1, value_dv0y_dy, NINE_OVER_8_DELTAY, ONE_OVER_24_DELTAY,NX,NY)
        
    call compute_decentered_dU_dx_in_i(v0y, i, Jm1, value_dv0y_dx_prec,ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
    call compute_decentered_dU_dz_in_i(v0x, i, j, value_dv0x_dy_prec,ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y)
    call compute_decentered_dU_dx_in_i(v0y, i, j, value_dv0y_dx_next,ONE_OVER_SIX_DELTAX,NX,NY,v0x_half_x)
    call compute_decentered_dU_dz_in_i(v0x, Ip1, j, value_dv0x_dy_next,ONE_OVER_SIX_DELTAY,NX,NY,v0y_half_y)
    value_dv0y_dx = 0.5d0 * (value_dv0y_dx_prec + value_dv0y_dx_next)  
    value_dv0x_dy = 0.5d0 * (value_dv0x_dy_prec + value_dv0x_dy_next)
    
    value_vax_v0x_dv0x_dx = vax_half_x * v0x_half_x * value_dv0x_dx
    value_vax_v0y_dv0x_dy = vax_half_x * v0y_half_y * value_dv0x_dy
    value_vay_v0x_dv0y_dx = vay_half_y * v0x_half_x * value_dv0y_dx
    value_vay_v0y_dv0y_dy = vay_half_y * v0y_half_y * value_dv0y_dy

    
    ! compute rhoa and pa
    !rhoa(i,j) = rhoa(i,j) -  aux_rhoa(i,j) * (value_dv0x_dx + value_dv0y_dy) * DELTAT
    !rhoa(i,j) = rhoa(i,j) + (value_drhoav0x_dx + value_drhoav0y_dy)          * DELTAT  
    !rhoa(i,j) = rhoa(i,j) - (value_vax_v0x_dv0x_dx + value_vax_v0y_dv0x_dy)  * DELTAT
    !rhoa(i,j) = rhoa(i,j) - (value_vay_v0x_dv0y_dx + value_vay_v0y_dv0y_dy)  * DELTAT
    
    !pa(i,j) = pa(i,j) -  gamma_chimie * aux_pa(i,j) * (value_dv0x_dx + value_dv0y_dy) * DELTAT
    pa(i,j) = pa(i,j) + (value_dvax_dx + value_dvay_dy)                               * DELTAT
    !pa(i,j) = pa(i,j) + (value_dpav0x_dx + value_dpav0y_dy)                           * DELTAT  

    pa(i,j) = pa(i,j) + adjoint_source_term(i,j) * DELTAT

   enddo     
  enddo


endsubroutine compute_adjoint



subroutine compute_adjoint_source(adjoint_source_term, it_step)

 use parameters, only : obspressure, calcpressure, norm_obs, ix_rec, iy_rec, NREC, NX, NY, NSTEP, TINYVAL
 double precision, dimension(0:NX+1, 0:NY+1) :: adjoint_source_term
  
 integer :: it_step, irec
 integer :: i
 
  adjoint_source_term(:,:) = 0.0d0
  
   do irec=1,NREC
       adjoint_source_term(ix_rec(irec), iy_rec(irec)) = &
                (obspressure(NSTEP-it_step+1,irec) - calcpressure(NSTEP-it_step+1,irec)) 
                
      if (abs(adjoint_source_term(ix_rec(irec), iy_rec(irec)))  < TINYVAL) then
        adjoint_source_term(ix_rec(irec), iy_rec(irec)) = 0.0d0
      endif
                        
      adjoint_source_term(ix_rec(irec), iy_rec(irec)) = adjoint_source_term(ix_rec(irec), iy_rec(irec)) / norm_obs 
       !adjoint_source_term(i, iy_rec(irec)) = (obspressure(NSTEP-it_step+1,irec) - calcpressure(NSTEP-it_step+1,irec)) / norm_obs 
       !adjoint_source_term(ix_rec(irec), i) = (obspressure(NSTEP-it_step+1,irec) - calcpressure(NSTEP-it_step+1,irec)) / norm_obs   
   enddo
   
endsubroutine compute_adjoint_source
