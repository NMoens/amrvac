!> This is a template for a new user problem
module mod_usr

! Include a physics module
use mod_rhd
use mod_fld
use mod_global_parameters

implicit none

  integer :: int_rho, int_r_e, int_p, int_m1, int_m2, int_t

  integer, parameter :: FW_cells = 66
  double precision, parameter :: M_sun = 1.99d33
  double precision, parameter :: R_sun = 6.96d10
  double precision, parameter :: year = 365.25*24*60*60


  double precision :: y_FW(1:FW_cells)
  double precision :: rho_FW(1:FW_cells)
  double precision :: v_FW(1:FW_cells)
  double precision :: pg_FW(1:FW_cells)
  double precision :: Tg_FW(1:FW_cells)
  double precision :: a2_FW(1:FW_cells)

  double precision :: er_FW(1:FW_cells)
  double precision :: tr_FW(1:FW_cells)
  double precision :: Heff_FW(1:FW_cells)

  double precision :: rstar
  double precision :: mstar

  double precision :: mg_val4

contains

!> This routine should set user methods, and activate the physics module
subroutine usr_init()
  use mod_global_parameters
  use mod_usr_methods
  use mod_constants

  call set_coordinate_system("Cartesian_2D")

  ! Initialize units
  usr_set_parameters => initglobaldata_usr

  ! A routine for initial conditions is always required
  usr_init_one_grid => initial_conditions

  ! Routine for setting special boundary conditions
  usr_special_bc => boundary_conditions
  usr_special_mg_bc => mg_boundary_conditions


  ! Graviatational field
  usr_gravity => set_gravitation_field

  ! Get time integrated values
  usr_modify_output => time_average_values

  ! Output routines
  usr_aux_output    => specialvar_output
  usr_add_aux_names => specialvarnames_output

  ! Active the physics module
  call rhd_activate()

  int_rho = var_set_extravar('int_rho','int_rho')
  int_r_e = var_set_extravar('int_r_e','int_r_e')
  int_p = var_set_extravar('int_p','int_p')
  int_m1 = var_set_extravar('int_m1','int_m1')
  int_m2 = var_set_extravar('int_m2','int_m2')
  int_t = var_set_extravar('int_t','int_t')

end subroutine usr_init

!==========================================================================================

subroutine initglobaldata_usr
  use mod_global_parameters

  integer :: i,j
  double precision :: dum1, dum2
  double precision :: geff(FW_cells)

  !> Reading in parameters
  OPEN(1,FILE='FastWind_profile.txt')
    READ(1,*)
    READ(1,*)
    READ(1,*)
    READ(1,*)
    do i = 1,FW_cells
        j = FW_cells - i + 1
        READ(1,*) y_FW(j), rho_FW(j), v_FW(j), pg_FW(j), Tg_FW(j), dum1,&
            a2_FW(j)
    enddo
  CLOSE(1)

  mstar = 66.8*M_sun
  rstar = y_FW(1)

  geff = const_G*mstar/rstar**2*(0.5)

  Heff_FW = a2_FW/geff
  er_FW = const_rad_a*tg_FW**4

  !> Define units
  unit_numberdensity = rho_FW(1)/((1.d0+4.d0*He_abundance)*mp_cgs)
  unit_temperature = Tg_FW(1)
  unit_length = Heff_FW(1)

  !> Remaining units
  unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
  unit_pressure=(2.d0+3.d0*He_abundance&
     )*unit_numberdensity*kB_cgs*unit_temperature
  unit_velocity=dsqrt(unit_pressure/unit_density)
  unit_time=unit_length/unit_velocity
  unit_radflux = unit_velocity*unit_pressure
  unit_opacity = one/(unit_density*unit_length)

  !> Make input dimensionless:
  y_FW = y_FW/unit_length
  rho_FW = rho_FW/unit_density
  v_FW = v_FW/unit_velocity
  tg_FW = tg_FW/unit_temperature
  pg_FW = pg_FW/unit_pressure
  er_FW = er_FW/unit_pressure

  y_FW = y_FW - y_FW(1)
  tr_FW  = (er_FW/const_rad_a*unit_pressure/unit_temperature**4)**(1.0/4.0)

  ! pg_FW = kB_cgs/(mp_cgs*fld_mu)*tr_FW*rho_FW*(unit_temperature*unit_density)/unit_pressure

  if (mype .eq. 0) then
    print*, 'unit_length', unit_length
    print*, 'unit_density', unit_density
    print*, 'unit_pressure', unit_pressure
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_radflux', unit_radflux
    print*, 'unit_opacity', unit_opacity
    print*, 'unit_time', unit_time
    print*, 'unit_velocity', unit_velocity

    print*, minval(y_FW), xprobmin2
    print*, maxval(y_FW), xprobmax2

    ! if (minval(y_FW) .gt. xprobmin2) call mpistop("Simulation space not covered")
    ! if (maxval(y_FW) .lt. xprobmax2) call mpistop("Simulation space not covered")
  endif

end subroutine initglobaldata_usr

!==========================================================================================

!> A routine for specifying initial conditions
subroutine initial_conditions(ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
   ixmax1,ixmax2, w, x)
  use mod_global_parameters
  use mod_constants
  use mod_variables

  integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
     ixmin2,ixmax1,ixmax2
  double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2, ndim)
  double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, nw)

  double precision :: x_vac(ixGmin2:ixGmax2)
  double precision :: rho_vac(ixGmin2:ixGmax2)
  double precision :: v_vac(ixGmin2:ixGmax2)
  double precision :: pg_vac(ixGmin2:ixGmax2)
  double precision :: er_vac(ixGmin2:ixGmax2)

  double precision :: pert(ixGmin1:ixGmax1,ixGmin2:ixGmax2), amplitude
  integer :: i
  double precision :: rbs, xc1,xc2

  x_vac(ixGmin2:ixGmax2) = x(nghostcells+1,ixGmin2:ixGmax2,2)
  call Interpolate(ixGmin2, ixGmax2, y_FW, x_vac, rho_FW, rho_vac,.false.)
  call Interpolate(ixGmin2, ixGmax2, y_FW, x_vac, v_FW, v_vac,.false.)
  call Interpolate(ixGmin2, ixGmax2, y_FW, x_vac, pg_FW, pg_vac,.false.)
  call Interpolate(ixGmin2, ixGmax2, y_FW, x_vac, er_FW, er_vac,.false.)


  ! print*, er_vac(1:5)


  do i = ixGmin1, ixGmax1
    w(i, :, rho_) = rho_vac(:)
    w(i, :, mom(1)) = zero
    w(i, :, mom(2)) = rho_vac(:)*v_vac(:)
    w(i, :, r_e) = er_vac(:)
    w(i, :, e_) = pg_vac(:)/(rhd_gamma-1.0) + &
       half*rho_vac(:)*v_vac(:)*v_vac(:)
  enddo


  ! !> perturb rho
  ! amplitude = 0.00d0
  ! call RANDOM_NUMBER(pert)
  !
  ! w(ixGmin1:ixGmax1, i, rho_) = w(ixGmin1:ixGmax1, i, rho_)&
  ! *(one + amplitude*pert(ixGmin1:ixGmax1, i))


  ! rbs = (xprobmin1+xprobmax1)*0.25d0
  ! xc1 = (xprobmin1+xprobmax1)*0.5d0
  ! xc2 = (xprobmin2+xprobmax2)*0.5d0
  ! where((x(ix^S,1)-xc1)**2+(x(ix^S,2)-xc2)**2<rbs**2)
  !   w(ix^S,rho_) = one
  ! endwhere

  call get_rad_extravars(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
     ixmax1,ixmax2)
  ! call set_mg_bounds()


end subroutine initial_conditions

!==========================================================================================

subroutine boundary_conditions(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
   ixBmin2,ixBmax1,ixBmax2,iB,w,x)
  use mod_global_parameters
  integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2, iB
  double precision, intent(in)    :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
     1:ndim)
  double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

  double precision :: x_vac(ixGmin2:ixGmax2)
  double precision :: rho_vac(ixGmin2:ixGmax2)
  double precision :: v_vac(ixGmin2:ixGmax2)
  double precision :: pg_vac(ixGmin2:ixGmax2)
  double precision :: er_vac(ixGmin2:ixGmax2)
  integer :: i,j

  x_vac(ixGmin2:ixGmax2) = x(nghostcells+1,ixGmin2:ixGmax2,2)
  call Interpolate(ixGmin2, ixGmax2, y_FW, x_vac, rho_FW, rho_vac,.false.)
  call Interpolate(ixGmin2, ixGmax2, y_FW, x_vac, v_FW, v_vac,.false.)
  call Interpolate(ixGmin2, ixGmax2, y_FW, x_vac, pg_FW, pg_vac,.false.)
  call Interpolate(ixGmin2, ixGmax2, y_FW, x_vac, er_FW, er_vac,.false.)

  select case (iB)

  case(3)
    do i = ixBmin2,ixBmax2
      w(ixGmin1:ixGmax1,i,rho_) = rho_vac(i)
      w(ixGmin1:ixGmax1,i,mom(1)) = w(ixGmin1:ixGmax1,ixBmax2+1,mom(1)) !rho_vac(i)*v_vac(i)
      w(ixGmin1:ixGmax1,i,mom(2)) = w(ixGmin1:ixGmax1,ixBmax2+1,mom(2))
      w(ixGmin1:ixGmax1,i,e_) = pg_vac(i)/(rhd_gamma-1.0) + &
         half*(w(ixGmin1:ixGmax1,i,mom(1))**2+w(ixGmin1:ixGmax1,i,&
         mom(1))**2)/rho_vac(i)
      w(ixGmin1:ixGmax1,i,r_e) = er_vac(i)
    enddo

  case(4)
    do i = ixBmin2,ixBmax2
      !> Conserve gradE/rho
      w(ixGmin1:ixGmax1,i,r_e) = w(ixGmin1:ixGmax1,i-1,rho_)/w(ixGmin1:ixGmax1,&
         i-2,rho_) *(w(ixGmin1:ixGmax1,i-1,r_e) - w(ixGmin1:ixGmax1,i-2,&
         r_e)) + w(ixGmin1:ixGmax1,i-1,r_e)
      ! do j = ixGmin1,ixGmax1
      !   w(j,i,r_e) = min(w(j,i,r_e),w(j,i-1,r_e))
      ! enddo
    enddo

    ! mg_val4 = sum(w(:,ixBmin2,r_e))/(ixGmax1 - ixGmin1)

  case default
    call mpistop('boundary not known')
  end select
end subroutine boundary_conditions

subroutine mg_boundary_conditions(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,iB,w,x)

  use mod_global_parameters

  integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, iB
  double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim)
  double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

  double precision :: grad4

  grad4 = sum((w(ixOmin1:ixOmax1,ixOmax2,r_e) - w(ixOmin1:ixOmax1,ixOmax2+1,&
      r_e)) / (x(ixOmin1:ixOmax1,ixOmax2,2) - x(ixOmin1:ixOmax1,ixOmax2+1,&
     2)))/(ixOmax1-ixOmin1)

  select case (iB)
    case (1)
       mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
    case (2)
       mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
    case (3)
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = 7.6447315544263788
    case (4)
      ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      ! mg%bc(iB, mg_iphi)%bc_value = 0.77780683570039721
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
      ! mg%bc(iB, mg_iphi)%bc_value = grad4 !min(grad4,zero)
      ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous

    case default
      print *, "Not a standard: ", trim(typeboundary(iw_r_e, iB))
      error stop "You have to set a user-defined boundary method"
  end select
end subroutine mg_boundary_conditions

!==========================================================================================

!> Calculate gravitational acceleration in each dimension
subroutine set_gravitation_field(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,wCT,x,gravity_field)
  use mod_global_parameters
  integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2
  double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
  double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
  double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,&
     ixImin2:ixImax2,ndim)

  gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,1) = zero
  gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
     2) = -const_G*mstar/rstar**2*(unit_time**2/unit_length)
end subroutine set_gravitation_field

!==========================================================================================

!> If defined, this routine is called before writing output, and it can
!> set/modify the variables in the w array.
subroutine time_average_values(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,qt,w,x)
  use mod_global_parameters
  use mod_physics, only: phys_get_pthermal

  integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2
  double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim)
  double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

  double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2)

  call phys_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,pth)

  if (global_time .eq. 1.d0) w(ixImin1:ixImax1,ixImin2:ixImax2,int_rho) = zero
  if (global_time .eq. 1.d0) w(ixImin1:ixImax1,ixImin2:ixImax2,int_r_e) = zero
  if (global_time .eq. 1.d0) w(ixImin1:ixImax1,ixImin2:ixImax2,int_p) = zero
  if (global_time .eq. 1.d0) w(ixImin1:ixImax1,ixImin2:ixImax2,int_m1) = zero
  if (global_time .eq. 1.d0) w(ixImin1:ixImax1,ixImin2:ixImax2,int_m2) = zero
  if (global_time .eq. 1.d0) w(ixImin1:ixImax1,ixImin2:ixImax2,int_t) = zero

  w(ixImin1:ixImax1,ixImin2:ixImax2,int_rho) = w(ixImin1:ixImax1,&
     ixImin2:ixImax2,int_rho) + w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)*dt
  w(ixImin1:ixImax1,ixImin2:ixImax2,int_r_e) = w(ixImin1:ixImax1,&
     ixImin2:ixImax2,int_r_e) + w(ixImin1:ixImax1,ixImin2:ixImax2,r_e)*dt
  w(ixImin1:ixImax1,ixImin2:ixImax2,int_p) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
     int_p) + pth(ixImin1:ixImax1,ixImin2:ixImax2)*dt
  w(ixImin1:ixImax1,ixImin2:ixImax2,int_m1) = w(ixImin1:ixImax1,&
     ixImin2:ixImax2,int_m1) + w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1))*dt
  w(ixImin1:ixImax1,ixImin2:ixImax2,int_m2) = w(ixImin1:ixImax1,&
     ixImin2:ixImax2,int_m2) + w(ixImin1:ixImax1,ixImin2:ixImax2,mom(2))*dt
  w(ixImin1:ixImax1,ixImin2:ixImax2,int_t) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
     int_t) + dt
end subroutine time_average_values




subroutine specialvar_output(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
  use mod_global_parameters
  use mod_physics
  use mod_fld

  integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2
  double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim)
  double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)
  double precision                   :: g_rad(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim), big_gamma(ixImin1:ixImax1,ixImin2:ixImax2)
  double precision                   :: Tgas(ixImin1:ixImax1,ixImin2:ixImax2),&
     Trad(ixImin1:ixImax1,ixImin2:ixImax2)
  integer                            :: idim

  do idim = 1,ndim
    g_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,i_op)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       i_flux(idim))/fld_speedofligt_0
  enddo
  big_gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = g_rad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,2)/(const_G*mstar/rstar**2*(unit_time**2/unit_length))

  call rhd_get_tgas(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, Tgas)
  call rhd_get_trad(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, Trad)

  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=Tgas(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)*unit_temperature
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=Trad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)*unit_temperature
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3)=big_gamma(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
end subroutine specialvar_output

subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
  use mod_global_parameters
  character(len=*) :: varnames

  varnames = 'Tgas Trad Gamma'
end subroutine specialvarnames_output


subroutine Interpolate(ixImin2, ixImax2, x_in, x_out, f_in, f_out, log)
  use mod_global_parameters
  use mod_constants
  use mod_variables

  integer, intent(in)             :: ixImin2,ixImax2
  double precision, intent(in)    :: x_in(FW_cells), f_in(FW_cells)
  double precision, intent(in)    :: x_out(ixImin2:ixImax2)
  logical, intent(in)             :: log
  double precision, intent(out)    :: f_out(ixImin2:ixImax2)

  double precision :: log_f_in_low, log_f_in_up, log_f_out
  double precision :: res(FW_cells), low_val
  integer :: up_i, low_i, i

  double precision :: x12, xx12, f12
  double precision :: x13, xx13, f13
  double precision :: a,b,c

  do i = ixImin2,ixImax2
    !> First look for the two surrounding cells in the in-grid
    res = x_in - x_out(i)

    low_val = maxval(res, MASK = res .le. 0) + x_out(i)
    low_i = minloc(abs(x_in - low_val),1)

    up_i = low_i + 1

    if (log) then
      !> Logarithmic Interpolation
      ! log_f_in_low = dlog(f_in(low_i))
      ! log_f_in_up = dlog(f_in(up_i))
      ! log_f_out = log_f_in_low + (x_out(i) - x_in(low_i))*(log_f_in_up - log_f_in_low)/(x_in(up_i) - x_in(low_i))
      ! f_out(i) = dexp(log_f_out)

      x12 = x_in(low_i) - x_in(up_i)
      x13 = x_in(low_i) - x_in(up_i + 1)

      xx12 = x_in(low_i)**2 - x_in(up_i)**2
      xx13 = x_in(low_i)**2 - x_in(up_i + 1)**2

      f12 = dlog(f_in(low_i)) - dlog(f_in(up_i))
      f13 = dlog(f_in(low_i)) - dlog(f_in(up_i + 1))

      b = (f13/xx13 - f12/xx12)/(x13/xx13 - x12/xx12)
      c = (f13/x13 - f12/x12)/(xx13/x13 - xx12/x12)
      a = f_in(low_i) - (b*x_in(low_i) + c*x_in(low_i)**2)

      f_out(i) = dexp(a + b*x_out(i) + c*x_out(i)**2)


    else
      ! !> Linear Interpolation
      ! f_out(i) = f_in(low_i) + (x_out(i) - x_in(low_i))*(f_in(up_i) - f_in(low_i))/(x_in(up_i) - x_in(low_i))

      x12 = x_in(low_i) - x_in(up_i)
      x13 = x_in(low_i) - x_in(up_i + 1)

      xx12 = x_in(low_i)**2 - x_in(up_i)**2
      xx13 = x_in(low_i)**2 - x_in(up_i + 1)**2

      f12 = f_in(low_i) - f_in(up_i)
      f13 = f_in(low_i) - f_in(up_i + 1)

      b = (f13/xx13 - f12/xx12)/(x13/xx13 - x12/xx12)
      c = (f13/x13 - f12/x12)/(xx13/x13 - xx12/x12)
      a = f_in(low_i) - (b*x_in(low_i) + c*x_in(low_i)**2)

      f_out(i) = a + b*x_out(i) + c*x_out(i)**2
    endif

    ! print*, x_in(low_i), x_in(up_i), x_out(i), f_in(low_i), f_in(up_i), f_out(i)
  enddo
end subroutine Interpolate


!==========================================================================================

end module mod_usr

!==========================================================================================
