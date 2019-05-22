!> This is a template for a new user problem
module mod_usr

! Include a physics module
use mod_rhd
use mod_fld
use mod_global_parameters

implicit none

  integer :: int_rho, int_r_e, int_p, int_m1, int_m2, int_t

  integer, parameter :: RK_cells = 256+2*4
  double precision, parameter :: M_sun = 1.99d33
  double precision, parameter :: R_sun = 6.96d10
  double precision, parameter :: year = 365.25*24*60*60

  double precision :: p_rk(RK_cells)
  double precision :: Er_rk(RK_cells)
  double precision :: rho_rk(RK_cells)

  double precision :: rstar, mstar

contains

!> This routine should set user methods, and activate the physics module
subroutine usr_init()
  use mod_global_parameters
  use mod_usr_methods
  use mod_constants

  call set_coordinate_system("Cartesian_2D")

  call initglobaldata_usr


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

  integer :: i
  double precision :: dy, p0, Er0
  double precision :: k1,k2,k3,k4
  double precision :: l1,l2,l3,l4

  !> Define units
  unit_numberdensity = 8.8837999999999995e-9/((1.d0+4.d0*He_abundance)*mp_cgs)
  unit_temperature = 106459.89999999999
  unit_length = 5075047950.3534832

  dy = 79297624.22436523
  p0 = 127551.50951160883
  Er0 = 967319.2978373199

  !> Remaining units
  unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
  unit_pressure=(2.d0+3.d0*He_abundance&
     )*unit_numberdensity*kB_cgs*unit_temperature
  unit_velocity=dsqrt(unit_pressure/unit_density)
  unit_time=unit_length/unit_velocity
  unit_radflux = unit_velocity*unit_pressure
  unit_opacity = one/(unit_density*unit_length)


  mstar = 66.8d0*M_sun
  rstar = 1252739648812.112


  if (mype .eq. 0) then
    print*, 'unit_length', unit_length
    print*, 'unit_density', unit_density
    print*, 'unit_pressure', unit_pressure
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_radflux', unit_radflux
    print*, 'unit_opacity', unit_opacity
    print*, 'unit_time', unit_time
    print*, 'unit_velocity', unit_velocity
  endif

  p_rk(1) = p0
  Er_rk(1) = Er0
  rho_rk(1) = get_rho(p0,Er0)

  do i = 1,RK_cells -1

    ! print*, i, p_rk(i), Er_rk(i), rho_rk(i)

    k1 = dy*f_pressure(p_rk(i),Er_rk(i))
    k2 = dy*f_pressure(p_rk(i)+k1/2.0,Er_rk(i))
    k3 = dy*f_pressure(p_rk(i)+k2/2.0,Er_rk(i))
    k4 = dy*f_pressure(p_rk(i)+k3,Er_rk(i))

    p_rk(i+1) = p_rk(i) + 1./6.*(k1+2*k2+2*k3+k4)

    l1 = dy*f_radiation(p_rk(i),Er_rk(i))
    l2 = dy*f_radiation(p_rk(i),Er_rk(i)+l1/2.0)
    l3 = dy*f_radiation(p_rk(i),Er_rk(i)+l2/2.0)
    l4 = dy*f_radiation(p_rk(i),Er_rk(i)+l3)

    Er_rk(i+1) = Er_rk(i) + 1./6.*(l1+2*l2+2*l3+l4)

    rho_rk(i+1) = get_rho(p_rk(i+1),Er_rk(i+1))
  enddo

  p_rk = p_rk/unit_pressure
  Er_rk = Er_rk/unit_pressure
  rho_rk = rho_rk/unit_density

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
  integer :: i

  do i = ixGmin2,ixGmax2
    w(:,i,rho_) = rho_rk(i)
    w(:,i,mom(:)) = zero
    w(:,i,e_) = p_rk(i)/(rhd_gamma - 1.d0)
    w(:,i,r_e) = Er_rk(i)
  enddo

  call get_rad_extravars(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
     ixmax1,ixmax2)

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
  double precision :: temp_vac(ixGmin2:ixGmax2)
  integer :: i,j

  select case (iB)

  case(3)
    do i = ixBmin2,ixBmax2
      w(ixGmin1:ixGmax1,i,rho_) = rho_rk(i)
      w(ixGmin1:ixGmax1,i,mom(1)) = w(ixGmin1:ixGmax1,ixBmax2+1,mom(1)) !rho_vac(i)*v_vac(i)
      w(ixGmin1:ixGmax1,i,mom(2)) = w(ixGmin1:ixGmax1,ixBmax2+1,mom(2))
      w(ixGmin1:ixGmax1,i,e_) = p_rk(i)/(rhd_gamma - 1.d0)
      w(ixGmin1:ixGmax1,i,r_e) = Er_rk(i)
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
      mg%bc(iB, mg_iphi)%bc_value = 7.3583819042386223 !7.6447315544263788
    case (4)
      ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      ! mg%bc(iB, mg_iphi)%bc_value = 0.77780683570039721
      ! mg%bc(iB, mg_iphi)%bc_value = grad4 !min(grad4,zero)
      ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous

      ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann

      if (sum(w(ixOmin1:ixOmax1,ixOmax2,r_e)) .le. sum(w(ixOmin1:ixOmax1,&
         ixOmax2+1, r_e))) then
         mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
         mg%bc(iB, mg_iphi)%bc_value = zero
      else
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
      endif

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

function get_rho(p_i,E_i) result(rho_o)
  use mod_constants
  use mod_fld

  double precision, intent(in) :: p_i, E_i
  double precision :: rho_o

  double precision :: Trad

  Trad = (E_i/const_rad_a)**0.25
  rho_o = p_i/Trad *const_mp*fld_mu/const_kB

end function

function get_kappa(p_i,E_i) result(kappa_o)
  use mod_constants

  double precision, intent(in) :: p_i, E_i
  double precision :: kappa_o

  double precision :: rho, a2
  double precision :: kappa_0, akram, bkram

  rho = get_rho(p_i,E_i)
  a2 = p_i/rho

  kappa_0 = 0.34d0
  akram = 13.1351597305
  bkram = -4.5182188206

  kappa_o = kappa_0*(1.d0+10.d0**akram*rho*(a2/1.d12)**bkram)
end function

function get_Gamma(p_i,E_i) result(Gamma_o)
  use mod_constants

  double precision, intent(in) :: p_i, E_i
  double precision :: Gamma_o

  double precision :: kappa, F_rk, grav

  grav = const_G*mstar/rstar**2
  F_rk = 194606471845852.0

  kappa = get_kappa(p_i,E_i)
  Gamma_o = kappa*F_rk/(grav*const_c)
end function

function f_pressure(p_i,E_i) result(f_o)
  use mod_constants

  double precision, intent(in) :: p_i, E_i
  double precision :: f_o

  double precision :: grav, rho, Gamma

  grav = const_G*mstar/rstar**2
  f_o = -get_rho(p_i,E_i)*grav*(1.d0-get_Gamma(p_i,E_i))
end function

function f_radiation(p_i,E_i) result(f_o)
  use mod_constants

  double precision, intent(in) :: p_i, E_i
  double precision :: f_o

  double precision :: F_rk

  F_rk = 194606471845852.0

  f_o = -3.d0*F_rk/const_c*get_kappa(p_i,E_i)*get_rho(p_i,E_i)
end function

!==========================================================================================

end module mod_usr

!==========================================================================================
