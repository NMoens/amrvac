
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

  ! Pseudo planar correction
  usr_source => PseudoPlanar

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
  unit_pressure=(2.d0+3.d0*He_abundance)*unit_numberdensity*kB_cgs*unit_temperature
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
subroutine initial_conditions(ixG^L, ix^L, w, x)
  use mod_global_parameters
  use mod_constants
  use mod_variables

  integer, intent(in)             :: ixG^L, ix^L
  double precision, intent(in)    :: x(ixG^S, ndim)
  double precision, intent(inout) :: w(ixG^S, nw)
  integer :: i

  do i = ixGmin2,ixGmax2
    w(:,i,rho_) = rho_rk(i)
    w(:,i,mom(:)) = zero
    w(:,i,e_) = p_rk(i)/(rhd_gamma - 1.d0)
    w(:,i,r_e) = Er_rk(i)
  enddo

  call get_rad_extravars(w, x, ixG^L, ix^L)

end subroutine initial_conditions

!==========================================================================================

subroutine boundary_conditions(qt,ixG^L,ixB^L,iB,w,x)
  use mod_global_parameters
  integer, intent(in)             :: ixG^L, ixB^L, iB
  double precision, intent(in)    :: qt, x(ixG^S,1:ndim)
  double precision, intent(inout) :: w(ixG^S,1:nw)

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
      w(ixGmin1:ixGmax1,i,r_e) = w(ixGmin1:ixGmax1,i-1,rho_)/w(ixGmin1:ixGmax1,i-2,rho_) &
      *(w(ixGmin1:ixGmax1,i-1,r_e) - w(ixGmin1:ixGmax1,i-2,r_e)) + w(ixGmin1:ixGmax1,i-1,r_e)
      do j = ixGmin1,ixGmax1
        w(j,i,r_e) = min(w(j,i,r_e),w(j,i-1,r_e))
      enddo
    enddo

  case default
    call mpistop('boundary not known')
  end select
end subroutine boundary_conditions

subroutine mg_boundary_conditions(qt,ixI^L,ixO^L,iB,w,x)

  use mod_global_parameters

  integer, intent(in)             :: ixI^L, ixO^L, iB
  double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
  double precision, intent(in)    :: w(ixI^S,1:nw)

  double precision :: grad4

  grad4 = sum((w(ixOmin1:ixOmax1,ixOmax2,r_e) - w(ixOmin1:ixOmax1,ixOmax2+1, r_e)) &
  / (x(ixOmin1:ixOmax1,ixOmax2,2) - x(ixOmin1:ixOmax1,ixOmax2+1,2)))/(ixOmax1-ixOmin1)

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

      if (sum(w(ixOmin1:ixOmax1,ixOmax2,r_e)) .le. &
         sum(w(ixOmin1:ixOmax1,ixOmax2+1, r_e))) then
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
subroutine set_gravitation_field(ixI^L,ixO^L,wCT,x,gravity_field)
  use mod_global_parameters
  integer, intent(in)             :: ixI^L, ixO^L
  double precision, intent(in)    :: x(ixI^S,1:ndim)
  double precision, intent(in)    :: wCT(ixI^S,1:nw)
  double precision, intent(out)   :: gravity_field(ixI^S,ndim)

  double precision :: radius(ixI^S)

  radius(ixI^S) = rstar+x(ixI^S,2)*unit_length

  gravity_field(ixI^S,1) = zero
  gravity_field(ixI^S,2) = -const_G*mstar/(radius(ixI^S))**2*(unit_time**2/unit_length)
end subroutine set_gravitation_field

!==========================================================================================

!> If defined, this routine is called before writing output, and it can
!> set/modify the variables in the w array.
subroutine time_average_values(ixI^L,ixO^L,qt,w,x)
  use mod_global_parameters
  use mod_physics, only: phys_get_pthermal

  integer, intent(in)             :: ixI^L,ixO^L
  double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
  double precision, intent(inout) :: w(ixI^S,1:nw)

  double precision :: pth(ixI^S)

  call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)

  if (global_time .eq. 1.d0) w(ixI^S,int_rho) = zero
  if (global_time .eq. 1.d0) w(ixI^S,int_r_e) = zero
  if (global_time .eq. 1.d0) w(ixI^S,int_p) = zero
  if (global_time .eq. 1.d0) w(ixI^S,int_m1) = zero
  if (global_time .eq. 1.d0) w(ixI^S,int_m2) = zero
  if (global_time .eq. 1.d0) w(ixI^S,int_t) = zero

  w(ixI^S,int_rho) = w(ixI^S,int_rho) + w(ixI^S,rho_)*dt
  w(ixI^S,int_r_e) = w(ixI^S,int_r_e) + w(ixI^S,r_e)*dt
  w(ixI^S,int_p) = w(ixI^S,int_p) + pth(ixI^S)*dt
  w(ixI^S,int_m1) = w(ixI^S,int_m1) + w(ixI^S,mom(1))*dt
  w(ixI^S,int_m2) = w(ixI^S,int_m2) + w(ixI^S,mom(2))*dt
  w(ixI^S,int_t) = w(ixI^S,int_t) + dt
end subroutine time_average_values




subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
  use mod_global_parameters
  use mod_physics
  use mod_fld

  integer, intent(in)                :: ixI^L,ixO^L
  double precision, intent(in)       :: x(ixI^S,1:ndim)
  double precision                   :: w(ixI^S,nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)
  double precision                   :: g_rad(ixI^S,1:ndim), big_gamma(ixI^S)
  double precision                   :: Tgas(ixI^S),Trad(ixI^S)
  integer                            :: idim

  do idim = 1,ndim
    g_rad(ixO^S,idim) = w(ixO^S,i_op)*w(ixO^S,i_flux(idim))/fld_speedofligt_0
  enddo
  big_gamma(ixO^S) = g_rad(ixO^S,2)/(const_G*mstar/rstar**2*(unit_time**2/unit_length))

  call rhd_get_tgas(w, x, ixI^L, ixO^L, Tgas)
  call rhd_get_trad(w, x, ixI^L, ixO^L, Trad)

  w(ixO^S,nw+1)=Tgas(ixO^S)*unit_temperature
  w(ixO^S,nw+2)=Trad(ixO^S)*unit_temperature
  w(ixO^S,nw+3)=big_gamma(ixO^S)
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

  kappa_o = kappa_0 !*(1.d0+10.d0**akram*rho*(a2/1.d12)**bkram)
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

!> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
!> iw=iwmin...iwmax.  wCT is at time qCT
subroutine PseudoPlanar(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
  use mod_global_parameters
  use mod_physics, only: phys_get_pthermal

  integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
  double precision, intent(in)    :: qdt, qtC, qt
  double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
  double precision, intent(inout) :: w(ixI^S,1:nw)

  double precision :: pth(ixI^S)
  double precision :: radius(ixI^S)
  integer :: rdir, pdir

  rdir = 2
  pdir = 1

  radius(ixI^S) = rstar/unit_length+x(ixI^S,2)

  !> Correction for spherical fluxes:
  !> drho/dt = -2 rho v_r/r
  w(ixI^S,rho_) = w(ixI^S,rho_) - qdt*two*wCT(ixI^S,mom(rdir))/radius(ixI^S)

  !> dm_r/dt = m_phi**2/r
  !> dm_phi/dt = - m_phi m_r/r
  w(ixI^S,mom(rdir)) = w(ixI^S,mom(rdir)) + qdt*(wCT(ixI^S,mom(pdir)))**two/(radius(ixI^S)*wCT(ixI^S,rho_))
  w(ixI^S,mom(pdir)) = w(ixI^S,mom(pdir)) - qdt*wCT(ixI^S,mom(rdir))*wCT(ixI^S,mom(pdir))/(radius(ixI^S)*wCT(ixI^S,rho_))

  !> de/dt = -2 (e+p)v_r/r
  call phys_get_pthermal(wCT,x,ixI^L,ixO^L,pth)
  w(ixI^S,e_) = w(ixI^S,e_) - qdt*two*(wCT(ixI^S,e_)+pth(ixI^S))*wCT(ixI^S,mom(rdir))/(wCT(ixI^S,rho_)*radius(ixI^S))

  !> dEr/dt = -2 (E v_r + F_r)/r
  if (rhd_radiation_diffusion) then
    call get_rad_extravars(w, x, ixI^L, ixO^L)
    w(ixI^S,r_e) = w(ixI^S,r_e) - qdt*two*(wCT(ixI^S,e_)*wCT(ixI^S,mom(rdir))/wCT(ixI^S,rho_) + w(ixI^S,i_flux(rdir)))/radius(ixI^S)
  else
    w(ixI^S,r_e) = w(ixI^S,r_e) - qdt*two*wCT(ixI^S,e_)*wCT(ixI^S,mom(rdir))/(wCT(ixI^S,rho_)*radius(ixI^S))
  endif

end subroutine PseudoPlanar

!==========================================================================================

end module mod_usr

!==========================================================================================
