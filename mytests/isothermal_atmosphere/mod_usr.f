!> This is a template for a new user problem
module mod_usr

! Include a physics module
use mod_rhd
use mod_fld

implicit none

double precision :: M_sun = 1.989d33
double precision :: L_sun = 3.827d33
double precision :: R_sun = 6.96d10

double precision :: M_star
double precision :: L_star
double precision :: T_star
double precision :: R_star

double precision :: Gamma, c_sound, Flux, g_eff, g_grav, H_eff, kappa, c_light

double precision :: Flux0, c_sound0, T_star0, kappa0
double precision :: c_light0, g0, geff0, heff0, Gamma_edd
double precision :: L_star0, R_star0, M_star0
double precision :: tau_bound,  P_bound, rho_bound

contains

!> This routine should set user methods, and activate the physics module
subroutine usr_init()
  use mod_global_parameters
  use mod_usr_methods
  use mod_constants

  call set_coordinate_system("Cartesian_2D")

  M_star = 150*M_sun
  L_star = (M_star/M_sun)**3.d0*L_sun !*unit_time/unit_pressure*unit_length**3.d0
  R_star = 30*R_sun
  T_star = (L_star/(4d0*dpi*R_star**2*5.67051d-5))**0.25d0
  tau_bound = 100.d0

  print*, M_star, L_star, R_star, T_star, tau_bound

  call initglobaldata_usr

  ! !Fix dimensionless stuff here
  ! unit_length        = R_star
  ! unit_numberdensity = 8.955d-8/((1.d0+4.d0*He_abundance)*mp_cgs)
  ! unit_temperature   = T_star

  ! Initialize units
  usr_set_parameters => initglobaldata_usr

  ! A routine for initial conditions is always required
  usr_init_one_grid => initial_conditions

  ! Special Boundary conditions
  usr_special_bc => special_bound
  usr_radiation_bc => radiation_bound

  ! Keep the internal energy constant with internal bound
  usr_internal_bc => constant_e

  ! Graviatational field
  usr_gravity => set_gravitation_field

  ! Output routines
  usr_aux_output    => specialvar_output
  usr_add_aux_names => specialvarnames_output

  ! Active the physics module
  call rhd_activate()

  print*, 'unit_time', unit_time
  print*, 'unit_temperature', unit_temperature
  print*, 'unit_length', unit_length
  print*, 'unit_density', unit_density
  print*, 'unit_numberdensity', unit_numberdensity
  print*, 'unit_velocity', unit_velocity
  print*, 'unit_pressure', unit_pressure
  print*, '================================================================'

end subroutine usr_init

!==========================================================================================

subroutine initglobaldata_usr
use mod_global_parameters

c_sound =  dsqrt((1.38d-16*T_star/(0.6*mp_cgs)))
g_grav = 6.67e-8*M_star/R_star**two

Flux =  L_star/(4*dpi*R_star**2)
kappa = 0.34d0
c_light = const_c
Gamma_edd = (kappa*Flux)/(c_light*g_grav)
g_eff = g_grav*(one - Gamma_edd)
H_eff = c_sound**2/g_eff

p_bound = g_eff*tau_bound/kappa
rho_bound = p_bound/c_sound**two

print*, "######################################################################"
print*, "######################################################################"
print*, c_sound, g_grav, Flux
print*, kappa, c_light, Gamma_edd
print*, g_eff, H_eff, p_bound
print*, rho_bound
print*, "######################################################################"
print*, "######################################################################"

unit_length        = H_eff
unit_numberdensity = rho_bound/((1.d0+4.d0*He_abundance)*mp_cgs)
unit_velocity      = c_sound

L_star0 = L_star*unit_time/(unit_pressure*unit_length**3.d0)
R_star0 = R_star/unit_length
M_star0 = M_star/(unit_density*unit_length**3.d0)

Flux0 = L_star0/(4*dpi*R_star0**2)
T_star0 = T_star/unit_temperature
c_sound0 = dsqrt((1.38d-16*T_star/(0.6*mp_cgs)))/unit_velocity

kappa0 = fld_kappa0
c_light0 = const_c/unit_velocity
g0 = 6.67e-8*M_star/R_star**2*(unit_density*unit_length/unit_pressure)
geff0 = g0*(one - (kappa0*Flux0)/(c_light0*g0))
heff0 = c_sound0**2/geff0
Gamma = (kappa0*Flux0)/(c_light0*g0)

P_bound = one
rho_bound = one

end subroutine initglobaldata_usr

!==========================================================================================

!> A routine for specifying initial conditions
subroutine initial_conditions(ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
   ixmax1,ixmax2, w, x)
  use mod_global_parameters
  use mod_constants
  use mod_variables
  use mod_physics, only: phys_get_pthermal

  integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
     ixmin2,ixmax1,ixmax2
  double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2, ndim)
  double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, nw)
  double precision :: density(ixGmin1:ixGmax1,ixGmin2:ixGmax2),&
      pressure(ixGmin1:ixGmax1,ixGmin2:ixGmax2), pert(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2), amplitude
  double precision :: rad_Flux(ixmin1:ixmax1,ixmin2:ixmax2, 1:ndim)
  double precision :: opacity(ixmin1:ixmax1,ixmin2:ixmax2),&
      Gamma_dep(ixmin1:ixmax1,ixmin2:ixmax2)
  integer :: i

  amplitude = 1.d-5 !3.d-2

  pressure(:,ixGmin2) = p_bound
  density(:,ixGmin2) = rho_bound

  do i=ixGmin2,ixGmax2
    pressure(:,i) = p_bound*dexp(-x(:,i,2)/heff0)
    density(:,i) = pressure(:,i)/c_sound0**two !rho_bound*dexp(-x(:,i,2)/heff0) !
  enddo

  ! Set initial values for w
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = density(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2)
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(:)) = zero
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, e_) = pressure(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2)/(rhd_gamma - one)
  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,r_e) = &
     3.d0*Gamma/(one-Gamma)*pressure(ixGmin1:ixGmax1,ixGmin2:ixGmax2) !> CHANGEd

  !> perturb rho
  call RANDOM_NUMBER(pert)

  w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = density(ixGmin1:ixGmax1,&
     ixGmin2:ixGmax2)*(one + amplitude*pert(ixGmin1:ixGmax1,ixGmin2:ixGmax2))

  call fld_get_opacity(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
     ixmax1,ixmax2)

  print*, "R_star", R_star0, L_star0
  print*, "R_star", R_star, L_star
  print*, "Flux", Flux0

  print*, "g0", g0
  print*, "geff0", geff0
  print*, "c_sound0", c_sound0
  print*, "Gamma", Gamma
  print*, "heff0", heff0
  print*, "Tstar0", T_star0
  print*, "Tstar", T_star

  ! print*, "density", w(5,3:10,rho_) *unit_density
  ! print*, "energy", w(5,3:10,e_) *unit_pressure
  ! print*, "rad_energy", w(5,3:10,r_e) *unit_pressure

  print*, rho_bound*unit_density, p_bound*unit_pressure
  print*, "factor", 3.d0*Gamma/(one-Gamma)

end subroutine initial_conditions

!==========================================================================================

! Extra routines can be placed here
! ...

!==========================================================================================

subroutine special_bound_old(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
   ixBmin2,ixBmax1,ixBmax2,iB,w,x)

  use mod_global_parameters
  use mod_variables
  use mod_physics, only: phys_get_pthermal

  integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixBmin1,ixBmin2,&
     ixBmax1,ixBmax2, iB
  double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
     1:ndim)
  double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
  double precision :: velocity(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndir),&
      pressure(ixGmin1:ixGmax1,ixGmin2:ixGmax2)
  double precision :: fld_R(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
  double precision :: Gamma_dep(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
  double precision :: fld_lambda(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
  double precision :: a(1:nw), b(1:nw), c(1:nw)
  integer :: i,j

  select case (iB)

  case(1)
    do i = ixGmin2, ixGmax2
      w(ixBmin1:ixBmax1,i,:) = w(ixGmax1-3:ixGmax1-2,i,:)
    enddo

  case(2)
    do i = ixGmin2, ixGmax2
      w(ixBmin1:ixBmax1,i,:) = w(ixGmin1+2:ixGmin1+3,i,:)
    enddo

  case(3)
    do i = ixBmin2,ixBmax2
      w(:,i, rho_) = p_bound*dexp(-x(:,i,2)/heff0)/c_sound0**2
      w(:,i, mom(1)) = w(:,i, mom(1))
      velocity(:,i,2) = two*w(:,i+1,mom(2))/w(:,i+1,rho_) - w(:,i+2,&
         mom(2))/w(:,i+2,rho_)
      w(:,i, mom(2)) = velocity(:,i,2)*w(:,i, rho_)
      w(:,i, e_) = p_bound*dexp(-x(:,i,2)/heff0)/(rhd_gamma-one)
      w(:,i, r_e) = 3.d0*Gamma/(one-Gamma)*p_bound*dexp(-x(:,i,2)/heff0)
    enddo

    !> Fixing the R_E boundary for correct gamma due to opacity fluctuation
    !-------------------------------------------------------------------------
    call fld_get_fluxlimiter(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixGmin1+2,&
       ixGmin2+2,ixGmax1-2,ixGmax2-2, fld_lambda, fld_R)

    Gamma_dep(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2) = one/(one + &
       (3.d0*p_bound*dexp(-x(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2,&
       2)/heff0)/w(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2, r_e)))

    do i = ixBmin2,ixBmax2
      w(:,i, r_e) = 3.d0*Gamma_dep(:,ixGmin2+2)/(one-Gamma_dep(:,&
         ixGmin2+2))*p_bound*dexp(-x(:,i,2)/heff0)
    enddo

    !> Corners?
    w(ixGmin1,ixBmin2:ixBmax2,r_e) = w(ixGmax1-3,ixBmin2:ixBmax2,r_e)
    w(ixGmin1+1,ixBmin2:ixBmax2,r_e) = w(ixGmax1-2,ixBmin2:ixBmax2,r_e)
    w(ixGmax1,ixBmin2:ixBmax2,r_e) = w(ixGmin1+3,ixBmin2:ixBmax2,r_e)
    w(ixGmax1-1,ixBmin2:ixBmax2,r_e) = w(ixGmin1+2,ixBmin2:ixBmax2,r_e)

    !-------------------------------------------------------------------------

  case(4)
    !> Linear interpolation
    do i = ixGmin1,ixGmax1
      w(i, ixBmin2, r_e) = -(w(i, ixBmin2-2, r_e) - w(i, ixBmin2-1, r_e))/(x(i,&
         ixBmin2-2,2) - x(i,ixBmin2-1,2)) * (x(i,ixBmin2-1,2) - x(i,ixBmin2,&
         2)) + w(i, ixBmin2-1, r_e)
      w(i, ixBmax2, r_e) = -(w(i, ixBmax2-2, r_e) - w(i, ixBmax2-1, r_e))/(x(i,&
         ixBmax2-2,2) - x(i,ixBmax2-1,2)) * (x(i,ixBmax2-1,2) - x(i,ixBmax2,&
         2)) + w(i, ixBmax2-1, r_e)

      do j = ixBmin2, ixBmax2
        w(i,j,r_e) = min(w(i,j,r_e),w(i,j-1,r_e))
        w(i,j,r_e) = max(w(i,j,r_e), zero)
      enddo
    enddo

    !> Corners?
    w(ixGmin1,ixBmin2:ixBmax2,r_e) = w(ixGmax1-3,ixBmin2:ixBmax2,r_e)
    w(ixGmin1+1,ixBmin2:ixBmax2,r_e) = w(ixGmax1-2,ixBmin2:ixBmax2,r_e)
    w(ixGmax1,ixBmin2:ixBmax2,r_e) = w(ixGmin1+3,ixBmin2:ixBmax2,r_e)
    w(ixGmax1-1,ixBmin2:ixBmax2,r_e) = w(ixGmin1+2,ixBmin2:ixBmax2,r_e)

  case default
    call mpistop("BC not specified")
  end select
end subroutine special_bound_old

!==========================================================================================

subroutine special_bound(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,ixBmin2,&
   ixBmax1,ixBmax2,iB,w,x)
  use mod_global_parameters
  integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2, iB
  double precision, intent(in)    :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
     1:ndim)
  double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
  double precision                :: w_rad(ixGmin1:ixGmax1,ixGmin2:ixGmax2)

  integer j

  select case (iB)
  case(3)
    do j = ixBmin2,ixBmax2
      w(ixGmin1:ixGmax1,j,rho_) = w(ixGmin1:ixGmax1,ixBmax2+1,rho_)
      w(ixGmin1:ixGmax1,j,mom(1)) = w(ixGmin1:ixGmax1,ixBmax2+1,mom(1))
      w(ixGmin1:ixGmax1,j,mom(2)) = w(ixGmin1:ixGmax1,ixBmax2+1,mom(2))
      w(ixGmin1:ixGmax1,j,e_) = w(ixGmin1:ixGmax1,ixBmax2+1,e_)
    enddo

    call radiation_bound(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,iB,w,w_rad,x)
    do j = ixBmin2,ixBmax2
      w(ixGmin1:ixGmax1,j,r_e) = w_rad(ixGmin1:ixGmax1,j)
    enddo

  case(4)
    do j = ixBmin2,ixBmax2
      w(ixGmin1:ixGmax1,j,rho_) = w(ixGmin1:ixGmax1,ixBmin2-1,rho_)
      w(ixGmin1:ixGmax1,j,mom(1)) = w(ixGmin1:ixGmax1,ixBmin2-1,mom(1))
      w(ixGmin1:ixGmax1,j,mom(2)) = w(ixGmin1:ixGmax1,ixBmin2-1,mom(2))
      w(ixGmin1:ixGmax1,j,e_) = w(ixGmin1:ixGmax1,ixBmin2-1,e_)
    enddo

    call radiation_bound(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,iB,w,w_rad,x)

    do j = ixBmin2,ixBmax2
      w(ixGmin1:ixGmax1,j,r_e) = w_rad(ixGmin1:ixGmax1,j)
    enddo

  case default
    call mpistop('boundary not known')
  end select

end subroutine special_bound

!==========================================================================================

subroutine radiation_bound(qt,ixImin1,ixImin2,ixImax1,ixImax2,iB,w,w_rad,x)
  use mod_global_parameters
  integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iB
  double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim)
  double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
  double precision, intent(out)   :: w_rad(ixImin1:ixImax1,ixImin2:ixImax2)

  integer i

  select case (iB)
  case(3)
    do i=ixImin1,ixImax1
      w_rad(i,2) = 2.d0*w(i,3,r_e)-w(i,4,r_e)
      w_rad(i,1) = 2.d0*w(i,2,r_e)-w(i,3,r_e)
    enddo

  case(4)
    do i=ixImin1,ixImax1
      w_rad(i,ixImax2-1) = 2.d0*w(i,ixImax2-2,r_e)-w(i,ixImax2-3,r_e)
      w_rad(i,ixImax2) = 2.d0*w(i,ixImax2-1,r_e)-w(i,ixImax2-2,r_e)
    enddo

  case default
    call mpistop('boundary not known')
  end select

end subroutine radiation_bound

!==========================================================================================

  !> internal boundary, user defined

  !> This subroutine can be used to artificially overwrite ALL conservative
  !> variables in a user-selected region of the mesh, and thereby act as
  !> an internal boundary region. It is called just before external (ghost cell)
  !> boundary regions will be set by the BC selection. Here, you could e.g.
  !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
  !> which can be used to identify the internal boundary region location.
  !> Its effect should always be local as it acts on the mesh.

  subroutine constant_e(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)

    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision :: pressure(ixImin1:ixImax1,ixImin2:ixImax2)

    pressure(ixImin1:ixImax1,ixImin2:ixImax2) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_)*c_sound0**2
    w(ixImin1:ixImax1,ixImin2:ixImax2, e_) = pressure(ixImin1:ixImax1,&
       ixImin2:ixImax2)/(rhd_gamma - one) + half*(w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(1))**two + w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(2))**two)/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)

  end subroutine constant_e

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
     2) = -6.67e-8*M_star/R_star**2*(unit_time**2/unit_length)
end subroutine set_gravitation_field


!==========================================================================================

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

  integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2
  double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim)
  double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)
  double precision                   :: rad_flux(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1:ndim), fld_lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
      fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
  double precision                   :: g_rad(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim), big_gamma(ixImin1:ixImax1,ixImin2:ixImax2), D(ixImin1:ixImax1,&
     ixImin2:ixImax2,1:ndim)
  integer                            :: idim

  call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, rad_flux)
  call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, fld_lambda, fld_R)
  call fld_get_diffcoef(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, D)

  do idim = 1,ndim
    g_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,i_op)*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idim)/c_light0
  enddo
  big_gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = g_rad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,2)/(6.67e-8*M_star/R_star**2*(unit_time**2/unit_length))

  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=rad_flux(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1)*(unit_pressure*unit_velocity)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=rad_flux(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,2)*(unit_pressure*unit_velocity)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3)=fld_lambda(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+4)=fld_R(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+5)=g_rad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1)*unit_length/(unit_time**2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+6)=g_rad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,2)*unit_length/(unit_time**2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+7)=big_gamma(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+8)=D(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+9)=D(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)

end subroutine specialvar_output

subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
  use mod_global_parameters
  character(len=*) :: varnames

  varnames = 'F1 F2 lamnda fld_R ar1 ar2 Gamma D1 D2'

end subroutine specialvarnames_output

!==========================================================================================

end module mod_usr

!==========================================================================================
