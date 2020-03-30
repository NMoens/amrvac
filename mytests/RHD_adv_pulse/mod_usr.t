!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: rho1 = 1.2d0
  double precision :: v1 = 1.d6
  double precision :: T1 = 1.d7
  double precision :: T2 = 2.d7
  double precision :: wdth = 24.d0

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    unit_velocity = v1 !r_arr(nghostcells) ! cm
    unit_numberdensity = rho1/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = wdth

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

  end subroutine initglobaldata_usr


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: temp(ixI^S), pth(ixI^S)
    double precision :: kappa(ixO^S), fld_R(ixO^S), lambda(ixO^S)

    temp(ixI^S) = T1 + (T2-T1)*dexp(-(x(ixI^S,1)*unit_length)**2.d0/(2*wdth**2))
    w(ixI^S,rho_) = rho1*T1/temp(ixI^S) &
    + const_rad_a*fld_mu*const_mp/(3.d0*const_kB) &
    * (T1**4.d0/temp(ixI^S) - temp(ixI^S)**3.d0)
    w(ixI^S,mom(1)) = w(ixI^S,rho_)*v1
    w(ixI^S,mom(2)) = 0.d0
    pth(ixI^S) = const_kB*temp(ixI^S)*w(ixI^S,rho_) &
    /(const_mp*fld_mu)
    w(ixI^S,e_) = pth(ixI^S)/(rhd_gamma-1.d0) + half*w(ixI^S,rho_)*v1**2
    w(ixI^S,r_e) = const_rad_a*temp(ixI^S)**4.d0

    w(ixI^S,rho_) = w(ixI^S,rho_)/unit_density
    w(ixI^S,mom(:)) = w(ixI^S,mom(:))/(unit_density*unit_velocity)
    w(ixI^S,e_) = w(ixI^S,e_)/unit_pressure
    w(ixI^S,r_e) = w(ixI^S,r_e)/unit_pressure

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))

  end subroutine initial_conditions

end module mod_usr
