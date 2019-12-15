!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd

  implicit none

  double precision :: rho0 = 3.216d-9
  double precision :: eg0 = 26.020d3
  double precision :: Er0! = 17.340d3

  double precision :: T0, a0, p0, ampl

  double precision :: wvl, frequency, tau_wave, wavenumber
  double precision :: Boltzmann_number, energy_ratio, L_damp, r_Bo

  double precision :: A_rho, A_v, A_p, A_e, A_Er

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 3D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_3D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Drive the wave using an internal boundary
    usr_internal_bc => Initialize_Wave

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! ...

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    p0 = eg0*(rhd_gamma - one)
    ! a0 = dsqrt(rhd_gamma*p0/rho0)
    a0 = dsqrt(p0/rho0)


    T0 = const_mp*fld_mu/const_kB*(p0/rho0)
    Er0 = const_rad_a*T0**4

    tau_wave = 1.d3

    wvl = tau_wave/(rho0*fld_kappa0)
    frequency = 2.d0*dpi*a0/wvl
    wavenumber = 2.d0*dpi/wvl

    Boltzmann_number = 4*rhd_gamma*a0*eg0/(const_c*Er0)
    r_Bo = a0/(const_c*Boltzmann_number)

    ! Choose independent normalization units if using dimensionless variables.
    unit_length = wvl ! cm
    unit_velocity   = a0 ! cm/s
    unit_numberdensity = rho0/((1.d0+4.d0*He_abundance)*mp_cgs) ! cm^-3

    unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kb)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    rho0 = rho0/unit_density
    a0 = a0/unit_velocity
    p0 = p0/unit_pressure
    eg0 = eg0/unit_pressure
    T0 = T0/unit_temperature
    Er0 = Er0/unit_pressure


    wvl = wvl/unit_length
    frequency = frequency*unit_time
    wavenumber = wavenumber*unit_length

    L_damp = const_c/unit_velocity*fld_kappa0/unit_opacity*rho0/frequency

    ampl = 1.d-5

    if (mype .eq. 0) then
      print*, 'unit_length', unit_length
      print*, 'unit_density', unit_density
      print*, 'unit_pressure', unit_pressure
      print*, 'unit_temperature', unit_temperature
      print*, 'unit_radflux', unit_radflux
      print*, 'unit_opacity', unit_opacity
      print*, 'unit_time', unit_time
      print*, 'unit_velocity', unit_velocity
      print*, '-----------------------------'
      print*, 'angular frequency', frequency
      print*, 'wvl', wvl
      print*, 'opt tickness 1 wvl', tau_wave
      print*, 'wave number', wavenumber
      print*, 'amplitude', ampl
      print*, 'Bo', Boltzmann_number
      print*, 'r_Bo', r_Bo
      print*, 'L_damp', L_damp
    endif

    A_rho = ampl
    A_v = frequency/(wavenumber*rho0)*A_rho
    A_p = frequency**2/wavenumber**2*A_rho
    A_e = rhd_gamma/(rhd_gamma-one)*p0/rho0*A_rho
    A_Er = Er0/rho0*A_rho

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: press(ixO^S), temp(ixO^S)
    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)

    ! Set initial values for w
    w(ixO^S, rho_) = rho0 + A_rho*dsin(wavenumber*x(ixO^S,1))
    w(ixO^S, mom(1)) = w(ixO^S, rho_)*A_v*dsin(wavenumber*x(ixO^S,1))
    w(ixO^S, mom(2)) = zero
    w(ixO^S, mom(3)) = zero
    w(ixO^S, e_) = eg0 + A_e*dsin(wavenumber*x(ixO^S,1))

    press(ixO^S) = p0 + A_p*dsin(wavenumber*x(ixO^S,1))
    temp(ixO^S) = (const_mp*fld_mu/const_kb)*press(ixO^S)/w(ixO^S, rho_)&
    *(unit_pressure/unit_density)/unit_temperature

    w(ixO^S, r_e) = const_rad_a*(temp(ixO^S)*unit_temperature)**4.d0/unit_pressure

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))

  end subroutine initial_conditions


  subroutine Initialize_Wave(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_fld
    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision :: press(ixO^S), temp(ixO^S)

    where (x(ixO^S,1) .lt. one)
      w(ixO^S, rho_) = rho0 + A_rho*dsin(wavenumber*x(ixO^S,1)-frequency*global_time)
      w(ixO^S, mom(1)) = w(ixO^S, rho_)*A_v*dsin(wavenumber*x(ixO^S,1)-frequency*global_time)
      w(ixO^S, mom(2)) = zero
      w(ixO^S, mom(3)) = zero
      w(ixO^S, e_) = eg0 + A_e*dsin(wavenumber*x(ixO^S,1)-frequency*global_time)

      press(ixO^S) = p0 + A_p*dsin(wavenumber*x(ixO^S,1)-frequency*global_time)
      temp(ixO^S) = (const_mp*fld_mu/const_kb)*press(ixO^S)/w(ixO^S, rho_)&
      *(unit_pressure/unit_density)/unit_temperature

      w(ixO^S, r_e) = const_rad_a*(temp(ixO^S)*unit_temperature)**4.d0/unit_pressure
    endwhere

  end subroutine Initialize_Wave

end module mod_usr
