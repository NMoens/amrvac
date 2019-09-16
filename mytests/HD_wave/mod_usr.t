!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd

  implicit none

  double precision :: rho0 = 3.216d-9
  double precision :: eg0 = 26.020d3
  double precision :: Er0! = 17.340d3

  double precision :: T0, a0, p0, ampl

  double precision :: wavelength, frequency, tau_wave, wavenumber
  double precision :: Boltzmann_number, energy_ratio, L_damp, r_Bo

  double precision :: A_rho, A_v, A_p, A_e


contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Drive the wave using an internal boundary
    usr_internal_bc => Initialize_Wave

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! ...

    ! Active the physics module
    call hd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters

    p0 = eg0*(hd_gamma - one)
    a0 = dsqrt(p0/rho0)

    tau_wave = 1.d9

    wavelength = tau_wave/(rho0*0.4d0)
    frequency = 2.d0*dpi*a0/wavelength
    wavenumber = 2.d0*dpi/wavelength

    ! Choose independent normalization units if using dimensionless variables.
    unit_length = wavelength ! cm
    unit_velocity   = a0 ! K
    unit_numberdensity = rho0/((1.d0+4.d0*He_abundance)*mp_cgs) ! cm^-3

    unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kb)
    unit_time=unit_length/unit_velocity

    rho0 = rho0/unit_density
    a0 = a0/unit_velocity
    p0 = p0/unit_pressure
    eg0 = eg0/unit_pressure


    wavelength = wavelength/unit_length
    frequency = frequency*unit_time
    wavenumber = wavenumber*unit_length

    ampl = 1.d-5

    if (mype .eq. 0) then
      print*, 'unit_length', unit_length
      print*, 'unit_density', unit_density
      print*, 'unit_pressure', unit_pressure
      print*, 'unit_temperature', unit_temperature
      print*, 'unit_radflux', unit_radflux
      print*, 'unit_time', unit_time
      print*, 'unit_velocity', unit_velocity
      print*, '-----------------------------'
      print*, 'angular frequency', frequency
      print*, 'wavelength', wavelength
      print*, 'opt tickness 1 wvl', tau_wave
      print*, 'wave number', wavenumber
      print*, 'amplitude', ampl
    endif

    A_rho = ampl
    A_v = frequency/(wavenumber*rho0)*A_rho
    A_p = frequency**2/wavenumber**2*A_rho
    A_e = hd_gamma/(hd_gamma-one)*p0/rho0*A_rho

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: press(ixI^S)

    ! Set initial values for w
    w(ixI^S, rho_) = rho0 + A_rho*dsin(wavenumber*x(ixI^S,1))
    w(ixI^S, mom(1)) = w(ixI^S, rho_)*A_v*dsin(wavenumber*x(ixI^S,1))
    w(ixI^S, mom(2)) = zero
    w(ixI^S, e_) = eg0 + A_e*dsin(wavenumber*x(ixI^S,1))

  end subroutine initial_conditions


  subroutine Initialize_Wave(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision :: press(ixI^S), temp(ixI^S)

    !> HD_WAVE

    where (x(ixI^S,1) .lt. one)
      w(ixI^S, rho_) = rho0 + A_rho*dsin(wavenumber*x(ixI^S,1)-frequency*global_time)
      w(ixI^S, mom(1)) = w(ixI^S, rho_)*A_v*dsin(wavenumber*x(ixI^S,1)-frequency*global_time)
      w(ixI^S, mom(2)) = zero
      w(ixI^S, e_) = eg0 + A_e*dsin(wavenumber*x(ixI^S,1)-frequency*global_time)
    endwhere

  end subroutine Initialize_Wave

end module mod_usr
