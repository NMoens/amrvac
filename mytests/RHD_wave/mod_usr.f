!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd

  implicit none

  double precision :: rho0 = 3.216d-9
  double precision :: eg0 = 26.020d3
  double precision :: Er0 = 17.340d3

  double precision :: T0, a0, p0

  double precision :: wavelength, frequency, tau_wave
  double precision :: Boltzmann_numer, energy_ratio

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
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters

    p0 = eg0*(rhd_gamma - one)
    a0 = dsqrt(p0/rho0)


    tau_wave = 1.d3
    wavelength = tau_wave/(rho0*0.4d0)
    frequency = 2.d0*dpi*a0/wavelength

    ! Choose independent normalization units if using dimensionless variables.
    unit_length = wavelength ! cm
    unit_velocity   = a0 ! K
    unit_numberdensity = rho0/((1.d0+4.d0*He_abundance)*mp_cgs) ! cm-3,cm-3

    unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+&
       3.d0*He_abundance)*unit_numberdensity*const_kb)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    rho0 = rho0/unit_density
    a0 = a0/unit_velocity
    p0 = p0/unit_pressure
    eg0 = eg0/unit_pressure
    Er0 = Er0/unit_pressure

    wavelength = wavelength/unit_length
    frequency = frequency*unit_time

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
      print*, 'wavelength', wavelength
      print*, 'opt tickness 1 wvl', tau_wave

    endif

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    ! Set initial values for w
    w(ixImin1:ixImax1,ixImin2:ixImax2, rho_) = rho0 !*(one + 0.1*dcos(2*dpi*x(ixImin1:ixImax1,ixImin2:ixImax2,1)/wavelength))
    w(ixImin1:ixImax1,ixImin2:ixImax2, mom(:)) = zero
    w(ixImin1:ixImax1,ixImin2:ixImax2, e_) = eg0 !*(one + 0.1*dcos(2*dpi*x(ixImin1:ixImax1,ixImin2:ixImax2,1)/wavelength))
    w(ixImin1:ixImax1,ixImin2:ixImax2, r_e) = Er0 !*(one + 0.1*dcos(2*dpi*x(ixImin1:ixImax1,ixImin2:ixImax2,1)/wavelength))

    call get_rad_extravars(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2)

  end subroutine initial_conditions


  subroutine Initialize_Wave(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    double precision :: ampl

    ampl = 1.d-5*p0

    where (x(ixImin1:ixImax1,ixImin2:ixImax2,1) .lt. one)
      w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = rho0 + &
         ampl*(2*dpi/(wavelength*frequency))**2 &
         *dsin(frequency*global_time)*dsin(2*dpi*x(ixImin1:ixImax1,&
         ixImin2:ixImax2,1)/wavelength)
      w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)) = &
         ampl*2*dpi/(wavelength*frequency) &
         *dcos(frequency*global_time)*dcos(2*dpi*x(ixImin1:ixImax1,&
         ixImin2:ixImax2,1)/wavelength)
      w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = eg0 + &
         ampl*(2*dpi/(wavelength*frequency))**2 &
         *((eg0+p0)/rho0)*dsin(frequency*global_time)*dsin(2*dpi*x(&
         ixImin1:ixImax1,ixImin2:ixImax2,1)/wavelength)
      w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = Er0 + &
         ampl*(2*dpi/(wavelength*frequency))**2 &
         *(Er0/rho0)*dsin(frequency*global_time)*dsin(2*dpi*x(ixImin1:ixImax1,&
         ixImin2:ixImax2,1)/wavelength)
    endwhere

  end subroutine Initialize_Wave


end module mod_usr
