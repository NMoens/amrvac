!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd

  implicit none

  double precision :: rho0 = 3.216d-9
  double precision :: eg0 = 26.020d3
  double precision :: Er0! = 17.340d3

  double precision :: T0, a0, p0, ampl

  double precision :: wavelength, frequency, tau_wave, wavenumber
  double precision :: Boltzmann_number, energy_ratio, L_damp, r_Bo

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
    use mod_fld

    p0 = eg0*(rhd_gamma - one)
    a0 = dsqrt(rhd_gamma*p0/rho0)

    T0 = const_mp*fld_mu/const_kB*(p0/rho0)
    Er0 = const_rad_a*T0**4

    tau_wave = 1.d9
    wavelength = tau_wave/(rho0*fld_kappa0)
    frequency = 2.d0*dpi*a0/wavelength
    wavenumber = 2.d0*dpi/wavelength

    Boltzmann_number = 4*rhd_gamma*a0*eg0/(const_c*Er0)
    r_Bo = a0/const_c*Boltzmann_number

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
    T0 = T0/unit_temperature
    Er0 = Er0/unit_pressure

    wavelength = wavelength/unit_length
    frequency = frequency*unit_time
    wavenumber = wavenumber*unit_length
    ampl = 1.d-2*p0

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
      print*, 'wave number', wavenumber
      print*, 'amplitude', ampl
      print*, 'Bo', Boltzmann_number
      print*, 'r_Bo', r_Bo

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
    w(ixImin1:ixImax1,ixImin2:ixImax2, rho_) = rho0 !*(one + 0.1*dsin(2*dpi*x(ixImin1:ixImax1,ixImin2:ixImax2,1)/wavelength))
    w(ixImin1:ixImax1,ixImin2:ixImax2, mom(:)) = zero
    w(ixImin1:ixImax1,ixImin2:ixImax2, e_) = eg0 !*(one + 0.1*dsin(2*dpi*x(ixImin1:ixImax1,ixImin2:ixImax2,1)/wavelength))
    w(ixImin1:ixImax1,ixImin2:ixImax2, r_e) = Er0 !*(one + 0.1*dsin(2*dpi*x(ixImin1:ixImax1,ixImin2:ixImax2,1)/wavelength))

    call get_rad_extravars(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2)

  end subroutine initial_conditions


  subroutine Initialize_Wave(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    use mod_fld
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    double precision :: press(ixImin1:ixImax1,ixImin2:ixImax2),&
        temp(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision :: a2

    a2 = p0/rho0

    where (x(ixImin1:ixImax1,ixImin2:ixImax2,1) .lt. one)
      w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = rho0 + &
         ampl*wavenumber**2/frequency*dsin(wavenumber*x(ixImin1:ixImax1,&
         ixImin2:ixImax2,1)-frequency*global_time)

      w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)) = &
         wavenumber*ampl*dsin(wavenumber*x(ixImin1:ixImax1,ixImin2:ixImax2,&
         1)-frequency*global_time)*w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
      w(ixImin1:ixImax1,ixImin2:ixImax2,mom(2)) = zero

      press(ixImin1:ixImax1,ixImin2:ixImax2) = p0 + &
         rho0*frequency*ampl*dsin(wavenumber*x(ixImin1:ixImax1,ixImin2:ixImax2,&
         1)-frequency*global_time)*w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)

      w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = press(ixImin1:ixImax1,&
         ixImin2:ixImax2)/(rhd_gamma-one) + half*(w(ixImin1:ixImax1,&
         ixImin2:ixImax2,mom(1))**2/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_))

      temp(ixImin1:ixImax1,ixImin2:ixImax2) = T0*(one + &
         rho0*frequency*ampl*dsin(wavenumber*x(ixImin1:ixImax1,ixImin2:ixImax2,&
         1)-frequency*global_time)*w(ixImin1:ixImax1,ixImin2:ixImax2,&
         rho_)/p0 - ampl*wavenumber**2/frequency*dsin(wavenumber*x(&
         ixImin1:ixImax1,ixImin2:ixImax2,1)-frequency*global_time)/rho0)

      w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = &
         const_rad_a*(temp(ixImin1:ixImax1,&
         ixImin2:ixImax2)*unit_temperature)**4.d0/unit_pressure
    endwhere

    ! Print out some information on Q:
    ! print*, 'Q',w(5,5,i_op)


  end subroutine Initialize_Wave


  ! function WaveSolution(ixI^L,ixO^L,x,time,my_w) result(my_w)
  !
  !   use mod_global_parameters
  !   use mod_fld
  !   integer, intent(in)             :: ixI^L,ixO^L
  !   double precision, intent(in)    :: time
  !   double precision, intent(inout) :: my_w(ixI^S,1:nw)
  !   double precision, intent(in)    :: x(ixI^S,1:ndim)
  !
  !   double precision :: press(ixI^S), temp(ixI^S)
  !
  !   double precision :: ampl, a2
  !
  !   my_w(ixI^S,rho_) = rho0 + ampl*wavenumber**2/frequency*dsin(wavenumber*x(ixI^S,1)-frequency*time)
  !
  !   my_w(ixI^S,mom(1)) = wavenumber*ampl*dsin(wavenumber*x(ixI^S,1)-frequency*time)*my_w(ixI^S,rho_)
  !   my_w(ixI^S,mom(2)) = zero
  !
  !   press(ixI^S) = p0 + rho0*frequency*ampl*dsin(wavenumber*x(ixI^S,1)-frequency*time)*my_w(ixI^S,rho_)
  !
  !   my_w(ixI^S,e_) = press(ixI^S)/(rhd_gamma-one) + half*(my_w(ixI^S,mom(1))**2/w(ixI^S,rho_))
  !
  !   temp(ixI^S) = T0*(one + rho0*frequency*ampl*dsin(wavenumber*x(ixI^S,1)-frequency*time)*my_w(ixI^S,rho_)/p0 &
  !   - ampl*wavenumber**2/frequency*dsin(wavenumber*x(ixI^S,1)-frequency*time)/rho0)
  !
  !   my_w(ixI^S,r_e) = const_rad_a*(temp(ixI^S)*unit_temperature)**4.d0/unit_pressure
  !
  !   where (x(ixI^S,1) .lt. one)
  !     my_w(ixI^S,rho_) = my_w(ixI^S,rho_)*dexp(-(x(ixI^S,1)-one)/L_damp)
  !
  !
  !   endwhere
  !
  ! end function WaveSolution


end module mod_usr
