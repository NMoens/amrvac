!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd

  implicit none

  double precision :: rho0
  double precision :: eg0
  double precision :: tau_wave
  double precision :: ampl

  double precision :: T0, a0, p0, Er0

  double precision :: wvl, frequency, wavenumber, tau_c, tau_a
  double precision :: Bo, energy_ratio, L_damp, r_Bo, ca

  double precision :: L_0, L_thin, L_thick, L_3

  double precision :: A_rho, A_v, A_p, A_e, A_Er

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

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    call params_read(par_files)


    p0 = eg0*(rhd_gamma - one)
    ca = dsqrt(rhd_gamma*p0/rho0)
    a0 = dsqrt(p0/rho0)


    T0 = const_mp*fld_mu/const_kB*(p0/rho0)
    ! Er0 = const_rad_a*T0**4

    wvl = tau_wave/(rho0*fld_kappa0)
    frequency = 2.d0*dpi*a0/wvl
    wavenumber = 2.d0*dpi/wvl

    Bo = 4*rhd_gamma*ca*eg0/(const_c*Er0)
    r_Bo = Er0/(4*rhd_gamma*eg0)

    !-------------------
    tau_c = const_c*fld_kappa0*rho0/frequency
    tau_a = a0*fld_kappa0*rho0/frequency

    L_thin = 1.d0/(2*dpi*tau_c)
    L_thick = 16.d0/(3*dpi*(rhd_gamma-1.d0))*a0**2/(const_c**2*Bo)*tau_a
    L_3    = 1.d0/(dsqrt(3.d0)*fld_kappa0*rho0)

    print*, L_thin/wvl
    print*, L_thick/wvl
    print*, L_3/wvl
    print*, 1.d0/(2*dpi*tau_a)/wvl
    !-------------------

    ! stop

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
      print*, 'Bo', Bo
      print*, 'r_Bo', r_Bo
      print*, 'L_damp', L_damp
    endif

    A_rho = ampl
    A_v = frequency/(wavenumber*rho0)*A_rho
    A_p = frequency**2/wavenumber**2*A_rho
    A_e = rhd_gamma/(rhd_gamma-one)*p0/rho0*A_rho
    A_Er = Er0/rho0*A_rho

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wave_list/ rho0, eg0, Er0, tau_wave, ampl

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wave_list, end=113)
113    close(unitpar)
    end do

  end subroutine params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: press(ixI^S), temp(ixI^S)
    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)

    ! Set initial values for w
    w(ixI^S, rho_) = rho0! + A_rho*dsin(wavenumber*x(ixI^S,1))
    w(ixI^S, mom(1)) = 0.d0 !w(ixI^S, rho_)*A_v*dsin(wavenumber*x(ixI^S,1))
    w(ixI^S, mom(2)) = zero
    w(ixI^S, e_) = eg0 !+ A_e*dsin(wavenumber*x(ixI^S,1))

    press(ixI^S) = p0 !+ A_p*dsin(wavenumber*x(ixI^S,1))
    temp(ixI^S) = (const_mp*fld_mu/const_kb)*press(ixI^S)/w(ixI^S, rho_)&
    *(unit_pressure/unit_density)/unit_temperature

    w(ixI^S, r_e) = const_rad_a*(temp(ixI^S)*unit_temperature)**4.d0/unit_pressure

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

    double precision :: press(ixI^S), temp(ixI^S)

    where (x(ixI^S,1) .lt. one)
      w(ixI^S, rho_) = rho0 + A_rho*dsin(wavenumber*x(ixI^S,1)-frequency*global_time)
      w(ixI^S, mom(1)) = w(ixI^S, rho_)*A_v*dsin(wavenumber*x(ixI^S,1)-frequency*global_time)
      w(ixI^S, mom(2)) = zero
      w(ixI^S, e_) = eg0 + A_e*dsin(wavenumber*x(ixI^S,1)-frequency*global_time)

      press(ixI^S) = p0 + A_p*dsin(wavenumber*x(ixI^S,1)-frequency*global_time)
      temp(ixI^S) = (const_mp*fld_mu/const_kb)*press(ixI^S)/w(ixI^S, rho_)&
      *(unit_pressure/unit_density)/unit_temperature

      w(ixI^S, r_e) = const_rad_a*(temp(ixI^S)*unit_temperature)**4.d0/unit_pressure
    endwhere

  end subroutine Initialize_Wave


  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: sol1(ixO^S), sol2(ixO^S), sol3(ixO^S)

    w(ixO^S,nw+1) = 1.d0/A_rho*(w(ixO^S,rho_)/rho0 - 1.d0)
    w(ixO^S,nw+2) = 1.d0/A_e*((w(ixO^S,e_)-half*w(ixO^S,mom(1))**2/w(ixO^S,rho_))/eg0 - 1.d0)
    w(ixO^S,nw+3) = 1.d0/A_er*(w(ixO^S,r_e)/Er0 - 1.d0)
    w(ixO^S,nw+4) = 1.d0/A_v*w(ixO^S,mom(1))/w(ixO^S,rho_)

    sol1(ixO^S) = dsin(wavenumber*x(ixO^S,1)-frequency*global_time)
    sol2(ixO^S) = dsin(wavenumber*x(ixO^S,1)-frequency*global_time)
    sol3(ixO^S) = dsin(wavenumber*x(ixO^S,1)-frequency*global_time)
    where (x(ixO^S,1) .gt. one)
      sol1(ixO^S) = sol1(ixO^S)*dexp(-(x(ixO^S,1)-one)/L_thin)
      sol2(ixO^S) = sol2(ixO^S)*dexp(-(x(ixO^S,1)-one)/L_thick)
      sol3(ixO^S) = sol3(ixO^S)*dexp(-(x(ixO^S,1)-one)/L_3)
    endwhere
    w(ixO^S,nw+5) = sol1(ixO^S)
    w(ixO^S,nw+6) = sol2(ixO^S)
    w(ixO^S,nw+7) = sol3(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'delta_rho delta_eg delta_Er v1 s1 s2 s3'
  end subroutine specialvarnames_output
end module mod_usr
