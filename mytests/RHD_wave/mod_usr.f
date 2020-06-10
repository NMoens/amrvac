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

  double precision :: wvl, omega, wavenumber, tau_c, tau_a
  double precision :: Bo, energy_ratio, L_damp, r_Bo, ca

  double precision :: L_0, L_thin, L_thick, L_3

  double precision :: A_rho, A_v, A_p, A_e, A_Er

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as n-D Cartesian
    
    
     call set_coordinate_system("Cartesian_3D")

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
    omega = 2.d0*dpi*a0/wvl
    wavenumber = 2.d0*dpi/wvl

    Bo = 4*rhd_gamma*ca*eg0/(const_c*Er0)
    r_Bo = Er0/(4*rhd_gamma*eg0)

    !-------------------
    tau_c = const_c*fld_kappa0*rho0/omega
    tau_a = a0*fld_kappa0*rho0/omega

    L_thin = 1.d0/(2*dpi*tau_c)
    L_thick = 16.d0/(3*dpi*(rhd_gamma-1.d0))*a0**2/(const_c**2*Bo)*tau_a
    L_3    = 1.d0/(dsqrt(3.d0)*fld_kappa0*rho0)

    print*, L_thin/wvl
    print*, L_thick/wvl
    print*, L_3/wvl
    print*, 1.d0/(2*dpi*tau_a)/wvl
    !-------------------

    ! Choose independent normalization units if using dimensionless variables.
    unit_length = wvl ! cm
    unit_velocity   = a0 ! cm/s
    unit_numberdensity = rho0/((1.d0+4.d0*He_abundance)*mp_cgs) !cm-3,cm-3,cm-3

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

    wvl = wvl/unit_length
    omega = omega*unit_time
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
      print*, 'angular omega', omega
      print*, 'wvl', wvl
      print*, 'opt tickness 1 wvl', tau_wave
      print*, 'wave number', wavenumber
      print*, 'amplitude', ampl
      print*, 'Bo', Bo
      print*, 'r_Bo', r_Bo
      print*, 'L_damp', L_damp
    endif

    A_rho = ampl
    A_v = omega/(wavenumber*rho0)*A_rho
    A_p = omega**2/wavenumber**2*A_rho
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
  subroutine initial_conditions(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: press(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        temp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    ! Set initial values for w
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, rho_) = rho0
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, mom(:)) = 0.d0
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, mom(1)) = 0.d0
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, e_) = eg0

    press(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = p0 !+ A_p*dsin(wavenumber*x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1))
    temp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = (const_mp*fld_mu/const_kb)*press(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)/w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, rho_)*(unit_pressure/unit_density)/unit_temperature

    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
        r_e) = const_rad_a*(temp(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)*unit_temperature)**4.d0/unit_pressure

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, lambda,&
        fld_R)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_diff_mg) = (const_c/unit_velocity)*lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)/(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_))

  end subroutine initial_conditions


  subroutine Initialize_Wave(level,qt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
    use mod_global_parameters
    use mod_fld
    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)

    double precision :: press(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        temp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    where (x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1) .lt. one)
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          rho_) = rho0 + A_rho*dsin(wavenumber*x(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1)-omega*global_time)
       w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, mom(2)) = zero
       w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3, mom(3)) = zero
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          mom(1)) = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          rho_)*A_v*dsin(wavenumber*x(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1)-omega*global_time)
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          e_) = eg0 + A_e*dsin(wavenumber*x(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3,1)-omega*global_time)

      press(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) = p0 + A_p*dsin(wavenumber*x(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3,1)-omega*global_time)
      temp(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3) = (const_mp*fld_mu/const_kb)*press(ixImin1:ixImax1,&
         ixImin2:ixImax2,ixImin3:ixImax3)/w(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3, rho_)*(unit_pressure/unit_density)/unit_temperature

      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
          r_e) = const_rad_a*(temp(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)*unit_temperature)**4.d0/unit_pressure
    endwhere

  end subroutine Initialize_Wave


  subroutine specialvar_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
       ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: sol1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        sol2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        sol3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+1) = 1.d0/A_rho*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)/rho0 - 1.d0)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+2) = 1.d0/A_e*((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       e_)-half*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(1))**2/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_))/eg0 - 1.d0)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+3) = 1.d0/A_er*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       r_e)/Er0 - 1.d0)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+4) = 1.d0/A_v*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(1))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)

    sol1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = dsin(wavenumber*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)-omega*global_time)
    sol2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = dsin(wavenumber*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)-omega*global_time)
    sol3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = dsin(wavenumber*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)-omega*global_time)
    where (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1) .gt. one)
      sol1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sol1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)-one)/L_thin)
      sol2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sol2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)-one)/L_thick)
      sol3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3) = sol3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)*dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,1)-one)/L_3)
    endwhere
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+5) = sol1(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+6) = sol2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       nw+7) = sol3(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'delta_rho delta_eg delta_Er v1 s1 s2 s3'
  end subroutine specialvarnames_output
end module mod_usr
