!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: rho0
  double precision :: eg0
  double precision :: Er0

  double precision :: eg_equib
  double precision :: Er_equib

  double precision :: Tgas = 0.d0
  double precision :: Trad = 0.d0

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Keep the radiative energy constant with internal bound
    usr_internal_bc => constant_r_e

    ! Active the physics module
    call rhd_activate()

    fld_kappa0 = fld_kappa0/unit_opacity

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

  double precision :: eg_eq, Er_eq
  double precision :: c0, c1, a0
  double precision :: total_energy

  call params_read(par_files)

  total_energy = Er0 + eg0
  a0 = const_rad_a*(const_mp*fld_mu/const_kB)**4

  c0 = total_energy*rho0**4/(a0*(rhd_gamma-1.d0)**4)
  c1 = rho0**4/(a0*(rhd_gamma-1.d0)**4)

  call loc_Halley_method(eg0, eg_equib, c0, c1)
  Er_equib = total_energy - eg_equib

  ! unit_density=rho0
  ! unit_temperature=Tgas
  ! unit_length=1.d0/(fld_kappa0*rho0)
  !
  ! unit_numberdensity=unit_density/((1.d0+4.d0*He_abundance)*const_mp)
  ! unit_pressure=(2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB*unit_temperature
  ! unit_velocity=dsqrt(unit_pressure/unit_density)
  ! unit_time=unit_length/unit_velocity

  fld_kappa0 = fld_kappa0*unit_opacity


  print*, 'Did halley work? :'
  print*, eg_equib**4 + c1*eg_equib - c0
  print*, 'Initial conditions'
  print*,  Er0, eg0
  print*, 'Equilibrium value'
  print*, Er_equib, eg_equib
  print*, 'Is energy conserved?:'
  print*, Er0 + eg0, Er_equib + eg_equib

  print*, 'timescale:'
  print*, 1/(const_c*fld_kappa0*rho0)

  ! rho0 = rho0/unit_density
  ! eg0 = eg0/unit_pressure
  ! Er0 = Er0/unit_pressure
  ! eg_equib = eg_equib/unit_pressure
  ! Er_equib = Er_equib/unit_pressure
  ! Tgas = Tgas/unit_temperature
  ! Trad = Trad/unit_temperature

  print*, fld_kappa0, fld_mu, rhd_gamma

  print*, 'Dimless:'
  print*, 'Er', Er0, 'eg', eg0
  print*, 'c k r E', const_c*fld_kappa0*rho0*Er0
  print*, '4 p k r B', fld_kappa0*rho0*const_c*const_rad_a*Tgas**4

end subroutine initglobaldata_usr

!==========================================================================================

!> Read parameters from a file
subroutine params_read(files)
  use mod_global_parameters, only: unitpar
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /my_list/ rho0, eg0, Er0, Trad, Tgas

  do n = 1, size(files)
     open(unitpar, file=trim(files(n)), status="old")
     rewind(unitpar)
     read(unitpar, my_list, end=113)
113    close(unitpar)
  end do

  if (Tgas .ne. 0.d0) then
    eg0 = rho0/(rhd_gamma - 1.d0) * const_kB*Tgas/(const_mp*fld_mu)
  else
    Tgas = const_mp*fld_mu/const_kB*(rhd_gamma-1)*eg0/rho0
  endif

  if (Trad .ne. 0.d0) then
    Er0 = const_rad_a*Trad**4.d0
  else
    Trad = (Er0/const_rad_a)**0.25d0
  endif

end subroutine params_read

!==========================================================================================

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    use mod_constants
    use mod_rhd_phys, only: rhd_get_pthermal

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)

    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2), fld_R(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    ! Set initial values for w
    w(ixImin1:ixImax1,ixImin2:ixImax2, rho_) = rho0
    w(ixImin1:ixImax1,ixImin2:ixImax2, mom(:)) = 0.d0
    w(ixImin1:ixImax1,ixImin2:ixImax2, e_) = eg0
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = Er0
    ! w(ixI^S,r_e) = Er_equib

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, lambda, fld_R)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_diff_mg) = &
       (const_c/unit_velocity)*lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/(kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

  end subroutine initial_conditions

!==========================================================================================

  ! Extra routines can be placed here
  ! ...

!==========================================================================================

!> internal boundary, user defined
  !
  !> This subroutine can be used to artificially overwrite ALL conservative
  !> variables in a user-selected region of the mesh, and thereby act as
  !> an internal boundary region. It is called just before external (ghost cell)
  !> boundary regions will be set by the BC selection. Here, you could e.g.
  !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
  !> which can be used to identify the internal boundary region location.
  !> Its effect should always be local as it acts on the mesh.

  subroutine constant_r_e(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = rho0
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = 0.d0
    ! w(ixI^S,r_e) = Er_equib

    Tgas = const_mp*fld_mu/const_kB*(rhd_gamma-1)*w(3,3,e_)/rho0
    Trad = (w(3,3,r_e)/const_rad_a)**0.25d0

    if (it .eq. 0) open(1,file='Halley1_1.d2')
    if (mod(it,100) .eq. 0.d0) write(1,*) global_time, Er_equib, eg_equib, w(3,&
       3,r_e), w(3,3,e_), Trad, Tgas, w(3,3,r_e) + w(3,3,e_)
    if (global_time .ge. time_max - dt) close(1)

  end subroutine constant_r_e

!###############################################################################
!###############################################################################
!###############################################################################
!### These are used for initial conditions
!###############################################################################
!###############################################################################
!###############################################################################

  !> Find the root of the 4th degree polynomial using the Halley method
  subroutine loc_Halley_method(eg_0, eg_eq, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: eg_0
    double precision, intent(out)   :: eg_eq

    double precision :: xval, yval, der, dder, deltax

    integer :: ii

    yval = bigdouble
    xval = eg_0
    der = one
    dder = one
    deltax = one

    ii = 0
    !> Compare error with dx = dx/dy dy
    do while (abs(deltax) .gt. 1.d-8)
      yval = loc_Polynomial_Bisection(xval, c0, c1)
      der = loc_dPolynomial_Bisection(xval, c0, c1)
      dder = loc_ddPolynomial_Bisection(xval, c0, c1)
      deltax = -two*yval*der/(two*der**2 - yval*dder)
      xval = xval + deltax
      ii = ii + 1
    enddo

    eg_eq = xval
  end subroutine loc_Halley_method

  !> Evaluate polynomial at argument e_gas
  function loc_Polynomial_Bisection(e_gas, c0, c1) result(val)
    use mod_global_parameters

    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: val

    val = e_gas**4.d0 + c1*e_gas - c0
  end function loc_Polynomial_Bisection

  !> Evaluate first derivative of polynomial at argument e_gas
  function loc_dPolynomial_Bisection(e_gas, c0, c1) result(der)
    use mod_global_parameters

    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: der

    der = 4.d0*e_gas**3.d0 + c1
  end function loc_dPolynomial_Bisection

  !> Evaluate second derivative of polynomial at argument e_gas
  function loc_ddPolynomial_Bisection(e_gas, c0, c1) result(dder)
    use mod_global_parameters

    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: dder

    dder = 4.d0*3.d0*e_gas**2.d0
  end function loc_ddPolynomial_Bisection

end module mod_usr

!==========================================================================================
