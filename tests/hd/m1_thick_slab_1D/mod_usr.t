!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd
  use mod_m1

  implicit none

  ! Custom global variables can be defined here
  ! User input parameters
  double precision :: mstar, rstar, Tbound, rhobound, beta
  double precision :: asound, Fbound, Ebound

  double precision :: lstar, vesc, gamma_e, logg, logge, heff
  double precision :: fred_bound

  integer :: i_prad, i_pht, i_fcmf, i_Tgas, i_Trad, i_chi, i_ff, i_lum
  integer :: i_SdivF, i_SdivP, i_SFdv, i_SgR

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see src/mod_usr_methods.t
    usr_special_bc => special_bound

    ! Fix negative radiation flux
    usr_set_F0 => stellar_flux

    !> Additional variables
    usr_process_grid => update_extravars

    ! Choose independent normalization units if using dimensionless variables.
    unit_length      = rstar    ! cm
    unit_temperature = Tbound    ! K
    unit_density     = rhobound ! g cm^-3

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

    !> Extra variables
    i_Prad = var_set_extravar('P_rr','P_rr')
    i_pht = var_set_extravar('pht','pht')
    i_fcmf = var_set_extravar('Fcmf','Fcmf')
    i_Tgas = var_set_extravar('Tgas','Tgas')
    i_Trad = var_set_extravar('Trad','Trad')
    i_ff = var_set_extravar('ff','ff')
    i_chi = var_set_extravar('chi','chi')
    i_lum = var_set_extravar('lum','lum')

    i_SdivF = var_set_extravar('Div_Fv','Div_Fv')
    i_SdivP = var_set_extravar('Div_c2P','Div_c2P')
    i_SFdv = var_set_extravar('DDv','DDv')
    i_SgR = var_set_extravar('gR','gR')

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /star_list/ mstar, rstar, Tbound, rhobound, gamma_e

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    enddo

    ! Scale to cgs units
    mstar = mstar * const_MSun
    rstar = rstar * const_RSun

  end subroutine usr_params_read

  subroutine initglobaldata_usr

    use mod_constants
    use mod_hd_phys

    ! Stellar structure
    lstar  = gamma_e*(4 *dpi * const_G * mstar * const_c)/(m1_kappa_e*unit_opacity)
    Fbound = lstar/(4*dpi*rstar**2)
    logg   = log10(const_G * mstar/rstar**2)
    logge  = logg + log10(1.0d0 - gamma_e)
    asound = sqrt(Tbound * kB_cgs/(hd_mu * mp_cgs))
    heff   = asound**2.0d0 / 10.0d0**logge
    Ebound = const_rad_a*Tbound**4

    if (mype == 0) then
      print*, '======================'
      print*, '   Unity quantities   '
      print*, '======================'
      print*, 'unit length        = ', unit_length
      print*, 'unit density       = ', unit_density
      print*, 'unit velocity      = ', unit_velocity
      print*, 'unit numberdensity = ', unit_numberdensity
      print*, 'unit pressure      = ', unit_pressure
      print*, 'unit temperature   = ', unit_temperature
      print*, 'unit time          = ', unit_time
      print*, '==============================================='
      print*, '   Stellar and wind parameters in CGS units    '
      print*, '==============================================='
      print*, 'L/Lsun                 = ', lstar/const_LSun
      print*, 'M/Msun                 = ', mstar/const_MSun
      print*, 'R/Rsun                 = ', rstar/const_RSun
      print*, 'Tbound                 = ', Tbound
      print*, 'Mean molecular weight  = ', hd_mu
      print*, 'log(g)                 = ', logg
      print*, 'eff. log(g)            = ', logge
      print*, 'heff/Rstar             = ', heff/rstar
      print*, 'Eddington gamma        = ', gamma_e
      print*, 'adiabatic gamma        = ', hd_gamma
      print*, 'isothermal asound      = ', asound
    endif

    ! Make some variables needed for computations dimensionless
    mstar    = mstar/(unit_density * unit_length**3.0d0)
    rstar    = rstar/unit_length
    asound   = asound/unit_velocity
    rhobound = rhobound/unit_density
    Fbound   = Fbound/unit_radflux
    Ebound   = Ebound/unit_pressure

    if (mype == 0) then
      print*, '========================================'
      print*, '  Dimensionless computation quantities  '
      print*, '========================================'
      print*, 'Mstar    = ', mstar
      print*, 'Rstar    = ', rstar
      print*, 'rhobound = ', rhobound
      print*, 'asound   = ', asound
    endif

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixI^S, rho_) = rhobound
    w(ixI^S, mom(:)) = 0.d0
    w(ixI^S, p_) = asound**2*rhobound
    w(ixI^S, Er_) = Ebound-(x(ixI^S,1)-xprobmin1)*3*w(ixI^S, rho_)*m1_kappa_e*Fbound/m1_c
    w(ixI^S, Fr(1)) = Fbound

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions


  subroutine special_bound(qt,ixI^L,ixB^L,iB,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: ir, ix1, ixOmin1
    double precision :: prad(ixI^S)
    double precision :: gradE
    double precision :: ff(ixI^S), chi(ixI^S)

    select case (iB)
    case(1)
      w(ixB^S, Er_) = Ebound-(x(ixB^S,1)-xprobmin1)*3*w(ixB^S, rho_)*m1_kappa_e*Fbound/m1_c
      ! w(ixB^S, Fr(1)) = Fbound

    case(2)
      w(ixB^S, Er_) = Ebound-(x(ixB^S,1)-xprobmin1)*3*w(ixB^S, rho_)*m1_kappa_e*Fbound/m1_c
      ! w(ixB^S, Fr(1)) = Fbound

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound


  subroutine stellar_flux(ixI^L,ixO^L,w,x,F0)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: w(ixI^S,1:ndim)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out) :: F0(ixI^S,1:ndir)

    F0(ixI^S,:) = 0.d0

    ! Only in radial direction
    F0(ixI^S,1) = Fbound

  end subroutine stellar_flux


  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    use mod_hd

    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: grav(ixI^S), grad(ixI^S), prad(ixI^S), v(ixI^S), gradv(ixI^S), pth(ixI^S)
    double precision :: T(ixI^S), divT(ixI^S)

    !> Photon tiring + Radflux CMF correction
    v(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, prad)
    call gradient_str(v, x, ixI^L, ixO^L, 1, gradv)

    w(ixI^S,i_Prad) = prad(ixI^S)
    w(ixO^S,i_pht) = - m1_RSOL*gradv(ixO^S)*prad(ixO^S)
    w(ixO^S,i_fcmf) = - m1_RSOL*w(ixO^S,Fr(1))*gradv(ixO^S)

    call hd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,i_Tgas) = (pth(ixO^S)*unit_pressure*const_mp*hd_mu) &
                   /(w(ixO^S,rho_)*unit_density*const_kB)

    w(ixO^S,i_Trad) = (w(ixO^S,Er_)*unit_pressure/const_rad_a)**0.25d0

    w(ixI^S,i_ff) = sum(w(ixI^S, Fr(:))**2, dim=ndim+1)/(m1_c*w(ixI^S,Er_))**2
    w(ixI^S,i_chi) = (3.d0+4*w(ixI^S,i_ff))/(5.d0+2*dsqrt(4.d0-3*w(ixI^S,i_ff)))
    w(ixI^S,i_lum) = 4*dpi*(x(ixO^S,1)*unit_length)**2 * w(ixO^S,Fr(1))*unit_radflux / lstar

    !> Flux equation terms:
    !> Directional derivative
    w(ixO^S,i_SFdv) = -w(ixO^S,Fr(1))*gradv(ixO^S)

    !> Shadowing radiation force source term
    w(ixI^S,i_SgR) = -m1_c*w(ixI^S,kF_)*w(ixI^S,rho_)*w(ixI^S,Fr(1))

    !> advected flux term
    T(ixI^S) = w(ixI^S,Fr(1))*w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradient_str(T, x, ixI^L, ixO^L, 1, divT)
    w(ixO^S,i_SdivF) = -divT(ixO^S)

    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, T)
    T = m1_c**2 * T
    call gradient_str(T, x, ixI^L, ixO^L, 1, divT)
    w(ixO^S,i_SdivP) = -divT(ixO^S)

    if (x(1,1) .lt. xprobmin1) then

    print*, 'iteration', it
    print*, '*-----------------------'
    print*, 'ff', w(1:5,i_ff)
    print*, 'chi', w(1:5,i_chi)
    print*, '*-----------------------'
    print*, 'S', w(1:5,i_SgR)
    print*, '*-----------------------'
    print*, 'E', w(1:5,Er_)/Ebound
    print*, 'F', w(1:5,Fr(1))/Fbound
    print*, 'P', 3*w(1:5,i_Prad)/Ebound
    print*, 'P/E', w(1:5,i_Prad)/w(1:5,Er_)
    print*, 'Div P', 'xxxxx', 'xxxxx', w(3:7,i_SdivP)

    end if

  end subroutine update_extravars

end module mod_usr
