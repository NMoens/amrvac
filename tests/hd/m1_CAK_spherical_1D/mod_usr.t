!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd

  implicit none

  ! Custom global variables can be defined here
  ! User input parameters
  double precision :: mstar, rstar, teff, rhobound, beta
  double precision :: mdot, vinf, asound, Ggrav

  double precision :: lstar, mumol, vesc, gammae, logg, logge, heff
  double precision :: cak_alpha, gayley_qbar, fred_bound

  integer :: i_Gamma, i_Mf, i_prad, i_pht, i_fcmf, i_Tgas, i_Trad, i_chi, i_ff

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("spherical")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see src/mod_usr_methods.t
    usr_special_bc     => special_bound
    usr_gravity        => stellar_gravity

    ! Fix negative radiation flux
    usr_internal_bc => Fix_radflux

    !> Additional variables
    usr_process_grid => update_extravars

    ! Choose independent normalization units if using dimensionless variables.
    unit_length      = rstar    ! cm
    unit_temperature = teff    ! K
    unit_density     = rhobound ! g cm^-3

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

    !> Extra variables
    i_Gamma = var_set_extravar('Gamma','Gamma')
    i_Mf = var_set_extravar('M_F','M_F')
    i_Prad = var_set_extravar('P_rr','P_rr')
    i_pht = var_set_extravar('pht','pht')
    i_fcmf = var_set_extravar('Fcmf','Fcmf')
    i_Tgas = var_set_extravar('Tgas','Tgas')
    i_Trad = var_set_extravar('Trad','Trad')
    i_ff = var_set_extravar('ff','ff')
    i_chi = var_set_extravar('chi','chi')

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /star_list/ mstar, rstar, teff, rhobound, beta, fred_bound

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

    !> initial conditions
    cak_alpha = 0.66d0
    gayley_qbar = 2000.d0

    ! Stellar structure
    lstar  = 4.0d0*dpi * rstar**2.0d0 * const_sigma * teff**4.0d0
    gammae = const_kappae * lstar/(4.0d0*dpi * const_G * mstar * const_c)
    logg   = log10(const_G * mstar/rstar**2.0d0)
    logge  = logg + log10(1.0d0 - gammae)
    mumol  = (1.0d0 + 4.0d0*He_abundance)/(2.0d0 + 3.0d0*He_abundance)
    asound = sqrt(teff * kB_cgs/(mumol * mp_cgs))
    heff   = asound**2.0d0 / 10.0d0**logge

    ! Wind quantities in CAK theory
    vesc   = sqrt(2.0d0 * const_G * mstar * (1.0d0 - gammae)/rstar)
    vinf   = vesc * sqrt(cak_alpha/(1.0d0 - cak_alpha))
    mdot   = lstar/const_c**2.0d0 * cak_alpha/(1.0d0 - cak_alpha) &
               * (gayley_qbar * gammae/(1.0d0 - gammae))**((1.0d0 - cak_alpha)/cak_alpha)

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
      print*, 'Twind                  = ', teff
      print*, 'Mean molecular weight  = ', mumol
      print*, 'log(g)                 = ', logg
      print*, 'eff. log(g)            = ', logge
      print*, 'heff/Rstar             = ', heff/rstar
      print*, 'Eddington gamma        = ', gammae
      print*, 'adiabatic gamma        = ', hd_gamma
      print*, 'isothermal asound      = ', asound
      print*, 'eff. vesc              = ', vesc
      print*, 'CAK vinf               = ', vinf
      print*, 'FD vinf                = ', 3.0d0 * vesc
      print*, 'analytic Mdot CAK      = ', mdot * (const_years/const_MSun)
      print*, '... with FD correction = ', mdot/(1.0d0 + cak_alpha)**(1.0d0/cak_alpha) * (const_years/const_MSun)
    endif

    ! Make some variables needed for computations dimensionless
    mstar    = mstar/(unit_density * unit_length**3.0d0)
    rstar    = rstar/unit_length
    asound   = asound/unit_velocity
    rhobound = rhobound/unit_density
    vinf     = vinf/unit_velocity
    mdot     = mdot * unit_time/(unit_density * unit_length**3.0d0)
    Ggrav    = const_G * unit_density * unit_time**2.0d0

    if (mype == 0) then
      print*, '========================================'
      print*, '  Dimensionless computation quantities  '
      print*, '========================================'
      print*, 'Mstar    = ', mstar
      print*, 'Rstar    = ', rstar
      print*, 'rhobound = ', rhobound
      print*, 'Mdot CAK = ', mdot
      print*, 'asound   = ', asound
      print*, 'vinf     = ', vinf
      print*, 'Ggrav    = ', Ggrav
    endif

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! Local variables
    double precision :: sfac

    ! Small offset (~vtherm/vinf) to avoid starting at terminal wind speed
    ! sfac = 1.0d0 - 1.0d-3**(1.0d0/beta)
    sfac = 1.d0

    !> Gas subsystem
    w(ixI^S,mom(1)) = vinf * ( 1.0d0 - sfac * rstar / x(ixI^S,1) )**beta
    w(ixI^S,rho_)   = mdot / (4.0d0*dpi * x(ixI^S,1)**2.0d0 * w(ixI^S,mom(1)))

    if (hd_energy) &
      w(ixI^S,p_) = asound**2.0 * rhobound * (w(ixI^S,rho_)/rhobound)**hd_gamma

    !> Radiation subsystem
    ! w(ixO^S,Fr(1)) = w(ixO^S,Er_)*const_c/unit_velocity
    w(ixI^S,Fr(1)) = lstar/(4*dpi*(x(ixI^S,1)*unit_length)**2)/unit_radflux
    w(ixI^S,Er_) = w(ixI^S,Fr(1))/m1_c
    ! w(ixI^S,Er_) = const_rad_a*teff**4/unit_pressure

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

  subroutine special_bound(qt,ixI^L,ixB^L,iB,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: ir, ix1
    double precision :: prad(ixI^S)
    double precision :: gradE(ixI^S)
    double precision :: ff(ixI^S), chi(ixI^S)

    select case (iB)
    case(1)
      call hd_to_primitive(ixI^L,ixI^L,w,x)

      ! Mass density fixed at boundary mass density
      w(ixBmax1,rho_) = rhobound
      do ix1 = ixBmax1-1, ixBmin1, -1
        w(ix1,rho_) = 10.d0**max(2 * dlog10(w(ix1+1,rho_)) - dlog10(w(ix1+2,rho_)), dlog10(w(ix1+1,rho_)))
      end do

      do ix1 = ixBmax1,ixBmin1,-1
        ! w(ix1,mom(1)) = max(2 * w(ix1+1,mom(1)) - w(ix1+2,mom(1)), w(ix1+1,mom(1)))
        w(ix1,mom(1)) = 2 * w(ix1+1,mom(1)) - w(ix1+2,mom(1))
      end do

      ! Avoid supersonic ghost cells, also avoid overloading too much
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)), asound)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), 0.d0)

      ! Gas pressure
      if (hd_energy) then
        ! w(ixB^S,p_) = asound**2.0 * rhobound * (w(ixB^S,rho_)/rhobound)**hd_gamma
        do ix1 = ixBmax1,ixBmin1,-1
          w(ix1,p_) = 2 * w(ix1+1,p_) - w(ix1+2,p_)
          ! w(ix1,p_) = max(w(ix1,p_), w(ix1+1,p_))
        end do
      end if

      ! Radiation flux
      w(ixB^S,Fr(1)) = lstar/(4*dpi*(x(ixB^S,1)*unit_length)**2)/unit_radflux

      ! Radiation energy
      !> Manual ff
      ff(ixB^S) = fred_bound**2
      w(ixB^S,Er_) = w(ixB^S,Fr(1))/(m1_c*dsqrt(ff(ixB^S) ))


      call hd_to_conserved(ixI^L,ixI^L,w,x)

    case (2)

      do ix1 = ixBmin1,ixBmax1
        w(ix1,mom(1)) = w(ix1-1,mom(1))!*(x(ix1,1)/x(ix1-1,1))**2
        w(ix1,Fr(1)) = w(ix1-1,Fr(1))*(x(ix1,1)/x(ix1-1,1))**2

        w(ix1,mom(1)) = max(0.d0,w(ix1,mom(1)))
        w(ix1,Fr(1)) = max(0.d0,w(ix1,Fr(1)))
      end do

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

  subroutine Fix_radflux(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    ! Local variables
    double precision :: sfac

    ! ! Small offset (~vtherm/vinf) to avoid starting at terminal wind speed
    ! ! sfac = 1.0d0 - 1.0d-3**(1.0d0/beta)
    ! sfac = 1.d0
    !
    ! !> Gas subsystem
    ! w(ixI^S,mom(1)) = vinf * ( 1.0d0 - sfac * rstar / x(ixI^S,1) )**beta
    ! w(ixI^S,rho_)   = mdot / (4.0d0*dpi * x(ixI^S,1)**2.0d0 * w(ixI^S,mom(1)))
    !
    ! if (hd_energy) &
    !   w(ixI^S,p_) = asound**2.0 * rhobound * (w(ixI^S,rho_)/rhobound)**hd_gamma
    ! call hd_to_conserved(ixI^L,ixO^L,w,x)

    !> Radiation subsystem
    ! w(ixI^S,Fr(1)) = lstar/(4*dpi*(x(ixI^S,1)*unit_length)**2)/unit_radflux
    ! w(ixI^S,Er_) = lstar/(4*dpi*(x(ixI^S,1)*unit_length)**2)/unit_radflux/m1_c

    ! where(w(ixI^S,mom(1)) .lt. 0.d0) w(ixI^S,mom(1)) = 0.d0
    ! where(w(ixI^S,Fr(1)) .lt. 0.d0) w(ixI^S,Fr(1)) = 0.d0

  end subroutine Fix_radflux

  subroutine stellar_gravity(ixI^L,ixO^L,wCT,x,gravity_field)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(out) :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:) = 0.0d0

    ! Only in radial direction
    gravity_field(ixO^S,1) = -Ggrav * mstar/x(ixO^S,1)**2.0d0

  end subroutine stellar_gravity

  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    use mod_hd

    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: grav(ixI^S), grad(ixI^S), prad(ixI^S), v(ixI^S), gradv(ixI^S), pth(ixI^S)

    grav(ixI^S) = w(ixI^S,rho_) * Ggrav * mstar/x(ixO^S,1)**2.0d0
    grad(ixI^S) = w(ixI^S,rho_) * w(ixI^S,kF_) * w(ixI^S,Fr(1)) / m1_c

    w(ixI^S,i_Gamma) = grad(ixI^S)/grav(ixI^S)
    w(ixI^S,i_Mf) = w(ixI^S,kF_)/(0.34d0/unit_opacity)

    !> Photon tiring + Radflux CMF correction
    v(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, prad)
    call gradient_str(v, x, ixI^L, ixO^L, 1, gradv)

    w(ixO^S,i_Prad) = prad(ixO^S)
    w(ixO^S,i_pht) = - m1_RSOL*gradv(ixO^S)*prad(ixO^S)
    w(ixO^S,i_fcmf) = - m1_RSOL*w(ixO^S,Fr(1))*gradv(ixO^S)

    call hd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,i_Tgas) = (pth(ixO^S)*unit_pressure*const_mp*hd_mu) &
                   /(w(ixO^S,rho_)*unit_density*const_kB)

    w(ixO^S,i_Trad) = (w(ixO^S,Er_)*unit_pressure/const_rad_a)**0.25d0

    w(ixO^S,i_ff) = sum(w(ixO^S, Fr(:))**2, dim=ndim+1)/(m1_c*w(ixO^S,Er_))**2
    w(ixO^S,i_chi) = (3.d0+4*w(ixO^S,i_ff))/(5.d0+2*dsqrt(4.d0-3*w(ixO^S,i_ff)))

  end subroutine update_extravars

end module mod_usr
