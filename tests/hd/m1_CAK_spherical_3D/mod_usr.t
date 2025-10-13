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
  double precision :: cak_alpha, gayley_qbar

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

    ! Choose independent normalization units if using dimensionless variables.
    unit_length      = rstar    ! cm
    unit_temperature = teff    ! K
    unit_density     = rhobound ! g cm^-3

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /star_list/ mstar, rstar, teff, rhobound, beta

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

    !> Initiate non-radial comps
    w(ixI^S,:) = 0.d0

    ! Small offset (~vtherm/vinf) to avoid starting at terminal wind speed
    sfac = 1.0d0 - 1.0d-3**(1.0d0/beta)

    !> Gas subsystem
    w(ixI^S,mom(1)) = vinf * ( 1.0d0 - sfac * rstar / x(ixI^S,1) )**beta
    w(ixI^S,rho_)   = mdot / (4.0d0*dpi * x(ixI^S,1)**2.0d0 * w(ixI^S,mom(1)))
    w(ixI^S,p_) = asound**2.0 * rhobound * (w(ixI^S,rho_)/rhobound)**hd_gamma

    !> Radiation subsystem
    ! w(ixO^S,Fr(1)) = w(ixO^S,Er_)*const_c/unit_velocity
    w(ixI^S,Fr(1)) = lstar/(4*dpi*(x(ixI^S,1)*unit_length)**2)/unit_radflux
    w(ixI^S,Er_) = w(ixI^S,Fr(1))/const_c*unit_velocity
    ! w(ixI^S,Er_) = const_rad_a*teff**4/unit_pressure

    ! if (mype .eq. 0) print*, 'F ic', w(1:5,Fr(1))
    ! if (mype .eq. 0) print*, 'E ic', w(1:5,Er_)

    call hd_to_conserved(ixI^L,ixO^L,w,x)

    ! if (any(w(ixO^S,:) .ne. w(ixO^S,:))) then
    !   print*, 'after ic'
    !   print*, w(ixImin1:5,5,5,:)
    !   print*, x(ixImin1:5,5,5,:)
    !   stop
    ! end if

  end subroutine initial_conditions

  subroutine special_bound(qt,ixI^L,ixB^L,iB,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: ir, ix1

    select case (iB)
    case(1)

      call hd_to_primitive(ixI^L,ixI^L,w,x)

      ! Mass density fixed at boundary mass density
      w(ixB^S,rho_) = rhobound

      ! Radial velocity field (constant slope extrapolation)
      ! do ix1 = ixBmax1,ixBmin1,-1
      !   if (ix1 == ixBmax1) then
      !     w(ix1,:,:,mom(1)) = w(ix1+1,:,:,mom(1)) - (w(ix1+2,:,:,mom(1)) - w(ix1+1,:,:,mom(1))) &
      !                       * (x(ix1+1,:,:,1) - x(ix1,:,:,1))/(x(ix1+2,:,:,1) - x(ix1+1,:,:,1))
      !   else
      !     w(ix1,:,:,mom(1)) = w(ix1+1,:,:,mom(1))
      !   endif
      ! enddo
      do ix1 = ixBmax1,ixBmin1,-1
        ! w(ix1,mom(1)) = max(2 * w(ix1+1,mom(1)) - w(ix1+2,mom(1)), w(ix1+1,mom(1)))
        w(ix1,:,:,mom(1)) = 2 * w(ix1+1,:,:,mom(1)) - w(ix1+2,:,:,mom(1))
      end do

      ! Avoid supersonic ghost cells, also avoid overloading too much
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)), 0.5d0*asound)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), -asound)

      w(ixB^S,mom(2)) = 0.d0
      w(ixB^S,mom(3)) = 0.d0

      ! Gas pressure
      ! w(ixB^S,p_) = asound**2.0 * rhobound * (w(ixB^S,rho_)/rhobound)**hd_gamma
      do ix1 = ixBmax1,ixBmin1,-1
        ! w(ix1,:,:,p_) = w(ix1+1,:,:,p_)
        w(ix1,:,:,p_) = max(2 * w(ix1+1,:,:,p_) - w(ix1+2,:,:,p_), w(ix1+1,:,:,p_))
      end do

      ! Radiation flux
      w(ixB^S,Fr(1)) = lstar/(4*dpi*(x(ixB^S,1)*unit_length)**2)/unit_radflux

      w(ixB^S,Fr(2)) = 0.d0
      w(ixB^S,Fr(3)) = 0.d0

      ! Radiation energy density
      ! w(ixB^S,Er_) = w(ixB^S,Fr(1))/const_c*unit_velocity
      do ix1 = ixBmax1,ixBmin1,-1
        ! w(ix1,Er_) = w(ix1+1,Er_)
        w(ix1,:,:,Er_) = max(2 * w(ix1+1,:,:,Er_) - w(ix1+2,:,:,Er_), w(ix1+1,:,:,Er_))
      end do

      call hd_to_conserved(ixI^L,ixI^L,w,x)

      ! if (any(w(ixB^S,:) .ne. w(ixB^S,:))) then
      !   print*, 'after bc'
      !   print*, w(ixBmin1:5,5,5,rho_)
      !   print*, w(ixBmin1:5,5,5,mom(1))
      !   print*, w(ixBmin1:5,5,5,mom(2))
      !   print*, w(ixBmin1:5,5,5,mom(3))
      !   print*, w(ixBmin1:5,5,5,e_)
      !   print*, w(ixBmin1:5,5,5,Er_)
      !   print*, w(ixBmin1:5,5,5,Fr(1))
      !   print*, w(ixBmin1:5,5,5,Fr(2))
      !   print*, w(ixBmin1:5,5,5,Fr(3))
      !   print*, x(ixBmin1:5,5,5,1)
      !   print*, '================'
      !   print*, w(ixB^S,:)
      !   stop
      ! end if


    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

  subroutine stellar_gravity(ixI^L,ixO^L,wCT,x,gravity_field)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(out) :: gravity_field(ixI^S,ndim)

    gravity_field(ixO^S,:) = 0.0d0

    ! Only in radial direction
    gravity_field(ixO^S,1) = -Ggrav * mstar/x(ixO^S,1)**2.0d0

  end subroutine stellar_gravity

end module mod_usr
