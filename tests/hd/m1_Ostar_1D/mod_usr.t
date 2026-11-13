!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd
  use mod_m1

  implicit none

  ! Custom global variables can be defined here
  ! User input parameters
  double precision :: mstar, rstar, Tbound, rhobound, beta
  double precision :: asound, Ggrav

  double precision :: lstar, vesc, gamma_e, logg, logge, heff
  double precision :: fred_bound

  integer :: i_Gamma, i_Mf, i_prad, i_pht, i_fcmf, i_Tgas, i_Trad, i_chi, i_ff, i_lum, i_lumg, i_lumr
  integer :: i_SdivF, i_SdivP, i_SFdv, i_SgR, i_ScdvP, i_tau_cell, i_LFLD, i_theta

  logical :: fixed_lum = .true.

  integer :: len_ic = 3070
  double precision, allocatable :: r_init(:), rho_init(:), vr_init(:), p_init(:), Er_init(:)
  character(len=32) :: ic_file = '1D_stable.txt'

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
    usr_set_F0 => fld_flux

    !> Additional variables
    usr_process_grid => update_extravars

    ! Choose independent normalization units if using dimensionless variables.
    unit_length      = rstar    ! cm
    unit_temperature = Tbound    ! K
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
    i_lum = var_set_extravar('lum','lum')

    i_lumg = var_set_extravar('lum_gas','lum_gas')
    i_lumr = var_set_extravar('lum_rad','lum_rad')


    i_LFLD = var_set_extravar('lum_fld','lum_fld')

    i_theta = var_set_extravar('theta','theta')

    i_tau_cell = var_set_extravar('tau_cell','tau_cell')

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /star_list/ mstar, rstar, Tbound, rhobound, beta, fred_bound, gamma_e, ic_file, len_ic, fixed_lum

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
    logg   = log10(const_G * mstar/rstar**2)
    logge  = logg + log10(1.0d0 - gamma_e)
    asound = sqrt(Tbound * kB_cgs/(hd_mu * mp_cgs))
    heff   = asound**2.0d0 / 10.0d0**logge

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
    Ggrav    = const_G * unit_density * unit_time**2.0d0

    if (mype == 0) then
      print*, '========================================'
      print*, '  Dimensionless computation quantities  '
      print*, '========================================'
      print*, 'Mstar    = ', mstar
      print*, 'Rstar    = ', rstar
      print*, 'rhobound = ', rhobound
      print*, 'asound   = ', asound
      print*, 'Ggrav    = ', Ggrav
    endif

    !> Allocate initial conditions arrays
    allocate(r_init(len_ic))
    allocate(rho_init(len_ic))
    allocate(vr_init(len_ic))
    allocate(p_init(len_ic))
    allocate(Er_init(len_ic))

    call read_init_conds()

  end subroutine initglobaldata_usr


  subroutine read_init_conds
    integer :: ii

    open(unit=1, file=trim(ic_file))
    read(1,*) !> first line of data
    do ii = 1,len_ic-2
      read(1,*) r_init(ii), rho_init(ii), vr_init(ii), p_init(ii), Er_init(ii)
    enddo
    CLOSE(1)

  end subroutine read_init_conds

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: ix1, nearest_i

    do ix1 = ixImin1, ixImax1
      !> Find nearest line in init
      nearest_i = minloc(abs(r_init - x(ix1,1)),DIM=1)

      !> Read closest values
      w(ix1, rho_) = rho_init(nearest_i)/unit_density
      w(ix1, mom(1)) = vr_init(nearest_i)/unit_velocity
      w(ix1, p_) = p_init(nearest_i)/unit_pressure
      w(ix1, Er_) = Er_init(nearest_i)/unit_pressure
    end do

    call hd_to_conserved(ixI^L,ixO^L,w,x)

    w(ixI^S,Fr(1)) = lstar/(4*dpi*(x(ixI^S,1)*unit_length)**2)/unit_radflux

  end subroutine initial_conditions

  subroutine special_bound(qt,ixI^L,ixB^L,iB,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: ir, ix1, ixOmin1
    double precision :: prad(ixI^S)
    double precision :: ff(ixI^S), chi(ixI^S)
    double precision :: F_adv(ixB^S)

    double precision :: gradE, D
    double precision :: fs, bs, cs

    select case (iB)
    case(1)

      call hd_to_primitive(ixI^L,ixI^L,w,x)

      ! Mass density fixed at boundary mass density
      w(ixBmax1,rho_) = rhobound
      do ix1 = ixBmax1-1, ixBmin1, -1
        w(ix1,rho_) = 10.d0**max(2 * dlog10(w(ix1+1,rho_)) - dlog10(w(ix1+2,rho_)), dlog10(w(ix1+1,rho_)))
        w(ix1,rho_) = max(w(ix1,rho_), rhobound)
      end do

      do ix1 = ixBmax1,ixBmin1,-1
        ! w(ix1,mom(1)) = max(2 * w(ix1+1,mom(1)) - w(ix1+2,mom(1)), w(ix1+1,mom(1)))
        w(ix1,mom(1)) = 2 * w(ix1+1,mom(1)) - w(ix1+2,mom(1))
      end do

      ! Avoid supersonic ghost cells, also avoid overloading too much
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)), asound)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), -0.d0*asound)

      !!! OVERWRITE VELOCITY TO 0
      w(ixB^S,mom(1)) = -w(ixBmax1+1,mom(1))

      ! ! Gas pressure.
      ! if (hd_energy) then
      !   !w(ixB^S,p_) = asound**2.0 * rhobound * (w(ixB^S,rho_)/rhobound)**hd_gamma
      !   do ix1 = ixBmax1, ixBmin1, -1
      !     w(ix1,p_) = 10.d0**(2 * dlog10(w(ix1+1,p_)) - dlog10(w(ix1+2,p_)))
      !   end do
      ! end if

      if (fixed_lum) then

        ! advected flux !TODO: Kinetic energy??
        F_adv(ixB^S) = (4.d0/3.d0*w(ixB^S,Er_) + hd_gamma*w(ixB^S,p_) + 0.5d0*w(ixB^S,mom(1))**2/w(ixB^S,rho_))*w(ixB^S,mom(1))

        ! Radiation flux
        w(ixB^S,Fr(1)) = lstar/(4*dpi*(x(ixB^S,1)*unit_length)**2)/unit_radflux &
                       - F_adv(ixB^S) !w(ixB^S,Er_)*w(ixB^S,mom(1))

        ixOmin1 = ixBmax1 + 1

        do ix1 = ixOmin1,ixBmin1+1,-1
          gradE = lstar/(4*dpi*(x(ix1,1)*unit_length)**2)/unit_radflux &
                - 4.d0/3*w(ix1,Er_)*w(ix1,mom(1))
          ! gradE = -3*w(ixOmin1,kF_)*w(ix1,rho_)*w(ix1,Fr(1))/m1_c
          gradE = -3*w(ixOmin1,kF_)*w(ix1,rho_)*gradE/m1_c

          !> Forward, backward and central derivs, needed for stretched grid gradient
          fs = (x(ix1,1) - x(ix1-1,1)) &
          / ((x(ix1+1,1) - x(ix1,1)) * (x(ix1+1,1) - x(ix1-1,1)))

          bs = -(x(ix1+1,1) - x(ix1,1)) &
          / ((x(ix1,1) - x(ix1-1,1)) * (x(ix1+1,1) - x(ix1-1,1)))

          cs = (x(ix1+1,1) + x(ix1-1,1) - 2*x(ix1,1)) &
          / ((x(ix1,1) - x(ix1-1,1)) * (x(ix1+1,1) - x(ix1,1)))

          w(ix1-1,Er_) = (gradE - fs*w(ix1+1,Er_) - cs*w(ix1,Er_))/bs

        end do

        ! Gas pressure
        if (hd_energy) then
          w(ixB^S,p_) = (w(ixB^S,Er_)*unit_pressure/const_rad_a)**0.25d0
          w(ixB^S,p_) = const_kB/(const_mp*hd_mu) * w(ixB^S,p_) * w(ixB^S,rho_)*unit_density
          w(ixB^S,p_) = w(ixB^S,p_)/unit_pressure
          !> limits
          w(ixB^S,p_) = min(w(ixB^S,p_),1.d1)
          w(ixB^S,p_) = max(w(ixB^S,p_),1.d-1)
        end if

        ! print*, it, '|' ,w(ixB^S,Er_), '|' ,w(ixB^S,Fr(1))

        ! print*, it, '------------------', w(ixBmax1,Fr(1))*x(ixBmax1,1)**2
        ! print*, 'kappa', w(1:5,kF_)/m1_kappa_e
        ! print*, 'E', w(1:5,Er_)
        ! print*, 'F', (w(1:5,Fr(1))+F_adv(1:5))*x(1:5,1)**2
        ! print*, 'vE', F_adv(1:2)/w(1:2,Fr(1))
        !
        ! stop


      else
        !> Gas pressure
        w(ixB^S,p_) = const_kB*Tbound*w(ixB^S,rho_)*unit_density/(const_mp*hd_mu)/unit_pressure
        w(ixB^S,Er_) = const_rad_a*Tbound**4/unit_pressure

        ! !> Radiation energy density (log extrapolated)
        ! do ix1 = ixBmax1-1, ixBmin1, -1
        !   w(ix1,Er_) = 10.d0**max(2 * dlog10(w(ix1+1,Er_)) - dlog10(w(ix1+2,Er_)), dlog10(w(ix1+1,Er_)))
        !   w(ix1,Er_) = max(w(ix1,Er_), w(ixBmax1,Er_))
        ! end do

        !> Radiation flux (const lum)
        w(ixB^S,Fr(:)) = 0.d0
        do ix1 = ixBmax1,ixBmin1,-1
          w(ix1,Fr(1)) = (w(ix1+1,Fr(1))-4.d0/3*w(ix1+1,mom(1))*w(ix1+1,Er_) )* (x(ix1+1,1)/x(ix1,1))**2
          ! w(ix1,Fr(1)) = 0.333d0*(w(ix1+1,Fr(1)) * (x(ix1+1,1)/x(ix1,1))**2 &
          !                       + w(ix1+2,Fr(1)) * (x(ix1+2,1)/x(ix1,1))**2 &
          !                       + w(ix1+3,Fr(1)) * (x(ix1+3,1)/x(ix1,1))**2)
        end do

        !> Radiation energy density (FLD extrapolated)
        do ix1 = ixBmax1-1, ixBmin1, -1
          ! Evaluate F and kappa*rho at i+1/2 (hence 1/2*1/2)
          D = 3.d0/4*(w(ix1,Fr(1)) + w(ix1+1,Fr(1)))*(w(ix1,kF_)*w(ix1,rho_) + w(ix1+1,kF_)*w(ix1+1,rho_))/m1_c
          w(ix1,Er_) = w(ix1,Er_) + D*(x(ix1+1,1) - x(ix1,1))
        end do

        !> Radiation flux (FLD extrapolated)
        ! do ix1 = ixBmax1,ixBmin1,-1
        !   !> Linear derivative
        !   w(ix1,Fr(1)) = -m1_c/(3*w(ix1,kF_)*w(ix1,rho_)) * (w(ix1+1,Er_) - w(ix1,Er_))/(x(ix1+1,1) - x(ix1,1))
        !   !> Log derivative
        !   ! w(ix1,Fr(1)) = -m1_c/(3*w(ix1,kF_)*w(ix1,rho_)) * (log(w(ix1+1,Er_)) - log(w(ix1,Er_)))/(x(ix1+1,1) - x(ix1,1))* w(ix1,Er_)
        ! end do

        ! print*, it, '------------------', w(ixBmax1,Fr(1))*x(ixBmax1,1)**2
        ! print*, 'kappa', w(1:5,kF_)/m1_kappa_e
        ! print*, 'E', w(1:5,Er_)/w(ixBmax1,Er_)
        ! print*, 'F', w(1:5,Fr(1))*x(1:5,1)**2 / (w(ixBmax1,Fr(1))*x(ixBmax1,1)**2)
        ! print*, 'vE', 4.d0/3*w(1:5,mom(1))*w(1:5,Er_)/w(1:5,Fr(1))

      end if

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

    double precision :: fld_flux(ixI^S), tau_cell(ixI^S)
    integer :: idir

    integer :: ix1, nearest_i

    do ix1 = ixImin1, ixImax1
      !> Find nearest line in init
      nearest_i = minloc(abs(r_init - x(ix1,1)),DIM=1)

      !> Read closest values
      w(ix1, rho_) = rho_init(nearest_i)/unit_density
      w(ix1, mom(1)) = vr_init(nearest_i)/unit_velocity!*w(ix1, rho_)
      w(ix1, p_) = p_init(nearest_i)/unit_pressure
      ! w(ix1, Er_) = Er_init(nearest_i)/unit_pressure
    end do

    call hd_to_conserved(ixI^L,ixO^L,w,x)


    !> Radiation subsystem
    ! w(ixI^S,Fr(1)) = lstar/(4*dpi*(x(ixI^S,1)*unit_length)**2)/unit_radflux
    ! w(ixI^S,Er_) = lstar/(4*dpi*(x(ixI^S,1)*unit_length)**2)/unit_radflux/m1_c

    ! where(w(ixI^.S,mom(1)) .lt. 0.d0) w(ixI^S,mom(1)) = 0.d0
    ! where(w(ixI^S,Fr(1)) .lt. 0.d0) w(ixI^S,Fr(1)) = 0.d0
    ! where(x(ixO^S,1) .lt. 1.2d0) w(ixO^S,Fr(1)) = lstar/(4*dpi*(x(ixO^S,1)*unit_length)**2)/unit_radflux&
    !                                       -w(ixO^S,mom(1))*w(ixO^S,Er_)/w(ixO^S,rho_)


  end subroutine Fix_radflux

  subroutine fld_flux(ixI^L,ixO^L,w,x,F0)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: w(ixI^S,1:ndim)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out) :: F0(ixI^S,1:ndir)

    double precision :: gradE(ixI^S), lambda(ixI^S), R(ixI^S)

    integer :: idir

    !> FLD approximation
    !> Get magnitude of grad E
    gradE = 0.d0
    do idir = 1, ndim
      call gradient_str(w(ixI^S,Er_),x,ixI^L,ixO^L,idir,lambda)
      gradE(ixI^S) = gradE(ixI^S) + lambda(ixI^S)**2
    end do
    gradE(ixI^S) = dsqrt(max(0.d0,gradE(ixI^S)))

    !> Get fld R
    R(ixI^S) = gradE(ixI^S)/(w(ixI^S,rho_)*w(ixI^S,kF_)*w(ixI^S,Er_))

    !> Get fld lambda
    lambda(ixI^S) = (2.d0+R(ixI^S))/(6.d0+3*R(ixI^S)+R(ixI^S)**2.d0)

    do idir = 1, ndim
      call gradient_str(w(ixI^S,Er_),x,ixI^L,ixO^L,idir,gradE)
      F0(ixO^S,idir) = -m1_c*lambda(ixO^S)/(w(ixO^S,rho_)*w(ixO^S,kF_))*gradE(ixO^S)
    end do

  end subroutine fld_flux

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
    double precision :: T(ixI^S), divT(ixI^S), F(ixI^S,1:ndim), F_adv_gas(ixI^S), F_adv_rad(ixI^S), theta(ixI^S)


    grav(ixI^S) = w(ixI^S,rho_) * Ggrav * mstar/x(ixI^S,1)**2.0d0
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
    w(ixO^S,i_lum) = 4*dpi*(x(ixO^S,1)*unit_length)**2 * w(ixO^S,Fr(1))*unit_radflux / lstar

    F_adv_gas(ixO^S) = (w(ixO^S,e_) + pth(ixO^S))*v(ixO^S)
    F_adv_rad(ixO^S) = (w(ixO^S,Er_) + prad(ixO^S))*v(ixO^S)

    w(ixO^S,i_lumg) = 4*dpi*(x(ixO^S,1)*unit_length)**2 * (F_adv_gas(ixO^S))*unit_radflux / lstar
    w(ixO^S,i_lumr) = 4*dpi*(x(ixO^S,1)*unit_length)**2 * (F_adv_rad(ixO^S))*unit_radflux / lstar


    call fld_flux(ixI^L,ixO^L,w,x,F)
    w(ixO^S,i_LFLD) = 4*dpi*F(ixO^S,1)*(x(ixO^S,1)*unit_length)**2*unit_radflux / lstar

    call m1_get_trapped_fraction(w, x, ixI^L, ixO^L, 1, theta)
    w(ixO^S, i_theta) = theta(ixO^S)

    w(ixO^S,i_tau_cell) = w(ixO^S,rho_)*w(ixO^S,kF_)*block%dx(ixO^S,1)

  end subroutine update_extravars

end module mod_usr
