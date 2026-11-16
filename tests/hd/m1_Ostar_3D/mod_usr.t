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

  integer :: i_Gamma, i_Mf, i_pht, i_fcmf, i_Tgas, i_Trad, i_chi, i_ff, i_lum, i_mdot, i_am
  integer :: i_SdivF, i_SdivP, i_SFdv, i_SgR, i_ScdvP, i_tau_cell
  integer :: i_Prr, i_Ptt, i_Ppp, i_Prt

  integer :: len_ic = 3000
  double precision, allocatable :: r_init(:), rho_init(:), vr_init(:), p_init(:), Er_init(:)
  character(len=32) :: ic_file = 'init_conds/Star89_init.txt'

  !> Equilibrium radiation energy density (Lucy 1971) on the level-1 radial grid
  double precision, allocatable :: E_eq(:), rc_eq(:)
  double precision :: Rph_eq

  logical :: fixed_lum = .true.

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
    i_Gamma = var_set_extravar('Gamma','Gamma')
    i_Mf = var_set_extravar('M_F','M_F')
    i_Prr = var_set_extravar('P_rr','P_rr')
    i_Ptt = var_set_extravar('P_tt','P_tt')
    i_Ppp = var_set_extravar('P_pp','P_pp')
    i_Prt = var_set_extravar('P_rt','P_rt')


    i_pht = var_set_extravar('pht','pht')
    i_fcmf = var_set_extravar('Fcmf','Fcmf')
    i_Tgas = var_set_extravar('Tgas','Tgas')
    i_Trad = var_set_extravar('Trad','Trad')
    i_ff = var_set_extravar('ff','ff')
    i_chi = var_set_extravar('chi','chi')
    i_lum = var_set_extravar('lum','lum')
    i_mdot = var_set_extravar('mdot','mdot')
    i_am = var_set_extravar('am','am')

    ! i_SdivF = var_set_extravar('Div_Fv','Div_Fv')
    ! i_SdivP = var_set_extravar('Div_c2P','Div_c2P')
    ! i_SFdv = var_set_extravar('DDv','DDv')
    ! i_SgR = var_set_extravar('gG','gG')
    ! i_ScdvP = var_set_extravar('corr_Div_c2P','corr_Div_c2P')

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

    !> Initial-condition tables and equilibrium integration are only needed for a new run:
    !> skip them when converting (-convert) or restarting from a snapshot (-if / -resume)
    if (.not. convert .and. restart_from_file == undefined) then

      !> Allocate initial conditions arrays
      allocate(r_init(len_ic))
      allocate(rho_init(len_ic))
      allocate(vr_init(len_ic))
      allocate(p_init(len_ic))
      allocate(Er_init(len_ic))

      if (mype .eq. 0) &
        call read_init_conds()

      call mpi_bcast(r_init, len_ic, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)
      call mpi_bcast(rho_init, len_ic, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)
      call mpi_bcast(vr_init, len_ic, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)
      call mpi_bcast(p_init, len_ic, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)
      call mpi_bcast(Er_init, len_ic, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)

      !> Equilibrium radiation field: computed on rank 0, broadcast to all
      allocate(E_eq(domain_nx1), rc_eq(domain_nx1))

      if (mype .eq. 0) call integrate_init_Erad()

      call mpi_bcast(E_eq,   domain_nx1, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)
      call mpi_bcast(rc_eq,  domain_nx1, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)
      call mpi_bcast(Rph_eq, 1,          MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)

    end if

  end subroutine initglobaldata_usr


  subroutine read_init_conds
    integer :: ii

    if (mype .eq. 0) print*, ic_file, len_ic, 'lines'

    open(unit=1, file=trim(ic_file))
    !read(1,*) !> first line of data
    do ii = 1,len_ic-2
      read(1,*) r_init(ii), rho_init(ii), vr_init(ii), p_init(ii), Er_init(ii)
    enddo
    CLOSE(1)

  end subroutine read_init_conds

  !> Radiation energy density from Lucy's (1971) unified atmosphere/wind approximation:
  !>   E(r) = L/(pi R_ph^2 c) W(r) + (3/c) int_r^rout rho kappa F dr'
  !> on the actual (possibly stretched) level-1 grid, with the same face-averaged
  !> rho*kappa as m1_flux_absorbtion. R_ph from tau_bar(R_ph) = 2/3, iterated.
  !> Called on rank 0 after the table broadcast; E_eq, rc_eq, Rph_eq are broadcast afterwards.
  subroutine integrate_init_Erad
    integer, parameter          :: max_pass = 30, max_it = 50
    double precision, parameter :: tol_R = 1.d-8, tol_E = 1.d-10, relax = 0.7d0

    integer          :: nx, i, k, it_k, ipass
    logical          :: found
    double precision :: q, dist, Fface, Lfac, Rph_old, Rph_new, Eold, g_i, g_ip1, wgt
    double precision :: Teff, Trad_ph, E_ph
    double precision, allocatable :: rf(:), rho_eq(:), rk(:), S(:)

    nx = domain_nx1
    allocate(rf(0:nx), rho_eq(nx), rk(nx), S(nx))

    !> Level-1 faces and centres, code units (same construction as mod_amr_solution_node)
    if (stretched_dim(1)) then
      if (stretch_type(1) /= stretch_uni) &
        call mpistop("integrate_init_Erad: only 'uni' stretching supported")
      q = qstretch(1,1)
      do i = 0, nx
        rf(i) = xprobmin1 + dxfirst(1,1)*(1.d0 - q**i)/(1.d0 - q)
      end do
    else
      do i = 0, nx
        rf(i) = xprobmin1 + i*(xprobmax1 - xprobmin1)/nx
      end do
    end if
    rf(nx) = xprobmax1
    rc_eq(1:nx) = 0.5d0*(rf(0:nx-1) + rf(1:nx))

    !> Density: same nearest-row lookup as initial_conditions (must stay identical!)
    do i = 1, nx
      k = minloc(abs(r_init(1:len_ic-2) - rc_eq(i)), DIM=1)
      rho_eq(i) = rho_init(k)                                   ! cgs
    end do

    Rph_eq = xprobmin1                                          ! first guess, code units

    do ipass = 1, max_pass
      Lfac = lstar/(dpi*(Rph_eq*unit_length)**2*const_c)        ! L/(pi R_ph^2 c), erg cm^-3

      !> Outer cell: pure dilution, S = 0
      S(nx)    = 0.d0
      E_eq(nx) = Lfac*lucy_W(rc_eq(nx), Rph_eq)
      rk(nx)   = rho_eq(nx)*kappa_cgs(rho_eq(nx), E_eq(nx))

      !> Inward integration of S with fixed point on kappa(T_rad)
      do i = nx-1, 1, -1
        dist  = (rc_eq(i+1) - rc_eq(i))*unit_length             ! centre-to-centre, cgs
        Fface = lstar/(4.d0*dpi*(rf(i)*unit_length)**2)         ! flux at shared face
        E_eq(i) = E_eq(i+1)
        do it_k = 1, max_it
          Eold    = E_eq(i)
          rk(i)   = rho_eq(i)*kappa_cgs(rho_eq(i), E_eq(i))
          S(i)    = S(i+1) + 3.d0*dist*Fface*0.5d0*(rk(i) + rk(i+1))/const_c
          E_eq(i) = Lfac*lucy_W(rc_eq(i), Rph_eq) + S(i)
          if (abs(E_eq(i) - Eold) <= tol_E*E_eq(i)) exit
        end do
      end do

      !> New photosphere: S(r) = L/(2 pi r^2 c)  <=>  tau_bar = 2/3
      found = .false.
      do i = nx-1, 1, -1
        g_i   = S(i)   - lstar/(2.d0*dpi*(rc_eq(i)  *unit_length)**2*const_c)
        g_ip1 = S(i+1) - lstar/(2.d0*dpi*(rc_eq(i+1)*unit_length)**2*const_c)
        if (g_i >= 0.d0 .and. g_ip1 < 0.d0) then
          wgt     = g_ip1/(g_ip1 - g_i)
          Rph_new = rc_eq(i+1) + wgt*(rc_eq(i) - rc_eq(i+1))
          found   = .true.
          exit
        end if
      end do
      if (.not. found) then
        if (S(nx-1) >= lstar/(2.d0*dpi*(rc_eq(nx-1)*unit_length)**2*const_c)) then
          Rph_new = rc_eq(nx)                                   ! outer edge still thick
          print*, 'WARNING integrate_init_Erad: outer boundary optically thick'
        else
          Rph_new = xprobmin1                                   ! whole domain thin
          print*, 'WARNING integrate_init_Erad: no photosphere in domain'
        end if
      end if

      Rph_old = Rph_eq
      Rph_eq  = Rph_old + relax*(Rph_new - Rph_old)             ! under-relaxed updat
      print '(a,i3,a,es14.7,a,es10.3)', ' Lucy pass ', ipass, &
                          '  R_ph = ', Rph_eq, '  dR = ', Rph_eq - Rph_old
      if (abs(Rph_new - Rph_old) <= tol_R*Rph_old) exit
    end do
    if (ipass > max_pass .and. mype == 0) &
      print*, 'WARNING integrate_init_Erad: R_ph not converged'

    !> Diagnostics: T_rad at R_ph should be ~ T_eff
    if (mype == 0) then
      Teff = (lstar/(4.d0*dpi*(Rph_eq*unit_length)**2*const_rad_a*const_c/4.d0))**0.25d0
      i = minloc(abs(rc_eq - Rph_eq), DIM=1)
      E_ph = E_eq(i)
      Trad_ph = (E_ph/const_rad_a)**0.25d0
      print*, 'Lucy: R_ph/R* = ', Rph_eq, '  T_eff = ', Teff, '  T_rad(R_ph) = ', Trad_ph
      print*, 'Lucy: E(base) = ', E_eq(1), '  E(top) = ', E_eq(nx), ' erg/cm^3'
      print*, 'Lucy: cE/F at top = ', E_eq(nx)*const_c/(lstar/(4.d0*dpi*(rc_eq(nx)*unit_length)**2))
    end if

    deallocate(rf, rho_eq, rk, S)
  end subroutine integrate_init_Erad

  !> Lucy dilution factor, W = 1/2 inside the photosphere
  pure double precision function lucy_W(r, Rph)
    double precision, intent(in) :: r, Rph
    if (r <= Rph) then
      lucy_W = 0.5d0
    else
      lucy_W = 0.5d0*(1.d0 - dsqrt(1.d0 - (Rph/r)**2))
    end if
  end function lucy_W

  !> OPAL opacity in cgs from cgs rho and E, same T floor as the m1 OPAL case
  double precision function kappa_cgs(rho_c, E_c)
    use mod_opal_opacity, only: set_opal_opacity
    double precision, intent(in) :: rho_c, E_c
    double precision :: Tr
    Tr = max((E_c/const_rad_a)**0.25d0, 1.d4)
    call set_opal_opacity(rho_c, Tr, kappa_cgs)
  end function kappa_cgs


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: wgt, Trad, rloc, Eloc
    integer :: ix1, nearest_i, ii

    do ix1 = ixImin1, ixImax1
      !> Find nearest line in init
      nearest_i = minloc(abs(r_init - x(ix1,ixOmin2,ixOmin3,1)),DIM=1)

      !> Read closest values
      w(ix1, :, :, rho_) = rho_init(nearest_i)/unit_density
      w(ix1, :, :, mom(1)) = vr_init(nearest_i)/unit_velocity
      w(ix1, :, :, mom(2)) = 0.d0
      w(ix1, :, :, mom(3)) = 0.d0
      w(ix1, :, :, p_) = p_init(nearest_i)/unit_pressure
      ! w(ix1, :, :, Er_) = Er_init(nearest_i)/unit_pressure

      !> Equilibrium E (Lucy), log-linear interpolation on the level-1 centres
      rloc = x(ix1,ixOmin2,ixOmin3,1)
      ii = minloc(abs(rc_eq - rloc), DIM=1)
      if (rloc < rc_eq(ii)) ii = ii - 1
      ii = min(max(ii, 1), domain_nx1-1)
      wgt = (rloc - rc_eq(ii))/(rc_eq(ii+1) - rc_eq(ii))
      wgt = min(max(wgt, 0.d0), 1.d0)                      ! ghost cells: clamp
      Eloc = dexp((1.d0-wgt)*dlog(E_eq(ii)) + wgt*dlog(E_eq(ii+1)))/unit_pressure
      w(ix1, :, :, Er_) = Eloc

      !> LTE: gas temperature equal to radiation temperature
      if (hd_energy) then
        Trad = (Eloc*unit_pressure/const_rad_a)**0.25d0
        w(ix1, :, :, p_) = const_kB/(const_mp*hd_mu)*Trad*w(ix1, :, :, rho_)*unit_density/unit_pressure
      end if
    end do

    !> Velocity perturbation
    ! w(ixI^S, mom(2)) = 1.d0*sin(x(ixI^S,2)*10)*sin((x(ixI^S,1)-1.d0)*10)

    call hd_to_conserved(ixI^L,ixO^L,w,x)

    w(ixI^S,Fr(1)) = lstar/(4*dpi*(x(ixI^S,1)*unit_length)**2)/unit_radflux
    w(ixI^S,Fr(2)) = 0.d0
    w(ixI^S,Fr(3)) = 0.d0

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
    double precision :: gradE(ixBmin2:ixBmax2,ixBmin3:ixBmax3)
    double precision :: fs(ixBmin2:ixBmax2,ixBmin3:ixBmax3), cs(ixBmin2:ixBmax2,ixBmin3:ixBmax3), bs(ixBmin2:ixBmax2,ixBmin3:ixBmax3)


    select case (iB)
    case(1)

      call hd_to_primitive(ixI^L,ixI^L,w,x)

      ! Mass density fixed at boundary mass density
      w(ixBmax1,:,:,rho_) = rhobound
      do ix1 = ixBmax1-1, ixBmin1, -1
        w(ix1,:,:,rho_) = 10.d0**max(2 * dlog10(w(ix1+1,:,:,rho_)) - dlog10(w(ix1+2,:,:,rho_)), dlog10(w(ix1+1,:,:,rho_)))
      end do

      w(ixB^S,mom(:)) = 0.d0
      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,:,:,mom(1)) = 2 * w(ix1+1,:,:,mom(1)) - w(ix1+2,:,:,mom(1))
      end do

      ! Avoid supersonic ghost cells, also avoid overloading too much
      w(ixB^S,mom(1)) = min(w(ixB^S,mom(1)), asound)
      w(ixB^S,mom(1)) = max(w(ixB^S,mom(1)), -0.d0*asound)


      if (fixed_lum) then

        ! Radiation flux
        w(ixB^S,Fr(:)) = 0.d0
        w(ixB^S,Fr(1)) = lstar/(4*dpi*(x(ixB^S,1)*unit_length)**2)/unit_radflux &
                       - w(ixB^S,Er_)*w(ixB^S,mom(1))

        ixOmin1 = ixBmax1 +1
        do ix1 = ixBmax1+1, ixBmin1+1, -1
          gradE(ixBmin2:ixBmax2,ixBmin3:ixBmax3) = lstar/(4*dpi*(x(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1)*unit_length)**2)/unit_radflux &
                - w(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,Er_)*w(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(1))
          gradE(ixBmin2:ixBmax2,ixBmin3:ixBmax3) = -3*w(ixBmax1 +1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,kF_)*w(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,rho_)*gradE(ixBmin2:ixBmax2,ixBmin3:ixBmax3)/m1_c

          !> Forward, backward and central derivs, needed for stretched grid gradient
          fs(ixBmin2:ixBmax2,ixBmin3:ixBmax3) = (x(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1) - x(ix1-1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1)) &
          / ((x(ix1+1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1) - x(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1)) * (x(ix1+1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1) - x(ix1-1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1)))

          bs(ixBmin2:ixBmax2,ixBmin3:ixBmax3) = -(x(ix1+1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1) - x(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1)) &
          / ((x(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1) - x(ix1-1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1)) * (x(ix1+1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1) - x(ix1-1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1)))

          cs(ixBmin2:ixBmax2,ixBmin3:ixBmax3) = (x(ix1+1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1) + x(ix1-1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1) - 2*x(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1)) &
          / ((x(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1) - x(ix1-1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1)) * (x(ix1+1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1) - x(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1)))

          w(ix1-1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,Er_) = (gradE(ixBmin2:ixBmax2,ixBmin3:ixBmax3) - fs(ixBmin2:ixBmax2,ixBmin3:ixBmax3)*w(ix1+1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,Er_) - cs(ixBmin2:ixBmax2,ixBmin3:ixBmax3)*w(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,Er_))/bs(ixBmin2:ixBmax2,ixBmin3:ixBmax3)
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

      else

        !> Gas pressure
        w(ixB^S,p_) = const_kB*Tbound*w(ixB^S,rho_)*unit_density/(const_mp*hd_mu)/unit_pressure

        !> Radiation energy density
        w(ixB^S,Er_) = const_rad_a*Tbound**4/unit_pressure

        !> Radiation flux
        w(ixB^S,Fr(:)) = 0.d0
        do ix1 = ixBmax1,ixBmin1,-1
          w(ix1,:,:,Fr(1)) = w(ix1+1,:,:,Fr(1)) * (x(ix1+1,:,:,1)/x(ix1,:,:,1))**2
        end do

      end if

      call hd_to_conserved(ixI^L,ixI^L,w,x)

    case (2)

      do ix1 = ixBmin1,ixBmax1
        w(ix1,:,:,mom(1)) = w(ix1-1,:,:,mom(1))*(x(ix1,:,:,1)/x(ix1-1,:,:,1))**2
        w(ix1,:,:,Fr(1)) = w(ix1-1,:,:,Fr(1))*(x(ix1,:,:,1)/x(ix1-1,:,:,1))**2

        w(ix1,:,:,mom(1)) = max(0.d0,w(ix1,:,:,mom(1)))
        w(ix1,:,:,Fr(1)) = max(0.d0,w(ix1,:,:,Fr(1)))
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
      nearest_i = minloc(abs(r_init - x(ix1,ixOmin2,ixOmin3,1)),DIM=1)

      !> Read closest values
      w(ix1, :, :, rho_) = rho_init(nearest_i)/unit_density
      w(ix1, :, :, mom(1)) = vr_init(nearest_i)/unit_velocity!*w(ix1, rho_)
      w(ix1, :, :, p_) = p_init(nearest_i)/unit_pressure
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

  subroutine stellar_flux(ixI^L,ixO^L,w,x,F0)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: w(ixI^S,1:ndim)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out) :: F0(ixI^S,1:ndir)

    double precision :: gradE(ixI^S), lambda(ixI^S), R(ixI^S)

    integer :: idir

    !> Fixed total luminosity
    ! F0(ixI^S,:) = 0.d0
    !
    ! ! Only in radial direction
    ! F0(ixI^S,1) = lstar/(4*dpi*(x(ixI^S,1)*unit_length)**2)/unit_radflux &
    ! -w(ixI^S,mom(1))*w(ixI^S,Er_)/w(ixI^S,rho_)

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
      F0(ixI^S,idir) = -m1_c*lambda(ixI^S)/(w(ixI^S,rho_)*w(ixI^S,kF_))*gradE(ixI^S)
    end do

  end subroutine stellar_flux

  subroutine stellar_gravity(ixI^L,ixO^L,wCT,x,gravity_field)

    ! Subroutine arguments
    integer, intent(in)  :: ixI^L, ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(out) :: gravity_field(ixI^S,ndim)

    gravity_field(ixI^S,:) = 0.0d0

    ! Only in radial direction
    gravity_field(ixI^S,1) = -Ggrav * mstar/x(ixI^S,1)**2

  end subroutine stellar_gravity

  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    use mod_hd

    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: grav(ixI^S), grad(ixI^S), prad(ixI^S), v(ixI^S), gradv(ixI^S), pth(ixI^S)
    double precision :: T(ixI^S), divT(ixI^S)

    grav(ixI^S) = w(ixI^S,rho_) * Ggrav * mstar/x(ixI^S,1)**2.0d0
    grad(ixI^S) = w(ixI^S,rho_) * w(ixI^S,kF_) * w(ixI^S,Fr(1)) / m1_c

    w(ixI^S,i_Gamma) = grad(ixI^S)/grav(ixI^S)
    w(ixI^S,i_Mf) = w(ixI^S,kF_)/(0.34d0/unit_opacity)

    !> Photon tiring + Radflux CMF correction
    v(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, prad)
    call gradient_str(v, x, ixI^L, ixO^L, 1, gradv)

    w(ixO^S,i_Prr) = prad(ixO^S)
    w(ixO^S,i_pht) = - m1_RSOL*gradv(ixO^S)*prad(ixO^S)
    w(ixO^S,i_fcmf) = - m1_RSOL*w(ixO^S,Fr(1))*gradv(ixO^S)

    call hd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,i_Tgas) = (pth(ixO^S)*unit_pressure*const_mp*hd_mu) &
                   /(w(ixO^S,rho_)*unit_density*const_kB)

    w(ixO^S,i_Trad) = (w(ixO^S,Er_)*unit_pressure/const_rad_a)**0.25d0

    w(ixO^S,i_ff) = sum(w(ixO^S, Fr(:))**2, dim=ndim+1)/(m1_c*w(ixO^S,Er_))**2
    w(ixO^S,i_chi) = (3.d0+4*w(ixO^S,i_ff))/(5.d0+2*dsqrt(4.d0-3*w(ixO^S,i_ff)))
    w(ixO^S,i_lum) = 4*dpi*(x(ixO^S,1)*unit_length)**2 * w(ixO^S,Fr(1))*unit_radflux / lstar

    w(ixO^S,i_mdot) = 4*dpi*w(ixO^S,mom(1))*x(ixO^S,1)**2 * unit_velocity*unit_density*unit_length**2/const_Msun*365*24*60*60
    w(ixO^S,i_am) = x(ixO^S,1)*w(ixO^S,mom(2))

    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, T)
    w(ixO^S,i_Prr) = T(ixO^S)/w(ixO^S,Er_)
    call m1_get_prad(w, x, ixI^L, ixO^L, 2, 2, T)
    w(ixO^S,i_Ptt) = T(ixO^S)/w(ixO^S,Er_)
    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 2, T)
    w(ixO^S,i_Prt) = T(ixO^S)/w(ixO^S,Er_)

    w(ixO^S,i_Ppp) = 1.d0 - w(ixO^S,i_Prr) - w(ixO^S,i_Ptt)

    ! !> Flux equation terms:
    ! !> Directional derivative
    ! w(ixO^S,i_SFdv) = -w(ixO^S,Fr(1))*gradv(ixO^S)
    !
    ! !> Shadowing radiation force source term
    ! w(ixO^S,i_SgR) = -m1_c*w(ixO^S,kF_)*w(ixO^S,rho_)*w(ixO^S,Fr(1))
    !
    ! !> advected flux term
    ! T(ixI^S) = w(ixI^S,Fr(1))*w(ixI^S,mom(1))/w(ixI^S,rho_)
    ! call gradient_str(T, x, ixI^L, ixO^L, 1, divT)
    ! w(ixO^S,i_SdivF) = -divT(ixO^S)
    !
    ! call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, T)
    ! T = m1_c**2 * T
    ! call gradient_str(T, x, ixI^L, ixO^L, 1, divT)
    ! w(ixO^S,i_SdivP) = -divT(ixO^S)
    ! w(ixO^S,i_ScdvP) = (1-T(ixO^S)/w(ixO^S,Er_))*w(ixO^S,Er_)/x(ixO^S,1)

    w(ixO^S,i_tau_cell) = w(ixO^S,rho_)*w(ixO^S,kF_)*m1_c*m1_RSOL

  end subroutine update_extravars

end module mod_usr
