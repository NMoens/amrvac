!> Nicolas Moens
!> Module for including flux limited diffusion in hydrodynamics simulations
!> Based on Turner and stone 2001
module mod_fld
    use mod_multigrid_coupling
    implicit none

    !> source split or not
    logical :: fld_split = .false.

    !> Opacity per unit of unit_density
    double precision, public :: fld_kappa0 = 0.34d0

    double precision, public :: fld_sigma_0

    !> mean particle mass
    double precision, public :: fld_mu = 0.6d0

    !> Maximum amount of pseudotimesteps before trying something else
    integer, public :: fld_maxdw = 100

    !> Tolerance for bisection method for Energy sourceterms
    !> This is a percentage of the minimum of gas- and radiation energy
    double precision, public :: fld_bisect_tol = 1.d-4

    !> Tolerance for adi method for radiative Energy diffusion
    double precision, public :: fld_diff_tol = 1.d-4

    !> Number for splitting the diffusion module
    double precision, public :: diff_crit

    double precision :: fld_max_fracdt = 50.d0

    !> Index for kappa
    integer, public :: i_op

    !> Index for flux limiter
    integer, public :: i_lambda

    !> Index for ratio of scaleheights R
    integer, public :: i_fld_R

    !> Index for testvariable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! DELETE WHEN DONE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, public :: i_test

    !> Index for Flux
    integer, allocatable, public :: i_flux(:)

    !> Indexes for Eddington Tensor
    integer, allocatable, public :: i_edd(:,:)

    !> Use constant Opacity?
    character(len=8) :: fld_opacity_law = 'const'

    !> Diffusion limit lambda = 0.33
    character(len=16) :: fld_fluxlimiter = 'Pomraning'

    !> diffusion coefficient for multigrid method
    integer :: i_diff_mg

    !> Which method to solve diffusion part
    character(len=8) :: fld_diff_scheme = 'adi'

    !> Which method to find the root for the energy interaction polynomial
    character(len=8) :: fld_interaction_method = 'Bisect'

    !> Set Diffusion coefficient to unity
    logical :: fld_diff_testcase = .false.

    !> Take running average for Diffusion coefficient
    logical :: diff_coef_filter = .false.
    integer :: size_D_filter = 1

    !> Use or don't use lineforce opacities
    logical :: Lineforce_opacities = .false.

    !> Index for Flux weighted opacities
    integer, allocatable, public :: i_opf(:)

    !> public methods
    !> these are called in mod_rhd_phys
    public :: get_fld_rad_force
    public :: get_fld_energy_interact
    public :: get_fld_diffusion
    public :: fld_init
    public :: fld_get_radflux
    public :: fld_get_radpress
    public :: fld_get_fluxlimiter
    public :: fld_get_opacity
    public :: get_rad_extravars
    public :: set_mg_bounds

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! GENERAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Reading in fld-list parameters from .par file
  subroutine fld_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /fld_list/ fld_kappa0, fld_split, fld_maxdw, &
    fld_bisect_tol, fld_diff_testcase, fld_diff_tol, fld_max_fracdt,&
    fld_opacity_law, fld_fluxlimiter, fld_diff_scheme, fld_interaction_method, &
    diff_coef_filter, size_D_filter, lineforce_opacities

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, fld_list, end=111)
       111    close(unitpar)
    end do
  end subroutine fld_params_read

  !> Initialising FLD-module:
  !> Read opacities
  !> Initialise Multigrid
  !> adimensionalise kappa
  !> Add extra variables to w-array, flux, kappa, eddington Tensor
  !> Lambda and R
  !> ...
  subroutine fld_init(He_abundance, rhd_radiation_diffusion)
    use mod_global_parameters
    use mod_variables
    use mod_physics, only: phys_global_source
    use mod_opacity, only: init_opal
    use mod_multigrid_coupling, only: mg_copy_boundary_conditions

    double precision, intent(in) :: He_abundance
    logical, intent(in) :: rhd_radiation_diffusion
    double precision :: sigma_thomson
    integer :: idir,jdir

    character(len=1) :: ind_1
    character(len=1) :: ind_2
    character(len=2) :: cmp_f
    character(len=5) :: cmp_e

    !> read par files
    call fld_params_read(par_files)

    !> Set radiative flux as variable
    allocate(i_flux(ndir))
    do idir = 1,ndir
      write(ind_1,'(I1)') idir
      cmp_f = 'F' // ind_1
      i_flux(idir) = var_set_extravar(cmp_f,cmp_f)
    enddo


    !> Set lineforce opacities as variable
    if (lineforce_opacities) then
      allocate(i_opf(ndir))
      do idir = 1,ndir
        write(ind_1,'(I1)') idir
        cmp_f = 'k' // ind_1
        i_opf(idir) = var_set_extravar(cmp_f,cmp_f)
      enddo
    endif

    !> Introduce opacity, lambda and R as global variables
    i_op = var_set_extravar("Kappa", "Kappa")
    i_lambda = var_set_extravar("lambda", "lambda")
    i_fld_R = var_set_extravar("fld_R", "fld_R")

    !> Introduce test variable globally
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! DELETE WHEN DONE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    i_test = var_set_extravar('test','test')

    allocate(i_edd(ndir,ndir))
    do idir = 1, ndir
      do jdir = 1, ndir
        write(ind_1,'(I1)') idir
        write(ind_2,'(I1)') jdir
        cmp_e = 'Edd' // ind_1 // ind_2
        i_edd(idir,jdir) = var_set_extravar(cmp_e, cmp_e)
      enddo
    enddo

    if (rhd_radiation_diffusion) then
      if (fld_diff_scheme .eq. 'mg') then

        use_multigrid = .true.

        if (rhd_radiation_diffusion) then
          phys_global_source => Diffuse_E_rad_mg
        endif

        mg_after_new_tree => set_mg_diffcoef

        mg%n_extra_vars = 1
        mg%operator_type = mg_vhelmholtz

        ! i_diff_mg = var_set_extravar("D", "D")
      endif
    endif
    i_diff_mg = var_set_extravar("D", "D")


    !> Check if fld_numdt is not 1
    if (fld_maxdw .lt. 2) call mpistop("fld_maxdw should be an integer larger than 1")

    !> Need mean molecular weight
    fld_mu = (1.+4*He_abundance)/(2.+3.*He_abundance)

    !> Dimensionless Boltzman constante sigma
    fld_sigma_0 = const_sigma*(unit_temperature**4.d0)/(unit_velocity*unit_pressure)

    !> Read in opacity table if necesary
    if (fld_opacity_law .eq. 'opal') call init_opal(He_abundance)
    if ((fld_opacity_law .eq. 'thomson') .or. (fld_opacity_law .eq. 'fastwind'))  then
      sigma_thomson = 6.6524585d-25
      fld_kappa0 = sigma_thomson/const_mp * (1.+2.*He_abundance)/(1.+4.*He_abundance)
    endif
  end subroutine fld_init

  !> Compute all extra variables in w-array:
  !> Flux, Eddington tensor, lambda, R, kappa
  subroutine get_rad_extravars(w, x, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    call fld_get_opacity(w, x, ixI^L, ixO^L)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L)
    call fld_get_radflux(w, x, ixI^L, ixO^L)
    call fld_get_eddington(w, x, ixI^L, ixO^L)

    if (fld_diff_scheme .eq. 'mg') then
      call fld_get_diffcoef_central(w, x, ixI^L, ixO^L)
      call set_mg_bounds(w, x, ixI^L, ixO^L)
      call get_diffusion_criterion(w, x, ixI^L, ixO^L)
    endif
  end subroutine get_rad_extravars

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  !> This subroutine handles the radiation force
  subroutine get_fld_rad_force(qdt,ixI^L,ixO^L,wCT,w,x,&
       energy,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods

    use mod_physics, only: phys_get_pthermal  !needed to get temp

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision                :: wCCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active
    double precision :: radiation_force(ixO^S,1:ndim)

    integer :: idir, i, jx^L

    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_split) then
      active = .true.

      wCCT = wCT
      call fld_get_radflux(wCCT, x, ixI^L, ixO^L)

      do idir = 1,ndir
        !> Radiation force = kappa*rho/c *Flux
        radiation_force(ixO^S,idir) = wCT(ixO^S,iw_rho)*wCT(ixO^S,i_op)*wCCT(ixO^S, i_flux(idir))/(const_c/unit_velocity)

        !> Momentum equation source term
        w(ixO^S,iw_mom(idir)) = w(ixO^S,iw_mom(idir)) &
            + qdt * radiation_force(ixO^S,idir)

        ! print*, it, 'Not Adding F_rad to kinetic energy'
        if (.not. block%e_is_internal) then
          !> Energy equation source term (kinetic energy)
          w(ixO^S,iw_e) = w(ixO^S,iw_e) &
              + qdt * radiation_force(ixO^S,idir) * wCT(ixO^S,iw_mom(idir))/wCT(ixO^S,iw_rho)
        endif
      enddo

    end if

  end subroutine get_fld_rad_force

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  !> This subroutine handles the energy exchange between gas and radiation
  subroutine get_fld_energy_interact(qdt,ixI^L,ixO^L,wCT,w,x,&
       energy,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods

    use mod_physics, only: phys_get_pthermal  !needed to get temp

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active

    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_split) then
      active = .true.
      !> Add energy sourceterms
      call Energy_interaction(w, x, ixI^L, ixO^L)
    end if
  end subroutine get_fld_energy_interact

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  !> This subroutine handles the diffusion of the radiation energy density,
  !> calling either a multigrid-method or an ADI-scheme (perhaps outdated? Need to check).
  !> To be added: 1D backward euler
  subroutine get_fld_diffusion(qdt,ixI^L,ixO^L,wCT,w,x,&
       energy,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods
    use mod_physics
    use mod_multigrid_coupling
    use m_diffusion
    use mpi

    use mod_physics, only: phys_get_pthermal  !needed to get temp

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: D_center(ixI^S)
    logical, intent(in) :: energy,qsourcesplit
    logical, intent(inout) :: active


    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_split) then
      active = .true.
      !> Begin by evolving the radiation energy field
      select case (fld_diff_scheme)
      case('adi')
        call Evolve_E_rad(w, x, ixI^L, ixO^L)
      case('mg')
        call fld_get_diffcoef_central(w, x, ixI^L, ixO^L)
        call set_mg_bounds(w, x, ixI^L, ixO^L)

        active = .true.

      case default
        call mpistop('Numerical diffusionscheme unknown, try adi or mg')
      end select
      end if
  end subroutine get_fld_diffusion

  !> Sets the opacity in the w-array
  !> by calling mod_opacity
  subroutine fld_get_opacity(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal
    use mod_physics, only: phys_get_tgas
    use mod_usr_methods
    use mod_opacity

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision :: fld_kappa(ixO^S)
    double precision :: Temp(ixI^S), pth(ixI^S), a2(ixO^S)
    double precision :: rho0,Temp0,n,sigma_b
    double precision :: akram, bkram
    double precision :: vth(ixO^S), gradv(ixI^S), eta(ixO^S), t(ixO^S)

    integer :: i,j,ix^D, idir

    select case (fld_opacity_law)
      case('const')
        fld_kappa = fld_kappa0/unit_opacity
      case('thomson')
        fld_kappa = fld_kappa0/unit_opacity
      case('kramers')
        rho0 = half !> Take lower value of rho in domain
        fld_kappa(ixO^S) = fld_kappa0/unit_opacity*((w(ixO^S,iw_rho)/rho0))
      case('bump')
        !> Opacity bump
        rho0 = 0.2d0 !0.5d-1
        n = 7.d0
        sigma_b = 2.d-2
        !fld_kappa(ixO^S) = fld_kappa0/unit_opacity*(one + n*dexp(-((rho0  - w(ixO^S,iw_rho))**two)/rho0))
        fld_kappa(ixO^S) = fld_kappa0/unit_opacity*(one + n*dexp(-one/sigma_b*(dlog(w(ixO^S,iw_rho)/rho0))**two))
      case('non_iso')
        call phys_get_pthermal(w,x,ixI^L,ixO^L,Temp)
        Temp(ixO^S)=Temp(ixO^S)/w(ixO^S,iw_rho)

        rho0 = 0.5d0 !> Take lower value of rho in domain
        Temp0 = one
        n = -7.d0/two
        fld_kappa(ixO^S) = fld_kappa0/unit_opacity*(w(ixO^S,iw_rho)/rho0)*(Temp(ixO^S)/Temp0)**n
      case('fastwind')
        call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)
        a2(ixO^S) = pth(ixO^S)/w(ixO^S,iw_rho)*unit_velocity**2.d0

        akram = 13.1351597305
        bkram = -4.5182188206

        fld_kappa(ixO^S) = fld_kappa0/unit_opacity &
        * (1.d0+10.d0**akram*w(ixO^S,iw_rho)*unit_density*(a2(ixO^S)/1.d12)**bkram)

        {do ix^D=ixOmin^D,ixOmax^D\ }
          !> Hard limit on kappa
          fld_kappa(ix^D) = min(fld_kappa(ix^D),2.3d0*fld_kappa0/unit_opacity)

          !> Limit kappa through T
          ! fld_kappa(ix^D) = fld_kappa0/unit_opacity &
          ! * (1.d0+10.d0**akram*w(ix^D,iw_rho)*unit_density &
          ! * (max(a2(ix^D),const_kB*5.9d4/(fld_mu*const_mp))/1.d12)**bkram)
        {enddo\ }

      case('opal')
        call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)
        {do ix^D=ixOmin^D,ixOmax^D\ }
            rho0 = w(ix^D,iw_rho)*unit_density
            Temp0 = Temp(ix^D)*unit_temperature
            call set_opal_opacity(rho0,Temp0,n)
            fld_kappa(ix^D) = n/unit_opacity
        {enddo\ }

      case('special')
        if (.not. associated(usr_special_opacity)) then
          call mpistop("special opacity not defined")
        endif
        call usr_special_opacity(ixI^L, ixO^L, w, x, fld_kappa)

      case default
        call mpistop("Doesn't know opacity law")
      end select

      w(ixO^S, i_op) = fld_kappa(ixO^S)
  end subroutine fld_get_opacity

  !> Set lineforce opacities
  subroutine fld_get_lineopacity(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    double precision :: vel(ixI^S), gradv(ixI^S), forceM(ixI^S)
    integer :: ix^D, idir

    double precision, parameter :: Qbar = 2000.d0
    double precision, parameter :: alpha = 0.6d0

    !> Set lineforce opacities
    !> Stevens & Kallman 1990
    !> Will this work on initialisation? F depends on kappa, kappa_f depends on F
    if (fld_opacity_law .ne. 'thomson') &
      call mpistop('When using line-opacities, you should use a thomson opacity law')

    if (lineforce_opacities) then

      !> Set t
      do idir = 1,ndir
        vel(ixI^S) = w(ixI^S,iw_mom(idir))/w(ixI^S,iw_rho)
        call gradient(vel,ixI^L,ixO^L,idir,gradv)
        forceM(ixI^S) = Qbar/(one-alpha) &
        *(gradv(ixI^S)/(w(ixI^S,iw_rho)*(const_c/unit_velocity)*Qbar*w(ixI^S,i_op)))**alpha
        w(ixO^S,i_opf(idir)) = w(ixO^S,i_op)*forceM(ixO^S)
      enddo
    else
      call mpistop("Lineforce opacities are not calculated")
    endif

  end subroutine fld_get_lineopacity

  !> Calculate fld flux limiter
  !> This subroutine calculates flux limiter lambda using the prescription
  !> stored in fld_fluxlimiter.
  !> It also calculates the ratio of radiation scaleheight and mean free path
  subroutine fld_get_fluxlimiter(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision :: fld_R(ixI^S), fld_lambda(ixI^S)
    double precision ::  normgrad2(ixI^S)
    double precision :: grad_r_e(ixI^S), rad_e(ixI^S)
    integer :: idir, i, j

    select case (fld_fluxlimiter)
    case('Diffusion')
      w(ixI^S,i_lambda) = one/3.d0
      w(ixI^S,i_fld_R) = zero

    case('FreeStream')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixI^S) = zero
      rad_e(ixI^S) = w(ixI^S, iw_r_e)
      do idir = 1,ndir
        !> gradient or gradientS ?!?!?!?!?!?
        call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        normgrad2(ixI^S) = normgrad2(ixI^S) + grad_r_e(ixI^S)**2
      end do

      fld_R(ixI^S) = dsqrt(normgrad2(ixI^S))/(w(ixI^S,i_op)*w(ixI^S,iw_rho)*w(ixI^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      fld_lambda(ixI^S) = one/fld_R(ixI^S)

      w(ixI^S,i_lambda) = fld_lambda(ixI^S)
      w(ixI^S,i_fld_R) = fld_R(ixI^S)

    case('Pomraning')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixI^S) = zero
      rad_e(ixI^S) = w(ixI^S, iw_r_e)
      do idir = 1,ndir
        !> gradient or gradientS ?!?!?!?!?!?
        call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        normgrad2(ixI^S) = normgrad2(ixI^S) + grad_r_e(ixI^S)**2
      end do

      fld_R(ixI^S) = dsqrt(normgrad2(ixI^S))/(w(ixI^S,i_op)*w(ixI^S,iw_rho)*w(ixI^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
      fld_lambda(ixI^S) = (2.d0+fld_R(ixI^S))/(6.d0+3*fld_R(ixI^S)+fld_R(ixI^S)**2.d0)

      w(ixI^S,i_lambda) = fld_lambda(ixI^S)
      w(ixI^S,i_fld_R) = fld_R(ixI^S)

    case('Pomraning2')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixI^S) = zero
      rad_e(ixI^S) = w(ixI^S, iw_r_e)
      do idir = 1,ndir
        !> gradient or gradientS ?!?!?!?!?!?
        call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        normgrad2(ixI^S) = normgrad2(ixI^S) + grad_r_e(ixI^S)**2
      end do

      fld_R(ixI^S) = dsqrt(normgrad2(ixI^S))/(w(ixI^S,i_op)*w(ixI^S,iw_rho)*w(ixI^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Levermore and Pomraning: lambda = 1/R(coth(R)-1/R)
      fld_lambda(ixI^S) = one/fld_R(ixI^S)*(one/dtanh(fld_R(ixI^S)) - one/fld_R(ixI^S))

      w(ixI^S,i_lambda) = fld_lambda(ixI^S)
      w(ixI^S,i_fld_R) = fld_R(ixI^S)


    case('Minerbo')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixI^S) = zero
      rad_e(ixI^S) = w(ixI^S, iw_r_e)
      do idir = 1,ndir
        !> gradient or gradientS ?!?!?!?!?!?
        call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        normgrad2(ixI^S) = normgrad2(ixI^S) + grad_r_e(ixI^S)**2
      end do

      fld_R(ixI^S) = dsqrt(normgrad2(ixI^S))/(w(ixI^S,i_op)*w(ixI^S,iw_rho)*w(ixI^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Minerbo:
      do i = ixImin1, ixImax1
        do j = ixImin2, ixImax2
          if (fld_R(i,j) .lt. 3.d0/2.d0) then
            fld_lambda(i,j) = 2.d0/(3.d0 + dsqrt(9.d0 + 12.d0*fld_R(i,j)**2.d0))
          else
            fld_lambda(i,j) = 1.d0/(1.d0 + fld_R(i,j) + dsqrt(1.d0 + 2.d0*fld_R(i,j)))
          endif
        enddo
      enddo

      w(ixI^S,i_lambda) = fld_lambda(ixI^S)
      w(ixI^S,i_fld_R) = fld_R(ixI^S)
    case('special')
      if (.not. associated(usr_special_fluxlimiter)) then
        call mpistop("special fluxlimiter not defined")
      endif
      call usr_special_fluxlimiter(ixI^L, ixO^L, w, x, fld_lambda, fld_R)
      w(ixI^S,i_lambda) = fld_lambda(ixI^S)
      w(ixI^S,i_fld_R) = fld_R(ixI^S)
    case default
      call mpistop('Fluxlimiter unknown')
    end select
  end subroutine fld_get_fluxlimiter

  !> Calculate Radiation Flux
  !> stores radiation flux in w-array
  subroutine fld_get_radflux(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision :: rad_flux(ixI^S, 1:ndim)
    double precision :: L_star, R_star
    double precision :: grad_r_e(ixI^S)
    double precision :: rad_e(ixI^S)
    integer :: ix^D, idir

    rad_e(ixI^S) = w(ixI^S, iw_r_e)

    !> Calculate the Flux using the fld closure relation
    !> F = -c*lambda/(kappa*rho) *grad E
    do idir = 1,ndir
      !> gradient or gradientS ?!?!?!?!?!?
      call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
      rad_flux(ixI^S, idir) = -(const_c/unit_velocity)*w(ixI^S,i_lambda)/(w(ixI^S,i_op)*w(ixI^S,iw_rho))*grad_r_e(ixI^S)
    end do

    w(ixI^S,i_flux(:)) = rad_flux(ixI^S,:)

    ! !>Cheaty bit:
    !   w(:,ixOmin2,i_flux(2)) = (x(:,ixOmin2+1,2)/x(:,ixOmin2,2))**2*w(:,ixOmin2+1,i_flux(2))
    !
    !   w(:,ixOmax2,i_flux(2)) = (x(:,ixOmax2,2)/x(:,ixOmax2-1,2))**2*w(:,ixOmax2-1,i_flux(2))

  end subroutine fld_get_radflux

  !> Calculate Eddington-tensor
  !> Stores Eddington-tensor in w-array
  subroutine fld_get_eddington(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision :: eddington_tensor(ixO^S,1:ndim,1:ndim)
    double precision :: tnsr2(ixO^S,1:ndim,1:ndim)
    double precision :: normgrad2(ixO^S), f(ixO^S)
    double precision :: grad_r_e(ixI^S, 1:ndim), rad_e(ixI^S)
    integer :: i,j, idir,jdir

    !> Calculate R everywhere
    !> |grad E|/(rho kappa E)
    normgrad2(ixO^S) = zero

    rad_e(ixI^S) = w(ixI^S, iw_r_e)
    grad_r_e(ixI^S,:) = zero
    do idir = 1,ndir
      !> gradient or gradientS ?!?!?!?!?!?
      call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e(ixI^S,idir))
      normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S,idir)**two
    end do

    !> Calculate radiation pressure
    !> P = (lambda + lambda^2 R^2)*E
    f(ixO^S) = w(ixO^S,i_lambda) + w(ixO^S, i_lambda)**two * w(ixO^S, i_fld_R)**two
    f(ixO^S) = one/two*(one-f(ixO^S)) + one/two*(3.d0*f(ixO^S) - one)

    do idir = 1,ndir
      eddington_tensor(ixO^S,idir,idir) = half*(one-f(ixO^S))
    enddo

    do idir = 1,ndir
      do jdir = 1,ndir
        if (idir .ne. jdir) eddington_tensor(ixO^S,idir,jdir) = zero
        tnsr2(ixO^S,idir,jdir) =  half*(3.d0*f(ixO^S) - 1)&
        *grad_r_e(ixO^S,idir)*grad_r_e(ixO^S,jdir)/normgrad2(ixO^S)
      enddo
    enddo

    do idir = 1,ndir
      do jdir = 1,ndir
        where ((tnsr2(ixO^S,idir,jdir) .eq. tnsr2(ixO^S,idir,jdir)) &
          .and. (normgrad2(ixO^S) .gt. smalldouble))
          eddington_tensor(ixO^S,idir,jdir) = eddington_tensor(ixO^S,idir,jdir) + tnsr2(ixO^S,idir,jdir)
        endwhere
      enddo
    enddo

    do idir = 1,ndir
      do jdir = 1,ndir
        w(ixO^S,i_edd(idir,jdir)) = eddington_tensor(ixO^S,idir,jdir)
      enddo
    enddo
  end subroutine fld_get_eddington

  !> Calculate Radiation Pressure
  !> Returns Radiation Pressure as tensor
  subroutine fld_get_radpress(w, x, ixI^L, ixO^L, rad_pressure)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision             :: eddington_tensor(ixO^S,1:ndim,1:ndim)
    double precision, intent(out):: rad_pressure(ixO^S,1:ndim,1:ndim)

    integer i,j

    do i=1,ndim
      do j=1,ndim
        rad_pressure(ixO^S,i,j) = w(ixO^S,i_edd(i,j))* w(ixO^S,iw_r_e)
      enddo
    enddo
  end subroutine fld_get_radpress

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! Multigrid diffusion
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calling all subroutines to perform the multigrid method
  !> Communicates rad_e and diff_coeff to multigrid library
  subroutine Diffuse_E_rad_mg(qdt, qt, active)
    use mod_global_parameters
    use mod_multigrid_coupling
    use m_diffusion

    double precision, intent(in) :: qdt, qt
    logical, intent(inout)       :: active
    double precision             :: max_res

    ! integer :: i, itdiff, Ndiff
    !
    ! call mg_copy_to_tree(iw_r_e, mg_iphi, .false., .false.)
    !
    ! if (diff_crit .lt. one) then
    !   call diffusion_solve_vcoeff(mg, qdt, 2, fld_diff_tol)
    ! else
    !   Ndiff = int(diff_crit)
    !   do itdiff = 1,Ndiff
    !     call diffusion_solve_vcoeff(mg, qdt/Ndiff, 2, fld_diff_tol)
    !   enddo
    ! endif
    !
    ! call mg_copy_from_tree(mg_iphi, iw_r_e)
    ! active = .true.

    print*, 'beginning diffusion', mype

    call mg_copy_to_tree(iw_r_e, mg_iphi, .false., .false.)
    call diffusion_solve_vcoeff(mg, qdt, 2, fld_diff_tol)
    call mg_copy_from_tree(mg_iphi, iw_r_e)
    active = .true.

    print*, 'ending diffusion', mype


  end subroutine Diffuse_E_rad_mg

  !> Calculates cell-centered diffusion coefficient to be used in multigrid
  subroutine fld_get_diffcoef_central(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    double precision :: max_D(ixI^S), grad_r_e(ixI^S), rad_e(ixI^S)
    integer :: idir,i,j, ix^D

    if (fld_diff_testcase) then

      w(ixO^S,i_diff_mg) = 1.d0

      ! w(ixO^S,i_diff_mg) = 1.d0/(1.d-3*(const_c/unit_velocity)*0.33333333d0/(w(ixO^S,i_op)*w(ixO^S,iw_rho)))
      !
      ! w(ixO^S,i_test) = (const_c/unit_velocity)*0.33333333d0/(w(ixO^S,i_op)*w(ixO^S,iw_rho))

    else
      !> calculate diffusion coefficient
      w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*w(ixO^S,i_lambda)/(w(ixO^S,i_op)*w(ixO^S,iw_rho))

      !> ghostcells
      do i = ixImin1,ixImax1
        w(i,ixImin2:ixImin2+nghostcells-1,i_diff_mg) = w(i,ixImin2+nghostcells,i_diff_mg)
        w(i,ixImax2-nghostcells+1:ixImin2,i_diff_mg) = w(i,ixImax2-nghostcells,i_diff_mg)
      enddo
      do i = ixImin2,ixImax2
        w(ixImin1:ixImin1+nghostcells-1,i,i_diff_mg) = w(ixImin1+nghostcells,i,i_diff_mg)
        w(ixImax1-nghostcells+1:ixImin1,i,i_diff_mg) = w(ixImax1-nghostcells,i,i_diff_mg)
      enddo

      ! !> Check if energy doesn't go faster than speed of light
      ! !for simplicity, only in direction 2
      ! rad_e(ixI^S) = w(ixI^S,iw_r_e)
      ! call gradient(rad_e,ixI^L,ixO^L,2,grad_r_e)
      ! max_D(ixO^S) = abs((const_c/unit_velocity)*rad_e(ixO^S)/grad_r_e(ixO^S))
      !
      ! {do ix^D = ixOmin^D,ixOmax^D\}
      !     ! call mpistop('You have reached maximal D, the D is too big')
      !     w(i,j,i_diff_mg) = min(w(i,j,i_diff_mg), max_D(i,j))
      ! {enddo\}

      if (diff_coef_filter) then
        !call mpistop('Hold your bloody horses, not implemented yet ')
        call fld_smooth_diffcoef(w, ixI^L, ixO^L)
      endif
    endif


    !> CHEATY
    ! w(:,ixOmax2,i_diff_mg) = w(:,ixOmax2-2,i_diff_mg)
    ! w(:,ixOmax2-1,i_diff_mg) = w(:,ixOmax2-2,i_diff_mg)

  end subroutine fld_get_diffcoef_central

  !> Use running average on Diffusion coefficient
  subroutine fld_smooth_diffcoef(w, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)

    double precision :: tmp_D(ixI^S), filtered_D(ixI^S)
    integer :: ix^D, filter, idim

    if (size_D_filter .lt. 1) call mpistop("D filter of size < 1 makes no sense")
    if (size_D_filter .gt. nghostcells) call mpistop("D filter of size > nghostcells makes no sense")

    tmp_D(ixI^S) = w(ixI^S,i_diff_mg)
    filtered_D(ixI^S) = zero

    do filter = 1,size_D_filter
      {do ix^D = ixOmin^D+size_D_filter,ixOmax^D-size_D_filter\}
        do idim = 1,ndim
          filtered_D(ix^D) = filtered_D(ix^D) &
                           + tmp_D(ix^D+filter*kr(idim,^D)) &
                           + tmp_D(ix^D-filter*kr(idim,^D))
        enddo
      {enddo\}
    enddo

    {do ix^D = ixOmin^D+size_D_filter,ixOmax^D-size_D_filter\}
      tmp_D(ix^D) = (tmp_D(ix^D)+filtered_D(ix^D))/(1+2*size_D_filter*ndim)
    {enddo\}

    w(ixO^S,i_diff_mg) = tmp_D(ixO^S)
  end subroutine fld_smooth_diffcoef


  !> Communicates diffusion coeff to multigrid library
  subroutine set_mg_diffcoef()
    call mg_copy_to_tree(i_diff_mg, mg_iveps, .true., .true.)
  end subroutine set_mg_diffcoef

  !> Sets boundary conditions for multigrid, based on hydro-bounds
  subroutine set_mg_bounds(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_usr_methods

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    integer :: iB

    do iB = 1,2*ndim
      select case (typeboundary(iw_r_e, iB))
      case ('symm')
         mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
         mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
      case ('asymm')
         mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
         mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
      case ('cont')
         mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
         mg%bc(iB, mg_iphi)%bc_value = 0.0_dp ! Not needed
      case ('periodic')
        !> Do nothing
      case ('special')

        if (.not. associated(usr_special_mg_bc)) call mpistop("Set special mg bound")
        call usr_special_mg_bc(global_time,ixI^L,ixO^L,iB,w,x)

      case default
         print *, "Not a standard: ", trim(typeboundary(iw_r_e, iB))
         error stop "You have to set a user-defined boundary method"
      end select
    enddo
  end subroutine set_mg_bounds

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! ADI
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Main loop for ADI scheme.
  !> Here, the splitting of the hydro-timestep, nr of pseudo-steps
  !> and error-controll are performed
  subroutine Evolve_E_rad(w, x, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision :: E_new(ixI^S), E_old(ixI^S), ADI_Error
    double precision :: frac_grid
    integer :: w_max, frac_dt
    logical :: converged

    integer :: i

    E_new(ixI^S) = w(ixI^S,iw_r_e)

    converged = .false.
    ADI_Error = bigdouble

    ! w_max = 1
    !> Trying out something new
    ! This should make sure that the amount of pseudotimesteps
    ! Goes down with one after completing a hydro step.
    w_max = max(1,w_max/2)

    frac_grid = two
    frac_dt = 1

    do while (converged .eqv. .false.)

      !> Check if solution converged
      if (ADI_Error .lt. fld_diff_tol) then
        !> If converged in former loop, break loop
        converged = .true.
      else
        !> Reset E_new
        E_old(ixI^S) = w(ixI^S,iw_r_e)
        E_new(ixI^S) = w(ixI^S,iw_r_e)

        !> If no convergence, adapt pseudostepping
        w_max = 2*w_max
        frac_grid = 2*frac_grid
      endif

      !> Evolve using ADI
      if (converged .eqv. .false.) then

        call Evolve_ADI(w, x, E_new, E_old, w_max, frac_grid, ixI^L, ixO^L)
        call Error_check_ADI(w, x, E_new, E_old, ixI^L, ixO^L, ADI_Error) !> SHOULD THIS BE DONE EVERY ITERATION???
        if (ADI_Error .lt. fld_diff_tol) then
          converged = .true.
        endif
      endif

      !> If adjusting pseudostep doesn't work, divide the actual timestep in smaller parts
      if (w_max .gt. fld_maxdw) then
        if (converged .eqv. .false.) then
          !> use a smaller timestep than the hydrodynamical one
          call half_timestep_ADI(w, x, E_new, E_old, ixI^L, ixO^L, converged)
          call Error_check_ADI(w, x, E_new, E_old, ixI^L, ixO^L, ADI_Error)
          if (ADI_Error .lt. fld_diff_tol) then
            converged = .true.
          endif
        endif
      endif
    enddo

    w(ixO^S,iw_r_e) = E_new(ixO^S)
  end subroutine Evolve_E_rad

  !> Perform ADI on half a hydro-timestep
  subroutine half_timestep_ADI(w, x, E_new, E_old, ixI^L, ixO^L, converged)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(in) :: E_old(ixI^S)
    double precision, intent(out) :: E_new(ixI^S)
    logical, intent(inout) :: converged
    double precision :: frac_grid
    double precision :: E_loc(ixI^S)
    double precision :: saved_dt, ADI_Error
    integer :: i,  w_max, frac_dt

    saved_dt = dt
    ADI_Error = bigdouble
    frac_dt = 1
    5231 frac_dt = 2*frac_dt
    w_max = 1
    frac_grid = two

    if (frac_dt .gt. fld_max_fracdt) call mpistop("No convergence after halving timestep N times")
    dt = dt/frac_dt

    E_loc = E_old

    do i = 1,frac_dt
      !---------------------------------------------------------------
      do while (converged .eqv. .false.)
        !> Check if solution converged
        if (ADI_Error .lt. fld_diff_tol) then
          !> If converged in former loop, break loop
          converged = .true.
          goto 7895
        else
          !> If no convergence, adapt pseudostepping
          w_max = 2*w_max
          frac_grid = 2*frac_grid
        endif

        !> Evolve using ADI
        call Evolve_ADI(w, x, E_new, E_loc, w_max, frac_grid, ixI^L, ixO^L)
        call Error_check_ADI(w, x, E_new, E_loc, ixI^L, ixO^L, ADI_Error) !> SHOULD THIS BE DONE EVERY ITERATION???

        !> If adjusting pseudostep doesn't work, divide the actual timestep in smaller parts
        if (w_max .gt. fld_maxdw) goto 5231

        if (ADI_Error .lt. fld_diff_tol) then
          converged = .true.
        endif

      enddo
      !---------------------------------------------------------------
      7895 E_loc = E_new
    enddo

    dt = saved_dt
  end subroutine half_timestep_ADI

  !> Calculate error after one ADI-timestep
  subroutine Error_check_ADI(w, x, E_new, E_old, ixI^L, ixO^L, ADI_Error)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S, 1:ndim), w(ixI^S, 1:nw)
    double precision, intent(in) :: E_new(ixI^S), E_old(ixI^S)
    double precision, intent(out) :: ADI_Error
    double precision :: LHS(ixO^S), RHS(ixO^S), D(ixI^S,1:ndim)
    integer :: jx1^L, hx1^L,jx2^L, hx2^L

    integer :: i

    jx1^L=ixO^L+kr(1,^D);
    hx1^L=ixO^L-kr(1,^D);
    jx2^L=ixO^L+kr(2,^D);
    hx2^L=ixO^L-kr(2,^D);

    call fld_get_diffcoef(w, x, ixI^L, ixO^L, D)

    !> LHS = dx^2/dt * (E_new - E_old)
    LHS(ixO^S) = (x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))*&
    (x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))/dt*&
    (E_new(ixO^S) - E_old(ixO^S))

    !> RHS = D1(E_+ - E) - D1(E - E_-) + D2(E_+ - E) - D2(E - E_-)
    RHS(ixO^S) = &
      D(jx1^S,1)*(E_new(jx1^S) - E_new(ixO^S)) &
    - D(ixO^S,1)*(E_new(ixO^S) - E_new(hx1^S)) &
    + D(jx2^S,2)*(E_new(jx2^S) - E_new(ixO^S)) &
    - D(ixO^S,2)*(E_new(ixO^S) - E_new(hx2^S))

    ADI_Error = maxval(abs((RHS-LHS)/(E_old/dt))) !> Try mean value or smtn
    !ADI_Error = sum(abs((RHS-LHS)/(E_old/dt)))/((ixOmax1-ixOmin1)*(ixOmax2-ixOmin2))
  end subroutine Error_check_ADI

  !> Do all pseudo-timesteps to advance one hydro timestep.
  !> This routine loops over the pseudo-steps, each time doing 2 matrix-inversions
  subroutine Evolve_ADI(w, x, E_new, E_old, w_max, frac_grid, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, w_max
    double precision, intent(in) :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim), frac_grid
    double precision, intent(in) :: E_old(ixI^S)
    double precision, intent(out):: E_new(ixI^S)
    double precision :: E_m(ixI^S), E_n(ixI^S)
    double precision :: diag1(ixImax1,ixImax2),sub1(ixImax1,ixImax2),sup1(ixImax1,ixImax2),bvec1(ixImax1,ixImax2)
    double precision :: diag2(ixImax2,ixImax1),sub2(ixImax2,ixImax1),sup2(ixImax2,ixImax1),bvec2(ixImax2,ixImax1)
    double precision :: Evec1(ixImin1:ixImax1), Evec2(ixImin2:ixImax2)
    double precision :: dw, w0, w1
    integer :: m, j, i

    w0 = (x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))&
    *(x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))/frac_grid
    w1 = (x(ixOmax1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))&
    *(x(ixOmin1,ixOmax2,2)-x(ixOmin1,ixOmin2,2))/frac_grid !4.d0

    E_m = E_old

    do m = 1,w_max
      E_n = E_old

      !> Set pseudotimestep
      dw = w0*(w1/w0)**((m-one)/(w_max-one))

      call ADI_boundary_conditions(ixI^L,ixO^L,E_m,w,x)

      !> Setup matrix and vector for sweeping in direction 1
      call make_matrix(x,w,dw,E_m,E_n,1,ixImax1,ixI^L, ixO^L,diag1,sub1,sup1,bvec1,diag2,sub2,sup2,bvec2)
      do j = ixImin2,ixImax2
        Evec1(ixImin1:ixImax1) = E_m(ixImin1:ixImax1,j)
        call solve_tridiag(ixOmin1,ixOmax1,ixImin1,ixImax1,diag1(:,j),sub1(:,j),sup1(:,j),bvec1(:,j),Evec1)
        !E_m(ixOmin1:ixOmax1,j) = Evec1(ixOmin1:ixOmax1)
        E_m(iximin1:ixImax1,j) = Evec1(ixImin1:ixImax1)
      enddo

      call ADI_boundary_conditions(ixI^L,ixO^L,E_m,w,x)

      !> Setup matrix and vector for sweeping in direction 2
      call make_matrix(x,w,dw,E_m,E_n,2,ixImax2,ixI^L,ixO^L,diag1,sub1,sup1,bvec1,diag2,sub2,sup2,bvec2)
      do j = ixImin1,ixImax1
        Evec2(ixImin2:ixImax2) = E_m(j,ixImin2:ixImax2)
        call solve_tridiag(ixOmin2,ixOmax2,ixImin2,ixImax2,diag2(:,j),sub2(:,j),sup2(:,j),bvec2(:,j),Evec2)
        !E_m(j,ixOmin2:ixOmax2) = Evec2(ixOmin2:ixOmax2)
        E_m(j,ixImin2:ixImax2) = Evec2(ixImin2:ixImax2)
      enddo

      call ADI_boundary_conditions(ixI^L,ixO^L,E_m,w,x)

    enddo
    E_new = E_m
  end subroutine Evolve_ADI

  !> Calculate cell-faced diffusion coefficient out 6 neighbouring cells
  subroutine fld_get_diffcoef(w, x, ixI^L, ixO^L, D)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: D(ixI^S,1:ndim)
    double precision :: D_center(ixI^S)
    integer :: idir,i,j

    if (fld_diff_testcase) then
      D(ixI^S,1:ndim) = one/(unit_length*unit_velocity)
    else
      !> calculate diffusion coefficient
      D_center(ixO^S) = (const_c/unit_velocity)*w(ixO^S,i_lambda)/(w(ixO^S,i_op)*w(ixO^S,iw_rho))

      !> Extrapolate lambda to ghostcells
      !> Edges
      !> To calculate the diffusion coefficient at the ghostcells, copy lambda from grid, but use correct kappa and rho
      do i = 0,nghostcells-1
        D_center(ixImin1+i,:) = D_center(ixImin1+nghostcells,:)
        D_center(ixImax1-i,:) = D_center(ixImax1-nghostcells,:)
        D_center(:,ixImin2+i) = D_center(:,ixImin2+nghostcells)
        D_center(:,ixImax2-i) = D_center(:,ixImax2-nghostcells)
      end do

      !call Diff_boundary_conditions(ixI^L,ixO^L,D)

      !> Corners
      do i = 0,nghostcells-1
        do j = 0, nghostcells-1
          D_center(ixImin1+i,ixImax2-j) = D_center(ixImin1+nghostcells,ixImax2-nghostcells)
          D_center(ixImax1-i,ixImax2-j) = D_center(ixImax1-nghostcells,ixImax2-nghostcells)
          D_center(ixImin1+i,ixImin2+j) = D_center(ixImin1+nghostcells,ixImin2+nghostcells)
          D_center(ixImax1-i,ixImin2+j) = D_center(ixImax1-nghostcells,ixImin2+nghostcells)
        end do
      end do

      !> Go from cell center to cell face
      do i = ixImin1+1, ixImax1
      do j = ixImin2+1, ixImax2
         ! D(i,j,1) = (D_center(i,j) + D_center(i-1,j))/two
         ! D(i,j,2) = (D_center(i,j) + D_center(i,j-1))/two
        D(i,j,1) = (2*D_center(i,j) + 2*D_center(i-1,j)&
                  + D_center(i,j+1) + D_center(i-1,j+1)&
                  + D_center(i,j-1) + D_center(i-1,j-1))/8.d0
        D(i,j,2) = (2*D_center(i,j) + 2*D_center(i,j-1)&
                  + D_center(i+1,j) + D_center(i+1,j-1)&
                  + D_center(i-1,j) + D_center(i-1,j-1))/8.d0
      enddo
      enddo
      D(ixImin1,:,1) = D_center(ixImin1,:)
      D(:,ixImin2,1) = D_center(:,ixImin2)
      D(ixImin1,:,2) = D_center(ixImin1,:)
      D(:,ixImin2,2) = D_center(:,ixImin2)

      !D(:,ixImax2-2,:) = D(:,ixImax2-3,:)

    endif
  end subroutine fld_get_diffcoef

  subroutine get_diffusion_criterion(w, x, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    double precision             :: Q(ixO^S,1:ndim),diff_crit_mype
    integer                      :: jxO^L, hxO^L, idir

    do idir = 1,ndir
      hxO^L=ixO^L-kr(idir,^D);
      jxO^L=ixO^L+kr(idir,^D);

      Q(ixO^S,idir) = w(ixO^S,i_diff_mg)*dt/((x(jxO^S,idir)-x(hxO^S,idir))/2)**2
    enddo
    diff_crit_mype = minval(Q(ixO^S,1:ndim))

    !> Communicate 'Q' over processors
    call MPI_ALLREDUCE(diff_crit_mype,diff_crit,1,MPI_DOUBLE_PRECISION,MPI_MIN, &
                       icomm,ierrmpi)


  end subroutine get_diffusion_criterion

  !> Construct the matrix out of rad-hydro variables
  subroutine make_matrix(x,w,dw,E_m,E_n,sweepdir,ixImax,ixI^L,ixO^L,diag1,sub1,sup1,bvec1,diag2,sub2,sup2,bvec2)
    use mod_global_parameters

    integer, intent(in) :: sweepdir, ixImax
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw), dw
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(in) :: E_n(ixI^S), E_m(ixI^S)
    double precision, intent(out):: diag1(ixImin1:ixImax1,ixImin2:ixImax2),sub1(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out):: sup1(ixImin1:ixImax1,ixImin2:ixImax2),bvec1(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out):: diag2(ixImin2:ixImax2,ixImin1:ixImax1),sub2(ixImin2:ixImax2,ixImin1:ixImax1)
    double precision, intent(out):: sup2(ixImin2:ixImax2,ixImin1:ixImax1),bvec2(ixImin2:ixImax2,ixImin1:ixImax1)
    double precision :: D(ixI^S,1:ndim), h, beta(0:ixImax), delta_x
    integer :: idir,i,j

    call fld_get_diffcoef(w, x, ixI^L, ixO^L, D)

    !calculate h
    if (sweepdir == 1) then
      delta_x = x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1)
      !delta_x = x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2)
    elseif (sweepdir == 2) then
      !delta_x = x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1)
      delta_x = x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2)
    endif
    h = dw/(two*delta_x**two)

    !> Matrix depends on sweepingdirection
    if (sweepdir == 1) then
      !calculate matrix for sweeping in 1-direction
      do j = ixImin2,ixImax2
       !calculate beta
       do i = ixImin1,ixImax1-1
         beta(i) = one + dw/(two*dt) + h*(D(i+1,j,1)+D(i,j,1))
       enddo

       do i = ixImin1,ixImax1-1
         diag1(i,j) = beta(i)
         sub1(i+1,j) = -h*D(i+1,j,1)
         sup1(i,j) = -h*D(i+1,j,1)
         bvec1(i,j) = (one - h*(D(i,j+1,2)+D(i,j,2)))*E_m(i,j) &
         + h*D(i,j+1,2)*E_m(i,j+1) + h*D(i,j,2)*E_m(i,j-1) + dw/(two*dt)*E_n(i,j)
       enddo

       !> Boundary conditions on matrix
       sub1(ixImin1,j) = zero
       sup1(ixImax1,j) = zero
       diag1(ixImin1,j) = beta(ixImin1) - h*D(ixImin1,j,1)
       diag1(ixImax1,j) = beta(ixImax1-1) - h*D(ixImax1,j,1)
       bvec1(ixImax1,j) = (one - h*(D(ixImax1,j+1,2)+D(ixImax1,j,2)))*E_m(ixImax1,j) &
       + h*D(ixImax1,j+1,2)*E_m(ixImax1,j+1) + h*D(ixImax1,j,2)*E_m(ixImax1,j-1) + dw/(two*dt)*E_n(ixImax1,j)

      enddo

    elseif ( sweepdir == 2 ) then
      !calculate matrix for sweeping in 2-direction
      do j = ixImin1,ixImax1
       !calculate beta
       do i = ixImin2,ixImax2-1
         beta(i) = one + dw/(two*dt) + h*(D(j,i+1,2)+D(j,i,2))
       enddo

       do i = ixImin2,ixImax2-1
         diag2(i,j) = beta(i)
         sub2(i+1,j) = -h*D(j,i+1,2)
         sup2(i,j) = -h*D(j,i+1,2)
         bvec2(i,j) = (one - h*(D(j+1,i,1)+D(j,i,1)))*E_m(j,i) &
         + h*D(j+1,i,1)*E_m(j+1,i) + h*D(j,i,1)*E_m(j-1,i) + dw/(two*dt)*E_n(j,i)
       enddo

       !> Boundary conditions on matrix
       sub2(ixImin2,j) = zero
       sup2(ixImax2,j) = zero
       diag2(ixImin2,j) = beta(ixImin2) - h*D(j,ixImin2,2)
       diag2(ixImax2,j) = beta(ixImax2-1) - h*D(j,ixImax2,2)
       bvec2(ixImax2,j) = (one - h*(D(j+1,ixImax2,1)+D(j,ixImax2,1)))*E_m(j,ixImax2) &
       + h*D(j+1,ixImax2,1)*E_m(j+1,ixImax2) + h*D(j,ixImax2,1)*E_m(j-1,ixImax2) + dw/(two*dt)*E_n(j,ixImax2)
      enddo

    else
      call mpistop("sweepdirection unknown")
    endif
  end subroutine make_matrix

  !> Invert a tridiagonal matrix using Thomas' method
  subroutine solve_tridiag(ixOmin,ixOmax,ixImin,ixImax,diag,sub,sup,bvec,Evec)
    use mod_global_parameters
    implicit none

    integer, intent(in) :: ixOmin,ixOmax,ixImin,ixImax
    double precision, intent(in) :: diag(ixImin:ixImax), bvec(ixImin:ixImax)
    double precision, intent(in) :: sub(ixImin:ixImax), sup(ixImin:ixImax)
    double precision, intent(out) :: Evec(ixImin:ixImax)
    double precision :: cp(ixImin:ixImax), dp(ixImin:ixImax)
    integer :: i

    ! initialize c-prime and d-prime
    cp(ixImin) = sup(ixImin)/diag(ixImin)
    dp(ixImin) = bvec(ixImin)/diag(ixImin)

    ! solve for vectors c-prime and d-prime
    do i = ixImin+1 ,ixImax-1
      cp(i) = sup(i)/(diag(i)-cp(i-1)*sub(i))
      dp(i) = (bvec(i)-dp(i-1)*sub(i))/(diag(i)-cp(i-1)*sub(i))
    enddo
    dp(ixImax) = (bvec(ixImax)-dp(ixImax-1)*sub(ixImax))/(diag(ixImax)-cp(ixImax-1)*sub(ixImax))

    ! initialize x
    Evec(ixImax-1) = dp(ixImax-1)

    ! solve for x from the vectors c-prime and d-prime
    do i = ixImax-2, ixImin+1, -1
      Evec(i) = dp(i)-cp(i)*Evec(i+1)
    end do
  end subroutine solve_tridiag

  !> Perform boundary conditions to the tridiagonal matrix
  subroutine ADI_boundary_conditions(ixI^L,ixO^L,E_m,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(inout) :: E_m(ixI^S)
    integer :: iB
    integer g, h


    !call mpistop('Check if typeboundary(:,4) is defined')

    !> Boundary conditions for bound 1 (left)
    select case (typeboundary(iw_r_e,1))
    case('periodic')
      E_m(ixImin1:ixImin1+nghostcells-1,:) = E_m(ixImax1-2*nghostcells+1:ixImax1-nghostcells,:)
    case('cont')
      do g=1,nghostcells
        E_m(ixOmin1-g,:) = 2.d0*E_m(ixOmin1-g+1,:) - E_m(ixOmin1-g+2,:)
      enddo
    case('noinflow')
      do g=1,nghostcells
        E_m(ixOmin1-g,:) = 2.d0*E_m(ixOmin1-g+1,:) - E_m(ixOmin1-g+2,:)
        do h=ixImin2,ixImax2
          E_m(g,h) = min(E_m(g,h),E_m(g+1,h))
        enddo
      enddo
    case('fixed')
      E_m(ixImin1:ixImin1+nghostcells-1,:) = w(ixImin1:ixImin1+nghostcells-1,:,iw_r_e)
    case('special')
      call mpistop("Special adi not existing anymore")
    case default
      call mpistop("ADI boundary not defined")
    end select

    !> Boundary conditions for bound 2 (right)
    select case (typeboundary(iw_r_e,2))
    case('periodic')
      E_m(ixImax1-nghostcells+1:ixImax1,:) = E_m(ixImin1+nghostcells:ixImin1+2*nghostcells-1,:)
    case('cont')
      do g=1,nghostcells
        E_m(ixOmax1+g,:) = 2.d0*E_m(ixOmax1+g-1,:) - E_m(ixOmax1+g-2,:)
      enddo
    case('fixed')
      E_m(ixOmax1+1:ixImax1,:) = w(ixOmax1+1:ixImax1,:,iw_r_e)
    case('special')
      call mpistop("Special adi not existing anymore")

    case default
      call mpistop("ADI boundary not defined")
    end select

    !> Boundary conditions for bound 3 (bottom)
    select case (typeboundary(iw_r_e,3))
    case('periodic')
      E_m(:,ixImin2:ixImin2+nghostcells-1) = E_m(:,ixImax2-2*nghostcells+1:ixImax2-nghostcells)
    case('cont')
      do g=1,nghostcells
        E_m(:,ixOmin2-g) = 2.d0*E_m(:,ixOmin2-g+1) - E_m(:,ixOmin2-g+2)
      enddo
    case('fixed')
      E_m(:,ixImin2:ixImin2+nghostcells-1) = w(:,ixImin2:ixImin2+nghostcells-1,iw_r_e)
    case('special')
      call mpistop("Special adi not existing anymore")

    case default
      call mpistop("ADI boundary not defined")
    end select

    !> Boundary conditions for bound 4 (top)
    select case (typeboundary(iw_r_e,4))
    case('periodic')
      E_m(:,ixImax2-nghostcells+1:ixImax2) = E_m(:,ixImin2+nghostcells:ixImin2+2*nghostcells-1)
    case('cont')
      do g = 1,nghostcells
        E_m(:,ixOmax2+g) = 2.d0*E_m(:,ixOmax2+g-1) - E_m(:,ixOmax2+g-2)
      enddo
    case('noinflow')
      do g = 1,nghostcells
        E_m(:,ixOmax2+g) = 2.d0*E_m(:,ixOmax2+g-1) - E_m(:,ixOmax2+g-2)
        do h=ixImin1,ixImax1
          E_m(h,g) = min(E_m(h,g),E_m(h,g-1))
        enddo
      enddo
    case('fixed')
      E_m(:,ixImax2-nghostcells+1:ixImax2) = w(:,ixImax2-nghostcells+1:ixImax2,iw_r_e)
    case('special')
      call mpistop("Special adi not existing anymore")

    case default
      call mpistop("ADI boundary not defined")
    end select

    !> Corners
    do g = 0,nghostcells-1
      do h = 0, nghostcells-1
        E_m(ixOmin1-g,ixOmin2-h) = w(ixOmin1,ixOmin2,iw_r_e)
        E_m(ixOmax1+g,ixOmin2-h) = w(ixOmax1,ixOmin2,iw_r_e)
        E_m(ixOmin1-g,ixOmax2+h) = w(ixOmin1,ixOmax2,iw_r_e)
        E_m(ixOmax1+g,ixOmax2+h) = w(ixOmax1,ixOmax2,iw_r_e)
      end do
    end do
  end subroutine ADI_boundary_conditions

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! Gas-Rad Energy interaction
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> This subroutine calculates the radiation heating, radiation cooling
  !> and photon tiring using an implicit scheme.
  !> These sourceterms are applied using the root-finding of a 4th order polynomial
  !> This routine loops over every cell in the domain
  !> and computes the coefficients of the polynomials in every cell
  subroutine Energy_interaction(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry
    use mod_physics

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: div_v(ixI^S,1:ndir,1:ndir)
    double precision :: divvP(ixO^S)
    double precision :: temperature(ixI^S), vel(ixI^S)
    double precision :: a1(ixO^S), a2(ixO^S), a3(ixO^S)
    double precision :: c0(ixO^S), c1(ixO^S)
    double precision :: e_gas(ixO^S), E_rad(ixO^S)
    double precision :: grad_v(ixI^S)

    integer :: i,j,idir,ix^D

    double precision ::     rhd_gamma = 1.6666667d0

    if (it == 0) print*, "ATTENTION RHD_GAMMA FAKE"

    !> calculate tensor div_v
    do i = 1,ndir
      do j = 1,ndim
        vel(ixI^S) = w(ixI^S,iw_mom(j))/w(ixI^S,iw_rho)
        call gradient(vel,ixI^L,ixO^L,i,grad_v)
        div_v(ixO^S,i,j) = grad_v(ixO^S)
      enddo
    enddo

    !> VARIABLE NAMES DIV ARE ACTUALLY GRADIENTS
    {^IFONED
    divvP(ixO^S) = div_v(ixO^S,1,1)*w(ixO^S,i_edd(1,1))  &
    }

    {^IFTWOD
    !>eq 34 Turner and stone (Only 2D)
    divvP(ixO^S) = div_v(ixO^S,1,1)*w(ixO^S,i_edd(1,1))  &
                 + div_v(ixO^S,1,2)*w(ixO^S,i_edd(1,2))  &
                 + div_v(ixO^S,2,1)*w(ixO^S,i_edd(2,1))  &
                 + div_v(ixO^S,2,2)*w(ixO^S,i_edd(2,2))
    }

    {^IFTHREED
    divvP(ixO^S) = div_v(ixO^S,1,1)*w(ixO^S,i_edd(1,1))  &
                 + div_v(ixO^S,1,2)*w(ixO^S,i_edd(1,2))  &
                 + div_v(ixO^S,1,3)*w(ixO^S,i_edd(1,3))  &
                 + div_v(ixO^S,2,1)*w(ixO^S,i_edd(2,1))  &
                 + div_v(ixO^S,2,2)*w(ixO^S,i_edd(2,2))  &
                 + div_v(ixO^S,2,3)*w(ixO^S,i_edd(2,3))  &
                 + div_v(ixO^S,3,1)*w(ixO^S,i_edd(3,1))  &
                 + div_v(ixO^S,3,2)*w(ixO^S,i_edd(3,2))  &
                 + div_v(ixO^S,3,3)*w(ixO^S,i_edd(3,3))
    }

    divvP(ixO^S) = divvP(ixO^S)*w(ixO^S,iw_r_e)

    !> e_gas is the INTERNAL ENERGY without KINETIC ENERGY
    e_gas(ixO^S) = w(ixO^S,iw_e) - half*sum(w(ixO^S, iw_mom(:))**2, dim=ndim+1)/w(ixO^S, iw_rho)
    E_rad(ixO^S) = w(ixO^S,iw_r_e)

    !> Calculate coefficients for polynomial
    a1(ixO^S) = 4*w(ixO^S,i_op)*w(ixO^S,iw_rho)*(const_sigma*(unit_temperature**4.d0)/(unit_velocity*unit_pressure))*((rhd_gamma-one)/w(ixO^S,iw_rho))**4.d0*dt
    a2(ixO^S) = (const_c/unit_velocity)*w(ixO^S,i_op)*w(ixO^S,iw_rho)*dt
    a3(ixO^S) = divvP(ixO^S)/E_rad(ixO^S)*dt

    c0(ixO^S) = ((one + a2(ixO^S) + a3(ixO^S))*e_gas(ixO^S) + a2(ixO^S)*E_rad(ixO^S))/(a1(ixO^S)*(one + a3(ixO^S)))
    c1(ixO^S) = (one + a2(ixO^S) + a3(ixO^S))/(a1(ixO^S)*(one + a3(ixO^S)))

    !> Loop over every cell for rootfinding method
    {do ix^D=ixOmin^D,ixOmax^D\ }
      select case(fld_interaction_method)
      case('Bisect')
        call Bisection_method(e_gas(ix^D), E_rad(ix^D), c0(ix^D), c1(ix^D))
      case('Newton')
        call Newton_method(e_gas(ix^D), E_rad(ix^D), c0(ix^D), c1(ix^D))
      case('Halley')
        call Halley_method(e_gas(ix^D), E_rad(ix^D), c0(ix^D), c1(ix^D))
      case default
        call mpistop('root-method not known')
      end select
    {enddo\}

    !> Update gas-energy in w, internal + kinetic
    w(ixO^S,iw_e) = e_gas(ixO^S) + half*sum(w(ixO^S, iw_mom(:))**2, dim=ndim+1)/w(ixO^S, iw_rho)

    !> advance E_rad
    E_rad(ixO^S) = (a1(ixO^S)*e_gas(ixO^S)**4.d0 + E_rad(ixO^S))/(one + a2(ixO^S) + a3(ixO^S))

    !> Update rad-energy in w
    w(ixO^S,iw_r_e) = E_rad(ixO^S)

  end subroutine Energy_interaction

  !> Find the root of the 4th degree polynomial using the bisection method
  subroutine Bisection_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas

    double precision :: bisect_a, bisect_b, bisect_c
    integer :: n, max_its

    n = 0
    max_its = 1d7

    bisect_a = zero
    bisect_b = 1.2d0*max(abs(c0/c1),abs(c0)**(1.d0/4.d0))

    ! do while (abs(Polynomial_Bisection(bisect_b, c0, c1)-Polynomial_Bisection(bisect_a, c0, c1))&
    !    .ge. fld_bisect_tol*min(e_gas,E_rad))
    do while (abs(bisect_b-bisect_a) .ge. fld_bisect_tol*min(e_gas,E_rad))
      bisect_c = (bisect_a + bisect_b)/two

      n = n +1
      if (n .gt. max_its) then
        goto 2435
        call mpistop('No convergece in bisection scheme')
      endif

      if (Polynomial_Bisection(bisect_a, c0, c1)*&
      Polynomial_Bisection(bisect_b, c0, c1) .lt. zero) then

        if (Polynomial_Bisection(bisect_a, c0, c1)*&
        Polynomial_Bisection(bisect_c, c0, c1) .lt. zero) then
          bisect_b = bisect_c
        elseif (Polynomial_Bisection(bisect_b, c0, c1)*&
        Polynomial_Bisection(bisect_c, c0, c1) .lt. zero) then
          bisect_a = bisect_c
        elseif (Polynomial_Bisection(bisect_a, c0, c1) .eq. zero) then
          bisect_b = bisect_a
          bisect_c = bisect_a
          goto 2435
        elseif (Polynomial_Bisection(bisect_b, c0, c1) .eq. zero) then
          bisect_a = bisect_b
          bisect_c = bisect_b
          goto 2435
        elseif (Polynomial_Bisection(bisect_c, c0, c1) .eq. zero) then
          bisect_a = bisect_c
          bisect_b = bisect_c
          goto 2435
        else
          call mpistop("Problem with fld bisection method")
        endif
      elseif (Polynomial_Bisection(bisect_a, c0, c1) &
        - Polynomial_Bisection(bisect_b, c0, c1) .lt. fld_bisect_tol*Polynomial_Bisection(bisect_a, c0, c1)) then
        goto 2435
      else
        bisect_a = e_gas
        bisect_b = e_gas
        print*, "IGNORING GAS-RAD ENERGY EXCHANGE ", c0, c1

        print*, Polynomial_Bisection(bisect_a, c0, c1), Polynomial_Bisection(bisect_b, c0, c1)

        if (Polynomial_Bisection(bisect_a, c0, c1) .le. smalldouble) then
          bisect_b = bisect_a
        elseif (Polynomial_Bisection(bisect_a, c0, c1) .le. smalldouble) then
          bisect_a = bisect_b
        endif

        goto 2435

      endif
    enddo

      2435 e_gas = (bisect_a + bisect_b)/two
  end subroutine Bisection_method

  !> Find the root of the 4th degree polynomial using the Newton-Ralphson method
  subroutine Newton_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas

    double precision :: xval, yval, der, deltax

    integer :: ii

    yval = bigdouble
    xval = e_gas
    der = one
    deltax = one

    ii = 0
    !> Compare error with dx = dx/dy dy
    do while (abs(deltax) .gt. fld_bisect_tol)
      yval = Polynomial_Bisection(xval, c0, c1)
      der = dPolynomial_Bisection(xval, c0, c1)
      deltax = -yval/der
      xval = xval + deltax
      ii = ii + 1
      if (ii .gt. 1d2) then
        call Bisection_method(e_gas, E_rad, c0, c1)
        return
      endif
    enddo

    e_gas = xval
  end subroutine Newton_method

  !> Find the root of the 4th degree polynomial using the Halley method
  subroutine Halley_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas

    double precision :: xval, yval, der, dder, deltax

    integer :: ii

    yval = bigdouble
    xval = e_gas
    der = one
    dder = one
    deltax = one

    ii = 0
    !> Compare error with dx = dx/dy dy
    do while (abs(deltax) .gt. fld_bisect_tol)
      yval = Polynomial_Bisection(xval, c0, c1)
      der = dPolynomial_Bisection(xval, c0, c1)
      dder = ddPolynomial_Bisection(xval, c0, c1)
      deltax = -two*yval*der/(two*der**2 - yval*dder)
      xval = xval + deltax
      ii = ii + 1
      if (ii .gt. 1d2) then
        call mpistop('Halley did not convergggge')
        call Newton_method(e_gas, E_rad, c0, c1)
        return
      endif
    enddo

    e_gas = xval
  end subroutine Halley_method

  !> Evaluate polynomial at argument e_gas
  function Polynomial_Bisection(e_gas, c0, c1) result(val)
    use mod_global_parameters

    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: val

    val = e_gas**4.d0 + c1*e_gas - c0
  end function Polynomial_Bisection

  !> Evaluate first derivative of polynomial at argument e_gas
  function dPolynomial_Bisection(e_gas, c0, c1) result(der)
    use mod_global_parameters

    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: der

    der = 4.d0*e_gas**3.d0 + c1
  end function dPolynomial_Bisection

  !> Evaluate second derivative of polynomial at argument e_gas
  function ddPolynomial_Bisection(e_gas, c0, c1) result(dder)
    use mod_global_parameters

    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: dder

    dder = 4.d0*3.d0*e_gas**2.d0
  end function ddPolynomial_Bisection

end module mod_fld
