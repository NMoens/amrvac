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
    i_lambda = var_set_extravar("Lambda", "Lambda")
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
      ! case('adi')
      !   call Evolve_E_rad(w, x, ixI^L, ixO^L)
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
        forceM(ixO^S) = Qbar/(one-alpha) &
        *(gradv(ixO^S)/(w(ixO^S,iw_rho)*(const_c/unit_velocity)*Qbar*w(ixO^S,i_op)))**alpha
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
    integer :: idir, i, j, ix^D

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
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      fld_R(ixI^S) = dsqrt(normgrad2(ixO^S))/(w(ixO^S,i_op)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      fld_lambda(ixO^S) = one/fld_R(ixO^S)

      w(ixO^S,i_lambda) = fld_lambda(ixO^S)
      w(ixO^S,i_fld_R) = fld_R(ixO^S)

    case('Pomraning')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixI^S) = zero
      rad_e(ixI^S) = w(ixI^S, iw_r_e)
      do idir = 1,ndir
        !> gradient or gradientS ?!?!?!?!?!?
        call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(w(ixO^S,i_op)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
      fld_lambda(ixO^S) = (2.d0+fld_R(ixO^S))/(6.d0+3*fld_R(ixO^S)+fld_R(ixO^S)**2.d0)

      w(ixO^S,i_lambda) = fld_lambda(ixO^S)
      w(ixO^S,i_fld_R) = fld_R(ixO^S)

    case('Pomraning2')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixI^S) = zero
      rad_e(ixI^S) = w(ixI^S, iw_r_e)
      do idir = 1,ndir
        !> gradient or gradientS ?!?!?!?!?!?
        call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(w(ixO^S,i_op)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Levermore and Pomraning: lambda = 1/R(coth(R)-1/R)
      fld_lambda(ixO^S) = one/fld_R(ixO^S)*(one/dtanh(fld_R(ixO^S)) - one/fld_R(ixO^S))

      w(ixO^S,i_lambda) = fld_lambda(ixO^S)
      w(ixO^S,i_fld_R) = fld_R(ixO^S)


    case('Minerbo')
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixI^S) = zero
      rad_e(ixI^S) = w(ixI^S, iw_r_e)
      do idir = 1,ndir
        !> gradient or gradientS ?!?!?!?!?!?
        call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(w(ixO^S,i_op)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Minerbo:
      {do ix^D=ixOmin^D,ixOmax^D\ }
          if (fld_R(ix^D) .lt. 3.d0/2.d0) then
            fld_lambda(ix^D) = 2.d0/(3.d0 + dsqrt(9.d0 + 12.d0*fld_R(ix^D)**2.d0))
          else
            fld_lambda(ix^D) = 1.d0/(1.d0 + fld_R(ix^D) + dsqrt(1.d0 + 2.d0*fld_R(ix^D)))
          endif
      {enddo\}

      w(ixO^S,i_lambda) = fld_lambda(ixO^S)
      w(ixO^S,i_fld_R) = fld_R(ixO^S)
    case('special')
      if (.not. associated(usr_special_fluxlimiter)) then
        call mpistop("special fluxlimiter not defined")
      endif
      call usr_special_fluxlimiter(ixI^L, ixO^L, w, x, fld_lambda, fld_R)
      w(ixO^S,i_lambda) = fld_lambda(ixO^S)
      w(ixO^S,i_fld_R) = fld_R(ixO^S)
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
      rad_flux(ixO^S, idir) = -(const_c/unit_velocity)*w(ixO^S,i_lambda)/(w(ixO^S,i_op)*w(ixO^S,iw_rho))*grad_r_e(ixO^S)
    end do

    w(ixO^S,i_flux(:)) = rad_flux(ixO^S,:)

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

    call mg_copy_to_tree(iw_r_e, mg_iphi, .false., .false.)
    call diffusion_solve_vcoeff(mg, qdt, 2, fld_diff_tol)
    call mg_copy_from_tree(mg_iphi, iw_r_e)
    active = .true.

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

    else

      !> calculate diffusion coefficient
      w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*w(ixO^S,i_lambda)/(w(ixO^S,i_op)*w(ixO^S,iw_rho))

      if (diff_coef_filter) then
        !call mpistop('Hold your bloody horses, not implemented yet ')
        call fld_smooth_diffcoef(w, ixI^L, ixO^L)
      endif
    endif

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
    call mg_copy_to_tree(i_diff_mg, mg_iveps, .false., .false.)
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
