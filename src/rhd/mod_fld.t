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

    !> mean particle mass
    double precision, public :: fld_mu = 0.6d0

    !> Dimensionless Boltzman constante sigma
    double precision, public :: fld_sigma_0

    !> Dimensionless speed of light
    double precision, public :: fld_speedofligt_0

    !> Maximum amount of pseudotimesteps before trying something else
    integer, public :: fld_maxdw = 100

    !> Tolerance for bisection method for Energy sourceterms
    !> This is a percentage of the minimum of gas- and radiation energy
    double precision, public :: fld_bisect_tol = 1.d-5

    !> Tolerance for adi method for radiative Energy diffusion
    double precision, public :: fld_adi_tol = 1.d-2

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
    logical :: fld_complete_diffusion_limit = .false.

    !> diffusion coefficient for multigrid method
    integer :: i_diff_mg

    !> Whicht method to solve diffusion part
    character(len=8) :: fld_diff_scheme = 'adi'

    !> Set Diffusion coefficient to unity
    logical :: fld_diff_testcase = .false.

    !> public methods
    !> these are called in mod_hd_phys
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

  subroutine fld_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /fld_list/ fld_kappa0, fld_split, fld_maxdw, &
    fld_bisect_tol, fld_diff_testcase, fld_adi_tol, fld_max_fracdt,&
    fld_opacity_law, fld_complete_diffusion_limit, fld_diff_scheme

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, fld_list, end=111)
       111    close(unitpar)
    end do
  end subroutine fld_params_read

  subroutine fld_init(He_abundance)
    use mod_global_parameters
    use mod_variables
    use mod_physics, only: phys_global_source
    use mod_opacity, only: init_opal
    use mod_multigrid_coupling, only: mg_copy_boundary_conditions

    double precision, intent(in) :: He_abundance
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

    if (fld_diff_scheme .eq. 'mg') then

      use_multigrid = .true.

      phys_global_source => Diffuse_E_rad_mg
      mg_after_new_tree => set_mg_diffcoef

      mg%n_extra_vars = 1
      mg%operator_type = mg_vhelmholtz

      i_diff_mg = var_set_extravar("D", "D")
    endif

    !> Check if fld_numdt is not 1
    if (fld_maxdw .lt. 2) call mpistop("fld_maxdw should be an integer larger than 1")

    !> Need mean molecular weight
    ! fld_mu = (1.d0+4.d0*He_abundance)/two
    fld_mu = (1.+4*He_abundance)/(2.+3.*He_abundance)


    !> Dimensionless speed of light
    fld_speedofligt_0 = const_c/unit_velocity

    !> Dimensionless Boltzman constante sigma
    fld_sigma_0 = const_sigma*(unit_temperature**4.d0)/(unit_velocity*unit_pressure)

    !> Make kappa dimensionless !!!STILL NEED TO MULTIPLY W RHO
    fld_kappa0 = fld_kappa0/unit_opacity

    !> Read in opacity table if necesary
    if (fld_opacity_law .eq. 'opal') call init_opal(He_abundance)
  end subroutine fld_init

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
      call set_mg_bounds()
    endif
  end subroutine get_rad_extravars


  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine get_fld_rad_force(qdt,ixI^L,ixO^L,wCT,w,x,&
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
    double precision :: radiation_force(ixO^S,1:ndim)

    integer :: idir, i, jx^L

    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_split) then
      active = .true.

      !> Add momentum sourceterms
      jx^L=ixO^L+kr(idir,^D);


      do idir = 1,ndir
        !> Radiation force = kappa*rho/c *Flux
        radiation_force(ixO^S,idir) = w(ixO^S,i_op)*wCT(ixO^S,iw_rho)/fld_speedofligt_0*w(ixO^S, i_flux(idir))

        !> Momentum equation source term
        w(ixO^S,iw_mom(idir)) = w(ixO^S,iw_mom(idir)) &
            + qdt * radiation_force(ixO^S,idir)
            !+ qdt * half*(radiation_force(ixO^S,idir) + radiation_force(jx^S,idir))
            !> NOT SURE ON HOW TO AVERAGE OVER LEFTHANDSIDE AND RIGHTHANDSIDE FLUX EDGE
      enddo
    end if
  end subroutine get_fld_rad_force

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
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
  subroutine get_fld_diffusion(qdt,ixI^L,ixO^L,wCT,w,x,&
       energy,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods
    use mod_physics

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
        !> Do nothing OR CHECK WHAT IS ALREADY DONE BY POINTING
        call fld_get_diffcoef_central(w, x, ixI^L, ixO^L)
        call set_mg_diffcoef()
        ! call Diffuse_E_rad_mg(qdt, qt, active)
        call phys_global_source(dt, global_time, active)
      case default
        call mpistop('Numerical diffusionscheme unknown, try adi or mg')
      end select
      end if

      !> Set Diffcoef for next timestep?
      call fld_get_diffcoef_central(w, x, ixI^L, ixO^L)
  end subroutine get_fld_diffusion

  !> Sets the opacity in the w-array
  subroutine fld_get_opacity(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal
    use mod_physics, only: phys_get_tgas
    use mod_opacity

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision :: fld_kappa(ixO^S)
    double precision :: Temp(ixI^S)
    double precision :: rho0,Temp0,n,sigma_b

    integer :: i,j

    select case (fld_opacity_law)
      case('const')
        fld_kappa = fld_kappa0
      case('kramers')
        rho0 = half !> Take lower value of rho in domain
        fld_kappa(ixO^S) = fld_kappa0*((w(ixO^S,iw_rho)/rho0))
      case('bump')
        !> Opacity bump
        rho0 = 0.2d0 !0.5d-1
        n = 7.d0
        sigma_b = 2.d-2
        !fld_kappa(ixO^S) = fld_kappa0*(one + n*dexp(-((rho0  - w(ixO^S,iw_rho))**two)/rho0))
        fld_kappa(ixO^S) = fld_kappa0*(one + n*dexp(-one/sigma_b*(dlog(w(ixO^S,iw_rho)/rho0))**two))
      case('non_iso')
        call phys_get_pthermal(w,x,ixI^L,ixO^L,Temp)
        Temp(ixO^S)=Temp(ixO^S)/w(ixO^S,iw_rho)

        rho0 = 0.5d0 !> Take lower value of rho in domain
        Temp0 = one
        n = -7.d0/two
        fld_kappa(ixO^S) = fld_kappa0*(w(ixO^S,iw_rho)/rho0)*(Temp(ixO^S)/Temp0)**n
      case('opal')
        !call mpistop("Not implemented yet, hold your bloody horses")
        call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)
        do i = ixOmin1,ixOmax1
          do j= ixOmin2,ixOmax2
            rho0 = w(i,j,iw_rho)*unit_density
            Temp0 = Temp(i,j)*unit_temperature

            call set_opal_opacity(rho0,Temp0,n)

            fld_kappa(i,j) = n/unit_opacity
          enddo
        enddo
      case default
        call mpistop("Doesn't know opacity law")
      end select

      w(ixO^S, i_op) = fld_kappa(ixO^S)
  end subroutine fld_get_opacity

  !> Calculate fld flux limiter
  subroutine fld_get_fluxlimiter(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision :: fld_R(ixO^S), fld_lambda(ixO^S)
    double precision ::  normgrad2(ixO^S)
    double precision :: grad_r_e(ixI^S)
    integer :: idir

    if (fld_complete_diffusion_limit) then
      fld_lambda = one/3.d0

      w(ixO^S,i_lambda) = fld_lambda
      w(ixO^S,i_fld_R) = zero
    else
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixO^S) = zero
      do idir = 1,ndir
        call gradient(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,grad_r_e)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(w(ixO^S,i_op)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e))

      !> Calculate the flux limiter, lambda
      !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
      fld_lambda(ixO^S) = (2.d0+fld_R(ixO^S))/(6.d0+3*fld_R(ixO^S)+fld_R(ixO^S)**2.d0)

      w(ixO^S,i_lambda) = fld_lambda(ixO^S)
      w(ixO^S,i_fld_R) = fld_R(ixO^S)
    endif
  end subroutine fld_get_fluxlimiter

  !> Calculate Radiation Flux
  !> Returns Radiation flux and radiation pressure
  subroutine fld_get_radflux(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision :: rad_flux(ixO^S, 1:ndim)
    double precision :: L_star, R_star
    double precision :: normgrad2(ixO^S), f(ixO^S)
    double precision :: grad_r_e(ixI^S)
    double precision :: rad_e(ixI^S)
    integer :: idir

    rad_e(ixI^S) = w(ixI^S, iw_r_e)

    !> Calculate the Flux using the fld closure relation
    !> F = -c*lambda/(kappa*rho) *grad E
    do idir = 1,ndir
      call gradient(rad_e,ixI^L,ixO^L,idir,grad_r_e)
      rad_flux(ixO^S, idir) = -fld_speedofligt_0*w(ixO^S,i_lambda)/(w(ixO^S,i_op)*w(ixO^S,iw_rho))*grad_r_e(ixO^S)
    end do

    w(ixO^S,i_flux(:)) = rad_flux(ixO^S,:)
  end subroutine fld_get_radflux

  !> Calculate Radiation Pressure
  !> Returns Radiation Pressure
  subroutine fld_get_eddington(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision :: eddington_tensor(ixO^S,1:ndim,1:ndim)
    double precision :: tnsr2(ixO^S,1:ndim,1:ndim)
    double precision :: normgrad2(ixO^S), f(ixO^S)
    double precision :: grad_r_e(ixI^S, 1:ndim)
    integer :: i,j, idir,jdir

    !> Calculate R everywhere
    !> |grad E|/(rho kappa E)
    normgrad2(ixO^S) = zero

    do idir = 1,ndir
      call gradient(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,grad_r_e(ixI^S,idir))
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
  !> Returns Radiation Pressure
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

  subroutine Diffuse_E_rad_mg(qdt, qt, active)
    use mod_global_parameters
    use mod_multigrid_coupling
    use m_diffusion

    double precision, intent(in) :: qdt, qt
    logical, intent(inout)       :: active
    double precision             :: max_res

    call mg_copy_to_tree(iw_r_e, mg_iphi, .false., .false.)
    call diffusion_solve_vcoeff(mg, qdt, 1, 1.d-4)
    call mg_copy_from_tree(mg_iphi, iw_r_e)
    active = .true.
  end subroutine Diffuse_E_rad_mg

  subroutine fld_get_diffcoef_central(w, x, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    integer :: idir,i,j

    if (fld_diff_testcase) then
      !w(ixI^S,i_diff_mg) = one!*unit_length/unit_velocity

      w(ixO^S,i_diff_mg) = fld_speedofligt_0/w(ixO^S,i_op)

    else
      !> calculate diffusion coefficient
      w(ixO^S,i_diff_mg) = fld_speedofligt_0*w(ixO^S,i_lambda)/(w(ixO^S,i_op)*w(ixO^S,iw_rho))

      do i = ixImin1,ixImax1
        w(i,ixImin2:ixImin2+nghostcells-1,i_diff_mg) = w(i,ixImin2+nghostcells,i_diff_mg)
        w(i,ixImax2-nghostcells+1:ixImin2,i_diff_mg) = w(i,ixImax2-nghostcells,i_diff_mg)
      enddo
      do i = ixImin2,ixImax2
        w(ixImin1:ixImin1+nghostcells-1,i,i_diff_mg) = w(ixImin1+nghostcells,i,i_diff_mg)
        w(ixImax1-nghostcells+1:ixImin1,i,i_diff_mg) = w(ixImax1-nghostcells,i,i_diff_mg)
      enddo
    endif
  end subroutine fld_get_diffcoef_central

  subroutine set_mg_diffcoef()
    call mg_copy_to_tree(i_diff_mg, mg_iveps, .true., .true.)
  end subroutine set_mg_diffcoef

  subroutine set_mg_bounds()
    use mod_global_parameters
    integer :: iB

    do iB = 1,4
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
        !call mpistop('Hold your bloody horses, not implemented yet.')
        select case (iB)
        case (1)
           mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
           mg%bc(iB, mg_iphi)%bc_value = 0.0_dp ! Not needed
        case (2)
           mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
           mg%bc(iB, mg_iphi)%bc_value = 0.0_dp ! Not needed
        case (3)
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
        case (4)
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
        case default
          print *, "Not a standard: ", trim(typeboundary(iw_r_e, iB))
          error stop "You have to set a user-defined boundary method"
        end select
      case default
         print *, "Not a standard: ", trim(typeboundary(iw_r_e, iB))
         error stop "You have to set a user-defined boundary method"
      end select
    enddo
  end subroutine set_mg_bounds

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!! ADI
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
      if (ADI_Error .lt. fld_adi_tol) then
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
        if (ADI_Error .lt. fld_adi_tol) then
          converged = .true.
        endif
      endif

      !> If adjusting pseudostep doesn't work, divide the actual timestep in smaller parts
      if (w_max .gt. fld_maxdw) then
        if (converged .eqv. .false.) then
          !> use a smaller timestep than the hydrodynamical one
          call half_timestep_ADI(w, x, E_new, E_old, ixI^L, ixO^L, converged)
          call Error_check_ADI(w, x, E_new, E_old, ixI^L, ixO^L, ADI_Error)
          if (ADI_Error .lt. fld_adi_tol) then
            converged = .true.
          endif
        endif
      endif
    enddo

    w(ixO^S,iw_r_e) = E_new(ixO^S)
  end subroutine Evolve_E_rad

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
        if (ADI_Error .lt. fld_adi_tol) then
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

        if (ADI_Error .lt. fld_adi_tol) then
          converged = .true.
        endif

      enddo
      !---------------------------------------------------------------
      7895 E_loc = E_new
    enddo

    dt = saved_dt
  end subroutine half_timestep_ADI

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

  subroutine fld_get_diffcoef(w, x, ixI^L, ixO^L, D)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: D(ixI^S,1:ndim)
    double precision :: D_center(ixI^S)
    integer :: idir,i,j

    if (fld_diff_testcase) then
      !dE/dt + D dE/dx2
      ! [D] = x2/t = [v*l]

      D(ixI^S,1:ndim) = one/(unit_length*unit_velocity)


    else
      !> calculate diffusion coefficient
      D_center(ixO^S) = fld_speedofligt_0*w(ixO^S,i_lambda)/(w(ixO^S,i_op)*w(ixO^S,iw_rho))

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

  subroutine ADI_boundary_conditions(ixI^L,ixO^L,E_m,w,x)
    use mod_global_parameters
    use mod_usr_methods, only: usr_radiation_bc

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
      if (.not. associated(usr_radiation_bc)) then
        call mpistop("special ADI boundary not defined")
      endif
      call usr_radiation_bc(global_time,ixI^L,1,w,E_m,x)
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
      if (.not. associated(usr_radiation_bc)) then
        call mpistop("special ADI boundary not defined")
      endif
      call usr_radiation_bc(global_time,ixI^L,2,w,E_m,x)
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
      if (.not. associated(usr_radiation_bc)) then
        call mpistop("special ADI boundary not defined")
      endif
      call usr_radiation_bc(global_time,ixI^L,3,w,E_m,x)
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
      if (.not. associated(usr_radiation_bc)) then
        call mpistop("special ADI boundary not defined")
      endif
      call usr_radiation_bc(global_time,ixI^L,4,w,E_m,x)
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

  subroutine Energy_interaction(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry
    use mod_physics, only: phys_get_tgas

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

    integer :: i,j,idir

    !> calculate tensor div_v
    do i = 1,ndim
      do j = 1,ndim
        vel(ixI^S) = w(ixI^S,iw_mom(j))
        call gradient(vel,ixI^L,ixO^L,i,grad_v)
        div_v(ixO^S,i,j) = grad_v(ixO^S)
      enddo
    enddo

    !>eq 34 Turner and stone (Only 2D)
    divvP(ixO^S) = div_v(ixO^S,1,1)*w(ixO^S,i_edd(1,1))  &
                 + div_v(ixO^S,2,2)*w(ixO^S,i_edd(2,2))  &
                 + div_v(ixO^S,1,2)*w(ixO^S,i_edd(1,2))  &
                 + div_v(ixO^S,2,1)*w(ixO^S,i_edd(2,1))

    divvP(ixO^S) = divvP(ixO^S)*w(ixO^S,iw_r_e)

    !> Get temperature
    call phys_get_tgas(w,x,ixI^L,ixO^L,temperature)

    e_gas(ixO^S) = w(ixO^S,iw_e)
    E_rad(ixO^S) = w(ixO^S,iw_r_e)

    !> Calculate coefficients for polynomial
    a1(ixO^S) = 4*w(ixO^S,i_op)*w(ixO^S,iw_rho)*fld_sigma_0*(temperature(ixO^S)/e_gas(ixO^S))**4.d0*dt
    a2(ixO^S) = fld_speedofligt_0*w(ixO^S,i_op)*w(ixO^S,iw_rho)*dt
    a3(ixO^S) = divvP(ixO^S)/E_rad(ixO^S)*dt

    c0(ixO^S) = ((one + a1(ixO^S) + a3(ixO^S))*e_gas(ixO^S) + a2(ixO^S)*E_rad(ixO^S))/(a1(ixO^S)*(one + a3(ixO^S)))
    c1(ixO^S) = (one + a1(ixO^S) + a3(ixO^S))/(a1(ixO^S)*(one + a3(ixO^S)))

    w(ixO^S,i_test) = a1(ixO^S)

    !> Loop over every cell for bisection method
    do i = ixOmin1,ixOmax1
    do j =  ixOmin2,ixOmax2
          call Bisection_method(e_gas(i,j), E_rad(i,j), c0(i,j), c1(i,j))
    enddo
    enddo

    !> Update gas-energy in w
    w(ixO^S,iw_e) = e_gas(ixO^S)

    !> Calculate new radiation energy
    !> Get temperature
    call phys_get_tgas(w,x,ixI^L,ixO^L,temperature)

    !> Update a1
    a1(ixO^S) = 4*w(ixO^S,i_op)*w(ixO^S,iw_rho)*fld_sigma_0*(temperature(ixO^S)/e_gas(ixO^S))**4.d0*dt

    !> advance E_rad
    E_rad(ixO^S) = (a1*e_gas(ixO^S)**4.d0 + E_rad(ixO^S))/(one + a2 + a3)

    !> Update rad-energy in w
    w(ixO^S,iw_r_e) = E_rad(ixO^S)

    print*, w(4,4,iw_e)

  end subroutine Energy_interaction

  subroutine Bisection_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas

    double precision :: bisect_a, bisect_b, bisect_c

    bisect_a = zero
    bisect_b = max(abs(c0/c1),abs(c0)**(1.d0/4.d0))

    do while (abs(Polynomial_Bisection(bisect_b, c0, c1)-Polynomial_Bisection(bisect_a, c0, c1))&
       .ge. fld_bisect_tol*min(e_gas,E_rad))
      bisect_c = (bisect_a + bisect_b)/two

      ! print*, bisect_a, bisect_b, bisect_c
      ! print*, Polynomial_Bisection(bisect_a, c0, c1), Polynomial_Bisection(bisect_b, c0, c1), Polynomial_Bisection(bisect_c, c0, c1)


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
      else
        bisect_a = e_gas
        bisect_b = e_gas
        print*, "IGNORING ENERGY GAS-RAD EXCHANGE ", c0, c1
        goto 2435
      endif
    enddo

      2435 e_gas = (bisect_a + bisect_b)/two
  end subroutine Bisection_method

  function Polynomial_Bisection(e_gas, c0, c1) result(pol_result)
    use mod_global_parameters

    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: pol_result

    pol_result = e_gas**4.d0 + c1*e_gas - c0
  end function Polynomial_Bisection

end module mod_fld
