!> RadiatHydrodynamics physics module
module mod_rhd_phys

  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: rhd_energy = .true.

  !> Whether thermal conduction is added
  logical, public, protected              :: rhd_thermal_conduction = .false.

  !> Whether radiative cooling is added
  logical, public, protected              :: rhd_radiative_cooling = .false.

  !> Whether dust is added
  logical, public, protected              :: rhd_dust = .false.

  !> Whether viscosity is added
  logical, public, protected              :: rhd_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: rhd_gravity = .false.

  !> Whether particles module is added
  logical, public, protected              :: rhd_particles = .false.

  !> Number of tracer species
  integer, public, protected              :: rhd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  !> The adiabatic index
  !> Radiative gasses have a adiabatic index of 4/3, and not 5/3
  double precision, public, protected     :: rhd_gamma = 4.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: rhd_adiab = 1.0d0

  !> The smallest allowed energy
  double precision, protected             :: small_e

  !> The smallest allowed radiation energy
  double precision, protected             :: small_r_e = 0.d0

  !> Helium abundance over Hydrogen
  double precision, public, protected     :: He_abundance = 0.1d0

  !> Index of the radiation energy
  integer, public, protected              :: r_e

  !> Formalism to treat radiation
  character(len=8), public :: rhd_radiation_formalism = 'fld'

  !> Treat radiation fld_Rad_force
  logical, public, protected :: rhd_radiation_force = .true.

  !> Treat radiation-gas energy interaction
  logical, public, protected :: rhd_energy_interact = .true.

  !> Treat radiation energy diffusion
  logical, public, protected :: rhd_radiation_diffusion = .true.

  !> Treat radiation advection
  logical, public, protected :: rhd_radiation_advection = .true.

  !> Do a running mean over the radiation pressure when determining dt
  logical, protected :: radio_acoustic_filter = .false.
  integer, protected :: size_ra_filter = 1

  !> kb/(m_p mu)* 1/a_rad**4,
  double precision, public :: kbmpmua4

  !> Use the speed of light to calculate the timestep
  logical :: dt_c = .false.


  ! Public methods
  public :: rhd_phys_init
  public :: rhd_kin_en
  public :: rhd_get_pthermal
  public :: rhd_get_pradiation
  public :: rhd_get_ptot
  public :: rhd_to_conserved
  public :: rhd_to_primitive
  public :: rhd_get_tgas
  public :: rhd_get_trad
  {^NOONED
  public :: rhd_set_mg_bounds
  }

contains

  !> Read this module's parameters from a file
  subroutine rhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /rhd_list/ rhd_energy, rhd_n_tracer, rhd_gamma, rhd_adiab, &
    rhd_dust, rhd_thermal_conduction, rhd_radiative_cooling, rhd_viscosity, &
    rhd_gravity, He_abundance, SI_unit, rhd_particles, rhd_radiation_formalism,&
    rhd_radiation_force, rhd_energy_interact, rhd_radiation_diffusion, &
    rhd_radiation_advection, radio_acoustic_filter, size_ra_filter

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, rhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine rhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine rhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = rhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine rhd_write_info

  !> Add fluxes in an angular momentum conserving way
  subroutine rhd_angmomfix(fC,x,wnew,ixI^L,ixO^L,idim)
    use mod_global_parameters
    use mod_dust, only: dust_n_species, dust_mom
    use mod_geometry
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision, intent(inout)    :: fC(ixI^S,1:nwflux,1:ndim),  wnew(ixI^S,1:nw)
    integer, intent(in)                :: ixI^L, ixO^L
    integer, intent(in)                :: idim
    integer                            :: hxO^L, kxC^L, iw
    double precision                   :: inv_volume(ixI^S)

    logical isangmom

    ! shifted indexes
    hxO^L=ixO^L-kr(idim,^D);
    ! all the indexes
    kxCmin^D=hxOmin^D;
    kxCmax^D=ixOmax^D;

    inv_volume(ixO^S) = 1.0d0/block%dvolume(ixO^S)

    select case(coordinate)
    case (cylindrical)
       do iw=1,nwflux
        isangmom = (iw==iw_mom(phi_))
        if (rhd_dust) &
             isangmom = (isangmom .or. any(dust_mom(phi_,1:dust_n_species) == iw))
        if (idim==r_ .and. isangmom) then
          fC(kxC^S,iw,idim)= fC(kxC^S,iw,idim)*(x(kxC^S,r_)+half*block%dx(kxC^S,idim))
          wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
               (inv_volume(ixO^S)/x(ixO^S,idim))
        else
          wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
                inv_volume(ixO^S)
        endif
      enddo
     case (spherical)
      if (rhd_dust) &
        call mpistop("Error: rhd_angmomfix is not implemented &\\
        &with dust and coordinate==sperical")
      do iw=1,nwflux
        if     (idim==r_ .and. (iw==iw_mom(2) .or. iw==iw_mom(phi_))) then
          fC(kxC^S,iw,idim)= fC(kxC^S,iw,idim)*(x(kxC^S,idim)+half*block%dx(kxC^S,idim))
          wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
               (inv_volume(ixO^S)/x(ixO^S,idim))
        elseif (idim==2  .and. iw==iw_mom(phi_)) then
          fC(kxC^S,iw,idim)=fC(kxC^S,iw,idim)*sin(x(kxC^S,idim)+half*block%dx(kxC^S,idim)) ! (x(4,3,1)-x(3,3,1)))
          wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
               (inv_volume(ixO^S)/sin(x(ixO^S,idim)))
        else
          wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idim)-fC(hxO^S,iw,idim)) * &
                inv_volume(ixO^S)
        endif
      enddo

    end select

  end subroutine rhd_angmomfix

  !> Initialize the module
  subroutine rhd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_dust, only: dust_init
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init
    use mod_fld
    use mod_physics

    integer :: itr, idir

    call rhd_read_params(par_files)

    physics_type = "rhd"
    phys_energy  = rhd_energy
    use_particles = rhd_particles

    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (rhd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if

    !> set radiation energy
    r_e = var_set_radiation_energy()

    phys_get_dt              => rhd_get_dt
    phys_get_cmax            => rhd_get_cmax
    phys_get_cbounds         => rhd_get_cbounds
    phys_get_flux            => rhd_get_flux
    phys_get_v_idim          => rhd_get_v
    phys_add_source_geom     => rhd_add_source_geom
    phys_add_source          => rhd_add_source
    phys_to_conserved        => rhd_to_conserved
    phys_to_primitive        => rhd_to_primitive
    phys_check_params        => rhd_check_params
    phys_check_w             => rhd_check_w
    phys_get_pthermal        => rhd_get_pthermal
    phys_get_tgas            => rhd_get_tgas
    phys_get_trad            => rhd_get_trad
    phys_write_info          => rhd_write_info
    phys_handle_small_values => rhd_handle_small_values
    phys_angmomfix           => rhd_angmomfix
    {^NOONED
    phys_set_mg_bounds       => rhd_set_mg_bounds
    }

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .true.

    ! derive units from basic units
    call rhd_physical_units()

    if (rhd_dust) call dust_init(rho_, mom(:), e_)

    select case (rhd_radiation_formalism)
    case('fld')
      call fld_init(He_abundance, rhd_radiation_diffusion, rhd_gamma)
    case default
      call mpistop('Radiation formalism unknown')
    end select

    allocate(tracer(rhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, rhd_n_tracer
       tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! initialize thermal conduction module
    if (rhd_thermal_conduction) then
      if (.not. rhd_energy) &
           call mpistop("thermal conduction needs rhd_energy=T")
      call thermal_conduction_init(rhd_gamma)
    end if

    ! Initialize radiative cooling module
    if (rhd_radiative_cooling) then
      if (.not. rhd_energy) &
           call mpistop("radiative cooling needs rhd_energy=T")
      call radiative_cooling_init(rhd_gamma,He_abundance)
    end if

    ! Initialize viscosity module
    if (rhd_viscosity) call viscosity_init(phys_wider_stencil,phys_req_diagonal)

    ! Initialize gravity module
    if (rhd_gravity) call gravity_init()

    ! Initialize particles module
    if (rhd_particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1

    kbmpmua4 = unit_pressure**(-3./4)*unit_density*const_kB/(const_mp*fld_mu)*const_rad_a**(-1.d0/4)

    if (.not. rhd_energy .and. rhd_energy_interact) &
      call mpistop('Energy interact. not possible without gas energy eq.')

  end subroutine rhd_phys_init

  subroutine rhd_check_params
    use mod_global_parameters
    use mod_dust, only: dust_check_params

    if (.not. rhd_energy) then
       if (rhd_gamma <= 0.0d0) call mpistop ("Error: rhd_gamma <= 0")
       if (rhd_adiab < 0.0d0) call mpistop  ("Error: rhd_adiab < 0")
       small_pressure= rhd_adiab*small_density**rhd_gamma
    else
       if (rhd_gamma <= 0.0d0 .or. rhd_gamma == 1.0d0) &
            call mpistop ("Error: rhd_gamma <= 0 or rhd_gamma == 1.0")
       small_e = small_pressure/(rhd_gamma - 1.0d0)
    end if

    small_r_e = small_e

    if (rhd_dust) call dust_check_params()

    {^NOONED
    if (use_multigrid) call rhd_set_mg_bounds()
    }
  end subroutine rhd_check_params

  {^NOONED
  !> Set the boundaries for the diffusion of E
  subroutine rhd_set_mg_bounds
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_usr_methods

    integer :: iB

    ! Set boundary conditions for the multigrid solver
    do iB = 1, 2*ndim
       select case (typeboundary(r_e, iB))
       case ('symm')
          ! d/dx u = 0
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
       case ('asymm')
          ! u = 0
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
       case ('cont')
          ! d/dx u = 0
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
          mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
       case ('periodic')
          ! Nothing to do here
       case ('noinflow')
          call usr_special_mg_bc(iB)
       case ('special')
          call usr_special_mg_bc(iB)
       case default
          call mpistop("divE_multigrid warning: unknown b.c.: ", &
               trim(typeboundary(r_e, iB)))
       end select
    end do
  end subroutine rhd_set_mg_bounds
  }

  subroutine rhd_physical_units
    use mod_global_parameters
    double precision :: mp,kB
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
    else
      mp=mp_cgs
      kB=kB_cgs
    end if
    if(unit_velocity==0) then
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=(2.d0+3.d0*He_abundance)*unit_numberdensity*kB*unit_temperature
      unit_velocity=dsqrt(unit_pressure/unit_density)
      unit_time=unit_length/unit_velocity
    else
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*kB)
      unit_time=unit_length/unit_velocity
    end if

      unit_radflux = unit_velocity*unit_pressure
      unit_opacity = one/(unit_density*unit_length)
  end subroutine rhd_physical_units

  !> Returns 0 in argument flag where values are ok
  subroutine rhd_check_w(primitive, ixI^L, ixO^L, w, flag)
    use mod_global_parameters

    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    integer, intent(inout)       :: flag(ixI^S)
    double precision             :: tmp(ixI^S)

    flag(ixO^S) = 0
    where(w(ixO^S, rho_) < small_density) flag(ixO^S) = rho_

    if (rhd_energy) then
       if (primitive) then
          where(w(ixO^S, e_) < small_pressure) flag(ixO^S) = e_
       else
          tmp(ixO^S) = (rhd_gamma - 1.0d0)*(w(ixO^S, e_) - &
               rhd_kin_en(w, ixI^L, ixO^L))
          where(tmp(ixO^S) < small_pressure) flag(ixO^S) = e_
       endif
    end if

    where(w(ixO^S, r_e) < small_r_e) flag(ixO^S) = r_e

  end subroutine rhd_check_w

  !> Transform primitive variables into conservative ones
  subroutine rhd_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_dust, only: dust_to_conserved
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision                :: invgam
    integer                         :: idir, itr

    if (check_small_values .and. small_values_use_primitive) then
      call rhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'rhd_to_conserved')
    end if

    if (rhd_energy) then
       invgam = 1.d0/(rhd_gamma - 1.0d0)
       ! Calculate total energy from pressure and kinetic energy
       w(ixO^S, e_) = w(ixO^S, e_) * invgam + &
            0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) * w(ixO^S, rho_)
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, rho_) * w(ixO^S, mom(idir))
    end do

    if (rhd_dust) then
      call dust_to_conserved(ixI^L, ixO^L, w, x)
    end if

    if (check_small_values .and. .not. small_values_use_primitive) then
      call rhd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'rhd_to_conserved')
    end if
  end subroutine rhd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine rhd_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_dust, only: dust_to_primitive
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: itr, idir
    double precision                :: inv_rho(ixO^S)

    if (check_small_values .and. .not. small_values_use_primitive) then
      call rhd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'rhd_to_primitive')
    end if

    inv_rho = 1.0d0 / w(ixO^S, rho_)

    if (rhd_energy) then
       ! Compute pressure
       w(ixO^S, e_) = (rhd_gamma - 1.0d0) * (w(ixO^S, e_) - &
            rhd_kin_en(w, ixI^L, ixO^L, inv_rho))
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, mom(idir)) * inv_rho
    end do

    ! Convert dust momentum to dust velocity
    if (rhd_dust) then
      call dust_to_primitive(ixI^L, ixO^L, w, x)
    end if

    if (check_small_values .and. small_values_use_primitive) then
      call rhd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'rhd_to_primitive')
    end if

  end subroutine rhd_to_primitive

  subroutine e_to_rhos(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    if (rhd_energy) then
       w(ixO^S, e_) = (rhd_gamma - 1.0d0) * w(ixO^S, rho_)**(1.0d0 - rhd_gamma) * &
            (w(ixO^S, e_) - rhd_kin_en(w, ixI^L, ixO^L))
    else
       call mpistop("energy from entropy can not be used with -eos = iso !")
    end if
  end subroutine e_to_rhos

  subroutine rhos_to_e(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    if (rhd_energy) then
       w(ixO^S, e_) = w(ixO^S, rho_)**(rhd_gamma - 1.0d0) * w(ixO^S, e_) &
            / (rhd_gamma - 1.0d0) + rhd_kin_en(w, ixI^L, ixO^L)
    else
       call mpistop("entropy from energy can not be used with -eos = iso !")
    end if
  end subroutine rhos_to_e

  !> Calculate v_i = m_i / rho within ixO^L
  subroutine rhd_get_v(w, x, ixI^L, ixO^L, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:ndim)
    double precision, intent(out) :: v(ixI^S)

    v(ixO^S) = w(ixO^S, mom(idim)) / w(ixO^S, rho_)
  end subroutine rhd_get_v

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine rhd_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    double precision                          :: csound(ixI^S)
    double precision                          :: v(ixI^S)

    call rhd_get_v(w, x, ixI^L, ixO^L, idim, v)
    call rhd_get_csound2(w,x,ixI^L,ixO^L,csound)
    csound(ixO^S) = sqrt(csound(ixO^S))
    cmax(ixO^S) = abs(v(ixO^S))+csound(ixO^S)

    if (rhd_dust) then
      call dust_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
    end if
  end subroutine rhd_get_cmax

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine rhd_get_cbounds(wLC, wRC, wLp, wRp, x, ixI^L, ixO^L, idim, cmax, cmin)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative left and right status
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    ! primitive left and right status
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    double precision :: wmean(ixI^S,nw)
    double precision, dimension(ixI^S) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3

    if (typeboundspeed/='cmaxmean') then
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.

      tmp1(ixO^S)=sqrt(wLp(ixO^S,rho_))
      tmp2(ixO^S)=sqrt(wRp(ixO^S,rho_))
      tmp3(ixO^S)=1.d0/(sqrt(wLp(ixO^S,rho_))+sqrt(wRp(ixO^S,rho_)))
      umean(ixO^S)=(wLp(ixO^S,mom(idim))*tmp1(ixO^S)+wRp(ixO^S,mom(idim))*tmp2(ixO^S))*tmp3(ixO^S)

      if(rhd_energy) then
        csoundL(ixO^S)=rhd_gamma*wLp(ixO^S,p_)/wLp(ixO^S,rho_)
        csoundR(ixO^S)=rhd_gamma*wRp(ixO^S,p_)/wRp(ixO^S,rho_)
      else
        csoundL(ixO^S)=rhd_gamma*kbmpmua4*wLp(ixO^S,r_e)**(1.d0/4)
        csoundR(ixO^S)=rhd_gamma*kbmpmua4*wRp(ixO^S,r_e)**(1.d0/4)
      end if

      dmean(ixO^S) = (tmp1(ixO^S)*csoundL(ixO^S)+tmp2(ixO^S)*csoundR(ixO^S)) * &
           tmp3(ixO^S) + 0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2 * &
           (wRp(ixO^S,mom(idim))-wLp(ixO^S,mom(idim)))**2

      dmean(ixO^S)=sqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S)=umean(ixO^S)+dmean(ixO^S)
      else
        cmax(ixO^S)=abs(umean(ixO^S))+dmean(ixO^S)
      end if

      if (rhd_dust) then
        wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
        call dust_get_cmax(wmean, x, ixI^L, ixO^L, idim, cmax, cmin)
      end if

    else

      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
      tmp1(ixO^S)=wmean(ixO^S,mom(idim))/wmean(ixO^S,rho_)
      call rhd_get_csound2(wmean,x,ixI^L,ixO^L,csoundR)
      csoundR(ixO^S) = sqrt(csoundR(ixO^S))

      if(present(cmin)) then
        cmax(ixO^S)=max(tmp1(ixO^S)+csoundR(ixO^S),zero)
        cmin(ixO^S)=min(tmp1(ixO^S)-csoundR(ixO^S),zero)
      else
        cmax(ixO^S)=abs(tmp1(ixO^S))+csoundR(ixO^S)
      end if

      if (rhd_dust) then
        call dust_get_cmax(wmean, x, ixI^L, ixO^L, idim, cmax, cmin)
      end if
    end if

  end subroutine rhd_get_cbounds

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine rhd_get_csound2(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)


    if(rhd_energy) then
      call rhd_get_ptot(w,x,ixI^L,ixO^L,csound2)
      csound2(ixO^S)=max(rhd_gamma,4.d0/3.d0)*csound2(ixO^S)/w(ixO^S,rho_)
    else
      call rhd_get_ptot(w,x,ixI^L,ixO^L,csound2)
      csound2(ixO^S)=max(rhd_gamma,4.d0/3.d0)*csound2(ixO^S)/w(ixO^S,rho_)
    end if
  end subroutine rhd_get_csound2

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho) within ixO^L
  subroutine rhd_get_pthermal(w, x, ixI^L, ixO^L, pth)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: pth(ixI^S)

    if (rhd_energy) then
       pth(ixO^S) = (rhd_gamma - 1.0d0) * (w(ixO^S, e_) - &
            rhd_kin_en(w, ixI^L, ixO^L))
    else
       pth(ixI^S) = kbmpmua4*w(ixI^S, r_e)**(1.d0/4)/w(ixI^S, rho_)
    end if

  end subroutine rhd_get_pthermal

  !> Calculate radiation pressure within ixO^L
  subroutine rhd_get_pradiation(w, x, ixI^L, ixO^L, prad)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: prad(ixO^S, 1:ndim, 1:ndim)

    select case (rhd_radiation_formalism)
    case('fld')
      call fld_get_radpress(w, x, ixI^L, ixO^L, prad)
    case default
      call mpistop('Radiation formalism unknown')
    end select
  end subroutine rhd_get_pradiation

  !> calculates the sum of the gas pressure and max Prad tensor element
  subroutine rhd_get_ptot(w, x, ixI^L, ixO^L, ptot)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision             :: pth(ixI^S)
    double precision             :: prad_tensor(ixO^S, 1:ndim, 1:ndim)
    double precision             :: prad_max(ixO^S)
    double precision, intent(out):: ptot(ixI^S)
    integer :: ix^D

    call rhd_get_pthermal(w, x, ixI^L, ixO^L, pth)
    call rhd_get_pradiation(w, x, ixI^L, ixO^L, prad_tensor)

    {do ix^D = ixOmin^D,ixOmax^D\}
      prad_max(ix^D) = maxval(prad_tensor(ix^D,:,:))
    {enddo\}

    !> filter cmax
    if (radio_acoustic_filter) then
      call rhd_radio_acoustic_filter(x, ixI^L, ixO^L, prad_max)
    endif

    ptot(ixO^S) = pth(ixO^S) + prad_max(ixO^S)

  end subroutine rhd_get_ptot

  !> Filter peaks in cmax due to radiation energy density
  subroutine rhd_radio_acoustic_filter(x, ixI^L, ixO^L, prad_max)
    use mod_global_parameters

    integer, intent(in)                       :: ixI^L, ixO^L
    double precision, intent(in)              :: x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: prad_max(ixO^S)

    double precision :: tmp_prad(ixI^S)
    integer :: ix^D, filter, idim

    if (size_ra_filter .lt. 1) call mpistop("ra filter of size < 1 makes no sense")
    if (size_ra_filter .gt. nghostcells) call mpistop("ra filter of size < nghostcells makes no sense")

    tmp_prad(ixI^S) = zero
    tmp_prad(ixO^S) = prad_max(ixO^S)

    do filter = 1,size_ra_filter
      do idim = 1,ndim
        ! {do ix^D = ixOmin^D+filter,ixOmax^D-filter\}
        {do ix^D = ixOmin^D,ixOmax^D\}
            prad_max(ix^D) = min(tmp_prad(ix^D),tmp_prad(ix^D+filter*kr(idim,^D)))
            prad_max(ix^D) = min(tmp_prad(ix^D),tmp_prad(ix^D-filter*kr(idim,^D)))
        {enddo\}
      enddo
    enddo
  end subroutine rhd_radio_acoustic_filter

  !> Calculates gas temperature
  subroutine rhd_get_tgas(w, x, ixI^L, ixO^L, tgas)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision             :: pth(ixI^S)
    double precision, intent(out):: tgas(ixI^S)

    double precision :: mu

    call rhd_get_pthermal(w, x, ixI^L, ixO^L, pth)

    mu = (1.d0+4.d0*He_abundance)/(2.d0+3.d0*He_abundance)

    tgas(ixI^S) = pth(ixI^S)/w(ixI^S,rho_)*const_mp*mu/const_kB &
    *unit_pressure/(unit_density*unit_temperature)

  end subroutine rhd_get_tgas

  !> Calculates radiation temperature
  subroutine rhd_get_trad(w, x, ixI^L, ixO^L, trad)
    use mod_global_parameters
    use mod_constants

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: trad(ixI^S)

    trad(ixI^S) = (w(ixI^S,r_e)&
    /(const_rad_a*unit_temperature**4.d0/unit_pressure))**(1.d0/4.d0)

  end subroutine rhd_get_trad

  ! Calculate flux f_idim[iw]
  subroutine rhd_get_flux_cons(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    use mod_dust, only: dust_get_flux

    integer, intent(in)             :: ixI^L, ixO^L, idim
    double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(out)   :: f(ixI^S, nwflux)
    double precision                :: pth(ixI^S), v(ixI^S)
    integer                         :: idir, itr

    call rhd_get_pthermal(w, x, ixI^L, ixO^L, pth)
    call rhd_get_v(w, x, ixI^L, ixO^L, idim, v)

    f(ixO^S, rho_) = v(ixO^S) * w(ixO^S, rho_)

    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
      f(ixO^S, mom(idir)) = v(ixO^S) * w(ixO^S, mom(idir))
    end do

    f(ixO^S, mom(idim)) = f(ixO^S, mom(idim)) + pth(ixO^S)

    if(rhd_energy) then
      ! Energy flux is v_i*e + v*p ! Check? m_i/rho*p
      f(ixO^S, e_) = v(ixO^S) * (w(ixO^S, e_) + pth(ixO^S))
    end if

    if (rhd_radiation_advection) then
      f(ixO^S, r_e) = v(ixO^S) * w(ixO^S,r_e)
    else
      f(ixO^S, r_e) = zero
    endif

    do itr = 1, rhd_n_tracer
       f(ixO^S, tracer(itr)) = v(ixO^S) * w(ixO^S, tracer(itr))
    end do

    ! Dust fluxes
    if (rhd_dust) then
      call dust_get_flux(w, x, ixI^L, ixO^L, idim, f)
    end if

  end subroutine rhd_get_flux_cons

  ! Calculate flux f_idim[iw]
  subroutine rhd_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    use mod_dust, only: dust_get_flux_prim
    use mod_viscosity, only: visc_get_flux_prim ! viscInDiv

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in)    :: wC(ixI^S, 1:nw)
    ! primitive w
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(out)   :: f(ixI^S, nwflux)
    double precision                :: pth(ixI^S)
    integer                         :: idir, itr

    if (rhd_energy) then
       pth(ixO^S) = w(ixO^S,p_)
    else
       call rhd_get_pthermal(w, x, ixI^L, ixO^L, pth)
    end if

    f(ixO^S, rho_) = w(ixO^S,mom(idim)) * w(ixO^S, rho_)

    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
      f(ixO^S, mom(idir)) = w(ixO^S,mom(idim)) * wC(ixO^S, mom(idir))
    end do

    f(ixO^S, mom(idim)) = f(ixO^S, mom(idim)) + pth(ixO^S)

    if(rhd_energy) then
      ! Energy flux is v_i*e + v*p ! Check? m_i/rho*p
      f(ixO^S, e_) = w(ixO^S,mom(idim)) * (wC(ixO^S, e_) + w(ixO^S,p_))
    end if

    if (rhd_radiation_advection) then
      f(ixO^S, r_e) = w(ixO^S,mom(idim)) * wC(ixO^S,r_e)
    else
      f(ixO^S, r_e) = zero
    endif

    do itr = 1, rhd_n_tracer
       f(ixO^S, tracer(itr)) = w(ixO^S,mom(idim)) * w(ixO^S, tracer(itr))
    end do

    ! Dust fluxes
    if (rhd_dust) then
      call dust_get_flux_prim(w, x, ixI^L, ixO^L, idim, f)
    end if

    ! Viscosity fluxes - viscInDiv
    if (rhd_viscosity) then
      call visc_get_flux_prim(w, x, ixI^L, ixO^L, idim, f, rhd_energy)
    endif

  end subroutine rhd_get_flux

  !> Add geometrical source terms to w
  !>
  !> Notice that the expressions of the geometrical terms depend only on ndir,
  !> not ndim. Eg, they are the same in 2.5D and in 3D, for any geometry.
  !>
  !> Ileyk : to do :
  !>     - address the source term for the dust in case (typeaxial == 'spherical')
  subroutine rhd_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters
    use mod_viscosity, only: visc_add_source_geom ! viscInDiv
    use mod_dust, only: dust_n_species, dust_mom, dust_rho, dust_small_to_zero, set_dusttozero, dust_min_rho
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    ! to change and to set as a parameter in the parfile once the possibility to
    ! solve the equations in an angular momentum conserving form has been
    ! implemented (change tvdlf.t eg)
    double precision :: pth(ixI^S), source(ixI^S), minrho
    integer                         :: iw,idir, h1x^L{^NOONED, h2x^L}
    integer :: mr_,mphi_ ! Polar var. names
    integer :: irho, ifluid, n_fluids

    if (rhd_dust) then
       n_fluids = 1 + dust_n_species
    else
       n_fluids = 1
    end if

    select case (coordinate)
    case (cylindrical)
       do ifluid = 0, n_fluids-1
          ! s[mr]=(pthermal+mphi**2/rho)/radius
          if (ifluid == 0) then
             ! gas
             irho  = rho_
             mr_   = mom(r_)
             if(phi_>0) mphi_ = mom(phi_)
             call rhd_get_pthermal(wCT, x, ixI^L, ixO^L, source)
             minrho = 0.0d0
          else
             ! dust : no pressure
             irho  = dust_rho(ifluid)
             mr_   = dust_mom(r_, ifluid)
             mphi_ = dust_mom(phi_, ifluid)
             source(ixI^S) = zero
             minrho = dust_min_rho
          end if
          if (phi_ > 0) then
             where (wCT(ixO^S, irho) > minrho)
                source(ixO^S) = source(ixO^S) + wCT(ixO^S, mphi_)**2 / wCT(ixO^S, irho)
                w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * source(ixO^S) / x(ixO^S, r_)
             end where
             ! s[mphi]=(-mphi*mr/rho)/radius
             if(.not. angmomfix) then
                where (wCT(ixO^S, irho) > minrho)
                   source(ixO^S) = -wCT(ixO^S, mphi_) * wCT(ixO^S, mr_) / wCT(ixO^S, irho)
                   w(ixO^S, mphi_) = w(ixO^S, mphi_) + qdt * source(ixO^S) / x(ixO^S, r_)
                end where
             end if
          else
             ! s[mr]=2pthermal/radius
             w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * source(ixO^S) / x(ixO^S, r_)
          end if
       end do
    case (spherical)
       if (rhd_dust) then
          call mpistop("Dust geom source terms not implemented yet with spherical geometries")
       end if
       mr_   = mom(r_)
       if(phi_>0) mphi_ = mom(phi_)
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       ! s[mr]=((mtheta**2+mphi**2)/rho+2*p)/r
       call rhd_get_pthermal(wCT, x, ixI^L, ixO^L, pth)
       source(ixO^S) = pth(ixO^S) * x(ixO^S, 1) &
            *(block%surfaceC(ixO^S, 1) - block%surfaceC(h1x^S, 1)) &
            /block%dvolume(ixO^S)
       if (ndir > 1) then
         do idir = 2, ndir
           source(ixO^S) = source(ixO^S) + wCT(ixO^S, mom(idir))**2 / wCT(ixO^S, rho_)
         end do
       end if
       w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * source(ixO^S) / x(ixO^S, 1)

       {^NOONED
       ! s[mtheta]=-(mr*mtheta/rho)/r+cot(theta)*(mphi**2/rho+p)/r
       source(ixO^S) = pth(ixO^S) * x(ixO^S, 1) &
            * (block%surfaceC(ixO^S, 2) - block%surfaceC(h2x^S, 2)) &
            / block%dvolume(ixO^S)
       if (ndir == 3) then
          source(ixO^S) = source(ixO^S) + (wCT(ixO^S, mom(3))**2 / wCT(ixO^S, rho_)) / tan(x(ixO^S, 2))
       end if
       if (.not. angmomfix) then
          source(ixO^S) = source(ixO^S) - (wCT(ixO^S, mom(2)) * wCT(ixO^S, mr_)) / wCT(ixO^S, rho_)
       end if
       w(ixO^S, mom(2)) = w(ixO^S, mom(2)) + qdt * source(ixO^S) / x(ixO^S, 1)

       if (ndir == 3) then
         ! s[mphi]=-(mphi*mr/rho)/r-cot(theta)*(mtheta*mphi/rho)/r
         if (.not. angmomfix) then
           source(ixO^S) = -(wCT(ixO^S, mom(3)) * wCT(ixO^S, mr_)) / wCT(ixO^S, rho_)&
                      - (wCT(ixO^S, mom(2)) * wCT(ixO^S, mom(3))) / wCT(ixO^S, rho_) / tan(x(ixO^S, 2))
           w(ixO^S, mom(3)) = w(ixO^S, mom(3)) + qdt * source(ixO^S) / x(ixO^S, 1)
         end if
       end if
       }
    end select

    if (rhd_dust .and. dust_small_to_zero) then
       call set_dusttozero(qdt, ixI^L, ixO^L,  wCT,  w, x)
    end if

    if (rhd_viscosity) call visc_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)

  end subroutine rhd_add_source_geom

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine rhd_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_dust, only: dust_add_source, dust_mom, dust_rho, dust_n_species
    use mod_viscosity, only: viscosity_add_source
    use mod_usr_methods, only: usr_gravity
    use mod_gravity, only: gravity_add_source, grav_split
    use mod_dust, only: dust_small_to_zero, set_dusttozero


    use mod_fld


    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    double precision :: gravity_field(ixI^S, 1:ndim)
    integer :: idust, idim
    integer :: i

    if(rhd_dust) then
      call dust_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    end if

    if(rhd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
           qsourcesplit,active)
    end if

    if(rhd_viscosity) then
      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
           rhd_energy,qsourcesplit,active)
    end if

    if (rhd_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
           rhd_energy,qsourcesplit,active)

      if (rhd_dust .and. qsourcesplit .eqv. grav_split) then
         active = .true.

         call usr_gravity(ixI^L, ixO^L, wCT, x, gravity_field)
         do idust = 1, dust_n_species
            do idim = 1, ndim
               w(ixO^S, dust_mom(idim, idust)) = w(ixO^S, dust_mom(idim, idust)) &
                    + qdt * gravity_field(ixO^S, idim) * wCT(ixO^S, dust_rho(idust))
            end do
         end do
         if (dust_small_to_zero) then
            call set_dusttozero(qdt, ixI^L, ixO^L,  wCT,  w, x)
         end if
      end if
    end if

    call rhd_add_radiation_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)

  end subroutine rhd_add_source

  subroutine rhd_add_radiation_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods
    use mod_fld

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in) :: qsourcesplit
    logical, intent(inout) :: active
    double precision :: cmax(ixI^S)

    if (fld_diff_scheme .eq. 'mg') call fld_get_diffcoef_central(w, wCT, x, ixI^L, ixO^L)
    ! if (fld_diff_scheme .eq. 'mg') call set_mg_bounds(wCT, x, ixI^L, ixO^L)

    select case(rhd_radiation_formalism)
    case('fld')
      !> radiation force
      ! print*, it, 'Doing radforce stuff'
      if (rhd_radiation_force) call get_fld_rad_force(qdt,ixI^L,ixO^L,wCT,w,x,&
        rhd_energy,qsourcesplit,active)
      !> photon tiring, heating and cooling
      ! print*, it, 'Doing bisection stuff'
      if (rhd_energy) then
      if (rhd_energy_interact) call get_fld_energy_interact(qdt,ixI^L,ixO^L,wCT,w,x,&
        rhd_energy,qsourcesplit,active)
      endif
    case default
      call mpistop('Radiation formalism unknown')
    end select

    !>  NOT necessary for calculation, just want to know the grid-dependent-timestep
    call rhd_get_cmax(w, x, ixI^L, ixO^L, 2, cmax)
    w(ixI^S,i_test) = cmax(ixI^S)

  end subroutine rhd_add_radiation_source

  subroutine rhd_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    use mod_dust, only: dust_get_dt
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt
    use mod_fld, only: fld_radforce_get_dt

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

    if (.not. dt_c) then

      if(rhd_dust) then
        call dust_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
      end if

      if(rhd_radiation_force) then
        call fld_radforce_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
      endif

      if(rhd_radiative_cooling) then
        call cooling_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
      end if

      if(rhd_viscosity) then
        call viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
      end if

      if(rhd_gravity) then
        call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
      end if
    else
      dtnew = min(dx^D*unit_velocity/const_c)
    endif

  end subroutine rhd_get_dt

  function rhd_kin_en(w, ixI^L, ixO^L, inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)                    :: ixI^L, ixO^L
    double precision, intent(in)           :: w(ixI^S, nw)
    double precision                       :: ke(ixO^S)
    double precision, intent(in), optional :: inv_rho(ixO^S)

    if (present(inv_rho)) then
       ke = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) * inv_rho
    else
       ke = 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) / w(ixO^S, rho_)
    end if
  end function rhd_kin_en

  function rhd_inv_rho(w, ixI^L, ixO^L) result(inv_rho)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: inv_rho(ixO^S)

    ! Can make this more robust
    inv_rho = 1.0d0 / w(ixO^S, rho_)
  end function rhd_inv_rho

  subroutine rhd_handle_small_values(primitive, w, x, ixI^L, ixO^L, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    double precision :: smallone
    integer :: idir, flag(ixI^S)

    if (small_values_method == "ignore") return

    call rhd_check_w(primitive, ixI^L, ixO^L, w, flag)

    if (any(flag(ixO^S) /= 0)) then
      select case (small_values_method)
      case ("replace")
        if (small_values_fix_iw(rho_)) then
          where(flag(ixO^S) == rho_) w(ixO^S,rho_) = small_density
        end if

        do idir = 1, ndir
          if (small_values_fix_iw(mom(idir))) then
            where(flag(ixO^S) == rho_) w(ixO^S, mom(idir)) = 0.0d0
          end if
        end do

        if (rhd_energy) then
          if (small_values_fix_iw(e_)) then
            if(primitive) then
              where(flag(ixO^S) /= 0) w(ixO^S,e_) = small_pressure
            else
              where(flag(ixO^S) /= 0)
                ! Add kinetic energy
                w(ixO^S,e_) = small_e + 0.5d0 * &
                     sum(w(ixO^S, mom(:))**2, dim=ndim+1) / w(ixO^S, rho_)
              end where
            end if
          end if
        end if

        if (small_values_fix_iw(r_e)) then
          where(flag(ixO^S) /= 0) w(ixO^S,r_e) = small_r_e
        end if
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag)
      case default
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if
  end subroutine rhd_handle_small_values

end module mod_rhd_phys
