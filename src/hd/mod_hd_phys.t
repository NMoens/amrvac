!> Hydrodynamics physics module
module mod_hd_phys
  use mod_thermal_conduction, only: tc_fluid
  use mod_radiative_cooling, only: rc_fluid
  use mod_thermal_emission, only: te_fluid
  use mod_physics
  use mod_comm_lib, only: mpistop
  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: hd_energy = .true.
  logical, public, protected              :: hd_m1 = .true.

  !> Whether thermal conduction is added
  logical, public, protected              :: hd_thermal_conduction = .false.
  type(tc_fluid), allocatable, public :: tc_fl
  type(te_fluid), allocatable, public :: te_fl_hd

  !> Whether radiative cooling is added
  logical, public, protected              :: hd_radiative_cooling = .false.
  type(rc_fluid), allocatable, public :: rc_fl

  !> Whether dust is added
  logical, public, protected              :: hd_dust = .false.

  !> Whether viscosity is added
  logical, public, protected              :: hd_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: hd_gravity = .false.

  !> Whether particles module is added
  logical, public, protected              :: hd_particles = .false.

  !> Whether rotating frame is activated
  logical, public, protected              :: hd_rotating_frame = .false.

  !> Whether CAK radiation line force is activated
  logical, public, protected              :: hd_cak_force = .false.

  !> Number of tracer species
  integer, public, protected              :: hd_n_tracer = 0

  !> Whether plasma is partially ionized
  logical, public, protected              :: hd_partial_ionization = .false.

  !> Compute m1 radiation force
  logical, public, protected              :: m1_force = .true.

  !> Compute m1 radiation heating and cooling
  logical, public, protected              :: m1_heatcool = .true.

  !> Force split m1 force
  logical, public, protected              :: m1_force_split = .false.

  !> Compute m1 pressure in diffusion limit
  logical, public, protected              :: m1_diffusion = .false.

  !> Reduced speed of light factor
  double precision, public                :: m1_RSOL = 1.d0

  !> Theta parameter for stiff sourceterms
  double precision, public                :: m1_theta_stiff = 0.51d0

  !> Electron scattering opacity
  double precision, public                :: m1_kappa_e = 0.34d0

  !> Use well-balanced scheme for flux force term
  logical, public                         :: m1_well_balanced = .false.

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

  !> Index of the radiation energy density (-1 if not present)
  integer, public, protected              :: Er_

  !> Indices of the radiation flux (-1 if not present)
  integer, allocatable, public, protected :: Fr(:)

  !> (Reduced?) speed of light used in m1 subsystem
  double precision, public, protected     :: m1_c
  double precision, public, protected     :: m1_sb

  !> Opacity prescription
  character(len=8)                        :: m1_opacity_type = 'kappa_e'

  !> Indices of opacities
  integer, public, protected              :: kF_, kE_, kB_

  !> Indices of temperature
  integer, public, protected              :: Te_

  !> Index of the cutoff temperature for the TRAC method
  integer, public, protected              :: Tcoff_

  !> The adiabatic index
  double precision, public                :: hd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: hd_adiab = 1.0d0

  !> The mean molecular weight
  double precision, public                :: hd_mu = 0.6d0

  !> The small_est allowed energy
  double precision, protected             :: small_e

  !> Whether TRAC method is used
  logical, public, protected              :: hd_trac = .false.
  integer, public, protected              :: hd_trac_type = 1

  !> Allows overruling default corner filling (for debug mode, since otherwise corner primitives fail)
  logical, public, protected              :: hd_force_diagonal = .false.

  !> Helium abundance over Hydrogen
  double precision, public, protected     :: He_abundance=0.1d0
  !> Ionization fraction of H
  !> H_ion_fr = H+/(H+ + H)
  double precision, public, protected     :: H_ion_fr=1d0
  !> Ionization fraction of He
  !> He_ion_fr = (He2+ + He+)/(He2+ + He+ + He)
  double precision, public, protected     :: He_ion_fr=1d0
  !> Ratio of number He2+ / number He+ + He2+
  !> He_ion_fr2 = He2+/(He2+ + He+)
  double precision, public, protected     :: He_ion_fr2=1d0
  ! used for eq of state when it is not defined by units,
  ! the units do not contain terms related to ionization fraction
  ! and it is p = RR * rho * T
  double precision, public, protected     :: RR=1d0
  ! remove the below flag  and assume default value = .false.
  ! when eq state properly implemented everywhere
  ! and not anymore through units
  logical, public, protected              :: eq_state_units = .true.

  !>temporal local defenition of line driving parameters:
  double precision, public, protected :: a_F = 0.7d0
  double precision, public, protected :: a_E = 0.66d0
  double precision, public, protected :: a_B = 0.66d0
  double precision, public, protected :: Qb_F = 2.d3
  double precision, public, protected :: Qb_E = 2.d2
  double precision, public, protected :: Qb_B = 2.d2
  double precision, public, protected :: Q0_F = 2.d3
  double precision, public, protected :: Q0_E = 2.d2
  double precision, public, protected :: Q0_B = 2.d2

  procedure(sub_get_pthermal), pointer :: hd_get_Rfactor   => null()
  ! Public methods
  public :: hd_phys_init
  public :: hd_kin_en
  public :: hd_get_pthermal
  public :: hd_get_csound2
  public :: hd_to_conserved
  public :: hd_to_primitive
  public :: hd_check_params
  public :: hd_check_w
  public :: m1_get_prad

contains

  !> Read this module's parameters from a file
  subroutine hd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /hd_list/ hd_energy, hd_m1, m1_force, hd_n_tracer, hd_gamma, hd_adiab, &
    hd_dust, hd_thermal_conduction, hd_radiative_cooling, hd_viscosity, &
    hd_gravity, He_abundance,H_ion_fr, He_ion_fr, He_ion_fr2, eq_state_units, &
    SI_unit, hd_particles, hd_rotating_frame, hd_trac, &
    hd_force_diagonal, hd_trac_type, hd_cak_force, hd_partial_ionization , &
    m1_opacity_type, m1_diffusion, m1_RSOL, m1_force_split, m1_theta_stiff, &
    m1_heatcool, m1_well_balanced, m1_kappa_e

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, hd_list, end=111)
111    close(unitpar)
    end do

  end subroutine hd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine hd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = hd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine hd_write_info

  !> Initialize the module
  !TODO: Clean up m1 initialisation to seperate subroutine
  subroutine hd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_dust, only: dust_init
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init
    use mod_rotating_frame, only:rotating_frame_init
    use mod_cak_force, only: cak_init
    use mod_opal_opacity, only: init_opal_table
    use mod_cak_opacity, only: init_cak
    use mod_windsim_tools, only: init_windsim_tools, wst_get_photosphere
    use mod_supertimestepping, only: sts_init, add_sts_method,&
            set_conversion_methods_to_head, set_error_handling_to_head
    use mod_ionization_degree
    use mod_usr_methods, only: usr_Rfactor

    integer :: itr, idir

    call hd_read_params(par_files)

    physics_type = "hd"
    phys_energy  = hd_energy
    phys_total_energy  = hd_energy
    phys_internal_e = .false.
    phys_gamma = hd_gamma
    phys_partial_ionization=hd_partial_ionization

    phys_trac=hd_trac
    if(phys_trac) then
      if(ndim .eq. 1) then
        if(hd_trac_type .gt. 2) then
          hd_trac_type=1
          if(mype==0) write(*,*) 'WARNING: set hd_trac_type=1'
        end if
        phys_trac_type=hd_trac_type
      else
        phys_trac=.false.
        if(mype==0) write(*,*) 'WARNING: set hd_trac=F when ndim>=2'
      end if
    end if

    ! set default gamma for polytropic/isothermal process
    if(.not.hd_energy) then
      if(hd_thermal_conduction) then
        hd_thermal_conduction=.false.
        if(mype==0) write(*,*) 'WARNING: set hd_thermal_conduction=F when hd_energy=F'
      end if
      if(hd_radiative_cooling) then
        hd_radiative_cooling=.false.
        if(mype==0) write(*,*) 'WARNING: set hd_radiative_cooling=F when hd_energy=F'
      end if
    end if
    if(.not.eq_state_units) then
      if(hd_partial_ionization) then
        hd_partial_ionization=.false.
        if(mype==0) write(*,*) 'WARNING: set hd_partial_ionization=F when eq_state_units=F'
      end if
    end if
    use_particles = hd_particles

    allocate(start_indices(number_species),stop_indices(number_species))

    ! set the index of the first flux variable for species 1
    start_indices(1)=1
    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (hd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if

    !> Set index of M1 radiation variables
    if (hd_m1) then
      er_ = var_set_radiation_energy()
      allocate(Fr(ndir))
      Fr(:) = var_set_radiation_flux(ndir)
      kF_ = var_set_extravar("kappa_F", "kappa_F")
      kE_ = var_set_extravar("kappa_E", "kappa_E")
      kB_ = var_set_extravar("kappa_B", "kappa_B")
    else
      er_ = -1
      Fr(:) = -1
      kF_ = -1
      kE_ = -1
      kB_ = -1
    end if

    phys_get_dt              => hd_get_dt
    phys_get_cmax            => hd_get_cmax
    phys_get_a2max           => hd_get_a2max
    phys_get_tcutoff         => hd_get_tcutoff
    phys_get_cbounds         => hd_get_cbounds
    phys_get_flux            => hd_get_flux
    phys_add_source_geom     => hd_add_source_geom
    phys_add_source          => hd_add_source
    phys_to_conserved        => hd_to_conserved
    phys_to_primitive        => hd_to_primitive
    phys_check_params        => hd_check_params
    phys_check_w             => hd_check_w
    phys_get_pthermal        => hd_get_pthermal
    phys_get_v               => hd_get_v
    phys_write_info          => hd_write_info
    phys_handle_small_values => hd_handle_small_values

    if (hd_m1) then
      !> Heating and cooling sourceterms
      phys_implicit_update => m1_implicit_update
      phys_evaluate_implicit => m1_evaluate_implicit
      !> Setting opacities
      phys_get_auxiliary => m1_set_opacities
    end if

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    ! derive units from basic units
    call hd_physical_units()

    if (hd_m1) then
      m1_c = const_c/unit_velocity
      m1_sb = const_sigma/unit_radflux*unit_temperature**4
    end if

    if (hd_dust) then
        call dust_init(rho_, mom(:), e_)
    endif

    if (hd_force_diagonal) then
       ! ensure corners are filled, otherwise divide by zero when getting primitives
       !  --> only for debug purposes
       phys_req_diagonal = .true.
    endif

    allocate(tracer(hd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, hd_n_tracer
       tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! set number of variables which need update ghostcells
    nwgc=nwflux

    ! set the index of the last flux variable for species 1
    stop_indices(1)=nwflux

    !  set temperature as an auxiliary variable to get ionization degree
    if(hd_partial_ionization) then
      Te_ = var_set_auxvar('Te','Te')
    else
      Te_ = -1
    end if

    if(hd_trac) then
      Tcoff_ = var_set_wextra()
      iw_tcoff=Tcoff_
    else
      Tcoff_ = -1
    end if

    ! choose Rfactor in ideal gas law
    if(hd_partial_ionization) then
      hd_get_Rfactor=>Rfactor_from_temperature_ionization
      phys_update_temperature => hd_update_temperature
    else if(associated(usr_Rfactor)) then
      hd_get_Rfactor=>usr_Rfactor
    else
      hd_get_Rfactor=>Rfactor_from_constant_ionization
    end if

    ! initialize thermal conduction module
    if (hd_thermal_conduction) then
      if (.not. hd_energy) &
           call mpistop("thermal conduction needs hd_energy=T")
      phys_req_diagonal = .true.

      call sts_init()
      call tc_init_params(hd_gamma)

      allocate(tc_fl)
      call tc_get_hd_params(tc_fl,tc_params_read_hd)
      call add_sts_method(hd_get_tc_dt_hd,hd_sts_set_source_tc_hd,e_,1,e_,1,.false.)
      call set_conversion_methods_to_head(hd_e_to_ei, hd_ei_to_e)
      call set_error_handling_to_head(hd_tc_handle_small_e)
      tc_fl%get_temperature_from_conserved => hd_get_temperature_from_etot
      tc_fl%get_temperature_from_eint => hd_get_temperature_from_eint
      tc_fl%get_rho => hd_get_rho
      tc_fl%e_ = e_
    end if

    ! Initialize radiative cooling module
    if (hd_radiative_cooling) then
      if (.not. hd_energy) &
           call mpistop("radiative cooling needs hd_energy=T")
      call radiative_cooling_init_params(hd_gamma,He_abundance)
      allocate(rc_fl)
      call radiative_cooling_init(rc_fl,rc_params_read)
      rc_fl%get_rho => hd_get_rho
      rc_fl%get_pthermal => hd_get_pthermal
      rc_fl%get_var_Rfactor => hd_get_Rfactor
      rc_fl%e_ = e_
      rc_fl%Tcoff_ = Tcoff_
    end if
    allocate(te_fl_hd)
    te_fl_hd%get_rho=> hd_get_rho
    te_fl_hd%get_pthermal=> hd_get_pthermal
    te_fl_hd%get_var_Rfactor => hd_get_Rfactor
{^IFTHREED
    phys_te_images => hd_te_images
}
    ! Initialize viscosity module
    if (hd_viscosity) call viscosity_init(phys_wider_stencil,phys_req_diagonal)

    ! Initialize gravity module
    if (hd_gravity) call gravity_init()

    ! Initialize rotating_frame module
    if (hd_rotating_frame) call rotating_frame_init()

    ! Initialize CAK radiation force module
    if (hd_cak_force) call cak_init(hd_gamma)

    ! Initialise opal module if necessary
    if (hd_m1 )then
      if (m1_opacity_type .eq. 'hybrid1') call init_opal_table('Y02800')

      if (m1_opacity_type .eq. 'hybrid2') then
        call init_opal_table('Y02800')
        call init_cak('Y02800')

        !> Photosphere needed
        call init_windsim_tools(rho_,m1_kappa_e)
        phys_process_global => wst_get_photosphere
      end if
    end if

    ! Initialize particles module
    if (hd_particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

    !> Source-splitting for m1 forces
    if (m1_force_split) any_source_split=.true.

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    if (hd_m1) then
      nvector      = 2 ! No. vector vars
      allocate(iw_vector(nvector))
      iw_vector(1) = mom(1) - 1
      iw_vector(2) = Fr(1) - 1
    else
      nvector      = 1 ! No. vector vars
      allocate(iw_vector(nvector))
      iw_vector(1) = mom(1) - 1
    end if


    ! initialize ionization degree table
    if(hd_partial_ionization) call ionization_degree_init()

  end subroutine hd_phys_init

{^IFTHREED
  subroutine hd_te_images
    use mod_global_parameters
    use mod_thermal_emission
    select case(convert_type)
      case('EIvtiCCmpi','EIvtuCCmpi')
        call get_EUV_image(unitconvert,te_fl_hd)
      case('ESvtiCCmpi','ESvtuCCmpi')
        call get_EUV_spectrum(unitconvert,te_fl_hd)
      case('SIvtiCCmpi','SIvtuCCmpi')
        call get_SXR_image(unitconvert,te_fl_hd)
      case('WIvtiCCmpi','WIvtuCCmpi')
        call get_whitelight_image(unitconvert,te_fl_hd)
      case default
        call mpistop("Error in synthesize emission: Unknown convert_type")
      end select
  end subroutine hd_te_images
}
!!start th cond
  ! wrappers for STS functions in thermal_conductivity module
  ! which take as argument the tc_fluid (defined in the physics module)
  subroutine  hd_sts_set_source_tc_hd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux)
    use mod_global_parameters
    use mod_fix_conserve
    use mod_thermal_conduction, only: sts_set_source_tc_hd
    integer, intent(in) :: ixI^L, ixO^L, igrid, nflux
    double precision, intent(in) ::  x(ixI^S,1:ndim)
    double precision, intent(inout) ::  wres(ixI^S,1:nw), w(ixI^S,1:nw)
    double precision, intent(in) :: my_dt
    logical, intent(in) :: fix_conserve_at_step
    call sts_set_source_tc_hd(ixI^L,ixO^L,w,x,wres,fix_conserve_at_step,my_dt,igrid,nflux,tc_fl)
  end subroutine hd_sts_set_source_tc_hd

  function hd_get_tc_dt_hd(w,ixI^L,ixO^L,dx^D,x) result(dtnew)
    !Check diffusion time limit dt < dx_i**2/((gamma-1)*tc_k_para_i/rho)
    !where                      tc_k_para_i=tc_k_para*B_i**2/B**2
    !and                        T=p/rho
    use mod_global_parameters
    use mod_thermal_conduction, only: get_tc_dt_hd

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision :: dtnew

    dtnew=get_tc_dt_hd(w,ixI^L,ixO^L,dx^D,x,tc_fl)
  end function hd_get_tc_dt_hd

  subroutine hd_tc_handle_small_e(w, x, ixI^L, ixO^L, step)
    ! move this in a different  routine as in mhd if needed in more places
    use mod_global_parameters
    use mod_small_values

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer, intent(in)    :: step

    integer :: idir
    logical :: flag(ixI^S,1:nw)
    character(len=140) :: error_msg

    flag=.false.
    where(w(ixO^S,e_)<small_e) flag(ixO^S,e_)=.true.
    if(any(flag(ixO^S,e_))) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixO^S,e_)) w(ixO^S,e_)=small_e
      case ("average")
        call small_values_average(ixI^L, ixO^L, w, x, flag, e_)
      case default
        ! small values error shows primitive variables
        w(ixO^S,e_)=w(ixO^S,e_)*(hd_gamma - 1.0d0)
        do idir = 1, ndir
           w(ixO^S, iw_mom(idir)) = w(ixO^S, iw_mom(idir))/w(ixO^S,rho_)
        end do
        write(error_msg,"(a,i3)") "Thermal conduction step ", step
        call small_values_error(w, x, ixI^L, ixO^L, flag, error_msg)
      end select
    end if
  end subroutine hd_tc_handle_small_e

    ! fill in tc_fluid fields from namelist
    subroutine tc_params_read_hd(fl)
      use mod_global_parameters, only: unitpar,par_files
      type(tc_fluid), intent(inout) :: fl
      integer                      :: n
      logical :: tc_saturate=.false.
      double precision :: tc_k_para=0d0

      namelist /tc_list/ tc_saturate, tc_k_para

      do n = 1, size(par_files)
         open(unitpar, file=trim(par_files(n)), status="old")
         read(unitpar, tc_list, end=111)
111      close(unitpar)
      end do
      fl%tc_saturate = tc_saturate
      fl%tc_k_para = tc_k_para

    end subroutine tc_params_read_hd

  subroutine hd_get_rho(w,x,ixI^L,ixO^L,rho)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: rho(ixI^S)

    rho(ixO^S) = w(ixO^S,rho_)

  end subroutine hd_get_rho

!!end th cond
!!rad cool
    subroutine rc_params_read(fl)
      use mod_global_parameters, only: unitpar,par_files
      use mod_constants, only: bigdouble
      use mod_basic_types, only: std_len
      type(rc_fluid), intent(inout) :: fl
      integer                      :: n
      ! list parameters
      integer :: ncool = 4000
      double precision :: cfrac=0.1d0

      !> Name of cooling curve
      character(len=std_len)  :: coolcurve='JCcorona'

      !> Name of cooling method
      character(len=std_len)  :: coolmethod='exact'

      !> Fixed temperature not lower than tlow
      logical    :: Tfix=.false.

      !> Lower limit of temperature
      double precision   :: tlow=bigdouble

      !> Add cooling source in a split way (.true.) or un-split way (.false.)
      logical    :: rc_split=.false.


      namelist /rc_list/ coolcurve, coolmethod, ncool, cfrac, tlow, Tfix, rc_split

      do n = 1, size(par_files)
        open(unitpar, file=trim(par_files(n)), status="old")
        read(unitpar, rc_list, end=111)
111     close(unitpar)
      end do

      fl%ncool=ncool
      fl%coolcurve=coolcurve
      fl%coolmethod=coolmethod
      fl%tlow=tlow
      fl%Tfix=Tfix
      fl%rc_split=rc_split
      fl%cfrac=cfrac
    end subroutine rc_params_read
!! end rad cool

  subroutine hd_check_params
    use mod_global_parameters
    use mod_dust, only: dust_check_params, dust_implicit_update, dust_evaluate_implicit

    if (.not. hd_energy) then
       if (hd_gamma <= 0.0d0) call mpistop ("Error: hd_gamma <= 0")
       if (hd_adiab < 0.0d0) call mpistop  ("Error: hd_adiab < 0")
       small_pressure= hd_adiab*small_density**hd_gamma
    else
       if (hd_gamma <= 0.0d0 .or. hd_gamma == 1.0d0) &
            call mpistop ("Error: hd_gamma <= 0 or hd_gamma == 1.0")
       small_e = small_pressure/(hd_gamma - 1.0d0)
    end if

    if (hd_dust) call dust_check_params()
    if(use_imex_scheme .and. hd_dust) then
        ! implicit dust update
        phys_implicit_update => dust_implicit_update
        phys_evaluate_implicit => dust_evaluate_implicit
    endif

  end subroutine hd_check_params

  subroutine hd_physical_units
    use mod_global_parameters
    double precision :: mp,kB
    double precision :: a,b

    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
    else
      mp=mp_cgs
      kB=kB_cgs
    end if
    if(eq_state_units) then
      a = 1d0 + 4d0 * He_abundance
      if(hd_partial_ionization) then
        b = 2.d0+3.d0*He_abundance
      else
        b = 1d0 + H_ion_fr + He_abundance*(He_ion_fr*(He_ion_fr2 + 1d0)+1d0)
      end if
      RR = 1d0
    else
      a = 1d0
      b = 1d0
      RR = (1d0 + H_ion_fr + He_abundance*(He_ion_fr*(He_ion_fr2 + 1d0)+1d0))/(1d0 + 4d0 * He_abundance)
    end if
    if(unit_density/=1.d0) then
      unit_numberdensity=unit_density/(a*mp)
    else
      ! unit of numberdensity is independent by default
      unit_density=a*mp*unit_numberdensity
    end if
    if(unit_velocity/=1.d0) then
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
    else if(unit_pressure/=1.d0) then
      unit_temperature=unit_pressure/(b*unit_numberdensity*kB)
      unit_velocity=sqrt(unit_pressure/unit_density)
    else
      ! unit of temperature is independent by default
      unit_pressure=b*unit_numberdensity*kB*unit_temperature
      unit_velocity=sqrt(unit_pressure/unit_density)
    end if
    if(unit_time/=1.d0) then
      unit_length=unit_time*unit_velocity
    else
      ! unit of length is independent by default
      unit_time=unit_length/unit_velocity
    end if

    unit_mass = unit_density * unit_length**3
    unit_energy = unit_pressure * unit_length**3

    !> Radiation quantities
    unit_opacity = unit_length**2/unit_mass
    unit_radflux = unit_pressure*unit_velocity

    !> Mean molecular weight
    hd_mu = (1.0d0 + 4.0d0*He_abundance)/(2.0d0 + 3.0d0*He_abundance)
    m1_kappa_e = m1_kappa_e/unit_opacity

  end subroutine hd_physical_units

  !> Returns logical argument flag where values are ok
  subroutine hd_check_w(primitive, ixI^L, ixO^L, w, flag)
    use mod_global_parameters
    use mod_dust, only: dust_check_w

    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    logical, intent(inout)       :: flag(ixI^S,1:nw)
    double precision             :: tmp(ixI^S)

    flag=.false.

    if (hd_energy) then
       if (primitive) then
          where(w(ixO^S, e_) < small_pressure) flag(ixO^S,e_) = .true.
       else
          tmp(ixO^S) = (hd_gamma - 1.0d0)*(w(ixO^S, e_) - &
               hd_kin_en(w, ixI^L, ixO^L))
          where(tmp(ixO^S) < small_pressure) flag(ixO^S,e_) = .true.
       endif
    end if

    where(w(ixO^S, rho_) < small_density) flag(ixO^S,rho_) = .true.

    if (hd_dust) call dust_check_w(ixI^L,ixO^L,w,flag)

    if (hd_m1) call m1_check_w(ixI^L,ixO^L,w,flag)

  end subroutine hd_check_w

  !> Transform primitive variables into conservative ones
  subroutine hd_to_conserved(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_dust, only: dust_to_conserved
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision                :: invgam
    integer                         :: idir, itr

    if (fix_small_values) then
     call hd_handle_small_values(.true., w, x, ixI^L, ixO^L, 'hd_to_conserved')
    end if

    if (hd_energy) then
       invgam = 1.d0/(hd_gamma - 1.0d0)
       ! Calculate total energy from pressure and kinetic energy
       w(ixO^S, e_) = w(ixO^S, e_) * invgam + &
            0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1) * w(ixO^S, rho_)
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, rho_) * w(ixO^S, mom(idir))
    end do

    if (hd_dust) then
      call dust_to_conserved(ixI^L, ixO^L, w, x)
    end if

  end subroutine hd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine hd_to_primitive(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_dust, only: dust_to_primitive
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    integer                         :: itr, idir
    double precision                :: inv_rho(ixO^S)

    if (fix_small_values) then
      call hd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'hd_to_primitive')
    end if

    inv_rho = 1.0d0 / w(ixO^S, rho_)

    if (hd_energy) then
       ! Compute pressure
       w(ixO^S, e_) = (hd_gamma - 1.0d0) * (w(ixO^S, e_) - &
            hd_kin_en(w, ixI^L, ixO^L, inv_rho))
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixO^S, mom(idir)) = w(ixO^S, mom(idir)) * inv_rho
    end do

    ! Convert dust momentum to dust velocity
    if (hd_dust) then
      call dust_to_primitive(ixI^L, ixO^L, w, x)
    end if

  end subroutine hd_to_primitive

  !> Transform internal energy to total energy
  subroutine hd_ei_to_e(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    ! Calculate total energy from internal and kinetic energy
    w(ixO^S,e_)=w(ixO^S,e_)&
               +hd_kin_en(w,ixI^L,ixO^L)

  end subroutine hd_ei_to_e

  !> Transform total energy to internal energy
  subroutine hd_e_to_ei(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)

    ! Calculate ei = e - ek
    w(ixO^S,e_)=w(ixO^S,e_)&
                -hd_kin_en(w,ixI^L,ixO^L)

  end subroutine hd_e_to_ei

  subroutine e_to_rhos(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    if (hd_energy) then
       w(ixO^S, e_) = (hd_gamma - 1.0d0) * w(ixO^S, rho_)**(1.0d0 - hd_gamma) * &
            (w(ixO^S, e_) - hd_kin_en(w, ixI^L, ixO^L))
    else
       call mpistop("energy from entropy can not be used with -eos = iso !")
    end if
  end subroutine e_to_rhos

  subroutine rhos_to_e(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision             :: w(ixI^S, nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)

    if (hd_energy) then
       w(ixO^S, e_) = w(ixO^S, rho_)**(hd_gamma - 1.0d0) * w(ixO^S, e_) &
            / (hd_gamma - 1.0d0) + hd_kin_en(w, ixI^L, ixO^L)
    else
       call mpistop("entropy from energy can not be used with -eos = iso !")
    end if
  end subroutine rhos_to_e

  !> Calculate v_i = m_i / rho within ixO^L
  subroutine hd_get_v_idim(w, x, ixI^L, ixO^L, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S, nw), x(ixI^S, 1:ndim)
    double precision, intent(out) :: v(ixI^S)

    v(ixO^S) = w(ixO^S, mom(idim)) / w(ixO^S, rho_)
  end subroutine hd_get_v_idim

  !> Calculate velocity vector v_i = m_i / rho within ixO^L
  subroutine hd_get_v(w,x,ixI^L,ixO^L,v)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,nw), x(ixI^S,1:^ND)
    double precision, intent(out) :: v(ixI^S,1:ndir)

    integer :: idir

    do idir=1,ndir
      v(ixO^S,idir) = w(ixO^S, mom(idir)) / w(ixO^S, rho_)
    end do

  end subroutine hd_get_v

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine hd_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax

    integer, intent(in)                       :: ixI^L, ixO^L, idim
    double precision, intent(in)              :: w(ixI^S, nw), x(ixI^S, 1:ndim)
    double precision, intent(inout)           :: cmax(ixI^S)
    double precision                          :: csound(ixI^S)
    double precision                          :: v(ixI^S)

    double precision                          :: fi(ixI^S), fj(ixI^S)
    double precision                          :: ff(ixI^S), xi(ixI^S)

    if (hd_m1) then
      ! TODO: Implement actual eigenvalues
      !> e.g. Gonzalez+ 2007
      ! cmax(ixO^S) = m1_c*m1_RSOL

      ! print*, "GUTENTAG"

      !TODO: Check this sum, it does not produce what is expected
      ff(ixI^S) = sum(w(ixI^S, Fr(:))**2, dim=ndim+1)/(m1_c*w(ixI^S,Er_))**2

      ! print*, 'ff'
      ! print*, dsqrt(ff(:,5))

      !> Eddington scalars
      fi(ixI^S) = abs(w(ixI^S,Fr(idim))/(m1_c*w(ixI^S,Er_)))
      ! print*, 'fi'
      ! print*, fi(:, 5)

      fj(ixI^S) = dsqrt(max(ff(ixI^S) - fi(ixI^S)**2,0.d0))
      ! print*, 'fj'
      ! print*, fj(:, 5)

      xi(ixI^S) = dsqrt(4-3*ff(ixI^S))
      ! print*, 'xi'
      ! print*, xi(:, 5)

      cmax(ixI^S) = m1_c*m1_RSOL*(fi(ixI^S)/xi(ixI^S) &
                  + dsqrt(2*(xi(ixI^S)-1)*(xi(ixI^S)+2)*(2*(xi(ixI^S)-1)*(xi(ixI^S)+2) + 3*fj(ixI^S)**2)) &
                  /(dsqrt(3.d0)*xi(ixI^S)*(xi(ixI^S)+2)))

      ! print*, 'cmax'
      ! print*, cmax(:, 5)

      ! print*, 'eigenvalue'
      ! print*, (fi(ixI^S)/xi(ixI^S) &
      !         + dsqrt(2*(xi(ixI^S)-1)*(xi(ixI^S)+2)*(2*(xi(ixI^S)-1)*(xi(ixI^S)+2) + 3*fj(ixI^S)**2)) &
      !         /(dsqrt(3.d0)*xi(ixI^S)*(xi(ixI^S)+2)))


    else
      call hd_get_v_idim(w, x, ixI^L, ixO^L, idim, v)
      call hd_get_csound2(w,x,ixI^L,ixO^L,csound)
      csound(ixO^S) = dsqrt(csound(ixO^S))

      cmax(ixO^S) = dabs(v(ixO^S))+csound(ixO^S)

      if (hd_dust) then
        call dust_get_cmax(w, x, ixI^L, ixO^L, idim, cmax)
      end if

    end if

  end subroutine hd_get_cmax

  subroutine hd_get_a2max(w,x,ixI^L,ixO^L,a2max)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: a2max(ndim)
    double precision :: a2(ixI^S,ndim,nw)
    integer :: gxO^L,hxO^L,jxO^L,kxO^L,i,j

    a2=zero
    do i = 1,ndim
      !> 4th order
      hxO^L=ixO^L-kr(i,^D);
      gxO^L=hxO^L-kr(i,^D);
      jxO^L=ixO^L+kr(i,^D);
      kxO^L=jxO^L+kr(i,^D);
      a2(ixO^S,i,1:nw)=dabs(-w(kxO^S,1:nw)+16.d0*w(jxO^S,1:nw)&
        -30.d0*w(ixO^S,1:nw)+16.d0*w(hxO^S,1:nw)-w(gxO^S,1:nw))
      a2max(i)=maxval(a2(ixO^S,i,1:nw))/12.d0/dxlevel(i)**2
    end do
  end subroutine hd_get_a2max

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine hd_get_tcutoff(ixI^L,ixO^L,w,x,tco_local,Tmax_local)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(out) :: tco_local, Tmax_local

    double precision, parameter :: trac_delta=0.25d0
    double precision :: tmp1(ixI^S),Te(ixI^S),lts(ixI^S)
    double precision :: ltr(ixI^S),ltrc,ltrp,tcoff(ixI^S)
    integer :: jxO^L,hxO^L
    integer :: jxP^L,hxP^L,ixP^L
    logical :: lrlt(ixI^S)

    {^IFONED
    call hd_get_temperature_from_etot(w,x,ixI^L,ixI^L,Te)

    Tco_local=zero
    Tmax_local=maxval(Te(ixO^S))
    select case(hd_trac_type)
    case(0)
      w(ixI^S,Tcoff_)=3.d5/unit_temperature
    case(1)
      hxO^L=ixO^L-1;
      jxO^L=ixO^L+1;
      lts(ixO^S)=0.5d0*dabs(Te(jxO^S)-Te(hxO^S))/Te(ixO^S)
      lrlt=.false.
      where(lts(ixO^S) > trac_delta)
        lrlt(ixO^S)=.true.
      end where
      if(any(lrlt(ixO^S))) then
        Tco_local=maxval(Te(ixO^S), mask=lrlt(ixO^S))
      end if
    case(2)
      !> iijima et al. 2021, LTRAC method
      ltrc=1.5d0
      ltrp=2.5d0
      ixP^L=ixO^L^LADD1;
      hxO^L=ixO^L-1;
      jxO^L=ixO^L+1;
      hxP^L=ixP^L-1;
      jxP^L=ixP^L+1;
      lts(ixP^S)=0.5d0*abs(Te(jxP^S)-Te(hxP^S))/Te(ixP^S)
      ltr(ixP^S)=max(one, (exp(lts(ixP^S))/ltrc)**ltrp)
      w(ixO^S,Tcoff_)=Te(ixO^S)*&
        (0.25*(ltr(jxO^S)+two*ltr(ixO^S)+ltr(hxO^S)))**0.4d0
    case default
      call mpistop("mhd_trac_type not allowed for 1D simulation")
    end select
    }
  end subroutine hd_get_tcutoff

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine hd_get_cbounds(wLC, wRC, wLp, wRp, x, ixI^L, ixO^L, idim,Hspeed,cmax, cmin)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax
    use mod_variables

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative left and right status
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    ! primitive left and right status
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: cmax(ixI^S,1:number_species)
    double precision, intent(inout), optional :: cmin(ixI^S,1:number_species)
    double precision, intent(in)    :: Hspeed(ixI^S,1:number_species)

    double precision :: wmean(ixI^S,nw)
    double precision, dimension(ixI^S) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3

    double precision                          :: fi(ixI^S), fj(ixI^S)
    double precision                          :: ff(ixI^S), xi(ixI^S)

    double precision                          :: lam_L_max(ixI^S), lam_L_min(ixI^S)
    double precision                          :: lam_R_max(ixI^S), lam_R_min(ixI^S)

    integer :: ix^D

    if (hd_m1) then

      ! TODO: Implement actual eigenvalues
      !> e.g. Gonzalez+ 2007
      ! cmax(ixI^S,1) = m1_c*m1_RSOL

      !TODO: Check this sum, it does not produce what is expected
      !> Left state
      ff(ixI^S) = sum(wLC(ixI^S, Fr(:))**2, dim=ndim+1)/(m1_c*wLC(ixI^S,Er_))**2
      ff(ixI^S) = min(ff(ixI^S),1.d0)

      !> Eddington scalars
      fi(ixI^S) = abs(wLC(ixI^S,Fr(idim))/(m1_c*wLC(ixI^S,Er_)))
      fi(ixI^S) = min(fi(ixI^S),1.d0)

      fj(ixI^S) = dsqrt(max(ff(ixI^S) - fi(ixI^S)**2,0.d0))

      xi(ixI^S) = dsqrt(4-3*ff(ixI^S))

      lam_L_max(ixI^S) = m1_c*m1_RSOL*(fi(ixI^S)/xi(ixI^S) &
                  + dsqrt(2*(xi(ixI^S)-1)*(xi(ixI^S)+2)*(2*(xi(ixI^S)-1)*(xi(ixI^S)+2) + 3*fj(ixI^S)**2)) &
                  /(dsqrt(3.d0)*xi(ixI^S)*(xi(ixI^S)+2)))

      lam_L_min(ixI^S) = m1_c*m1_RSOL*(fi(ixI^S)/xi(ixI^S) &
                  - dsqrt(2*(xi(ixI^S)-1)*(xi(ixI^S)+2)*(2*(xi(ixI^S)-1)*(xi(ixI^S)+2) + 3*fj(ixI^S)**2)) &
                  /(dsqrt(3.d0)*xi(ixI^S)*(xi(ixI^S)+2)))

      !> Right state
      ff(ixI^S) = sum(wRC(ixI^S, Fr(:))**2, dim=ndim+1)/(m1_c*wRC(ixI^S,Er_))**2
      ff(ixI^S) = min(ff(ixI^S),1.d0)

      !> Eddington scalars
      fi(ixI^S) = abs(wRC(ixI^S,Fr(idim))/(m1_c*wRC(ixI^S,Er_)))
      fi(ixI^S) = min(fi(ixI^S),1.d0)

      fj(ixI^S) = dsqrt(max(ff(ixI^S) - fi(ixI^S)**2,0.d0))

      xi(ixI^S) = dsqrt(4-3*ff(ixI^S))

      lam_R_max(ixI^S) = m1_c*m1_RSOL*(fi(ixI^S)/xi(ixI^S) &
                  + dsqrt(2*(xi(ixI^S)-1)*(xi(ixI^S)+2)*(2*(xi(ixI^S)-1)*(xi(ixI^S)+2) + 3*fj(ixI^S)**2)) &
                  /(dsqrt(3.d0)*xi(ixI^S)*(xi(ixI^S)+2)))

      lam_R_min(ixI^S) = m1_c*m1_RSOL*(fi(ixI^S)/xi(ixI^S) &
                  - dsqrt(2*(xi(ixI^S)-1)*(xi(ixI^S)+2)*(2*(xi(ixI^S)-1)*(xi(ixI^S)+2) + 3*fj(ixI^S)**2)) &
                  /(dsqrt(3.d0)*xi(ixI^S)*(xi(ixI^S)+2)))

      !>Skinner Ostricker 2013, Eq 39
      !TODO: Is below correctly implemented?
      cmax(ixI^S,1) = max(lam_L_max(ixI^S),lam_R_max(ixI^S))
      cmax(ixI^S,1) = max(cmax(ixI^S,1),0.d0)

      if (present(cmin)) then
        cmin(ixI^S,1) = min(lam_L_min(ixI^S),lam_R_min(ixI^S))
        cmin(ixI^S,1) = min(cmin(ixI^S,1),0.d0)
      end if

    else

      select case(boundspeed)
      case (1)
        ! This implements formula (10.52) from "Riemann Solvers and Numerical
        ! Methods for Fluid Dynamics" by Toro.

        tmp1(ixO^S)=dsqrt(wLp(ixO^S,rho_))
        tmp2(ixO^S)=dsqrt(wRp(ixO^S,rho_))
        tmp3(ixO^S)=1.d0/(dsqrt(wLp(ixO^S,rho_))+dsqrt(wRp(ixO^S,rho_)))
        umean(ixO^S)=(wLp(ixO^S,mom(idim))*tmp1(ixO^S)+wRp(ixO^S,mom(idim))*tmp2(ixO^S))*tmp3(ixO^S)

        if(hd_energy) then
          csoundL(ixO^S)=hd_gamma*wLp(ixO^S,p_)/wLp(ixO^S,rho_)
          csoundR(ixO^S)=hd_gamma*wRp(ixO^S,p_)/wRp(ixO^S,rho_)
        else
          call hd_get_csound2(wLC,x,ixI^L,ixO^L,csoundL)
          call hd_get_csound2(wRC,x,ixI^L,ixO^L,csoundR)
        end if

        dmean(ixO^S) = (tmp1(ixO^S)*csoundL(ixO^S)+tmp2(ixO^S)*csoundR(ixO^S)) * &
             tmp3(ixO^S) + 0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2 * &
             (wRp(ixO^S,mom(idim))-wLp(ixO^S,mom(idim)))**2

        dmean(ixO^S)=dsqrt(dmean(ixO^S))
        if(present(cmin)) then
          cmin(ixO^S,1)=umean(ixO^S)-dmean(ixO^S)
          cmax(ixO^S,1)=umean(ixO^S)+dmean(ixO^S)
          if(H_correction) then
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
              cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
              cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
            {end do\}
          end if
        else
          cmax(ixO^S,1)=dabs(umean(ixO^S))+dmean(ixO^S)
        end if

        if (hd_dust) then
          wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
          call dust_get_cmax(wmean, x, ixI^L, ixO^L, idim, cmax, cmin)
        end if

      case (2)
        wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
        tmp1(ixO^S)=wmean(ixO^S,mom(idim))/wmean(ixO^S,rho_)
        call hd_get_csound2(wmean,x,ixI^L,ixO^L,csoundR)
        csoundR(ixO^S) = dsqrt(csoundR(ixO^S))

        if(present(cmin)) then
          cmax(ixO^S,1)=max(tmp1(ixO^S)+csoundR(ixO^S),zero)
          cmin(ixO^S,1)=min(tmp1(ixO^S)-csoundR(ixO^S),zero)
          if(H_correction) then
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
              cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
              cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
            {end do\}
          end if
        else
          cmax(ixO^S,1)=dabs(tmp1(ixO^S))+csoundR(ixO^S)
        end if

        if (hd_dust) then
          call dust_get_cmax(wmean, x, ixI^L, ixO^L, idim, cmax, cmin)
        end if
      case (3)
        ! Miyoshi 2005 JCP 208, 315 equation (67)
        if(hd_energy) then
          csoundL(ixO^S)=hd_gamma*wLp(ixO^S,p_)/wLp(ixO^S,rho_)
          csoundR(ixO^S)=hd_gamma*wRp(ixO^S,p_)/wRp(ixO^S,rho_)
        else
          call hd_get_csound2(wLC,x,ixI^L,ixO^L,csoundL)
          call hd_get_csound2(wRC,x,ixI^L,ixO^L,csoundR)
        end if
        csoundL(ixO^S)=max(dsqrt(csoundL(ixO^S)),dsqrt(csoundR(ixO^S)))
        if(present(cmin)) then
          cmin(ixO^S,1)=min(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))-csoundL(ixO^S)
          cmax(ixO^S,1)=max(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))+csoundL(ixO^S)
          if(H_correction) then
            {do ix^DB=ixOmin^DB,ixOmax^DB\}
              cmin(ix^D,1)=sign(one,cmin(ix^D,1))*max(abs(cmin(ix^D,1)),Hspeed(ix^D,1))
              cmax(ix^D,1)=sign(one,cmax(ix^D,1))*max(abs(cmax(ix^D,1)),Hspeed(ix^D,1))
            {end do\}
          end if
        else
          cmax(ixO^S,1)=max(wLp(ixO^S,mom(idim)),wRp(ixO^S,mom(idim)))+csoundL(ixO^S)
        end if
        if (hd_dust) then
          wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
          call dust_get_cmax(wmean, x, ixI^L, ixO^L, idim, cmax, cmin)
        end if
      end select

    end if

  end subroutine hd_get_cbounds

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine hd_get_csound2(w,x,ixI^L,ixO^L,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixI^S)

    call hd_get_pthermal(w,x,ixI^L,ixO^L,csound2)
    csound2(ixO^S)=hd_gamma*csound2(ixO^S)/w(ixO^S,rho_)

  end subroutine hd_get_csound2

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho) within ixO^L
  subroutine hd_get_pthermal(w, x, ixI^L, ixO^L, pth)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_pthermal
    use mod_small_values, only: trace_small_values

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: pth(ixI^S)
    integer                      :: iw, ix^D

    if (hd_energy) then
       pth(ixO^S) = (hd_gamma - 1.0d0) * (w(ixO^S, e_) - &
            hd_kin_en(w, ixI^L, ixO^L))
    else
       if (.not. associated(usr_set_pthermal)) then
          pth(ixO^S) = hd_adiab * w(ixO^S, rho_)**hd_gamma
       else
          call usr_set_pthermal(w,x,ixI^L,ixO^L,pth)
       end if
    end if

    if (fix_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(pth(ix^D)<small_pressure) then
            pth(ix^D)=small_pressure
         endif
      {enddo^D&\}
    else if (check_small_values) then
      {do ix^DB= ixO^LIM^DB\}
         if(pth(ix^D)<small_pressure) then
           write(*,*) "Error: small value of gas pressure",pth(ix^D),&
                " encountered when call hd_get_pthermal"
           write(*,*) "Iteration: ", it, " Time: ", global_time
           write(*,*) "Location: ", x(ix^D,:)
           write(*,*) "Cell number: ", ix^D
           do iw=1,nw
             write(*,*) trim(cons_wnames(iw)),": ",w(ix^D,iw)
           end do
           ! use erroneous arithmetic operation to crash the run
           if(trace_small_values) write(*,*) dsqrt(pth(ix^D)-bigdouble)
           write(*,*) "Saving status at the previous time step"
           crash=.true.
         end if
      {enddo^D&\}
    end if

  end subroutine hd_get_pthermal

  !> Calculate the radiation pressure tensor components
  !> using the M1 closure
  !> https://www.sciencedirect.com/science/article/pii/S0022407314001691
  subroutine m1_get_prad(w, x, ixI^L, ixO^L, idir, jdir, prad, str)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: prad(ixI^S)
    integer, intent(in)          :: idir, jdir

    double precision :: Erad(ixI^S)
    logical, optional, intent(in):: str
    logical :: streaming = .false.

    integer                      :: iw, ix^D
    double precision :: fi(ixI^S), fj(ixI^S), ff(ixI^S), chi(ixI^S)

    if (present(str)) streaming = str

    if (m1_diffusion) then

      if (idir == jdir) then
        prad(ixI^S) = 1.d0/3.d0
      else
        prad(ixI^S) = 0.d0
      end if

    else

      if (streaming) then
        call m1_get_streaming_fraction(w, x, ixI^L, ixO^L, idir, Erad)
        Erad(ixI^S) = w(ixI^S,Er_)*Erad(ixI^S)
      else
        Erad(ixI^S) = w(ixI^S,Er_)
      end if

      !> Eddington scalars
      fi(ixI^S) = w(ixI^S,Fr(idir))/(m1_c*Erad(ixI^S))
      fj(ixI^S) = w(ixI^S,Fr(jdir))/(m1_c*Erad(ixI^S))

      !TODO: Check this sum, it does not produce what is expected
      ff(ixI^S) = sum(w(ixI^S, Fr(:))**2, dim=ndim+1)/(m1_c*Erad(ixI^S))**2

      !> Cap ff, this should not be necessary, but the sum above does not behave properly
      fi(ixI^S) = min(fi(ixI^S),1.d0)
      fj(ixI^S) = min(fj(ixI^S),1.d0)
      ff(ixI^S) = min(ff(ixI^S),1.d0)

      !> Eddington factor
      chi(ixI^S) = (3.d0+4*ff(ixI^S))/(5.d0+2*dsqrt(4.d0-3*ff(ixI^S)))
      chi(ixI^S) = max(chi(ixI^S),1.d0/3)
      chi(ixI^S) = min(chi(ixI^S),1.d0)

      !> Temporary (??) fix
      where (chi(ixI^S) .ne. chi(ixI^S))
        chi(ixI^S) = 1.d0 !5.d0/3.d0
      end where

      !> Eddington tensor
      prad(ixI^S) = (3.d0*chi(ixI^S)-1)/2.d0 * fi(ixI^S)*fj(ixI^S)/ff(ixI^S)

      !> Nan when F == 0 as fi*fj / ff = 0/0
      !> Tensor term should be 0
      where (prad(ixI^S) .ne. prad(ixI^S))
        prad(ixI^S) = 0.d0
      end where

      if (idir == jdir) then
        prad(ixI^S) = max(0.d0, prad(ixI^S) + (1.d0-chi(ixI^S))/2.d0)
      end if

    end if

    !> Radiation pressure tensor
    prad(ixI^S) = prad(ixI^S) * Erad(ixI^S)

  end subroutine m1_get_prad

  subroutine m1_get_chi(w, x, ixI^L, ixO^L, chi)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: chi(ixI^S)

    !TODO: Check this sum, it does not produce what is expected
    chi(ixI^S) = sum(w(ixI^S, Fr(:))**2, dim=ndim+1)/(m1_c*w(ixI^S,Er_))**2

    !> Cap ff, this should not be necessary, but the sum above does not behave properly
    chi(ixI^S) = min(chi(ixI^S),1.d0)

    !> Eddington factor
    chi(ixI^S) = (3.d0+4*chi(ixI^S))/(5.d0+2*dsqrt(4.d0-3*chi(ixI^S)))
    chi(ixI^S) = max(chi(ixI^S),1.d0/3)
    chi(ixI^S) = min(chi(ixI^S),1.d0)

    !> Temporary (??) fix
    where (chi(ixI^S) .ne. chi(ixI^S))
      chi(ixI^S) = 1.d0 !5.d0/3.d0
    end where

    end subroutine m1_get_chi

    subroutine m1_get_trapped_fraction(w, x, ixI^L, ixO^L, idim, theta)
      use mod_global_parameters

      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, 1:nw)
      double precision, intent(in) :: x(ixI^S, 1:ndim)
      double precision, intent(out):: theta(ixI^S)
      integer, intent(in)          :: idim

      integer :: hxO^L,jxO^L

      double precision :: tau_c(ixI^S)

      ! !> Cell optical depth in idim
      ! ! tau_c(ixI^S) = block%dx(ixI^S,idim)*w(ixI^S,rho_)*w(ixI^S,kF_)
      tau_c(ixI^S) = block%dx(ixI^S,idim)*w(ixI^S,rho_)*m1_kappa_e
      !
      ! ! !> Minimum of its neighbours
      hxO^L=ixO^L-kr(idim,^D);
      jxO^L=ixO^L+kr(idim,^D);
      tau_c(ixO^S) = min(tau_c(ixO^S),tau_c(hxO^S),tau_c(jxO^S))
      ! !TODO: stencil minimum in first GC's

      ! tau_c(ixO^S) = min(w(hxO^S,rho_),w(ixO^S,rho_),w(jxO^S,rho_))*w(ixO^S,kF_)

      ! theta(ixI^S) = 2.d0/(2 + 3*tau_c(ixI^S))
      theta(ixI^S) = dexp(-2.d0/(3*tau_c(ixI^S)))

    end subroutine m1_get_trapped_fraction

    subroutine m1_get_streaming_fraction(w, x, ixI^L, ixO^L, idim, theta)
      use mod_global_parameters

      integer, intent(in)          :: ixI^L, ixO^L
      double precision, intent(in) :: w(ixI^S, 1:nw)
      double precision, intent(in) :: x(ixI^S, 1:ndim)
      double precision, intent(out):: theta(ixI^S)
      integer, intent(in)          :: idim

      call m1_get_trapped_fraction(w, x, ixI^L, ixO^L, idim, theta)
      theta(ixI^S) = 1.d0 - theta(ixI^S)

    end subroutine m1_get_streaming_fraction

  !> Calculate temperature=p/rho when in e_ the  total energy is stored
  subroutine hd_get_temperature_from_etot(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)

    double precision :: R(ixI^S)

    call hd_get_Rfactor(w,x,ixI^L,ixO^L,R)
    call hd_get_pthermal(w, x, ixI^L, ixO^L, res)
    res(ixO^S)=res(ixO^S)/(R(ixO^S)*w(ixO^S,rho_))
  end subroutine hd_get_temperature_from_etot

  !> Calculate temperature=p/rho when in e_ the  internal energy is stored
  subroutine hd_get_temperature_from_eint(w, x, ixI^L, ixO^L, res)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: res(ixI^S)

    double precision :: R(ixI^S)

    call hd_get_Rfactor(w,x,ixI^L,ixO^L,R)
    res(ixO^S) = (hd_gamma - 1.0d0) * w(ixO^S, e_)/(w(ixO^S,rho_)*R(ixO^S))
  end subroutine hd_get_temperature_from_eint

  ! ! Calculate flux f_idim[iw]
  ! subroutine hd_get_flux_cons(w, x, ixI^L, ixO^L, idim, f)
  !   use mod_global_parameters
  !   use mod_dust, only: dust_get_flux
  !
  !   integer, intent(in)             :: ixI^L, ixO^L, idim
  !   double precision, intent(in)    :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
  !   double precision, intent(out)   :: f(ixI^S, nwflux)
  !   double precision                :: prad_ii(ixI^S), pth(ixI^S), v(ixI^S),frame_vel(ixI^S)
  !   integer                         :: idir, itr
  !
  !   call hd_get_pthermal(w, x, ixI^L, ixO^L, pth)
  !   call hd_get_v_idim(w, x, ixI^L, ixO^L, idim, v)
  !
  !   f(ixO^S, rho_) = v(ixO^S) * w(ixO^S, rho_)
  !
  !   ! Momentum flux is v_i*m_i, +p in direction idim
  !   do idir = 1, ndir
  !      f(ixO^S, mom(idir)) = v(ixO^S) * w(ixO^S, mom(idir))
  !   end do
  !
  !   f(ixO^S, mom(idim)) = f(ixO^S, mom(idim)) + pth(ixO^S)
  !
  !   if(hd_energy) then
  !     ! Energy flux is v_i*(e + p)
  !     f(ixO^S, e_) = v(ixO^S) * (w(ixO^S, e_) + pth(ixO^S))
  !   end if
  !
  !   if (hd_m1) then
  !     !> Radiation energy flux is v_i E + F_i
  !     f(ixO^S, Er_) = v(ixO^S) * w(ixO^S, Er_) + w(ixO^S, Fr(idim))
  !
  !     !> Radiation flux flux is v_i F_i + c^2 P_ii
  !     do idir = 1,ndir
  !       call m1_get_prad(w, x, ixI^L, ixO^L, idir, idim, prad_ii)
  !       f(ixO^S, Fr(idir)) = v(ixO^S) * w(ixO^S,Fr(idir)) + m1_RSOL * m1_c**2 * prad_ii(ixO^S)
  !     end do
  !   end if
  !
  !   do itr = 1, hd_n_tracer
  !      f(ixO^S, tracer(itr)) = v(ixO^S) * w(ixO^S, tracer(itr))
  !   end do
  !
  !   ! Dust fluxes
  !   if (hd_dust) then
  !     call dust_get_flux(w, x, ixI^L, ixO^L, idim, f)
  !   end if
  !
  ! end subroutine hd_get_flux_cons

  ! Calculate flux f_idim[iw]
  subroutine hd_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
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
    double precision                :: prad_ii(ixI^S), pth(ixI^S),frame_vel(ixI^S)
    integer                         :: idir, itr

    if (hd_energy) then
       pth(ixO^S) = w(ixO^S,p_)
    else
       call hd_get_pthermal(wC, x, ixI^L, ixO^L, pth)
    end if

    f(ixO^S, rho_) = w(ixO^S,mom(idim)) * w(ixO^S, rho_)

    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
       f(ixO^S, mom(idir)) = w(ixO^S,mom(idim)) * wC(ixO^S, mom(idir))
    end do

    f(ixO^S, mom(idim)) = f(ixO^S, mom(idim)) + pth(ixO^S)

    if(hd_energy) then
      ! Energy flux is v_i*(e + p)
      f(ixO^S, e_) = w(ixO^S,mom(idim)) * (wC(ixO^S, e_) + w(ixO^S,p_))
    end if

    if (hd_m1) then
      !> Radiation energy flux is v_i E + F
      f(ixO^S,Er_) = m1_RSOL * (w(ixO^S,mom(idim)) * w(ixO^S, Er_) + w(ixO^S, Fr(idim)))

      !> Radiation flux flux is v_i F_i + c^2 P_ii
      do idir = 1,ndir
        call m1_get_prad(w, x, ixI^L, ixO^L, idir, idim, prad_ii, .true.)
        f(ixO^S,Fr(idir)) = m1_RSOL * (w(ixO^S,mom(idim)) * w(ixO^S,Fr(idir)) + m1_c**2 * prad_ii(ixO^S))
      end do
    end if

    do itr = 1, hd_n_tracer
       f(ixO^S, tracer(itr)) = w(ixO^S,mom(idim)) * w(ixO^S, tracer(itr))
    end do

    ! Dust fluxes
    if (hd_dust) then
      call dust_get_flux_prim(w, x, ixI^L, ixO^L, idim, f)
    end if

    ! Viscosity fluxes - viscInDiv
    if (hd_viscosity) then
      call visc_get_flux_prim(w, x, ixI^L, ixO^L, idim, f, hd_energy)
    endif

  end subroutine hd_get_flux

  !> Add geometrical source terms to w
  !>
  !> Notice that the expressions of the geometrical terms depend only on ndir,
  !> not ndim. Eg, they are the same in 2.5D and in 3D, for any geometry.
  !>
  !> Ileyk : to do :
  !>     - address the source term for the dust in case (coordinate == spherical)
  subroutine hd_add_source_geom(qdt, dtfactor, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters
    use mod_usr_methods, only: usr_set_surface
    use mod_viscosity, only: visc_add_source_geom ! viscInDiv
    use mod_rotating_frame, only: rotating_frame_add_source
    use mod_dust, only: dust_n_species, dust_mom, dust_rho
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)
    ! to change and to set as a parameter in the parfile once the possibility to
    ! solve the equations in an angular momentum conserving form has been
    ! implemented (change tvdlf.t eg)
    double precision :: pth(ixI^S), source(ixI^S), minrho, prad(ixI^S)
    integer                         :: iw,idir, h1x^L{^NOONED, h2x^L}
    integer :: mr_,mphi_,Fr_,Fphi_ ! Polar var. names
    integer :: irho, ifluid, n_fluids
    double precision :: exp_factor(ixI^S), del_exp_factor(ixI^S), exp_factor_primitive(ixI^S)

    if (hd_dust) then
       n_fluids = 1 + dust_n_species
    else
       n_fluids = 1
    end if

    select case (coordinate)

    case(Cartesian_expansion)
      !the user provides the functions of exp_factor and del_exp_factor
      if(associated(usr_set_surface)) call usr_set_surface(ixI^L,x,block%dx,exp_factor,del_exp_factor,exp_factor_primitive)
      call hd_get_pthermal(wCT, x, ixI^L, ixO^L, source)
      source(ixO^S) = source(ixO^S)*del_exp_factor(ixO^S)/exp_factor(ixO^S)
      w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt*source(ixO^S)

    case (cylindrical)
      if (hd_m1) then
         call mpistop("m1 geom source terms not implemented yet with cylindrical geometries")
      end if
       do ifluid = 0, n_fluids-1
          ! s[mr]=(pthermal+mphi**2/rho)/radius
          if (ifluid == 0) then
             ! gas
             irho  = rho_
             mr_   = mom(r_)
             if(phi_>0) mphi_ = mom(phi_)
             call hd_get_pthermal(wCT, x, ixI^L, ixO^L, source)
             minrho = 0.0d0
          else
             ! dust : no pressure
             irho  = dust_rho(ifluid)
             mr_   = dust_mom(r_, ifluid)
             if(phi_>0) mphi_ = dust_mom(phi_, ifluid)
             source(ixI^S) = zero
             minrho = 0.0d0
          end if
          if (phi_ > 0) then
             where (wCT(ixO^S, irho) > minrho)
                source(ixO^S) = source(ixO^S) + wCT(ixO^S, mphi_)**2 / wCT(ixO^S, irho)
                w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * source(ixO^S) / x(ixO^S, r_)
             end where
             ! s[mphi]=(-mphi*mr/rho)/radius
             where (wCT(ixO^S, irho) > minrho)
                source(ixO^S) = -wCT(ixO^S, mphi_) * wCT(ixO^S, mr_) / wCT(ixO^S, irho)
                w(ixO^S, mphi_) = w(ixO^S, mphi_) + qdt * source(ixO^S) / x(ixO^S, r_)
             end where
          else
             ! s[mr]=2pthermal/radius
             w(ixO^S, mr_) = w(ixO^S, mr_) + qdt * source(ixO^S) / x(ixO^S, r_)
          end if
       end do
    case (spherical)
       if (hd_dust) then
          call mpistop("Dust geom source terms not implemented yet with spherical geometries")
       end if
       !> momentum
       mr_   = mom(r_)
       if(phi_>0) mphi_ = mom(phi_)
       h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
       ! s[mr]=((mtheta**2+mphi**2)/rho+2*p)/r
       call hd_get_pthermal(wCT, x, ixI^L, ixO^L, pth)
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
       source(ixO^S) = source(ixO^S) - (wCT(ixO^S, mom(2)) * wCT(ixO^S, mr_)) / wCT(ixO^S, rho_)
       w(ixO^S, mom(2)) = w(ixO^S, mom(2)) + qdt * source(ixO^S) / x(ixO^S, 1)

       if (ndir == 3) then
         ! s[mphi]=-(mphi*mr/rho)/r-cot(theta)*(mtheta*mphi/rho)/r
         source(ixO^S) = -(wCT(ixO^S, mom(3)) * wCT(ixO^S, mr_)) / wCT(ixO^S, rho_)&
                        - (wCT(ixO^S, mom(2)) * wCT(ixO^S, mom(3))) / wCT(ixO^S, rho_) / tan(x(ixO^S, 2))
         w(ixO^S, mom(3)) = w(ixO^S, mom(3)) + qdt * source(ixO^S) / x(ixO^S, 1)
       end if
       }

    end select

    if (hd_m1) call m1_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)

    if (hd_viscosity) call visc_add_source_geom(qdt,ixI^L,ixO^L,wCT,w,x)

    if (hd_rotating_frame) then
       if (hd_dust) then
          call mpistop("Rotating frame not implemented yet with dust")
       else
          call rotating_frame_add_source(qdt,dtfactor,ixI^L,ixO^L,wCT,w,x)
       end if
    end if

  end subroutine hd_add_source_geom

  subroutine m1_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim)
    double precision, intent(inout) :: wCT(ixI^S, 1:nw), w(ixI^S, 1:nw)

    double precision                :: chi(ixI^S), ff(ixI^S), fi(ixI^S)
    double precision                :: source(ixI^S)
    double precision                :: prr(ixI^S), ptt(ixI^S), ppp(ixI^S)
    double precision                :: ptr(ixI^S)
    double precision                :: ppr(ixI^S), ppt(ixI^S)
    integer                         :: h1x^L{^NOONED, h2x^L}

    !> Radiation flux
    !TODO: Fix geometric sourceterms for non-Cartesian multi-D
    !SEE APPENDIX B IN GONZALEZ 2007 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !if (any(w(ixO^S,Er_) .lt. 0.d0)) then
    !  stop "Er zero before m1_geom_source"
    !end if

    !if (any(w(ixO^S,Er_) .ne. w(ixO^S,Er_))) then
    !  stop "Er NaN before m1_geom_source"
    !end if

    !if (any(w(ixO^S,Fr(1)) .ne. w(ixO^S,Fr(1)))) then
    !  stop "Fr(1) NaN before m1_geom_source"
    !end if


    select case (coordinate)
    case (cylindrical,Cartesian_expansion)
      call mpistop("m1 geom source terms not implemented yet with cylindrical geometries")
    case(spherical)
      h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
      select case (ndir)
      case (1)
       !> Get relevant pressure tensor
       call m1_get_prad(wCT, x, ixI^L, ixO^L, 1, 1, prr, .true.)
       ptt = m1_c**2 * (wCT(ixI^S, Er_) - prr) /2
       ppp = ptt

       w(ixI^S, Fr(1)) = w(ixI^S, Fr(1)) + qdt * m1_RSOL * (ptt(ixI^S) + ppp(ixI^S))/ x(ixI^S, 1)

      {^NOONED
      case (2)
        ! if (ndim .eq. 1) call mpistop("1.x D not implemented yet")

        !> Get relevant pressure tensors:
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 1, 1, prr)
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 2, 2, ptt)
        ppp = wCT(ixI^S,Er_) - prr - ptt !> trace is one
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 2, 1, ptr)

        !> Add advection part
        ! prr = m1_c**2 * prr + wCT(ixI^S, Fr(1))*wCT(ixI^S, mom(1))/wCT(ixI^S, rho_)
        ptt = m1_c**2 * ptt + wCT(ixI^S, Fr(2))*wCT(ixI^S, mom(2))/wCT(ixI^S, rho_)
        ptr = m1_c**2 * ptr + wCT(ixI^S, Fr(2))*wCT(ixI^S, mom(1))/wCT(ixI^S, rho_) !< IS THIS ONE CORRECT?
        ! ptr = m1_c**2 * ptr + wCT(ixI^S, Fr(1))*wCT(ixI^S, mom(2))/wCT(ixI^S, rho_) !< IS THIS ONE CORRECT?
        ppp = m1_c**2 * ppp

        !> Perform corrections
        !> SWAP SIGNS AS ACCORDING TO WIKI: F = F - dt*(nabla T + S_corr)
        w(ixO^S, Fr(1)) = w(ixO^S, Fr(1)) + qdt * m1_RSOL * (ptt(ixO^S) + ppp(ixO^S))/ x(ixO^S, 1)
        w(ixO^S, Fr(2)) = w(ixO^S, Fr(2)) + qdt * m1_RSOL * (-ptr(ixO^S)/x(ixO^S, 1)  + ppp(ixO^S)/(x(ixO^S, 1)*dtan(x(ixO^S,2))))

      case(3)
        ! if ((ndim .eq. 1) .or. (ndim .eq. 2)) call mpistop("1.x D not implemented yet")

        !> Get relevant pressure tensors:
        ! call m1_get_prad(wCT, x, ixI^L, ixO^L, 1, 1, prr)
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 2, 2, ptt)
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 3, 3, ppp)
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 2, 1, ptr)
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 3, 2, ppt)
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 3, 1, ppr)

        !> Add advection part
        ! prr = m1_c**2 * prr + wCT(ixI^S, Fr(1))*wCT(ixI^S, mom(1))/wCT(ixI^S, rho_)
        ptt = m1_c**2 * ptt + wCT(ixI^S, Fr(2))*wCT(ixI^S, mom(2))/wCT(ixI^S, rho_)
        ppp = m1_c**2 * ppp + wCT(ixI^S, Fr(3))*wCT(ixI^S, mom(3))/wCT(ixI^S, rho_)

        ptr = m1_c**2 * ptr + wCT(ixI^S, Fr(2))*wCT(ixI^S, mom(1))/wCT(ixI^S, rho_) !< IS THIS ONE CORRECT?
        ppr = m1_c**2 * ppr + wCT(ixI^S, Fr(3))*wCT(ixI^S, mom(1))/wCT(ixI^S, rho_) !< IS THIS ONE CORRECT?
        ppt = m1_c**2 * ppt + wCT(ixI^S, Fr(3))*wCT(ixI^S, mom(2))/wCT(ixI^S, rho_) !< IS THIS ONE CORRECT?

        !> Perform corrections
        !> SWAP SIGNS AS ACCORDING TO WIKI: F = F - dt*nabla T
        w(ixO^S, Fr(1)) = w(ixO^S, Fr(1)) + qdt * m1_RSOL * (ptt(ixO^S) + ppp(ixO^S))/ x(ixO^S, 1)
        w(ixO^S, Fr(2)) = w(ixO^S, Fr(2)) + qdt * m1_RSOL * (-ptr(ixO^S)/x(ixO^S, 1)  + ppp(ixO^S)/(x(ixO^S, 1)*dtan(x(ixO^S,2))))

        !> 2 ptp BELOW: P is symmetric
        w(ixO^S, Fr(3)) = w(ixO^S, Fr(3)) + qdt * m1_RSOL * (-ppr(ixO^S)/x(ixO^S, 1) - 2*ppt(ixO^S)/(x(ixO^S, 1)*dtan(x(ixO^S,2))))
      }

      end select

    end select

!    if (any(w(ixO^S,Er_) .lt. 0.d0)) then
!      stop "Er zero after m1_geom_source"
!    end if
!
!    if (any(w(ixO^S,Er_) .ne. w(ixO^S,Er_))) then
!      stop "Er NaN after m1_geom_source"
!    end if
!
!    if (any(w(ixO^S,Fr(1)) .ne. w(ixO^S,Fr(1)))) then
!      print*, source(ixO^S)
!      stop "Fr(1) NaN after m1_geom_source"
!    end if


  end subroutine m1_add_source_geom

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine hd_add_source(qdt,dtfactor, ixI^L,ixO^L,wCT,wCTprim,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_dust, only: dust_add_source, dust_mom, dust_rho, dust_n_species
    use mod_viscosity, only: viscosity_add_source
    use mod_usr_methods, only: usr_gravity
    use mod_gravity, only: gravity_add_source, grav_split
    use mod_cak_force, only: cak_add_source

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, dtfactor
    double precision, intent(in)    :: wCT(ixI^S, 1:nw),wCTprim(ixI^S,1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    double precision :: gravity_field(ixI^S, 1:ndim)
    integer :: idust, idim

    if(hd_dust .and. .not. use_imex_scheme) then
      call dust_add_source(qdt,ixI^L,ixO^L,wCT,w,x,qsourcesplit,active)
    end if

    if(hd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,&
           qsourcesplit,active, rc_fl)
    end if

    if(hd_viscosity) then
      call viscosity_add_source(qdt,ixI^L,ixO^L,wCT,w,x,&
           hd_energy,qsourcesplit,active)
    end if

    if (hd_gravity) then
      call gravity_add_source(qdt,ixI^L,ixO^L,wCT,wCTprim,w,x,&
           hd_energy,.false.,qsourcesplit,active)

      if (hd_dust .and. qsourcesplit .eqv. grav_split) then
         active = .true.

         call usr_gravity(ixI^L, ixO^L, wCT, x, gravity_field)
         do idust = 1, dust_n_species
            do idim = 1, ndim
               w(ixO^S, dust_mom(idim, idust)) = w(ixO^S, dust_mom(idim, idust)) &
                    + qdt * gravity_field(ixO^S, idim) * wCT(ixO^S, dust_rho(idust))
            end do
         end do
      end if
    end if

    if (hd_cak_force) then
      call cak_add_source(qdt,ixI^L,ixO^L,wCT,w,x,hd_energy,qsourcesplit,active)
    end if

    if(hd_partial_ionization) then
      if(.not.qsourcesplit) then
        active = .true.
        call hd_update_temperature(ixI^L,ixO^L,wCT,w,x)
      end if
    end if

    if (hd_m1) then
      if (m1_force) call m1_hydro_source(qdt,ixI^L,ixO^L,wCT,w,x,hd_energy,qsourcesplit,active)
      call m1_rad_source(qdt,ixI^L,ixO^L,wCT,w,x,hd_energy,qsourcesplit,active)
    end if

  end subroutine hd_add_source

  subroutine hd_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    use mod_global_parameters
    use mod_dust, only: dust_get_dt
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt
    use mod_cak_force, only: cak_get_dt

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S, 1:^ND)
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

    if(hd_dust) then
      call dust_get_dt(w, ixI^L, ixO^L, dtnew, dx^D, x)
    end if

    if(hd_radiative_cooling) then
      call cooling_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x,rc_fl)
    end if

    if(hd_viscosity) then
      call viscosity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    end if

    if(hd_gravity) then
      call gravity_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
   end if

   if (hd_cak_force) then
     call cak_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
   end if

   if (hd_m1) then
     call m1_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
   end if

  end subroutine hd_get_dt

  function hd_kin_en(w, ixI^L, ixO^L, inv_rho) result(ke)
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
  end function hd_kin_en

  function hd_inv_rho(w, ixI^L, ixO^L) result(inv_rho)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S, nw)
    double precision              :: inv_rho(ixO^S)

    ! Can make this more robust
    inv_rho = 1.0d0 / w(ixO^S, rho_)
  end function hd_inv_rho

  subroutine hd_handle_small_values(primitive, w, x, ixI^L, ixO^L, subname)
    ! handles hydro (density,pressure,velocity) bootstrapping
    ! any negative dust density is flagged as well (and throws an error)
    ! small_values_method=replace also for dust
    use mod_global_parameters
    use mod_small_values
    use mod_dust, only: dust_n_species, dust_mom, dust_rho
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    double precision :: Fmag(ixI^S)

    integer :: n,idir
    logical :: flag(ixI^S,1:nw)

    call hd_check_w(primitive, ixI^L, ixO^L, w, flag)

    if (any(flag)) then
      select case (small_values_method)
      case ("replace")
        where(flag(ixO^S,rho_)) w(ixO^S,rho_) = small_density
        do idir = 1, ndir
          if(small_values_fix_iw(mom(idir))) then
            where(flag(ixO^S,rho_)) w(ixO^S, mom(idir)) = 0.0d0
          end if
        end do
        if(hd_energy)then
          if(small_values_fix_iw(e_)) then
            if(primitive) then
              where(flag(ixO^S,rho_)) w(ixO^S, p_) = small_pressure
            else
              where(flag(ixO^S,rho_)) w(ixO^S, e_) = small_e + hd_kin_en(w,ixI^L,ixO^L)
            endif
          end if
        endif

        if(hd_energy) then
          if(primitive) then
            where(flag(ixO^S,e_)) w(ixO^S,p_) = small_pressure
          else
            where(flag(ixO^S,e_))
              ! Add kinetic energy
              w(ixO^S,e_) = small_e + hd_kin_en(w,ixI^L,ixO^L)
            end where
          end if
        end if

        if(hd_dust)then
           do n=1,dust_n_species
              where(flag(ixO^S,dust_rho(n))) w(ixO^S,dust_rho(n)) = 0.0d0
              do idir = 1, ndir
                  where(flag(ixO^S,dust_rho(n))) w(ixO^S,dust_mom(idir,n)) = 0.0d0
              enddo
           enddo
        endif

        !> Radiation flux
        if (hd_m1) then

          where(flag(ixO^S,Er_)) w(ixO^S, Er_) = small_e

          if (any(flag(ixO^S,Fr(:)))) &
            Fmag(ixI^S) = dsqrt(sum(w(ixI^S,Fr(:))**2, dim=ndim+1))

          do idir = 1,ndir
            where(flag(ixI^S,Fr(idir))) w(ixI^S,Fr(idir)) = w(ixI^S,Fr(idir))*m1_c*w(ixI^S,Er_)/Fmag(ixI^S)
          end do
        end if


      case ("average")
        if(primitive)then
           ! averaging for all primitive fields, including dust
           call small_values_average(ixI^L, ixO^L, w, x, flag)
        else
           ! do averaging of density
           call small_values_average(ixI^L, ixO^L, w, x, flag, rho_)
           if(hd_energy) then
             ! do averaging of pressure
             w(ixI^S,p_)=(hd_gamma-1.d0)*(w(ixI^S,e_) &
              -0.5d0*sum(w(ixI^S, mom(:))**2, dim=ndim+1)/w(ixI^S,rho_))
             call small_values_average(ixI^L, ixO^L, w, x, flag, p_)
             w(ixI^S,e_)=w(ixI^S,p_)/(hd_gamma-1.d0) &
               +0.5d0*sum(w(ixI^S, mom(:))**2, dim=ndim+1)/w(ixI^S,rho_)
           end if
           if(hd_dust)then
              do n=1,dust_n_species
                 where(flag(ixO^S,dust_rho(n))) w(ixO^S,dust_rho(n)) = 0.0d0
                 do idir = 1, ndir
                    where(flag(ixO^S,dust_rho(n))) w(ixO^S,dust_mom(idir,n)) = 0.0d0
                 enddo
              enddo
          endif
        endif

      case default
        if(.not.primitive) then
          !convert w to primitive
          ! Calculate pressure = (gamma-1) * (e-ek)
          if(hd_energy) then
            w(ixO^S,p_)=(hd_gamma-1.d0)*(w(ixO^S,e_)-hd_kin_en(w,ixI^L,ixO^L))
          end if
          ! Convert gas momentum to velocity
          do idir = 1, ndir
            w(ixO^S, mom(idir)) = w(ixO^S, mom(idir))/w(ixO^S,rho_)
          end do
        end if
        ! NOTE: dust entries may still have conserved values here
        call small_values_error(w, x, ixI^L, ixO^L, flag, subname)
      end select
    end if

  end subroutine hd_handle_small_values

  subroutine Rfactor_from_temperature_ionization(w,x,ixI^L,ixO^L,Rfactor)
    use mod_global_parameters
    use mod_ionization_degree
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: Rfactor(ixI^S)

    double precision :: iz_H(ixO^S),iz_He(ixO^S)

    call ionization_degree_from_temperature(ixI^L,ixO^L,w(ixI^S,Te_),iz_H,iz_He)
    ! assume the first and second ionization of Helium have the same degree
    Rfactor(ixO^S)=(1.d0+iz_H(ixO^S)+0.1d0*(1.d0+iz_He(ixO^S)*(1.d0+iz_He(ixO^S))))/2.3d0

  end subroutine Rfactor_from_temperature_ionization

  subroutine Rfactor_from_constant_ionization(w,x,ixI^L,ixO^L,Rfactor)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: Rfactor(ixI^S)

    Rfactor(ixO^S)=RR

  end subroutine Rfactor_from_constant_ionization

  subroutine hd_update_temperature(ixI^L,ixO^L,wCT,w,x)
    use mod_global_parameters
    use mod_ionization_degree

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: iz_H(ixO^S),iz_He(ixO^S), pth(ixI^S)

    call ionization_degree_from_temperature(ixI^L,ixO^L,wCT(ixI^S,Te_),iz_H,iz_He)

    call hd_get_pthermal(w,x,ixI^L,ixO^L,pth)

    w(ixO^S,Te_)=(2.d0+3.d0*He_abundance)*pth(ixO^S)/(w(ixO^S,rho_)*(1.d0+iz_H(ixO^S)+&
     He_abundance*(iz_He(ixO^S)*(iz_He(ixO^S)+1.d0)+1.d0)))

  end subroutine hd_update_temperature

  subroutine m1_rad_source(qdt,ixI^L,ixO^L,wCT,w,x,energy,qsourcesplit,active)
    use mod_global_parameters
    use mod_geometry
    use mod_comm_lib, only: mpistop

    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)    :: energy, qsourcesplit
    logical, intent(inout) :: active

    integer :: idir, idim
    integer                         :: jxO^L, hxO^L

    double precision :: prad(ixI^S), gradv(ixI^S)
    double precision :: Fp(ixO^S), Fm(ixO^S), sm(ixO^S), sp(ixO^S)
    double precision :: geom_corr(ixO^S), vec_grad(ixO^S), vel(ixI^S,1:ndir)

!    if (any(w(ixO^S,Er_) .lt. 0.d0)) then
!      stop "Er zero before m1_rad_source"
!    end if
!
!    if (any(w(ixO^S,Er_) .ne. w(ixO^S,Er_))) then
!      stop "Er NaN before m1_rad_source"
!    end if
!
!    if (any(w(ixO^S,Fr(1)) .ne. w(ixO^S,Fr(1)))) then
!      stop "Fr(1) NaN before m1_rad_source"
!    end if
!

    if(qsourcesplit .eqv. m1_force_split) then
      active = .true.

      do idir = 1,ndir
        vel(ixI^S,idir) = wCT(ixI^S,mom(idir))/wCT(ixI^S,rho_)
      end do

      do idir = 1,ndir

        !> Photon tiring + Radflux CMF correction
        do idim = 1,ndim

          !TODO: implement vector gradient corrections in spherical!
          !TODO: implement directional derivative in spherical!
          !> Both are expressed at: https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
          call gradient_str(vel(ixI^S,idim), x, ixI^L, ixO^L, idir, gradv)
          call m1_get_prad(wCT, x, ixI^L, ixO^L, idir, idim, prad)

          !> Geometric correction terms to vector gradient
          vec_grad(ixO^S) = gradv(ixO^S)

          select case(coordinate)
          case(cylindrical)
            call mpistop("vector gradient not implemented for cylindrical coords")
          case(spherical)
            if (idim .eq. 1) then
            {^NOONED
              if (idir .eq. 2) then
                vec_grad(ixO^S) = vec_grad(ixO^S) - vel(ixO^S,2)/x(ixO^S,1)
              else if (idir .eq. 3) then
                vec_grad(ixO^S) = vec_grad(ixO^S) -  vel(ixO^S,3)/x(ixO^S,1)
              end if
            else if (idim .eq. 2) then
              if (idir .eq. 2) then
                vec_grad(ixO^S) = vec_grad(ixO^S) +  vel(ixO^S,1)/x(ixO^S,1)
              else if (idir .eq. 3) then
                vec_grad(ixO^S) = vec_grad(ixO^S) -  vel(ixO^S,3)/(x(ixO^S,1)*dtan(x(ixO^S,2)))
              end if
            }
            {^IFTHREED
            else if (idim .eq. 3) then
              if (idir .eq. 3) then
                vec_grad(ixO^S) = vec_grad(ixO^S) +  vel(ixO^S,1)/x(ixO^S,1)
                vec_grad(ixO^S) = vec_grad(ixO^S) +  vel(ixO^S,2)/(x(ixO^S,1)*dtan(x(ixO^S,2)))
              end if
            }
            end if
          end select

          !TODO: Well balanced source-terms here? Requires gradient in GC
          w(ixO^S,Er_) = w(ixO^S,Er_) - qdt*m1_RSOL*vec_grad(ixO^S)*prad(ixO^S)
          w(ixO^S,Fr(idim)) = w(ixO^S,Fr(idim)) - qdt*m1_RSOL*wCT(ixO^S,Fr(idir))*gradv(ixO^S)
        end do

        !> Geometric correction terms to directional derivative
        select case(coordinate)
        case(cylindrical)
          call mpistop("directional derivative not implemented for cylindrical coords")
        case(spherical)
          select case(idir)
          {^NOONED
          case(1)
            geom_corr(ixO^S) = -(wCT(ixO^S,Fr(2)) * vel(ixO^S,2))/x(ixO^S,1)
            if (ndir > 2) then
              geom_corr(ixO^S) = geom_corr(ixO^S) - (wCT(ixO^S,Fr(3))* vel(ixO^S,3))/x(ixO^S,1)
            end if
          case(2)
            geom_corr(ixO^S) = (wCT(ixO^S,Fr(2)) * vel(ixO^S,1))/x(ixO^S,1)
            if (ndir > 2) then
              geom_corr(ixO^S) = geom_corr(ixO^S) - (wCT(ixO^S,Fr(3))* vel(ixO^S,3))/(x(ixO^S,1)*dtan(x(ixO^S,2)))
            end if
          }
          {^IFTHREED
          case(3)
            geom_corr(ixO^S) = (wCT(ixO^S,Fr(3)) * vel(ixO^S,1))/x(ixO^S,1)
            geom_corr(ixO^S) = geom_corr(ixO^S) + (wCT(ixO^S,Fr(3)) * vel(ixO^S,2))/(x(ixO^S,1)*dtan(x(ixO^S,2)))
          }
          end select
          {^NOONED
          !> add correction
          w(ixO^S,Fr(idir)) = w(ixO^S,Fr(idir)) + qdt*m1_RSOL*geom_corr(ixO^S)
          }
        end select

      end do

    end if

!    if (any(w(ixO^S,Er_) .lt. 0.d0)) then
!      stop "Er zero in m1_rad_source"
!    end if
!
!    if (any(w(ixO^S,Er_) .ne. w(ixO^S,Er_))) then
!      stop "Er NaN in m1_rad_source"
!    end if
!
!    if (any(w(ixO^S,Fr(1)) .ne. w(ixO^S,Fr(1)))) then
!      stop "Fr(1) NaN after m1_rad_source"
!    end if

     if (fix_small_values) then
       call hd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'rad_source')
     end if


  end subroutine m1_rad_source

  subroutine m1_hydro_source(qdt,ixI^L,ixO^L,wCT,w,x,energy,qsourcesplit,active)
    use mod_global_parameters
    use mod_geometry
    use mod_comm_lib, only: mpistop

    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)    :: energy, qsourcesplit
    logical, intent(inout) :: active

    integer :: idir
    integer                         :: jxO^L, hxO^L
    double precision :: Fp(ixO^S), Fm(ixO^S), sm(ixO^S), sp(ixO^S)


!    if (any(w(ixO^S,e_) .lt. 0.d0)) then
!      stop "e zero before m1_hydro_source 1"
!    end if
!
!    if (any(w(ixO^S,e_) .ne. w(ixO^S,e_))) then
!      stop "e NaN before m1_hydro_source 1"
!    end if

    if(qsourcesplit .eqv. m1_force_split) then
      active = .true.

      do idir = 1,ndir
        !> Force in momentum
        ! TODO: Well-balanced scheme?!
        ! w(ixO^S,mom(idir)) = w(ixO^S,mom(idir)) + qdt*wCT(ixO^S,rho_)*wCT(ixO^S,kF_)*wCT(ixO^S,Fr(idir))/m1_c
        !> Well-balanced formulation, Bloch+ 2020
        ! hxO^L=ixO^L-kr(idir,^D);
        ! jxO^L=ixO^L+kr(idir,^D);

        !> Flux face averages
        ! sm(ixO^S) = 1.d0/2.d0*(wCT(hxO^S,rho_)*wCT(hxO^S,kF_) + wCT(ixO^S,rho_)*wCT(ixO^S,kF_))
        ! sp(ixO^S) = 1.d0/2.d0*(wCT(jxO^S,rho_)*wCT(jxO^S,kF_) + wCT(ixO^S,rho_)*wCT(ixO^S,kF_))
        ! Fm(ixO^S) = 1.d0/2.d0*(wCT(hxO^S,Fr(idir)) + wCT(ixO^S,Fr(idir)))
        ! Fp(ixO^S) = 1.d0/2.d0*(wCT(jxO^S,Fr(idir)) + wCT(ixO^S,Fr(idir)))

        !> Update momentum
        w(ixO^S,mom(idir)) = w(ixO^S,mom(idir)) + qdt*wCT(ixO^S,rho_)*wCT(ixO^S,kF_)*wCT(ixO^S,Fr(idir))/m1_c
        ! w(ixO^S,mom(idir)) = w(ixO^S,mom(idir)) + qdt/2*(sm(ixO^S)*Fm(ixO^S) + sp(ixO^S)*Fp(ixO^S))/m1_c

        !> Force work term
        if (energy) &
          w(ixO^S,e_) = w(ixO^S,e_) + qdt*wCT(ixO^S,mom(idir))*wCT(ixO^S,kF_)*wCT(ixO^S,Fr(idir))/m1_c

       ! !> Well balanced scheme
       ! ! if (energy) &
       !   ! w(ixO^S,e_) = w(ixO^S,e_) + qdt*wCT(ixO^S,mom(idir))/wCT(ixO^S,rho_)*(sm(ixO^S)*Fm(ixO^S) + sp(ixO^S)*Fp(ixO^S))/(2*m1_c)

       !   if (any(w(ixO^S,e_) .lt. 0.d0)) then
       !     ! print*, 'e', w(ixO^S,e_)
       !     ! print*, 'p', w(ixO^S,e_) - w(ixO^S,mom(idir))**2/(2*w(ixO^S,rho_))
       !     ! print*, 'S', qdt*wCT(ixO^S,mom(idir))*wCT(ixO^S,kF_)*wCT(ixO^S,Fr(idir))/m1_c
       !     ! print*, 'kappa', wCT(ixO^S,kF_)
       !     ! print*, 'F', wCT(ixO^S,Fr(idir))
       !     ! print*, 'mom', wCT(ixO^S,mom(idir))
       !     stop "e zero after m1_hydro_source 1"
       !   end if

       !   if (any(w(ixO^S,e_) .ne. w(ixO^S,e_))) then
       !     ! print*, 'e'
       !     ! print*, wCT(ixO^S,e_)
       !     ! print*, 'A', qdt*wCT(ixO^S,mom(idir))/wCT(ixO^S,rho_)/2
       !     ! print*, 'B',(sm(ixO^S)*Fm(ixO^S)+ sp(ixO^S)*Fp(ixO^S))/m1_c
       !     ! print*, 'B1', sm(ixO^S)*Fm(ixO^S)
       !     ! print*, 'B2', sp(ixO^S)*Fp(ixO^S)
       !     ! print*, 'F', wCT(ixO^S,Fr(idir))
       !     ! print*, 'k', wCT(ixO^S,kF_)
       !     ! print*, 'Se'
       !     ! print*,  qdt*wCT(ixO^S,mom(idir))/wCT(ixO^S,rho_)/2*(sm(ixO^S)*Fm(ixO^S) &
       !     ! + sp(ixO^S)*Fp(ixO^S))/m1_c

       !     stop "e NaN after m1_hydro_source 1"
       !   end if

      end do
    end if

   if (fix_small_values) then
     call hd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'Hydro_source')
   end if


  end subroutine m1_hydro_source

  subroutine m1_implicit_update(dtfactor,qdt,qtC,psa,psb)
   use mod_global_parameters
   type(state), target :: psa(max_blocks)
   type(state), target :: psb(max_blocks)
   double precision, intent(in) :: qdt
   double precision, intent(in) :: qtC
   double precision, intent(in) :: dtfactor

   integer :: ixO^L,iigrid,igrid

   ixO^L=ixG^LL^LSUBnghostcells;

   if (m1_heatcool .and. hd_energy) then
     !$OMP PARALLEL DO PRIVATE(igrid)
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        call m1_heating_cooling(psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L, dtfactor, qdt);
     end do
     !$OMP END PARALLEL DO
   end if

   !$OMP PARALLEL DO PRIVATE(igrid)
   do iigrid=1,igridstail; igrid=igrids(iigrid);
      call m1_flux_absorbtion(psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L, dtfactor, qdt);
   end do
   !$OMP END PARALLEL DO

 end subroutine m1_implicit_update

 subroutine m1_evaluate_implicit(qtC,psa)
  use mod_global_parameters
  type(state), target :: psa(max_blocks)
  double precision, intent(in) :: qtC

  integer :: ixO^L,iigrid,igrid

  ixO^L=ixG^LL^LSUBnghostcells;

  if (m1_heatcool .and. hd_energy) then
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       call m1_heating_cooling_evaluate(psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L);
    end do
    !$OMP END PARALLEL DO
  end if

  !$OMP PARALLEL DO PRIVATE(igrid)
  do iigrid=1,igridstail; igrid=igrids(iigrid);
     call m1_flux_absorbtion_evaluate(psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L);
  end do
  !$OMP END PARALLEL DO

end subroutine m1_evaluate_implicit

 subroutine m1_flux_absorbtion(w, x, ixI^L, ixO^L, dtfactor, qdt)
   use mod_global_parameters
   use mod_geometry
   use mod_usr_methods, only: usr_set_F0

   integer, intent(in)             :: ixI^L,ixO^L
   double precision, intent(in)    :: x(ixI^S,1:ndim)
   double precision, intent(in)    :: dtfactor, qdt
   double precision, intent(inout)    :: w(ixI^S,1:nw)

   double precision :: F0(ixI^S,1:ndir), lambda(ixI^S), tau(ixI^S)
   double precision :: gradE(ixI^S), R(ixI^S)

   double precision :: aa(ixI^S), bb(ixI^S)

   integer :: idir, idim

   ! !> Assymptotic preserving IMEX implicit step:
   ! if (associated(usr_set_F0)) then
   !   call usr_set_F0(ixI^L,ixO^L,w,x,F0)
   ! else
   !   !> FLD approximation
   !   !> Get magnitude of grad E
   !   gradE = 0.d0
   !   do idir = 1, ndim
   !     call gradient_str(w(ixI^S,Er_),x,ixI^L,ixO^L,idir,lambda)
   !     gradE(ixI^S) = gradE(ixI^S) + lambda(ixI^S)**2
   !   end do
   !   gradE(ixI^S) = dsqrt(max(0.d0,gradE(ixI^S)))
   !
   !   !> Get fld R
   !   R(ixI^S) = gradE(ixI^S)/(w(ixI^S,rho_)*w(ixI^S,kF_)*w(ixI^S,Er_))
   !
   !   !> Get fld lambda (Pomraning flux limiter)
   !   lambda(ixI^S) = (2.d0+R(ixI^S))/(6.d0+3*R(ixI^S)+R(ixI^S)**2.d0)
   !
   !   !> initialise for ndir != ndim
   !   F0(ixI^S,:) = 0.d0
   !   do idir = 1, ndim
   !     call gradient_str(w(ixI^S,Er_),x,ixI^L,ixO^L,idir,gradE)
   !     F0(ixI^S,idir) = -m1_c*lambda(ixI^S)/(w(ixI^S,rho_)*w(ixI^S,kF_))*gradE(ixI^S)
   !   end do
   ! end if
   !
   ! !> AP lambda
   ! lambda(ixI^S) = m1_RSOL*m1_c*w(ixI^S,rho_)*w(ixI^S,kF_)
   ! tau(ixI^S) = lambda(ixI^S)*dtfactor*dt
   !
   ! do idir = 1, ndir
   !   w(ixO^S,Fr(idir)) = F0(ixO^S,idir) + (w(ixO^S,Fr(idir)) - F0(ixO^S,idir))*dexp(-tau(ixO^S))
   ! end do


   ! > Exact solution: dF/dt = - a F - b
   !                  <=> F = C1 * exp(-a*t) - b/a
   !                  w. C1 = F0 + b/a

  do idim = 1, ndim
     aa(ixI^S) = m1_RSOL*m1_c*w(ixI^S,rho_)*w(ixI^S,kF_)
     call m1_get_trapped_fraction(w, x, ixI^L, ixO^L, idim, lambda)
     call gradient_str(w(ixI^S,Er_)*lambda(ixI^S),x,ixI^L,ixO^L,idim,gradE)
     bb(ixI^S) = m1_RSOL*m1_c**2*gradE(ixI^S)/3.d0
     w(ixO^S,Fr(idim)) = (w(ixO^S,Fr(idim)) + bb(ixO^S)/aa(ixO^S)) * dexp(-aa(ixO^S)*dtfactor*dt) - bb(ixO^S)/aa(ixO^S)
  end do

 end subroutine m1_flux_absorbtion

 subroutine m1_flux_absorbtion_evaluate(w, x, ixI^L, ixO^L)
   use mod_global_parameters
   use mod_geometry
   use mod_usr_methods, only: usr_set_F0

   integer, intent(in)             :: ixI^L,ixO^L
   double precision, intent(in)    :: x(ixI^S,1:ndim)
   double precision, intent(inout) :: w(ixI^S,1:nw)

   double precision :: F0(ixI^S,1:ndir), lambda(ixI^S) , tau(ixI^S)
   double precision :: gradE(ixI^S), R(ixI^S)
   integer :: idir
   !
   ! !> Fully implicit update
   ! do idir = 1, ndir
   !   w(ixO^S,Fr(idir)) = -m1_RSOL*m1_c*w(ixO^S,rho_)*w(ixO^S,kF_)*w(ixO^S,Fr(idir))
   ! end do

   !> Assymptotic preserving IMEX implicit evaluation:
   if (associated(usr_set_F0)) then
     call usr_set_F0(ixI^L,ixO^L,w,x,F0)
   else
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

     !> Get fld lambda (Pomraning flux limiter)
     lambda(ixI^S) = (2.d0+R(ixI^S))/(6.d0+3*R(ixI^S)+R(ixI^S)**2.d0)

     !> initialise for ndir != ndim
     F0(ixI^S,:) = 0.d0

     do idir = 1, ndim
       call gradient_str(w(ixI^S,Er_),x,ixI^L,ixO^L,idir,gradE)
       F0(ixI^S,idir) = -m1_c*lambda(ixI^S)/(w(ixI^S,rho_)*w(ixI^S,kF_))*gradE(ixI^S)
     end do
   end if

   lambda(ixI^S) = m1_RSOL*m1_c*w(ixI^S,rho_)*w(ixI^S,kF_)

   do idir = 1, ndir
     w(ixO^S,Fr(idir)) = -lambda(ixO^S)*(w(ixO^S,Fr(idir)) - F0(ixO^S,idir))
   end do

 end subroutine m1_flux_absorbtion_evaluate


 subroutine m1_heating_cooling(w, x, ixI^L, ixO^L, dtfactor, qdt)
   use mod_global_parameters
   use mod_constants

   integer, intent(in)             :: ixI^L,ixO^L
   double precision, intent(in)    :: x(ixI^S,1:ndim)
   double precision, intent(in)    :: dtfactor, qdt
   double precision, intent(inout) :: w(ixI^S,1:nw)

   double precision :: a4(ixO^S), a0(ixO^S), alpha
   double precision :: c0(ixO^S), c1(ixO^S)

   integer :: ix^D, idir

   !> Get gas internal energy
   call hd_e_to_ei(ixI^L,ixO^L,w,x)

   !> Set polynomial coefficients
   alpha = const_rad_a*m1_c*(unit_temperature**4/unit_pressure) !> mu, mp, kb hidden in units

   a4(ixO^S) = alpha*w(ixO^S,kB_)*(hd_gamma-1.d0)**4/w(ixO^S,rho_)**3*dtfactor*dt
   a0(ixO^S) = m1_c*w(ixO^S,kE_)*w(ixO^S,rho_)*dtfactor*qdt

   c1(ixO^S) = 1.d0/(a4(ixO^S) - m1_RSOL*a4(ixO^S)*a0(ixO^S)/(1.d0 + m1_RSOL*a0(ixO^S)))
   c0(ixO^S) = (a0(ixO^S)*w(ixO^S,Er_)/(1.d0 + m1_RSOL*a0(ixO^S)) + w(ixO^S,e_))*c1(ixO^S)

   !> Solve polynomial
   {do ix^D = ixOmin^D,ixOmax^D\}
    w(ix^D,e_) = max(w(ix^D,e_), small_e)
    call Bisection_method(w(ix^D,e_),w(ix^D,Er_),c0(ix^D),c1(ix^D))
   {enddo\}

   !> Set radiation energy density
   w(ixO^S,Er_) = (w(ixO^S,Er_) + m1_RSOL*a4(ixO^S)*w(ixO^S,e_)**4.d0)/(1.d0 + m1_RSOL*a0(ixO^S))

   !> Get gas total energy
   call hd_ei_to_e(ixI^L,ixO^L,w,x)

   if (fix_small_values) then
     call hd_handle_small_values(.false., w, x, ixI^L, ixO^L, 'Energy_interaction')
   end if

 end subroutine m1_heating_cooling

 subroutine m1_heating_cooling_evaluate(w, x, ixI^L, ixO^L)
   use mod_global_parameters

   integer, intent(in)             :: ixI^L,ixO^L
   double precision, intent(in)    :: x(ixI^S,1:ndim)
   double precision, intent(inout) :: w(ixI^S,1:nw)

   double precision :: Bb(ixI^S)

   !> Get pressure
   call hd_get_pthermal(w,x,ixI^L,ixO^L,Bb)

   !> Gas temperature
   Bb(ixI^S) = const_mp*hd_mu/const_kB*Bb(ixI^S)/w(ixO^S,rho_)
   Bb(ixI^S) = Bb(ixI^S)*unit_pressure/unit_density

   !> Integrated black body mean intensity
   Bb(ixI^S) = const_sigma*Bb(ixI^S)**4/dpi
   Bb(ixI^S) = Bb(ixI^S)*unit_time/unit_pressure

   !> Heating and cooling sourceterms
   w(ixO^S,e_) = m1_RSOL*m1_c*w(ixO^S,kE_)*w(ixO^S,rho_)*w(ixO^S,Er_) &
                - 4*dpi*m1_RSOL*w(ixO^S,kB_)*w(ixO^S,rho_)*Bb(ixO^S)

   w(ixO^S,Er_) = - m1_RSOL*m1_c*w(ixO^S,kE_)*w(ixO^S,rho_)*w(ixO^S,Er_) &
                  + 4*dpi*m1_RSOL*w(ixO^S,kB_)*w(ixO^S,rho_)**Bb(ixO^S)

 end subroutine m1_heating_cooling_evaluate

  subroutine m1_set_opacities(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)

    double precision :: tau_sob(ixI^S), vr(ixI^S), dvrdr(ixI^S)
    double precision :: fdfac(ixO^S), beta_fd(ixO^S)

    select case (m1_opacity_type)

    case('zero')
      w(ixI^S,kF_) = 0.d0
      w(ixI^S,kE_) = 0.d0
      w(ixI^S,kB_) = 0.d0

    case('kappa_e')
      w(ixI^S,kF_) = m1_kappa_e!/unit_opacity
      w(ixI^S,kE_) = m1_kappa_e!/unit_opacity
      w(ixI^S,kB_) = m1_kappa_e!/unit_opacity

    case('Gayley')

      !> Calculate velocity gradient
      vr(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
      call gradient_str(vr,x,ixI^L,ixO^L,1,dvrdr)

      !> Calculate sobolev optical depth parameter
      tau_sob(ixI^S) = dabs(m1_c*m1_kappa_e*w(ixI^S,rho_)/(dvrdr(ixI^S)+smalldouble))

      w(ixI^S,kF_) = m1_kappa_e * Qb_F/(1-a_F) * ((1.d0+Q0_F*tau_sob(ixI^S))**(1.d0-a_F) - 1.d0)/(Q0_F*tau_sob(ixI^S))
      w(ixI^S,kF_) = w(ixI^S,kF_) + m1_kappa_e

      ! w(ixI^S,kE_) = m1_kappa_e * Qb_E/(1-a_E)*((1+Q0_E*tau_sob(ixI^S))**(1-a_E) - 1)/(Q0_E*tau_sob(ixI^S))
      ! w(ixI^S,kB_) = m1_kappa_e * Qb_B/(1-a_B)*((1+Q0_B*tau_sob(ixI^S))**(1-a_B) - 1)/(Q0_B*tau_sob(ixI^S))
      w(ixI^S,kE_) = m1_kappa_e!*Qb_E
      w(ixI^S,kB_) = m1_kappa_e!*Qb_B

      !> Finite disc correction factor
      beta_fd(ixO^S) = ( 1.0d0 - vr(ixO^S)/(x(ixO^S,1) * dvrdr(ixO^S)) ) &
                        * (xprobmin1/x(ixO^S,1))**2

      where (beta_fd(ixO^S) >= 1.0d0)
        fdfac(ixO^S) = 1.0d0/(1.0d0 + a_F)
      elsewhere (beta_fd(ixO^S) < -1.0d10)
        fdfac(ixO^S) = abs(beta_fd(ixO^S))**a_F / (1.0d0 + a_F)
      elsewhere (abs(beta_fd(ixO^S)) > 1.0d-3)
        fdfac(ixO^S) = (1.0d0 - (1.0d0 - beta_fd(ixO^S))**(1.0d0 + a_F)) &
                       / (beta_fd(ixO^S)*(1.0d0 + a_F))
      elsewhere
        fdfac(ixO^S) = 1.0d0 - 0.5d0*a_F*beta_fd(ixO^S) &
                       * (1.0d0 + 1.0d0/3.0d0 * (1.0d0 - a_F)*beta_fd(ixO^S))
      end where

      w(ixO^S,kF_) = w(ixO^S,kF_)*fdfac(ixO^S)

      where (w(ixI^S,kF_) .gt. 15*m1_kappa_e) w(ixI^S,kF_) = 15*m1_kappa_e
      ! where (w(ixI^S,kE_) .gt. 20*m1_kappa_e) w(ixI^S,kE_) = 20*m1_kappa_e
      ! where (w(ixI^S,kB_) .gt. 20*m1_kappa_e) w(ixI^S,kB_) = 20*m1_kappa_e

      !> ghostcells on opacity for force-balancing
      {^IFONED
      w(ixOmin1-1,kF_) = w(ixOmin1,kF_)
      w(ixOmax1+1,kF_) = w(ixOmax1,kF_)
      }
      {^IFTWOD
      w(ixOmin1-1,:,kF_) = w(ixOmin1,:,kF_)
      w(ixOmax1+1,:,kF_) = w(ixOmax1,:,kF_)
      }
      {^IFTHREED
      w(ixOmin1-1,:,:,kF_) = w(ixOmin1,:,:,kF_)
      w(ixOmax1+1,:,:,kF_) = w(ixOmax1,:,:,kF_)
      }

    case('hybrid1')

      call hybrid1_opacity(ixI^L,ixO^L,w,x)

    case('hybrid2')

      call hybrid2_opacity(ixI^L,ixO^L,w,x)

    case('special')
      if (.not. associated(usr_special_opacity_m1)) then
        call mpistop("special opacity not defined")
      endif
      call usr_special_opacity_m1(ixI^L, ixO^L, w, x)

    case default
      call mpistop("Opacity type not known...")

    end select

  end subroutine m1_set_opacities

  subroutine hybrid1_opacity(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_opal_opacity
    use mod_geometry

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)

    double precision :: tau_sob(ixI^S), vr(ixI^S), dvrdr(ixI^S)

    double precision :: OPAL(ixO^S), rho, Trad, kappa
    integer :: ix^D

    !> Velocity stretched opacities
    !> Calculate velocity gradient
    vr(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradient_str(vr,x,ixI^L,ixO^L,1,dvrdr)

    !> Calculate sobolev optical depth parameter
    tau_sob(ixO^S) = dabs(m1_c*m1_kappa_e*w(ixO^S,rho_)/(dvrdr(ixO^S)+smalldouble))

    w(ixO^S,kF_) = m1_kappa_e * Qb_F/(1-a_F) * ((1.d0+Q0_F*tau_sob(ixO^S))**(1.d0-a_F) - 1.d0)/(Q0_F*tau_sob(ixO^S))
    w(ixO^S,kE_) = m1_kappa_e * Qb_E/(1-a_E)*((1+Q0_E*tau_sob(ixO^S))**(1-a_E) - 1)/(Q0_E*tau_sob(ixO^S))
    w(ixO^S,kB_) = m1_kappa_e * Qb_B/(1-a_B)*((1+Q0_B*tau_sob(ixO^S))**(1-a_B) - 1)/(Q0_B*tau_sob(ixO^S))

    !> Opal Opacities
    {do ix^D=ixOmin^D,ixOmax^D\ }
        rho = w(ix^D,rho_)*unit_density
        Trad = (w(ix^D,Er_)*unit_pressure/const_rad_a)**(0.25d0)
        Trad = max(Trad,1.d4)
        call set_opal_opacity(rho,Trad,kappa)
        OPAL(ix^D) = kappa/unit_opacity
    {enddo\ }

    !> Combine OPAL and CAK-like
    w(ixO^S,kF_) = min(w(ixO^S,kF_), 50.d0*m1_kappa_e) + OPAL(ixO^S)
    w(ixO^S,kE_) = w(ixO^S,kE_) + OPAL(ixO^S) - m1_kappa_e
    w(ixO^S,kB_) = w(ixO^S,kB_) + OPAL(ixO^S) - m1_kappa_e

    w(ixO^S,kF_) = max(w(ixO^S,kF_), 0.1d0*m1_kappa_e)
    w(ixO^S,kE_) = max(w(ixO^S,kE_), 1.d-5*m1_kappa_e)
    w(ixO^S,kB_) = max(w(ixO^S,kB_), 1.d-5*m1_kappa_e)

    if (any(w(ixO^S,kF_) .ne. w(ixO^S,kF_))) then
      stop "kappa_F NaN in hybrid opacity"
    end if

    if (any(w(ixO^S,kE_) .ne. w(ixO^S,kE_))) then
      stop "kappa_E NaN in hybrid opacity"
    end if

    if (any(w(ixO^S,kB_) .ne. w(ixO^S,kB_))) then
      stop "kappa_B NaN in hybrid opacity"
    end if

  end subroutine hybrid1_opacity

  subroutine hybrid2_opacity(ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_opal_opacity
    use mod_cak_opacity
    use mod_geometry
    use mod_windsim_tools, only: R_phot, i_phot

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)

    double precision :: tau_sob(ixI^S), vr(ixI^S), dvrdr(ixI^S)

    double precision :: OPAL(ixO^S), rho, Trad, kappa
    double precision :: CAK(ixO^S), alpha, Qbar, Q0, kappa_e_t, gradv, tau

    double precision :: fdfac, beta_fd

    integer :: ix^D

    !> Velocity stretched opacities
    !> Calculate velocity gradient
    vr(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradient_str(vr,x,ixI^L,ixO^L,1,dvrdr)

    !if ((modulo(it,1000) == 0) .and. (mype .eq. 0)) then
    !  print*, 'WST photosphere: ', mype, i_phot, R_phot
    !end if

    !> Opal Opacities
    {do ix^D=ixOmin^D,ixOmax^D\ }
        rho = w(ix^D,rho_)*unit_density
        Trad = (w(ix^D,Er_)*unit_pressure/const_rad_a)**(0.25d0)
        Trad = max(Trad,1.d4)
        gradv = abs(dvrdr(ix^D))*unit_velocity/unit_length

        call set_cak_opacity(rho, Trad, gradv ,alpha, Qbar, Q0, kappa_e_t)

        !> Finite disc correction factor
        beta_fd = ( 1.0d0 - vr(ix^D)/(x(ix^D,1) * dvrdr(ix^D)) ) &
        * (R_phot/x(ix^D,1))**2

        if (x(ix^D,1) .le. R_phot) then
          fdfac = 1.d0/(1.d0 + alpha)
        else if (beta_fd >= 1.0d0) then
          fdfac = 1.0d0/(1.0d0 + alpha)
        else if (beta_fd < -1.0d10) then
          fdfac = abs(beta_fd)**alpha / (1.0d0 + alpha)
        else if (abs(beta_fd) > 1.0d-3) then
          fdfac = (1.0d0 - (1.0d0 - beta_fd)**(1.0d0 + alpha)) &
          / (beta_fd*(1.0d0 + alpha))
        else
          fdfac = 1.0d0 - 0.5d0*alpha*beta_fd &
          * (1.0d0 + 1.0d0/3.0d0 * (1.0d0 - alpha)*beta_fd)
        end if

        tau = (m1_kappa_e*unit_opacity)*rho*m1_c*unit_velocity/gradv
        CAK(ix^D) = m1_kappa_e * fdfac * Qbar/(1-alpha)*((1+Q0*tau)**(1-alpha) - 1)/(Q0*tau)

        call set_opal_opacity(rho,Trad,kappa)
        OPAL(ix^D) = kappa/unit_opacity
    {enddo\ }

    !> Combine OPAL and CAK-like
    w(ixO^S,kF_) = min(CAK(ixO^S), 500.d0*m1_kappa_e) + OPAL(ixO^S)
    w(ixO^S,kE_) = 1.d-1 * CAK(ixO^S) + OPAL(ixO^S) - m1_kappa_e
    w(ixO^S,kB_) = 1.d-1 * CAK(ixO^S) + OPAL(ixO^S) - m1_kappa_e

 !   w(ixO^S,kE_) = OPAL(ixO^S) - m1_kappa_e
  !  w(ixO^S,kB_) = OPAL(ixO^S) - m1_kappa_e

    w(ixO^S,kF_) = max(w(ixO^S,kF_), 0.1d0*m1_kappa_e)
    w(ixO^S,kE_) = max(w(ixO^S,kE_), 1.d-1*m1_kappa_e)
    w(ixO^S,kB_) = max(w(ixO^S,kB_), 1.d-1*m1_kappa_e)

    !> Correction for line driving opacity at outer physical cells
    if (any(x(ixI^S,1) .gt. xprobmax1)) then
      {^IFONED w(ixOmax1,kF_) = max(w(ixOmax1-1,kF_), w(ixOmax1,kF_))}
      {^IFTWOD w(ixOmax1,:,kF_) = max(w(ixOmax1-1,:,kF_), w(ixOmax1,:,kF_))}
      {^IFTHREED w(ixOmax1,:,:,kF_) = max(w(ixOmax1-1,:,:,kF_), w(ixOmax1,:,:,kF_))}
    end if


    if (any(w(ixO^S,kF_) .ne. w(ixO^S,kF_))) then
      stop "kappa_F NaN in hybrid opacity"
    end if

    if (any(w(ixO^S,kE_) .ne. w(ixO^S,kE_))) then
      stop "kappa_E NaN in hybrid opacity"
    end if

    if (any(w(ixO^S,kB_) .ne. w(ixO^S,kB_))) then
      stop "kappa_B NaN in hybrid opacity"
    end if

  end subroutine hybrid2_opacity

  !> Check time step for total radiation contribution
  subroutine m1_get_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters

    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = courantpar*minval(block%dx(ixO^S,:)/(m1_c*m1_RSOL))

  end subroutine m1_get_dt

  !> Returns logical argument flag where values are ok
  subroutine m1_check_w(ixI^L, ixO^L, w, flag)
    use mod_global_parameters
    use mod_dust, only: dust_check_w

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, nw)
    logical, intent(inout)       :: flag(ixI^S,1:nw)
    double precision             :: Fmag2(ixI^S)

    integer :: idir

    !> Compute radflux magnitude
    Fmag2(ixI^S) = sum(w(ixI^S,Fr(:))**2, dim=ndim+1)

    !> Flag negative values of Erad
    where (w(ixO^S, Er_) < small_e) flag(ixO^S, Er_) = .true.

    !> Flag flux values where |F| > cE
    !> Flag all dirs
    do idir =1,ndir
      where (Fmag2(ixO^S) > (m1_c*w(ixO^S, Er_))**2)
        flag(ixO^S,Fr(idir)) = .true.
      end where
    end do

  end subroutine m1_check_w

  !> Find the root of the 4th degree polynomial using the bisection method
   subroutine Bisection_method(e_gas, E_rad, c0, c1)
     use mod_global_parameters
     double precision, intent(in)    :: c0, c1
     double precision, intent(in)    :: E_rad
     double precision, intent(inout) :: e_gas
     double precision :: e_a, e_b, e_c !> x-axis values (energy) to evaluat in bisection
     double precision :: p_a, p_b, p_c !> polynomial values
     integer :: n, max_its

     !> initialise counters and bisection limits
     n = 0
     max_its = 50
     e_a = zero
     e_b = min(abs(c0/c1),abs(c0)**(1.d0/4.d0))+smalldouble

     !> Start bisection algorithm
     do while((abs(e_b-e_a) .ge. 1.d-15*min(e_gas,E_rad)) .and. n < max_its)
       !> Bisection point
       e_c = (e_a + e_b)/two

       !> Evaluate polynomial
       p_a = Polynomial_Bisection(e_a, c0, c1)
       p_b = Polynomial_Bisection(e_b, c0, c1)
       p_c = Polynomial_Bisection(e_c, c0, c1)

       !> Sanity check NEEDS TO BE REMOVED AFTER TESTING
       if (p_a*p_b .gt. 0) then
         !print*, 'fail on it :', it
         !print*, 'p_a, p_b, n_hc',p_a, p_b, n
         !print*, 'e, E', e_gas, E_rad
         !print*, 'x**4 +',c1,'x -',c0
         !print*, 'between x=', e_a, 'and', e_b
         !call mpistop("Problems with bisection scheme")

        n = max_its
        e_a = small_e
        e_b = small_e
        e_c = small_e

       end if

       !> Bisection
       if (p_a*p_c .le. 0) then
         e_b = e_c
       else if (p_b*p_c .le. 0) then
         e_a = e_c
       else

         !print*, 'fail on it :', it
         !print*, 'e_g, E_r:',e_gas, E_rad
         !print*, 'c0, c1', c0, c1
         !print*, n, 'a, c, b'
         !print*, e_a, e_c, e_b
         !print*, p_a, p_c, p_b

         call mpistop('mehhh')

         e_a = e_c
         e_b = e_c

         print*, 'heatcool:', e_c

       end if

       !> counter
       n = n + 1

     end do

     !> Final bisection
     e_gas = (e_a + e_b)/two

   end subroutine Bisection_method

   function Polynomial_Bisection(e_gas, c0, c1) result(p)
      use mod_global_parameters
      double precision, intent(in) :: e_gas
      double precision, intent(in) :: c0, c1
      double precision :: p

      p = e_gas**4.d0 + c1*e_gas - c0
    end function Polynomial_Bisection

end module mod_hd_phys
