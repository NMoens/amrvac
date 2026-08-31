module mod_m1
  use mod_comm_lib, only: mpistop


  implicit none

  !> Compute m1 radiation force
  logical, public, protected              :: m1_force = .true.

  !> Compute m1 radiation heating and cooling
  logical, public                         :: m1_heatcool = .true.

  !> Force split m1 force
  logical, public, protected              :: m1_force_split = .false.

  !> Compute m1 pressure in diffusion limit
  logical, public, protected              :: m1_diffusion = .false.

  !> Split trapped vs streaming photons in Flux eq
  logical, public, protected              :: m1_split_trapped = .true.

  !> Have different values for kF, kE, kB opacities
  logical, public, protected              :: m1_split_opacities = .true.

  !> small_e set in hd module
  double precision, public                :: m1_small_e
  double precision, public                :: m1_gamma
  logical, public                         :: m1_req_diagonal = .false.

  !> Reduced speed of light factor
  double precision, public                :: m1_RSOL = 1.d0

  !> Electron scattering opacity
  double precision, public                :: m1_kappa_e = 0.34d0

  !> The mean molecular weight
  double precision, public                :: m1_mu = 0.6d0

  !> Speed of light used in m1 subsystem
  double precision, public, protected     :: m1_c
  double precision, public, protected     :: m1_sb
  double precision, public, protected     :: alpha_qdot, beta_qdot

  !> Opacity prescription
  character(len=8)                        :: m1_opacity_type = 'kappa_e'

  !> Indices of opacities
  integer, public, protected              :: kF_, kE_, kB_

  !> local defenition of line driving parameters:
  double precision, public, protected :: a_F = 0.7d0
  double precision, public, protected :: a_E = 0.66d0
  double precision, public, protected :: a_B = 0.66d0
  double precision, public, protected :: Qb_F = 2.d3
  double precision, public, protected :: Qb_E = 2.d2
  double precision, public, protected :: Qb_B = 2.d2
  double precision, public, protected :: Q0_F = 2.d3
  double precision, public, protected :: Q0_E = 2.d2
  double precision, public, protected :: Q0_B = 2.d2

contains

  subroutine m1_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /m1_list/ m1_force, m1_opacity_type, m1_diffusion, m1_RSOL, &
    m1_force_split, m1_split_trapped, m1_split_opacities, &
    m1_heatcool, m1_kappa_e

    do n = 1, size(files)
      open(unitpar, file=trim(files(n)), status="old")
      read(unitpar, m1_list, end=111)
      111   close(unitpar)
    end do

  end subroutine m1_read_params

  subroutine m1_init(phys_mu, small_e, He_abundance)
    use mod_global_parameters
    use mod_opal_opacity, only: init_opal_table
    use mod_cak_opacity, only: init_cak
    use mod_windsim_tools, only: init_windsim_tools, wst_get_photosphere
    use mod_physics
    use mod_constants

    double precision, intent(in) :: phys_mu, small_e, He_abundance

    double precision :: sigma_thom

    call m1_read_params(par_files)

    !> Set index of opacity variables
    kF_ = var_set_extravar("kappa_F", "kappa_F")
    if (m1_split_opacities) then
      kE_ = var_set_extravar("kappa_E", "kappa_E")
      kB_ = var_set_extravar("kappa_B", "kappa_B")
    else
      kE_ = kF_
      kB_ = kF_
    end if

    !> Heating and cooling, absorbtion sourceterms
    if (m1_opacity_type .ne. 'zero') then
      phys_implicit_update => m1_implicit_update
      phys_evaluate_implicit => m1_evaluate_implicit
    end if

    !> Setting opacities
    phys_get_auxiliary => m1_set_opacities

    !> Radiation units
    unit_opacity = unit_length**2/unit_mass
    unit_radflux = unit_pressure*unit_velocity

    !> Local copy of mean molecular weight
    m1_mu = phys_mu
    m1_gamma = phys_gamma
    m1_small_e = small_e

    if (SI_unit) then
      !> Thomson crossection
      sigma_thom = sigma_thom_SI

      !> dimensionless speed of light and stefan boltzmann
      m1_c = c_SI/unit_velocity
      m1_sb = sigma_SI/unit_radflux*unit_temperature**4

      !> Coefficient in heating-cooling sourceterm
      alpha_qdot = rad_a_SI*m1_c*(unit_temperature**4/unit_pressure) !> mu, mp, kb hidden in units
      beta_qdot = mp_SI*m1_mu/kB_SI*unit_pressure/unit_density

      !> Electron scattering opacity
      m1_kappa_e = sigma_thom_SI/mp_SI * (1 + 2*He_abundance)/(1 + 4*He_abundance)
      m1_kappa_e = m1_kappa_e/unit_opacity
    else
      !> Thomson crossection
      sigma_thom = sigma_thom_cgs

      !> dimensionless speed of light and stefan boltzmann
      m1_c = const_c/unit_velocity
      m1_sb = const_sigma/unit_radflux*unit_temperature**4

      !> Coefficient in heating-cooling sourceterm
      alpha_qdot = const_rad_a*m1_c*(unit_temperature**4/unit_pressure) !> mu, mp, kb hidden in units
      beta_qdot = const_mp*m1_mu/const_kB*unit_pressure/unit_density

      !> Electron scattering opacity
      m1_kappa_e = sigma_thom_cgs/mp_cgs * (1 + 2*He_abundance)/(1 + 4*He_abundance)
      m1_kappa_e = m1_kappa_e/unit_opacity
    end if

    !> Input from phys_init
    !> heatcool depends on energy eq
    if (.not. phys_energy) m1_heatcool = .false.

    !> diagonal requirements in gc
    m1_req_diagonal = phys_req_diagonal

    !> Initialise opacity modules if necessary
    if (m1_opacity_type .eq. 'hybrid1') call init_opal_table('Y02800')
    if (m1_opacity_type .eq. 'hybrid2') then
      call init_opal_table('Y02800')
      call init_cak('Y02800')

      !> Photosphere needed
      call init_windsim_tools(iw_rho,m1_kappa_e)
      phys_process_global => wst_get_photosphere
    end if

    !> Source-splitting for m1 forces
    if (m1_force_split) any_source_split=.true.

    !> Allocate vector index array
    nvector      = 2 ! No. vector vars
    deallocate(iw_vector)
    allocate(iw_vector(nvector))
    iw_vector(1) = iw_mom(1) - 1
    iw_vector(2) = iw_Fr(1) - 1

    ! set number of variables which need update ghostcells
    ! For M1, this includes the opacities (Needed for cbounds)
    if (m1_split_opacities) then
      nwgc=nwflux+3
    else
      nwgc=nwflux+1
    end if

  end subroutine m1_init


  subroutine m1_get_cbounds(wLC, wRC, wLp, wRp, x, ixI^L, ixO^L, idim, Hspeed, cmax, cmin)
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

    double precision                          :: m1_cmax(ixI^S), m1_cmin(ixI^S)

    double precision                          :: fi(ixI^S), fj(ixI^S)
    double precision                          :: ff(ixI^S), xi(ixI^S)

    double precision                          :: lam_L_max(ixI^S), lam_L_min(ixI^S)
    double precision                          :: lam_R_max(ixI^S), lam_R_min(ixI^S)

    double precision                          :: tau_c(ixI^S), kappa, mu(ixI^S)

    kappa = m1_kappa_e
    if (m1_opacity_type .eq. 'zero') kappa = 0.d0

    !TODO: Check this sum, it does not produce what is expected
    !> Left state
    ff(ixI^S) = sum(wLC(ixI^S, iw_Fr(:))**2, dim=ndim+1)/(m1_c*wLC(ixI^S,iw_r_e))**2
    ff(ixI^S) = min(ff(ixI^S),1.d0)

    !> Eddington scalars
    fi(ixI^S) = abs(wLC(ixI^S,iw_Fr(idim))/(m1_c*wLC(ixI^S,iw_r_e)))
    fi(ixI^S) = min(fi(ixI^S),1.d0)
    mu(ixI^S) = abs(wLC(ixI^S,iw_Fr(idim)))/dsqrt(sum(wLC(ixI^S, iw_Fr(:))**2, dim=ndim+1))
    mu(ixI^S) = min(mu(ixI^S),1.d0)
    mu(ixI^S) = max(mu(ixI^S),-1.d0)

    fj(ixI^S) = dsqrt(max(ff(ixI^S) - fi(ixI^S)**2,0.d0))

    xi(ixI^S) = dsqrt(4-3*ff(ixI^S))

    !TODO: Using same cell optical depth for both left and right state. Left state should use tau at i-1 i think?
    tau_c(ixI^S) = block%dx(ixI^S,idim)*wLC(ixI^S,iw_rho)*wLp(ixI^S,kF_)

    ! lam_L_max(ixI^S) = m1_c*m1_RSOL*(fi(ixI^S)/xi(ixI^S) &
    ! + dsqrt(2*(xi(ixI^S)-1)*(xi(ixI^S)+2)*(2*(xi(ixI^S)-1)*(xi(ixI^S)+2) + 3*fj(ixI^S)**2)) &
    ! /(dsqrt(3.d0)*xi(ixI^S)*(xi(ixI^S)+2)))
    lam_L_max(ixI^S) = m1_c*m1_RSOL*(fi(ixI^S)*dsqrt(ff(ixI^S)) &
    + dsqrt(2.d0/3*(4 - 3*ff(ixI^S) - xi(ixI^S)) &
    + 2*mu(ixI^S)**2*(2 - ff(ixI^S) - xi(ixI^S))))&
    /xi(ixI^S)


    ! if (block%igrid .eq. 1) then
    ! print*, 'max L'
    ! print*, lam_L_max(ixO^S)/m1_c
    ! print*, 'diff L'
    ! print*, 4.d0*m1_RSOL/(3*tau_c(ixO^S))
    ! end if

    ! Fuksman Melon 2019, Sadowski 20212
    ! lam_L_max(ixI^S) = min(lam_L_max(ixI^S), 2.d0/(3*tau_c(ixI^S)))
    lam_L_max(ixI^S) = min(lam_L_max(ixI^S), 4.d0*m1_c*m1_RSOL/(3*tau_c(ixI^S)))

    if (present(cmin)) then
      ! lam_L_min(ixI^S) = m1_c*m1_RSOL*(fi(ixI^S)/xi(ixI^S) &
      ! - dsqrt(2*(xi(ixI^S)-1)*(xi(ixI^S)+2)*(2*(xi(ixI^S)-1)*(xi(ixI^S)+2) + 3*fj(ixI^S)**2)) &
      ! /(dsqrt(3.d0)*xi(ixI^S)*(xi(ixI^S)+2)))
      lam_L_min(ixI^S) = m1_c*m1_RSOL*(mu(ixI^S)*dsqrt(ff(ixI^S)) &
      - dsqrt(2.d0/3*(4 - 3*ff(ixI^S) - xi(ixI^S)) &
      + 2*mu(ixI^S)**2*(2 - ff(ixI^S) - xi(ixI^S))))&
      /xi(ixI^S)

      ! Fuksman Melon 2019, Sadowski 20212
      ! lam_L_min(ixI^S) = max(lam_L_min(ixI^S), -2.d0/(3*tau_c(ixI^S)))
      lam_L_min(ixI^S) = max(lam_L_min(ixI^S), -4.d0*m1_c*m1_RSOL/(3*tau_c(ixI^S)))
    end if

    !> Right state
    ff(ixI^S) = sum(wRC(ixI^S, iw_Fr(:))**2, dim=ndim+1)/(m1_c*wRC(ixI^S,iw_r_e))**2
    ff(ixI^S) = min(ff(ixI^S),1.d0)

    !> Eddington scalars
    fi(ixI^S) = abs(wRC(ixI^S,iw_Fr(idim))/(m1_c*wRC(ixI^S,iw_r_e)))
    fi(ixI^S) = min(fi(ixI^S),1.d0)
    mu(ixI^S) = abs(wRC(ixI^S,iw_Fr(idim)))/dsqrt(sum(wRC(ixI^S, iw_Fr(:))**2, dim=ndim+1))
    mu(ixI^S) = min(mu(ixI^S),1.d0)
    mu(ixI^S) = max(mu(ixI^S),-1.d0)

    fj(ixI^S) = dsqrt(max(ff(ixI^S) - fi(ixI^S)**2,0.d0))

    xi(ixI^S) = dsqrt(4-3*ff(ixI^S))

    tau_c(ixI^S) = block%dx(ixI^S,idim)*wRC(ixI^S,iw_rho)*wRp(ixI^S,kF_)

    ! lam_R_max(ixI^S) = m1_c*m1_RSOL*(fi(ixI^S)/xi(ixI^S) &
    ! + dsqrt(2*(xi(ixI^S)-1)*(xi(ixI^S)+2)*(2*(xi(ixI^S)-1)*(xi(ixI^S)+2) + 3*fj(ixI^S)**2)) &
    ! /(dsqrt(3.d0)*xi(ixI^S)*(xi(ixI^S)+2)))
    lam_R_max(ixI^S) = m1_c*m1_RSOL*(mu(ixI^S)*dsqrt(ff(ixI^S)) &
    + dsqrt(2.d0/3*(4 - 3*ff(ixI^S) - xi(ixI^S)) &
    + 2*mu(ixI^S)**2*(2 - ff(ixI^S) - xi(ixI^S))))&
    /xi(ixI^S)

    ! if (block%igrid .eq. 1) then
    ! ! print*, 'max L'
    ! ! print*, lam_L_max(ixO^S)/m1_c
    ! print*, 'max R'
    ! print*, lam_R_max(ixO^S)/m1_c
    ! print*, 'diff R'
    ! print*,  4.d0*m1_RSOL/(3*tau_c(ixO^S))
    ! end if

    ! Fuksman Melon 2019, Sadowski 20212
    ! lam_R_max(ixI^S) = min(lam_R_max(ixI^S), 2.d0/(3*tau_c(ixI^S)))
    lam_R_max(ixI^S) = min(lam_R_max(ixI^S), 4.d0*m1_c*m1_RSOL/(3*tau_c(ixI^S)))

    if (present(cmin)) then
      ! lam_R_min(ixI^S) = m1_c*m1_RSOL*(fi(ixI^S)/xi(ixI^S) &
      ! - dsqrt(2*(xi(ixI^S)-1)*(xi(ixI^S)+2)*(2*(xi(ixI^S)-1)*(xi(ixI^S)+2) + 3*fj(ixI^S)**2)) &
      ! /(dsqrt(3.d0)*xi(ixI^S)*(xi(ixI^S)+2)))
      lam_R_min(ixI^S) = m1_c*m1_RSOL*(mu(ixI^S)*dsqrt(ff(ixI^S)) &
      - dsqrt(2.d0/3*(4 - 3*ff(ixI^S) - xi(ixI^S)) &
      + 2*mu(ixI^S)**2*(2 - ff(ixI^S) - xi(ixI^S))))&
      /xi(ixI^S)

      ! Fuksman Melon 2019, Sadowski 20212
      ! lam_R_min(ixI^S) = max(lam_R_min(ixI^S), -2.d0/(3*tau_c(ixI^S)))
      lam_R_min(ixI^S) = max(lam_R_min(ixI^S), -4.d0*m1_c*m1_RSOL/(3*tau_c(ixI^S)))
    end if

    ! Skinner Ostricker 2013, Eq 39,
    m1_cmax(ixI^S) = max(lam_L_max(ixI^S),lam_R_max(ixI^S))
    m1_cmax(ixI^S) = max(m1_cmax(ixI^S), 0.d0)
    ! if (block%igrid .eq. 1) then
    ! ! print*, 'max L'
    ! ! print*, lam_L_max(ixO^S)/m1_c
    ! ! print*, 'max R'
    ! ! print*, lam_R_max(ixO^S)/m1_c
    ! print*, 'cmax hd'
    ! print*, cmax(ixO^S,1)/m1_c
    ! print*, 'cmax m1'
    ! print*, m1_cmax(ixO^S)/m1_c
    ! print*, 'tau cell'
    ! print*, tau_c(ixO^S)
    ! print*, 'ratio'
    ! print*, m1_cmax(ixO^S)/cmax(ixO^S,1)
    ! end if
    cmax(ixI^S,1) = max(cmax(ixI^S,1), m1_cmax(ixI^S))

    if (present(cmin)) then
      m1_cmin(ixI^S) = min(lam_L_min(ixI^S),lam_R_min(ixI^S))
      m1_cmin(ixI^S) = min(m1_cmin(ixI^S), 0.d0)
      cmin(ixI^S,1) = min(cmin(ixI^S,1), m1_cmin(ixI^S))
    end if

  end subroutine m1_get_cbounds


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


    if (present(str)) then
      streaming = str
    else
      streaming = .false.
    end if

    if (streaming) then
      call m1_get_streaming_fraction(w, x, ixI^L, ixO^L, idir, Erad)
      Erad(ixI^S) = w(ixI^S,iw_r_e)*Erad(ixI^S)
    else
      Erad(ixI^S) = w(ixI^S,iw_r_e)
    end if

    if (m1_diffusion) then

      if (idir == jdir) then
        prad(ixI^S) = 1.d0/3.d0
      else
        prad(ixI^S) = 0.d0
      end if

    else

      !> Reduced fluxes
      fi(ixI^S) = dabs(w(ixI^S,iw_Fr(idir)))/(m1_c*Erad(ixI^S))
      fj(ixI^S) = dabs(w(ixI^S,iw_Fr(jdir)))/(m1_c*Erad(ixI^S))

      ff(ixI^S) = sum(w(ixI^S, iw_Fr(:))**2, dim=ndim+1)/(m1_c*Erad(ixI^S))**2

      !> Cap ff, fi, fj
      fi(ixI^S) = min(fi(ixI^S),1.d0)
      fj(ixI^S) = min(fj(ixI^S),1.d0)

      fi(ixI^S) = max(fi(ixI^S),smalldouble)
      fj(ixI^S) = max(fj(ixI^S),smalldouble)

      fi(ixI^S) = sign(fi(ixI^S),w(ixI^S,iw_Fr(idir)))
      fj(ixI^S) = sign(fj(ixI^S),w(ixI^S,iw_Fr(jdir)))

      ff(ixI^S) = min(ff(ixI^S),1.d0)
      ff(ixI^S) = max(ff(ixI^S),smalldouble)

      !> Eddington factor
      chi(ixI^S) = (3.d0+4*ff(ixI^S))/(5.d0+2*dsqrt(4.d0-3*ff(ixI^S)))
      chi(ixI^S) = max(chi(ixI^S),1.d0/3)
      chi(ixI^S) = min(chi(ixI^S),1.d0)

      !> Temporary (??) fix
      where (chi(ixI^S) .ne. chi(ixI^S))
        chi(ixI^S) = 1.d0 !3.d0/5.d0
      end where

      !> Refactor ff, fi, fj
      ! ff(ixI^S) = fi(ixI^S)*fj(ixI^S)/ff(ixI^S)
      ! !TODO IMPLEMENT CORRECT LIMIT!!! DEPENDS ON DIRECTION OF F OFC!!!!
      ! where(fi(ixI^S)*fj(ixI^S) .lt. smalldouble)
      !   ff(ixI^S) = 0.d0
      ! end where

      ! where (ff(ixI^S) .ne. ff(ixI^S))
      !   ff(ixI^S) = sign(1.d0,w(ixI^S,iw_Fr(idir))*w(ixI^S,iw_Fr(jdir)))
      ! end where

      !> Eddington tensor
      prad(ixI^S) = (3.d0*chi(ixI^S)-1)/2.d0 * fi(ixI^S)*fj(ixI^S)/ff(ixI^S)
      ! prad(ixI^S) = (3.d0*chi(ixI^S)-1)/2.d0 * ff(ixI^S)
      ! prad(ixI^S) = (3.d0*chi(ixI^S)-1)/2.d0 * sign(max(fi(ixI^S)*fj(ixI^S)/ff(ixI^S), smalldouble),w(ixI^S,iw_Fr(idir))*w(ixI^S,iw_Fr(idir)))

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
    chi(ixI^S) = sum(w(ixI^S, iw_Fr(:))**2, dim=ndim+1)/(m1_c*w(ixI^S,iw_r_e))**2

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

    double precision :: tau_c(ixI^S), fred(ixI^S)
    integer                         :: jxO^L, hxO^L

    hxO^L=ixO^L-kr(idim,^D);
    jxO^L=ixO^L+kr(idim,^D);

    !> Cell optical depth in idim
    tau_c(ixI^S) = block%dx(ixI^S,idim)*w(ixI^S,iw_rho)*w(ixI^S,kF_)

    !> Reduced flux
    fred(ixI^S) = dsqrt(sum(w(ixI^S, iw_Fr(:))**2,dim=ndim+1))/(m1_c*w(ixI^S,iw_r_e))

    ! theta(ixI^S) = 1.d0/(1.d0 + 2.d0/(3.d0*tau_c(ixI^S)))
    ! theta(ixI^S) = dexp(-2.d0/(3*tau_c(ixI^S)))
    ! theta(ixI^S) = 1.d0 - dexp(-tau_c(ixI^S))
    ! theta(ixI^S) = 1 - fred(ixI^S)
    ! theta(ixI^S) = (1-fred(ixI^S))/(1.d0 + 2.d0/(3.d0*tau_c(ixI^S)))
    theta(ixI^S) = (1-fred(ixI^S))*dexp(-2.d0/(3*tau_c(ixI^S))) !> THIS IS THE ONE
    ! theta(ixI^S) = (1-fred(ixI^S))*(1.d0 - dexp(-tau_c(ixI^S)))


    !> Neighbouring cells following
    theta(ixO^S) = min(theta(hxO^S),theta(ixO^S),theta(jxO^S))

    theta(ixI^S) = 1.d0

    !> Physical limits
    theta = max(0.d0, min(1.d0, theta))

    if (m1_opacity_type .eq. 'zero') theta = 0.d0
    if (.not. m1_split_trapped) theta = 0.d0

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

  ! Calculate flux f_idim[iw]
  subroutine m1_get_flux(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters
    ! use mod_geometry

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative w
    ! double precision, intent(in)    :: wC(ixI^S, 1:nw)
    ! primitive w
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout)   :: f(ixI^S, nwflux)
    double precision                :: prad_ii(ixI^S), ftrap(ixI^S), gradE(ixI^S)
    integer :: idir

    !> Radiation energy flux is v_i E + F_i
    f(ixO^S,iw_r_e) = m1_RSOL * (w(ixO^S,iw_mom(idim)) * w(ixO^S, iw_r_e) + w(ixO^S, iw_Fr(idim)))

    !> Radiation flux flux is v_i F_i + c^2 P_ii
    do idir = 1,ndir
      call m1_get_prad(w, x, ixI^L, ixO^L, idir, idim, prad_ii, .true.)
      f(ixO^S,iw_Fr(idir)) = m1_RSOL * (w(ixO^S,iw_mom(idim)) * w(ixO^S,iw_Fr(idir)) + m1_c**2 * prad_ii(ixO^S))
    end do

  end subroutine m1_get_flux

  subroutine m1_add_source_geom(qdt, ixI^L, ixO^L, wCT, w, x)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt, x(ixI^S, 1:ndim), wCT(ixI^S, 1:nw)
    double precision, intent(inout) :: w(ixI^S, 1:nw)

    double precision                :: chi(ixI^S), ff(ixI^S), fi(ixI^S), theta(ixI^S)
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

    select case (coordinate)
    case (cylindrical,Cartesian_expansion)
      call mpistop("m1 geom source terms not implemented yet with cylindrical geometries")
    case(spherical)
      h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
      select case (ndir)
      case (1)
        !> Get relevant pressure tensor elements
        call m1_get_streaming_fraction(wCT, x, ixI^L, ixO^L, 1, theta)
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 1, 1, prr, .true.)
        ptt(ixI^S) = m1_c**2 * (wCT(ixI^S, iw_r_e)*theta(ixI^S) - prr(ixI^S)) /2
        ppp(ixI^S) = ptt(ixI^S)

        w(ixO^S, iw_Fr(1)) = w(ixO^S, iw_Fr(1)) + qdt * m1_RSOL * (ptt(ixO^S) + ppp(ixO^S))/ x(ixO^S, 1)

        {^NOONED
      case (2)

        !> Get relevant pressure tensor elements
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 1, 1, prr, .true.)
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 2, 2, ptt, .true.)
        call m1_get_streaming_fraction(wCT, x, ixI^L, ixO^L, 1, theta)
        ppp = wCT(ixI^S,iw_r_e)*theta(ixI^S) - prr - ptt !> trace is one
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 2, 1, ptr, .true.)

        !> Add advection part
        ptt(ixI^S) = m1_c**2 * ptt(ixI^S) + wCT(ixI^S, iw_Fr(2))*wCT(ixI^S, iw_mom(2))/wCT(ixI^S, iw_rho)
        ptr(ixI^S) = m1_c**2 * ptr(ixI^S) + wCT(ixI^S, iw_Fr(2))*wCT(ixI^S, iw_mom(1))/wCT(ixI^S, iw_rho) !< IS THIS ONE CORRECT?
        ppp(ixI^S) = m1_c**2 * ppp(ixI^S)

        !> Perform corrections
        !> SWAP SIGNS AS ACCORDING TO WIKI: F = F - dt*(nabla T + S_corr)
        w(ixO^S, iw_Fr(1)) = w(ixO^S, iw_Fr(1)) + qdt * m1_RSOL * (ptt(ixO^S) + ppp(ixO^S))/ x(ixO^S, 1)
        w(ixO^S, iw_Fr(2)) = w(ixO^S, iw_Fr(2)) + qdt * m1_RSOL * (-ptr(ixO^S)/x(ixO^S, 1)  + ppp(ixO^S)/(x(ixO^S, 1)*dtan(x(ixO^S,2))))

      case(3)
        ! if ((ndim .eq. 1) .or. (ndim .eq. 2)) call mpistop("1.x D not implemented yet")

        !> Get relevant pressure tensor elements
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 2, 2, ptt, .true.)
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 3, 3, ppp, .true.)
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 2, 1, ptr, .true.)
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 3, 2, ppt, .true.)
        call m1_get_prad(wCT, x, ixI^L, ixO^L, 3, 1, ppr, .true.)

        !> Add advection part
        ptt = m1_c**2 * ptt + wCT(ixI^S, iw_Fr(2))*wCT(ixI^S, iw_mom(2))/wCT(ixI^S, iw_rho)
        ppp = m1_c**2 * ppp + wCT(ixI^S, iw_Fr(3))*wCT(ixI^S, iw_mom(3))/wCT(ixI^S, iw_rho)

        ptr = m1_c**2 * ptr + wCT(ixI^S, iw_Fr(2))*wCT(ixI^S, iw_mom(1))/wCT(ixI^S, iw_rho) !< IS THIS ONE CORRECT?
        ppr = m1_c**2 * ppr + wCT(ixI^S, iw_Fr(3))*wCT(ixI^S, iw_mom(1))/wCT(ixI^S, iw_rho) !< IS THIS ONE CORRECT?
        ppt = m1_c**2 * ppt + wCT(ixI^S, iw_Fr(3))*wCT(ixI^S, iw_mom(2))/wCT(ixI^S, iw_rho) !< IS THIS ONE CORRECT?

        !> Perform corrections
        !> SWAP SIGNS AS ACCORDING TO WIKI: F = F - dt*nabla T
        w(ixO^S, iw_Fr(1)) = w(ixO^S, iw_Fr(1)) + qdt * m1_RSOL * (ptt(ixO^S) + ppp(ixO^S))/ x(ixO^S, 1)
        w(ixO^S, iw_Fr(2)) = w(ixO^S, iw_Fr(2)) + qdt * m1_RSOL * (-ptr(ixO^S)/x(ixO^S, 1)  + ppp(ixO^S)/(x(ixO^S, 1)*dtan(x(ixO^S,2))))

        !> 2 ptp BELOW: P is symmetric
        w(ixO^S, iw_Fr(3)) = w(ixO^S, iw_Fr(3)) + qdt * m1_RSOL * (-ppr(ixO^S)/x(ixO^S, 1) - 2*ppt(ixO^S)/(x(ixO^S, 1)*dtan(x(ixO^S,2))))
        }

      end select

    end select

  end subroutine m1_add_source_geom

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

    double precision :: prad(ixI^S), gradv(ixI^S)
    double precision :: Fp(ixO^S), Fm(ixO^S), sm(ixO^S), sp(ixO^S)
    double precision :: geom_corr(ixO^S), vec_grad(ixO^S), vel(ixI^S,1:ndir)

    if(qsourcesplit .eqv. m1_force_split) then
      active = .true.

      do idir = 1,ndir
        vel(ixI^S,idir) = wCT(ixI^S,iw_mom(idir))/wCT(ixI^S,iw_rho)
      end do

      do idim = 1, ndim

        !> Photon tiring + Radflux CMF correction
        do idir = 1, ndir

          !TODO: implement vector gradient corrections in spherical!
          !TODO: implement directional derivative in spherical!
          !> Both are expressed at: https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
          ! call gradient_str(vel(ixI^S,idir), x, ixI^L, ixO^L, idim, gradv)
          call gradientS(vel(ixI^S,idir), ixI^L, ixO^L, idim, gradv)
          call m1_get_prad(wCT, x, ixI^L, ixO^L, idim, idir, prad)

          !> Geometric correction terms to vector gradient
          vec_grad(ixO^S) = gradv(ixO^S)

          select case(coordinate)
          case(cylindrical)
            call mpistop("vector gradient not implemented for cylindrical coords")
          case(spherical)
            if (idir .eq. 1) then
              {^NOONED
              if (idim .eq. 2) then
                vec_grad(ixO^S) = vec_grad(ixO^S) - vel(ixO^S,idim)/x(ixO^S,1)
              else if (idim .eq. 3) then
                vec_grad(ixO^S) = vec_grad(ixO^S) -  vel(ixO^S,idim)/x(ixO^S,1)
              end if
            else if (idir .eq. 2) then
              if (idim .eq. 2) then
                vec_grad(ixO^S) = vec_grad(ixO^S) +  vel(ixO^S,1)/x(ixO^S,1)
              else if (idim .eq. 3) then
                vec_grad(ixO^S) = vec_grad(ixO^S) -  vel(ixO^S,idim)/(x(ixO^S,1)*dtan(x(ixO^S,2)))
              end if
              }
              {^IFTHREED
            else if (idir .eq. 3) then
              if (idim .eq. 3) then
                vec_grad(ixO^S) = vec_grad(ixO^S) +  vel(ixO^S,1)/x(ixO^S,1)
                vec_grad(ixO^S) = vec_grad(ixO^S) +  vel(ixO^S,2)/(x(ixO^S,1)*dtan(x(ixO^S,2)))
              end if
              }
            end if
          end select

          w(ixO^S,iw_r_e) = w(ixO^S,iw_r_e) - qdt*m1_RSOL*vec_grad(ixO^S)*prad(ixO^S)
          w(ixO^S,iw_Fr(idir)) = w(ixO^S,iw_Fr(idir)) - qdt*m1_RSOL*wCT(ixO^S,iw_Fr(idim))*gradv(ixO^S)

        end do

        !> Geometric correction terms to directional derivative
        select case(coordinate)
        case(cylindrical)
          call mpistop("directional derivative not implemented for cylindrical coords")
        case(spherical)
          select case(idim)
            {^NOONED
          case(1)
            geom_corr(ixO^S) = -(wCT(ixO^S,iw_Fr(2)) * vel(ixO^S,2))/x(ixO^S,1)
            if (ndir > 2) then
              geom_corr(ixO^S) = geom_corr(ixO^S) - (wCT(ixO^S,iw_Fr(3))* vel(ixO^S,3))/x(ixO^S,1)
            end if
          case(2)
            geom_corr(ixO^S) = (wCT(ixO^S,iw_Fr(2)) * vel(ixO^S,1))/x(ixO^S,1)
            if (ndir > 2) then
              geom_corr(ixO^S) = geom_corr(ixO^S) - (wCT(ixO^S,iw_Fr(3))* vel(ixO^S,3))/(x(ixO^S,1)*dtan(x(ixO^S,2)))
            end if
            }
            {^IFTHREED
          case(3)
            geom_corr(ixO^S) = (wCT(ixO^S,iw_Fr(3)) * vel(ixO^S,1))/x(ixO^S,1)
            geom_corr(ixO^S) = geom_corr(ixO^S) + (wCT(ixO^S,iw_Fr(3)) * vel(ixO^S,2))/(x(ixO^S,1)*dtan(x(ixO^S,2)))
            }
          end select
          {^NOONED
          !> add correction
          w(ixO^S,iw_Fr(idim)) = w(ixO^S,iw_Fr(idim)) + qdt*m1_RSOL*geom_corr(ixO^S)
          }
        end select
      end do

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

    if(qsourcesplit .eqv. m1_force_split) then
      active = .true.

      do idir = 1,ndir
        !> Force in momentum
        ! TODO: Well-balanced scheme?!

        !> Update momentum
        w(ixO^S,iw_mom(idir)) = w(ixO^S,iw_mom(idir)) + qdt*wCT(ixO^S,iw_rho)*wCT(ixO^S,kF_)*wCT(ixO^S,iw_Fr(idir))/m1_c

        !> Force work term
        if (energy) &
        w(ixO^S,iw_e) = w(ixO^S,iw_e) + qdt*wCT(ixO^S,iw_mom(idir))*wCT(ixO^S,kF_)*wCT(ixO^S,iw_Fr(idir))/m1_c

      end do
    end if

  end subroutine m1_hydro_source

  subroutine m1_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_ghostcells_update
    type(state), target :: psa(max_blocks) !< Compute implicit part from this state and update it
    type(state), target :: psb(max_blocks) !< Will be unchanged, as on entry
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    integer :: ixO^L,iigrid,igrid

    ixO^L=ixG^LL^LSUBnghostcells;

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      call m1_flux_absorbtion(psa(igrid)%w, psb(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L, dtfactor, qdt);
    end do
    !$OMP END PARALLEL DO

    if (m1_heatcool) then
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        call m1_heating_cooling(psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L, dtfactor, qdt);
      end do
      !$OMP END PARALLEL DO
    end if

  end subroutine m1_implicit_update

  subroutine m1_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)
    double precision, intent(in) :: qtC

    integer :: ixO^L,iigrid,igrid

    call mpistop("mulit-step implicit IMEX requiring evaluate Not working yet")

    ixO^L=ixG^LL^LSUBnghostcells;

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)
      call m1_flux_absorbtion_evaluate(psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L);
    end do
    !$OMP END PARALLEL DO

    if (m1_heatcool) then
      !$OMP PARALLEL DO PRIVATE(igrid)
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        call m1_heating_cooling_evaluate(psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L);
      end do
      !$OMP END PARALLEL DO
    end if

    !> Set all other terms to zero
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      call m1_evaluate_rest(psa(igrid)%w, psa(igrid)%x, ixG^LL, ixO^L);
    end do
    !$OMP END PARALLEL DO


  end subroutine m1_evaluate_implicit

  subroutine m1_flux_absorbtion(w_out, w, x, ixI^L, ixO^L, dtfactor, qdt)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_set_F0

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: dtfactor, qdt
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: w_out(ixI^S,1:nw)

    double precision :: F0(ixI^S,1:ndir), theta(ixI^S), tau(ixI^S)
    double precision :: gradE(ixI^S)

    double precision :: aa(ixI^S), bb(ixI^S)

    integer :: idim
    integer                         :: jxO^L, hxO^L


    ! do idim = 1, ndim
    !
    !
    !   hxO^L=ixO^L-kr(idim,^D);
    !   jxO^L=ixO^L+kr(idim,^D);
    !
    !
    !   ! > Exact solution: dF/dt = - a F - b
    !   !                  <=> F = C1 * exp(-a*t) - b/a
    !   !                  w. C1 = F0 + b/a
    !
    !
    !   ! aa(ixI^S) = m1_RSOL*m1_c*w(ixI^S,iw_rho)*w(ixI^S,kF_)
    !   !> rho evaluated at cell face for well-balancing
    !   rho(ixO^S) = (w(hxO^S,iw_rho) + 2*w(ixO^S,iw_rho) + w(jxO^S,iw_rho))/4
    !   aa(ixO^S) = m1_RSOL*m1_c*rho(ixO^S)*w(ixO^S,kF_)

    ! > Exact solution: dF/dt = - a F - b
    !                  <=> F = C1 * exp(-a*t) - b/a
    !                  w. C1 = F0 + b/a

    ! do idim = 1, ndim
    !
    !
    !   hxO^L=ixO^L;
    !   jxO^L=ixO^L+kr(idim,^D);
    !
    !   aa(ixI^S) = m1_RSOL*m1_c*w(ixI^S,iw_rho)*w(ixI^S,kF_)
    !
    !   call gradientF(w(ixI^S,iw_r_e),x,ixI^L,ixO^L,idim,gradE)
    !
    !   !> Well-balanced flux
    !   aa(ixO^S) = 0.5d0*(aa(hxO^S) + aa(jxO^S))
    !   gradE(ixO^S) = 0.5d0*(gradE(hxO^S) + gradE(jxO^S))
    !
    !   call m1_get_trapped_fraction(w, x, ixI^L, ixO^L, idim, theta)
    !
    !   ! call gradientS(w(ixI^S,iw_r_e),ixI^L,ixO^L,idim,gradE)
    !   bb(ixI^S) = m1_RSOL*m1_c**2*gradE(ixI^S)*theta(ixI^S)/3.d0
    !
    !   w(ixO^S,iw_Fr(idim)) = (w(ixO^S,iw_Fr(idim)) + bb(ixO^S)/aa(ixO^S)) &
    !   * dexp(-aa(ixO^S)*dtfactor*dt) - bb(ixO^S)/aa(ixO^S)
    !
    ! end do

    w_out = w

    aa(ixI^S) = m1_c*w(ixI^S,iw_rho)*w_out(ixI^S,kF_)
    do idim = 1, ndim
      call m1_get_trapped_fraction(w, x, ixI^L, ixO^L, idim, theta)

      call gradientS(w(ixI^S,iw_r_e),ixI^L,ixO^L,idim,gradE)
      bb(ixI^S) = m1_c**2*gradE(ixI^S)*theta(ixI^S)/3.d0

      !> Analytic
      w_out(ixO^S,iw_Fr(idim)) = (w(ixO^S,iw_Fr(idim)) + bb(ixO^S)/aa(ixO^S)) &
      * dexp(-aa(ixO^S)*m1_RSOL*dtfactor*qdt) - bb(ixO^S)/aa(ixO^S)
    end do

  end subroutine m1_flux_absorbtion

  subroutine m1_flux_absorbtion_evaluate(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_set_F0

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: F0(ixI^S,1:ndir), theta(ixI^S) , tau(ixI^S)
    double precision :: gradE(ixI^S), R(ixI^S)

    double precision :: aa(ixI^S), bb(ixI^S)

    integer :: idim

    ! > Exact solution: dF/dt = - a F - b
    !                  <=> F = C1 * exp(-a*t) - b/a
    !                  w. C1 = F0 + b/a

    aa(ixI^S) = m1_RSOL*m1_c*w(ixI^S,iw_rho)*w(ixI^S,kF_)
    do idim = 1, ndim
      call m1_get_trapped_fraction(w, x, ixI^L, ixO^L, idim, theta)
      ! call gradient_str(w(ixI^S,iw_r_e)*lambda(ixI^S),x,ixI^L,ixO^L,idim,gradE)
      call gradientS(w(ixI^S,iw_r_e),ixI^L,ixO^L,idim,gradE)
      bb(ixI^S) = m1_RSOL*m1_c**2*gradE(ixI^S)*theta(ixI^S)/3.d0
      w(ixO^S,iw_Fr(idim)) = - aa(ixO^S)*w(ixO^S,iw_Fr(idim)) - bb(ixO^S)
    end do

  end subroutine m1_flux_absorbtion_evaluate


  subroutine m1_heating_cooling(w, x, ixI^L, ixO^L, dtfactor, qdt)
    use mod_global_parameters
    use mod_constants

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: dtfactor, qdt
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: ke(ixO^S)
    double precision :: a4(ixO^S), a0(ixO^S)
    double precision :: c0(ixO^S), c1(ixO^S)

    integer :: ix^D, idir

    !> kinetic energy
    ke(ixO^S) = 0.5d0 * sum(w(ixO^S, iw_mom(:))**2, dim=ndim+1) / w(ixO^S, iw_rho)

    !> Get gas internal energy
    w(ixO^S,iw_e) = w(ixO^S,iw_e) - ke(ixO^S)

    !> Set polynomial coefficients
    a4(ixO^S) = alpha_qdot*w(ixO^S,kB_)*(m1_gamma-1.d0)**4/w(ixO^S,iw_rho)**3*dtfactor*dt
    a0(ixO^S) = m1_c*w(ixO^S,kE_)*w(ixO^S,iw_rho)*dtfactor*qdt

    c1(ixO^S) = 1.d0/(a4(ixO^S) - m1_RSOL*a4(ixO^S)*a0(ixO^S)/(1.d0 + m1_RSOL*a0(ixO^S)))
    c0(ixO^S) = (a0(ixO^S)*w(ixO^S,iw_r_e)/(1.d0 + m1_RSOL*a0(ixO^S)) + w(ixO^S,iw_e))*c1(ixO^S)

    !> Solve polynomial
    {do ix^D = ixOmin^D,ixOmax^D\}
    w(ix^D,iw_e) = max(w(ix^D,iw_e), m1_small_e)
    call Bisection_method(w(ix^D,iw_e),w(ix^D,iw_r_e),c0(ix^D),c1(ix^D))
    {enddo\}

    !> Set radiation energy density
    w(ixO^S,iw_r_e) = (w(ixO^S,iw_r_e) + m1_RSOL*a4(ixO^S)*w(ixO^S,iw_e)**4.d0)/(1.d0 + m1_RSOL*a0(ixO^S))

    !> Get gas total energy
    w(ixO^S,iw_e) = w(ixO^S,iw_e) + ke(ixO^S)

  end subroutine m1_heating_cooling

  subroutine m1_heating_cooling_evaluate(w, x, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Bb(ixI^S)

    !> kinetic energy
    Bb(ixO^S) = 0.5d0 * sum(w(ixO^S, iw_mom(:))**2, dim=ndim+1) / w(ixO^S, iw_rho)

    !> pthermal pressure
    Bb(ixO^S) = (m1_gamma - 1.0d0) * (w(ixO^S, iw_e) - Bb(ixO^S))

    !> Gas temperature
    Bb(ixI^S) = beta_qdot*Bb(ixI^S)/w(ixO^S,iw_rho)

    !> Integrated black body mean intensity
    !TODO: Fix units here
    Bb(ixI^S) = const_sigma*Bb(ixI^S)**4/dpi
    Bb(ixI^S) = Bb(ixI^S)*unit_time/unit_pressure

    !> Heating and cooling sourceterms
    w(ixO^S,iw_e) = m1_c*w(ixO^S,kE_)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e) &
    - 4*dpi*w(ixO^S,kB_)*w(ixO^S,iw_rho)*Bb(ixO^S)

    w(ixO^S,iw_r_e) = - m1_RSOL*m1_c*w(ixO^S,kE_)*w(ixO^S,iw_rho)*w(ixO^S,iw_r_e) &
    + 4*dpi*m1_RSOL*w(ixO^S,kB_)*w(ixO^S,iw_rho)*Bb(ixO^S)

  end subroutine m1_heating_cooling_evaluate

  subroutine m1_evaluate_rest(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_geometry
    use mod_usr_methods, only: usr_set_F0

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout)    :: w(ixI^S,1:nw)

    !> Density and momentum don't have implicit sourceterms
    w(ixI^S,iw_rho) = 0.d0
    w(ixI^S,iw_mom(:)) = 0.d0

    !> in the case of no energy exchange
    if (.not. m1_heatcool) then
      w(ixI^S,iw_e) = 0.d0
      w(ixI^S,iw_r_e) = 0.d0
    end if

    !> In the case of ndim = ndir, final F-comp is not updated
    if (ndim .ne. ndir) then
      if (ndir .eq. 3) then !ndim = 2, ndir = 3
        w(ixI^S,iw_Fr(3)) = 0.d0
        if (ndim .eq. 1) w(ixI^S,iw_Fr(2)) = 0.d0 !ndim = 1, ndir = 3
      else if (ndir .eq. 2) then !ndim = 1, ndir = 2
        w(ixI^S,iw_Fr(2)) = 0.d0
      end if
    end if

  end subroutine m1_evaluate_rest

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
      !> GC's for cbounds
      w(ixI^S,kF_) = m1_kappa_e

      !> Calculate velocity gradient
      vr(ixI^S) = w(ixI^S,iw_mom(1))/w(ixI^S,iw_rho)
      call gradient_str(vr,x,ixI^L,ixO^L,1,dvrdr)

      !> Calculate sobolev optical depth parameter
      tau_sob(ixI^S) = dabs(m1_c*m1_kappa_e*w(ixI^S,iw_rho)/(dvrdr(ixI^S)+smalldouble))

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

      ! !> ghostcells on opacity for force-balancing
      ! {^IFONED
      ! w(ixOmin1-1,kF_) = w(ixOmin1,kF_)
      ! w(ixOmax1+1,kF_) = w(ixOmax1,kF_)
      ! }
      ! {^IFTWOD
      ! w(ixOmin1-1,:,kF_) = w(ixOmin1,:,kF_)
      ! w(ixOmax1+1,:,kF_) = w(ixOmax1,:,kF_)
      ! }
      ! {^IFTHREED
      ! w(ixOmin1-1,:,:,kF_) = w(ixOmin1,:,:,kF_)
      ! w(ixOmax1+1,:,:,kF_) = w(ixOmax1,:,:,kF_)
      ! }

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

    !> GC's for cbounds
    w(ixI^S,kF_) = m1_kappa_e

    !> Velocity stretched opacities
    !> Calculate velocity gradient
    vr(ixI^S) = w(ixI^S,iw_mom(1))/w(ixI^S,iw_rho)
    call gradient_str(vr,x,ixI^L,ixO^L,1,dvrdr)

    !> Calculate sobolev optical depth parameter
    tau_sob(ixO^S) = dabs(m1_c*m1_kappa_e*w(ixO^S,iw_rho)/(dvrdr(ixO^S)+smalldouble))

    w(ixO^S,kF_) = m1_kappa_e * Qb_F/(1-a_F) * ((1.d0+Q0_F*tau_sob(ixO^S))**(1.d0-a_F) - 1.d0)/(Q0_F*tau_sob(ixO^S))
    w(ixO^S,kE_) = m1_kappa_e * Qb_E/(1-a_E)*((1+Q0_E*tau_sob(ixO^S))**(1-a_E) - 1)/(Q0_E*tau_sob(ixO^S))
    w(ixO^S,kB_) = m1_kappa_e * Qb_B/(1-a_B)*((1+Q0_B*tau_sob(ixO^S))**(1-a_B) - 1)/(Q0_B*tau_sob(ixO^S))

    !> Opal Opacities
    {do ix^D=ixOmin^D,ixOmax^D\ }
    rho = w(ix^D,iw_rho)*unit_density
    Trad = (w(ix^D,iw_r_e)*unit_pressure/const_rad_a)**(0.25d0)
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

    !> Inner bound opacity
    if (any(x(ixI^S,1) .lt. xprobmin1)) then
      {^IFONED w(ixImin1:ixOmin1-1,kF_) = w(ixOmin1,kF_)}
      {^IFTWOD
          do ix2 = ixImin2, ixImax2
            w(ixImin1:ixOmin1-1,ix2,kF_) = w(ixOmin1,ix2,kF_)
          end do
      }
      {^IFTHREED
      do ix2 = ixImin2, ixImax2
        do ix3 = ixImin3, ixImax3
          w(ixImin1:ixOmin1-1,ix2,ix3,kF_) = w(ixOmin1,ix2,ix3,kF_)
        end do
      end do
      }
      ! {^IFTHREED w(ixOmax1,:,:,kF_) = max(w(ixOmax1-1,:,:,kF_), w(ixOmax1,:,:,kF_))}
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

    !> GC's for cbounds
    w(ixI^S,kF_) = m1_kappa_e

    !> Velocity stretched opacities
    !> Calculate velocity gradient
    vr(ixI^S) = w(ixI^S,iw_mom(1))/w(ixI^S,iw_rho)
    ! call gradient_str(vr,x,ixI^L,ixO^L,1,dvrdr)
    call gradientS(vr,ixI^L,ixO^L,1,dvrdr)

    !if ((modulo(it,1000) == 0) .and. (mype .eq. 0)) then
    !  print*, 'WST photosphere: ', mype, i_phot, R_phot
    !end if

    !> Opal Opacities
    {do ix^D=ixOmin^D,ixOmax^D\ }
      rho = w(ix^D,iw_rho)*unit_density
      Trad = (w(ix^D,iw_r_e)*unit_pressure/const_rad_a)**(0.25d0)
      Trad = max(Trad,1.d4)
      ! gradv = max(dvrdr(ix^D),smalldouble)*unit_velocity/unit_length
      gradv = max(dabs(dvrdr(ix^D)),smalldouble)*unit_velocity/unit_length

      call set_cak_opacity(rho, Trad, gradv ,alpha, Qbar, Q0, kappa_e_t)

      !> Finite disc correction factor
      beta_fd = ( 1.0d0 - vr(ix^D)/(x(ix^D,1) * max(abs(dvrdr(ix^D)),smalldouble)) ) &
      ! beta_fd = ( 1.0d0 - vr(ix^D)/(x(ix^D,1) * abs(dvrdr(ix^D))) ) &
      * (R_phot/x(ix^D,1))**2

      if (x(ix^D,1) .le. R_phot) then
        fdfac = 1.d0/(1.d0 + alpha)
      else if (beta_fd >= 1.0d0) then
        fdfac = 1.0d0/(1.0d0 + alpha)
      else if (beta_fd < -1.0d10) then
        fdfac = abs(beta_fd)**alpha / (1.0d0 + alpha)
        ! fdfac = 1.d10**alpha / (1.0d0 + alpha)
        ! fdfac = 1.d0
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

      ! print*, x(ix^D,1)/R_phot, fdfac, tau, CAK(ix^D)/fdfac/m1_kappa_e, beta_fd, alpha, Qbar, Q0

      ! print*, block%igrid, x(ix^D,1), alpha, Qbar, Q0, gradv, tau, fdfac
    {enddo\ }

    !> Combine OPAL and CAK-like
    w(ixO^S,kF_) = min(CAK(ixO^S), 50.d0*m1_kappa_e) + OPAL(ixO^S)
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

    if (any(x(ixI^S,1) .lt. xprobmin1)) then
      {^IFONED w(ixImin1:ixOmin1-1,kF_) = w(ixOmin1,kF_)}
      {^IFTWOD
          do ix2 = ixImin2, ixImax2
            w(ixImin1:ixOmin1-1,ix2,kF_) = w(ixOmin1,ix2,kF_)
          end do
      }
      {^IFTHREED
      do ix2 = ixImin2, ixImax2
        do ix3 = ixImin3, ixImax3
          w(ixImin1:ixOmin1-1,ix2,ix3,kF_) = w(ixOmin1,ix2,ix3,kF_)
        end do
      end do
      }
      ! {^IFTHREED w(ixOmax1,:,:,kF_) = max(w(ixOmax1-1,:,:,kF_), w(ixOmax1,:,:,kF_))}
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
    Fmag2(ixI^S) = sum(w(ixI^S,iw_Fr(:))**2, dim=ndim+1)

    !> Flag negative values of Erad
    where (w(ixO^S, iw_r_e) < m1_small_e) flag(ixO^S, iw_r_e) = .true.

    !> Flag flux values where |F| > cE
    !> Flag all dirs
    do idir =1,ndir
      where (Fmag2(ixO^S) > (m1_c*w(ixO^S, iw_r_e))**2)
        flag(ixO^S,iw_Fr(idir)) = .true.
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
        e_a = m1_small_e
        e_b = m1_small_e
        e_c = m1_small_e

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

        ! call mpistop('mehhh')

        e_a = e_c
        e_b = e_c

        ! print*, 'heatcool:', e_c

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

end module mod_m1
