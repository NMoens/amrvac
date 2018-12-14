!> Nicolas Moens
!> Module for including flux limited diffusion in hydrodynamics simulations
!> Based on Turner and stone 2001
module mod_fld
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
    double precision, public :: fld_bisect_tol = 1.d-3

    !> Tolerance for adi method for radiative Energy diffusion
    double precision, public :: fld_adi_tol = 1.d-2

    double precision :: fld_max_fracdt = 50.d0

    !> Switch different terms on/off
    !> Solve parabolic system using ADI (Diffusion)
    logical :: fld_Diffusion = .true.

    !> Use radiation force sourceterm
    logical :: fld_Rad_force = .true.

    !> Use Heating and Cooling sourceterms in e_
    logical :: fld_Energy_interact = .true.

    !> Let Vac advect radiative energy
    logical :: fld_Energy_advect = .true.

    !> Use constant Opacity?
    character(len=8) :: fld_opacity_law = 'const'

    !> Diffusion limit lambda = 0.33
    logical :: fld_complete_diffusion_limit = .false.

    !> Boundary conditions for radiative Energy in ADI.
    character(len=8) :: fld_bound_min1 = 'periodic'
    character(len=8) :: fld_bound_max1 = 'periodic'
    character(len=8) :: fld_bound_min2 = 'periodic'
    character(len=8) :: fld_bound_max2 = 'periodic'

    !> Set Diffusion coefficient to unity
    logical :: fld_diff_testcase = .false.

    !> public methods
    !> these are called in mod_hd_phys
    public :: fld_add_source
    public :: fld_get_flux
    public :: fld_get_csound2
    !> these are used in specialvar_output
    public :: fld_get_radflux
    public :: fld_get_radpress
    public :: fld_get_fluxlimiter
    public :: fld_get_opacity

  contains

  subroutine fld_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /fld_list/ fld_kappa0, fld_mu, fld_split, fld_maxdw, fld_Diffusion,&
    fld_Rad_force, fld_Energy_interact, fld_Energy_advect, fld_bisect_tol, fld_diff_testcase,&
    fld_bound_min1, fld_bound_max1, fld_bound_min2, fld_bound_max2, fld_adi_tol, fld_max_fracdt,&
    fld_opacity_law, fld_complete_diffusion_limit

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, fld_list, end=111)
       111    close(unitpar)
    end do
  end subroutine fld_params_read

  subroutine fld_init()
    use mod_global_parameters
    use mod_variables

    !> read par files
    call fld_params_read(par_files)
    !call params_read

    !> Check if fld_numdt is not 1
    if (fld_maxdw .lt. 2) call mpistop("fld_maxdw should be an integer larger than 1")

    !> Make kappa dimensionless !!!STILL NEED TO MULTIPLY W RHO
    fld_kappa0 = fld_kappa0*unit_time*unit_velocity*unit_density

    !> Dimensionless speed of light
    fld_speedofligt_0 = const_c/unit_velocity

    !> Dimensionless Boltzman constante sigma
    fld_sigma_0 = const_sigma*(unit_temperature**4.d0)/(unit_velocity*unit_pressure)
  end subroutine fld_init


  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine fld_rad_force(qdt,ixI^L,ixO^L,wCT,w,x,&
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

    double precision :: fld_kappa(ixO^S)

    double precision :: rad_flux(ixO^S,1:ndim)
    double precision :: radiation_force(ixO^S,1:ndim)

    integer :: idir, i, jx^L

    !> Calculate and add sourceterms
    if(qsourcesplit .eqv. fld_split) then
      active = .true.

      !> Add momentum sourceterms
      if (fld_Rad_force) then
        jx^L=ixO^L+kr(idir,^D);
        !> Calculate the radiative flux using the FLD Approximation
        call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
        call fld_get_radflux(wCT, x, ixI^L, ixO^L, rad_flux)

        do idir = 1,ndir
          !> Radiation force = kappa*rho/c *Flux
          radiation_force(ixO^S,idir) = fld_kappa(ixO^S)*wCT(ixO^S,iw_rho)/fld_speedofligt_0*rad_flux(ixO^S, idir)

          !> Momentum equation source term
          w(ixO^S,iw_mom(idir)) = w(ixO^S,iw_mom(idir)) &
              + qdt * half*(radiation_force(ixO^S,idir) + radiation_force(jx^S,idir))
              !> NOT SURE ON HOW TO AVERAGE OVER LEFTHANDSIDE AND RIGHTHANDSIDE FLUX EDGE
        enddo
      endif
    end if
  end subroutine fld_rad_force

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine fld_energy_interact(qdt,ixI^L,ixO^L,wCT,w,x,&
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
      if (fld_Energy_interact) then
        call Energy_interaction(w, x, ixI^L, ixO^L)
      endif

    end if
  end subroutine fld_energy_interact

  !> w[iw]=w[iw]+qdt*S[wCT,qtC,x] where S is the source based on wCT within ixO
  subroutine fld_radiation_diffusion(qdt,ixI^L,ixO^L,wCT,w,x,&
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

      !> Begin by evolving the radiation energy field
      if (fld_Diffusion) then
        call Evolve_E_rad(w, x, ixI^L, ixO^L)
      endif

    end if
  end subroutine fld_radiation_diffusion


  subroutine fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: fld_kappa(ixO^S)
    double precision :: Temp(ixI^S)
    double precision :: rho0,Temp0,n,sigma_b

    integer :: i

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
      case default
        call mpistop("Doesn't know opacity law")
      end select

  end subroutine fld_get_opacity


  !> Calculate fld flux limiter
  subroutine fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out) :: fld_R(ixO^S), fld_lambda(ixO^S)
    double precision :: fld_kappa(ixO^S)
    double precision ::  normgrad2(ixO^S)
    double precision :: grad_r_e(ixO^S)
    integer :: idir

    if (fld_complete_diffusion_limit) then
      fld_lambda = one/3.d0
    else
      !> Calculate R everywhere
      !> |grad E|/(rho kappa E)
      normgrad2(ixO^S) = zero
      do idir = 1,ndir
        call grad(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,x,grad_r_e)
        normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S)**2
      end do

      call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
      fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(fld_kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,r_e))

      !> Calculate the flux limiter, lambda
      !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
      fld_lambda(ixO^S) = (2+fld_R(ixO^S))/(6+3*fld_R(ixO^S)+fld_R(ixO^S)**2)
    endif
  end subroutine fld_get_fluxlimiter


  !> Calculate Radiation Flux
  !> Returns Radiation flux and radiation pressure
  subroutine fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: rad_flux(ixO^S, 1:ndim)
    double precision             :: L_star, R_star
    double precision :: fld_kappa(ixO^S)
    double precision :: fld_lambda(ixO^S), fld_R(ixO^S), normgrad2(ixO^S), f(ixO^S)
    double precision :: grad_r_e(ixO^S, 1:ndim)
    integer :: idir

    !> Calculate R everywhere
    !> |grad E|/(rho kappa E)
    normgrad2(ixO^S) = zero
    do idir = 1,ndir
      call grad(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,x,grad_r_e(ixO^S,idir))
      normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S,idir)**2
    end do

    call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
    fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(fld_kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,r_e))

    !> Calculate the flux limiter, lambda
    !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
    if (fld_complete_diffusion_limit) then
      fld_lambda = one/3.d0
    else
      fld_lambda(ixO^S) = (two+fld_R(ixO^S))/(6.d0+3.d0*fld_R(ixO^S)+fld_R(ixO^S)**two)
    endif

    !> Calculate the Flux using the fld closure relation
    !> F = -c*lambda/(kappa*rho) *grad E
    do idir = 1,ndir
      call gradE(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,x,grad_r_e(ixO^S,idir))
      rad_flux(ixO^S, idir) = -fld_speedofligt_0*fld_lambda(ixO^S)/(fld_kappa(ixO^S)*w(ixO^S,iw_rho)) *grad_r_e(ixO^S,idir)
    end do

    ! !> Cheaty
    ! L_star = 2724846166.4085770d0
    ! R_star = 252.29283539564574d0
    !
    ! rad_flux(:, ixOmin2, 2) = L_star/(4.d0*dpi*R_star**2.d0)
  end subroutine fld_get_radflux


  !> Calculate Radiation Pressure
  !> Returns Radiation Pressure
  subroutine fld_get_radpress(w, x, ixI^L, ixO^L, rad_pressure)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: rad_pressure(ixO^S)
    double precision :: fld_kappa(ixO^S)
    double precision :: fld_lambda(ixO^S), fld_R(ixO^S), normgrad2(ixO^S), f(ixO^S)
    double precision :: grad_r_e(ixO^S, 1:ndim)
    integer :: idir

    !> Calculate R everywhere
    !> |grad E|/(rho kappa E)
    normgrad2(ixO^S) = zero

    do idir = 1,ndir
      call grad(w(ixI^S, iw_r_e),ixI^L,ixO^L,idir,x,grad_r_e(ixO^S,idir))
      normgrad2(ixO^S) = normgrad2(ixO^S) + grad_r_e(ixO^S,idir)**two
    end do

    call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
    fld_R(ixO^S) = dsqrt(normgrad2(ixO^S))/(fld_kappa(ixO^S)*w(ixO^S,iw_rho)*w(ixO^S,r_e))

    !> Calculate the flux limiter, lambda
    !> Levermore and Pomraning: lambda = (2 + R)/(6 + 3R + R^2)
    if (fld_complete_diffusion_limit) then
      fld_lambda = one/3.d0
    else
      fld_lambda(ixO^S) = (two+fld_R(ixO^S))/(6.d0+3.d0*fld_R(ixO^S)+fld_R(ixO^S)**two)
    endif

    !> Calculate radiation pressure
    !> P = (lambda + lambda^2 R^2)*E
    f(ixO^S) = fld_lambda(ixO^S) + fld_lambda(ixO^S)**two * fld_R(ixO^S)**two
    f(ixO^S) = one/two*(one-f(ixO^S)) + one/two*(3.d0*f(ixO^S) - one)

    rad_pressure(ixO^S) = f(ixO^S) * w(ixO^S, iw_r_e)
  end subroutine fld_get_radpress


  subroutine fld_get_flux(w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, idim
    double precision, intent(in) :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(inout) :: f(ixI^S, nwflux)

    if (fld_Energy_advect) then
      f(ixO^S, iw_r_e) = w(ixO^S,iw_mom(idim)) * w(ixO^S, iw_r_e)
    else
      f(ixO^S, iw_r_e) = zero
    endif
  end subroutine fld_get_flux



  subroutine grad(q,ixI^L,ixO^L,idir,x,gradq)
    ! Compute the true gradient of a scalar q within ixL in direction idir ie :
    !  - in cylindrical : d_/dr , (1/r)*d_/dth , d_/dz
    !  - in spherical   : d_/dr , (1/r)*d_/dth , (1/rsinth)*d_/dphi
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, idir
    double precision, intent(in) :: q(ixI^S), x(ixI^S,1:ndim)
    double precision, intent(out) ::  gradq(ixO^S)

    integer :: jx^L, hx^L
    !-----------------------------------------------------------------------------
    jx^L=ixO^L+kr(idir,^D);
    hx^L=ixO^L-kr(idir,^D);

    gradq(ixO^S)=(q(jx^S)-q(hx^S))/(x(jx^S,idir)-x(hx^S,idir))

    select case (typeaxial)
    case('slab') ! nothing to do
    case('cylindrical')
      if (idir==phi_) gradq(ixO^S)=gradq(ixO^S)/ x(ixO^S,r_)
    case('spherical')
      if (idir==2   ) gradq(ixO^S)=gradq(ixO^S)/ x(ixO^S,r_)
      if (idir==phi_) gradq(ixO^S)=gradq(ixO^S)/(x(ixO^S,r_)*dsin(x(ixO^S,2)))
    case default
      call mpistop('Unknown geometry')
    end select
  end subroutine grad

  subroutine gradE(q,ixI^L,ixO^L,idir,x,gradq)
    ! Compute the true gradient of a scalar q within ixL in direction idir ie :
    !  - in cylindrical : d_/dr , (1/r)*d_/dth , d_/dz
    !  - in spherical   : d_/dr , (1/r)*d_/dth , (1/rsinth)*d_/dphi
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, idir
    double precision, intent(in) :: q(ixI^S), x(ixI^S,1:ndim)
    double precision, intent(out) ::  gradq(ixO^S)

    integer :: jx^L
    !-----------------------------------------------------------------------------
    jx^L=ixO^L+kr(idir,^D);

    gradq(ixO^S)=(q(jx^S)-q(ixO^S))/(x(jx^S,idir)-x(ixO^S,idir))
  end subroutine gradE

end module mod_fld
