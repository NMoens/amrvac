!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd

  implicit none

  ! Custom global variables can be defined here
  ! User input parameters
  double precision :: mstar, rstar, T_inf, rho_inf
  double precision :: mdot, asound, Ggrav

  double precision :: lstar, mumol, gammae, logg

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
    unit_temperature = T_inf    ! K
    unit_density     = rho_inf ! g cm^-3

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /star_list/ mstar, rstar, T_inf, rho_inf, gammae

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

    ! Stellar structure
    lstar  = gammae*4.0d0*dpi * const_G * mstar * const_c/0.34d0
    logg   = log10(const_G * mstar/rstar**2.0d0)
    mumol  = (1.0d0 + 4.0d0*He_abundance)/(2.0d0 + 3.0d0*He_abundance)
    asound = sqrt(T_inf * kB_cgs/(mumol * mp_cgs))

    ! Make some variables needed for computations dimensionless
    mstar    = mstar/(unit_density * unit_length**3.0d0)
    rstar    = rstar/unit_length
    asound   = asound/unit_velocity
    rho_inf = rho_inf/unit_density
    Ggrav    = const_G * unit_density * unit_time**2.0d0

    if (mype == 0) then
      print*, '========================================'
      print*, '  Dimensionless computation quantities  '
      print*, '========================================'
      print*, 'Mstar    = ', mstar
      print*, 'Rstar    = ', rstar
      print*, 'rhobound = ', rho_inf
      print*, 'asound   = ', asound
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

    !> Gas subsystem
    w(ixO^S,rho_) = rho_inf
    w(ixO^S,mom(1)) = 0.d0

    !> Radiation subsystem
    w(ixO^S,Fr(1)) = lstar/(4*dpi*(x(ixO^S,1)*unit_length)**2)/unit_radflux
    w(ixO^S,Er_) = w(ixO^S,Fr(1))/const_c*unit_velocity

    call hd_to_conserved(ixI^L,ixO^L,w,x)

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

      do ix1 = ixBmax1, ixBmin1, -1
        w(ix1,rho_) = w(ix1+1,rho_)
        w(ix1,mom(1)) = w(ix1+1,mom(1))

        w(ix1,Fr(1)) = lstar/(4*dpi*(x(ix1,1)*unit_length)**2)/unit_radflux
        w(ix1,Er_) = w(ix1,Fr(1))/const_c*unit_velocity
      end do

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
