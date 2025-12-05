!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd

  implicit none

  double precision :: L, rho, T_bg, T_beam, rho_clump

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    ! usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions
    usr_special_bc => boundary_conditions
    usr_special_opacity_m1 => opacity_law

    ! Choose independent normalization units if using dimensionless variables.
    ! unit_length      = 1.d0 ! cm
    ! unit_temperature = 1.d0 ! K
    ! unit_density     = 1.d0 ! g cm^-3

    unit_length      = L    ! cm
    unit_temperature = T_bg ! K
    unit_density     = rho ! g cm^-3


    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /stream_list/ L, rho, T_bg, T_beam, rho_clump

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, stream_list, end=111)
       111 close(unitpar)
    enddo

  end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Delta(ixI^S)

      !> clump structure
      Delta(ixI^S) = 10.d0*( (x(ixI^S,1)/0.1d0)**2.d0 - 1.d0)
      Delta = min(Delta,1.d2)

      !> densities
      w(ixI^S,rho_) = rho + (rho_clump - rho) /(1.d0 + dexp(Delta(ixI^S)))
      w(ixI^S,mom(:)) = 0.d0
      w(ixI^S,e_) = const_kB*T_bg*w(ixI^S,rho_)/(const_mp*hd_mu)/(hd_gamma - 1.d0)
      w(ixI^S,Er_) = const_rad_a*T_bg**4
      w(ixI^S,Fr(:)) = 0.d0

      !> Units
      w(ixI^S,rho_) = w(ixI^S,rho_)/unit_density
      w(ixI^S,e_) = w(ixI^S,e_)/unit_pressure
      w(ixI^S,Er_) = w(ixI^S,Er_)/unit_pressure

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    use mod_hd

    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: ii

    select case (iB)

    case(1)

      w(ixB^S,Er_) = const_rad_a*T_beam**4/unit_pressure
      w(ixB^S,Fr(1)) = const_rad_a*T_beam**4/unit_pressure*m1_c

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

subroutine opacity_law(ixI^L,ixO^L,w,x)
  use mod_global_parameters
  use mod_hd

  integer, intent(in)             :: ixI^L, ixO^L
  double precision, intent(in)    :: x(ixI^S,1:ndim)
  double precision, intent(inout) :: w(ixI^S,1:nw)

  double precision :: Tgas(ixI^S), kappa(ixI^S)

  Tgas(ixI^S) = w(ixI^S, e_)*(hd_gamma - 1)*unit_pressure
  Tgas(ixI^S) = Tgas(ixI^S)*const_mp*hd_mu/(w(ixI^S,rho_)*unit_density*const_kb)

  kappa(ixI^S) = 1.d2*(Tgas(ixI^S)/T_bg)**(-3.5d0)*w(ixI^S,rho_)

  w(ixI^S,kF_) = kappa(ixI^S)/unit_opacity
  w(ixI^S,kE_) = kappa(ixI^S)/unit_opacity
  w(ixI^S,kB_) = kappa(ixI^S)/unit_opacity

end subroutine opacity_law

end module mod_usr
