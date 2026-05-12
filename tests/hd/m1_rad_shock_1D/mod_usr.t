!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd
  use mod_m1

  implicit none

  double precision :: rho1, rho2, v1, v2, p1, p2, Er1, Er2
  double precision :: my_unit_density, mu

  integer :: i_massflux, i_momflux, i_eflux, i_trapped

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    !> Additional variables
    usr_process_grid => update_extravars

    ! Choose independent normalization units if using dimensionless variables.
    unit_time        = 1.d0   ! cm
    unit_velocity    = const_c !
    unit_density     = my_unit_density !100.d0 + smalldouble! g cm^-3

    call hd_activate()

    i_massflux = var_set_extravar('massflux','massflux')
    i_momflux = var_set_extravar('momflux','momflux')
    i_eflux = var_set_extravar('eflux','eflux')
    i_trapped = var_set_extravar('trapped','trapped')


  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /shock_list/ rho1, rho2, v1, v2, p1, p2, Er1, Er2

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, shock_list, end=111)
       111 close(unitpar)
    enddo

    !> Set unit_density such that Trad = Tgas
    mu = (1.0d0 + 4.0d0*He_abundance)/(2.0d0 + 3.0d0*He_abundance)
    my_unit_density = 4*const_sigma/Er1 * (const_mp*mu/const_kB)**4 * (p1/rho1)**4 * const_c**5

  end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_constants

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    where (x(ixI^S,1) .lt. 0.d0)
      w(ixI^S,rho_) = rho1
      w(ixI^S,mom(1)) = v1
      w(ixI^S,p_) = p1
      w(ixI^S,Er_) = Er1
      w(ixI^S,Fr(1)) = 0.d0 !1.d-2*Er1 !- v1*Er1
    elsewhere
      w(ixI^S,rho_) = rho2
      w(ixI^S,mom(1)) = v2
      w(ixI^S,p_) = p2
      w(ixI^S,Er_) = Er2
      w(ixI^S,Fr(1)) = 0.d0 !1.d-2*Er2 !- v2*Er2
    end where

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    use mod_hd

    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: prad(ixI^S), v(ixI^S), pth(ixI^S)
    double precision :: T(ixI^S), divT(ixI^S), F(ixI^S,1:ndir)
    double precision :: theta(ixI^S)

    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, prad)
    call hd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    v(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)

    w(ixO^S,i_massflux) = w(ixO^S,mom(1))
    w(ixO^S,i_momflux) = w(ixO^S,rho_)*v(ixO^S)**2 + pth(ixO^S) &
                       + prad(ixO^S) - v(ixO^S)*w(ixO^S,Fr(1))/m1_c**2
    w(ixO^S,i_eflux) = w(ixO^S,e_)*v(ixO^S) + pth(ixO^S)*v(ixO^S) &
                       + w(ixO^S,Fr(1)) - w(ixO^S,Er_)*v(ixO^S) - prad(ixO^S)*v(ixO^S)

    call m1_get_trapped_fraction(w, x, ixI^L, ixO^L, 1, theta)

    w(ixO^S,i_trapped) = theta(ixO^S)

  end subroutine update_extravars


end module mod_usr
