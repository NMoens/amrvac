!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd
  use mod_m1

  implicit none

  double precision :: L, rho, Er, Bg_E
  integer :: i_Pxx,i_Pxy,i_Pyx,i_Pyy
  integer :: i_r, i_sol, i_trap

  double precision :: r0 = 1.d0

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("spherical_2D")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    ! usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions

    !> Additional variables
    usr_process_grid => update_extravars

    ! Choose independent normalization units if using dimensionless variables.
    ! unit_length      = 1.d0 ! cm
    ! unit_temperature = 1.d0 ! K
    ! unit_density     = 1.d0 ! g cm^-3

    unit_length      = L    ! cm
    unit_temperature = (Bg_E/const_rad_a)**0.25d0  ! K
    unit_density     = rho ! g cm^-3


    rho = rho/unit_density
    Bg_E = Bg_E/unit_pressure

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

    i_Pxx = var_set_extravar('P_rr','P_rr')
    i_Pxy = var_set_extravar('P_rt','P_rt')
    i_Pyx = var_set_extravar('P_tr','P_tr')
    i_Pyy = var_set_extravar('P_tt','P_tt')
    i_r = var_set_extravar('r','r')

    i_trap = var_set_extravar('trap','trap')

    i_sol = var_set_extravar('sol','sol')

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /stream_list/ L, rho, Er, Bg_E

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
    double precision :: sol(ixI^S)

    call solution(ixI^L, ixO^L, x, 0.d0, sol)

    w(ixI^S,rho_) = rho
    w(ixI^S,mom(:)) = 0.d0
    w(ixI^S,Er_) = sol(ixI^S) !Bg_E*(1+Er*dexp(-distance2(ixI^S)/0.01d0))/unit_pressure
    w(ixI^S,Fr(:)) = 0.d0

  end subroutine initial_conditions

  subroutine special_bound(qt,ixI^L,ixB^L,iB,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixB^S,rho_) = rho/unit_density
    w(ixB^S,mom(:)) = 0.d0
    w(ixB^S,Er_) = Bg_E
    w(ixB^S,Fr(:)) = 0.d0

  end subroutine special_bound

  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    use mod_hd

    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: prad(ixI^S), r0, sol(ixI^S), theta(ixI^S)

    r0 = 3.d0
    w(ixI^S,i_r) = dsqrt(x(ixI^S,1)**2 + r0**2 - 2*x(ixI^S,1)*r0*cos(dpi/2-x(ixI^S,2)))

    ! w(ixI^S,i_ff) = sum(w(ixI^S, Fr(:))**2, dim=ndim+1)/(m1_c*w(ixI^S,Er_))**2
    ! w(ixI^S,i_chi) = (3.d0+4*w(ixI^S,i_ff))/(5.d0+2*dsqrt(4.d0-3*w(ixI^S,i_ff)))

    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, prad)
    w(ixO^S,i_Pxx) = prad(ixO^S)/w(ixO^S,Er_)
    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 2, prad)
    w(ixO^S,i_Pxy) = prad(ixO^S)/w(ixO^S,Er_)
    call m1_get_prad(w, x, ixI^L, ixO^L, 2, 1, prad)
    w(ixO^S,i_Pyx) = prad(ixO^S)/w(ixO^S,Er_)
    call m1_get_prad(w, x, ixI^L, ixO^L, 2, 2, prad)
    w(ixO^S,i_Pyy) = prad(ixO^S)/w(ixO^S,Er_)

    call m1_get_trapped_fraction(w, x, ixI^L, ixO^L, 1, theta)
    w(ixO^S,i_trap) = theta(ixO^S)

    call solution(ixI^L, ixO^L, x, global_time, sol)
    w(ixO^S,i_sol) = sol(ixO^S)

    w(ixI^S,i_r) = x(ixI^S,1)**2 + r0**2 - 2*x(ixI^S,1)*r0*cos(dpi/2-x(ixI^S,2))

  end subroutine update_extravars

  subroutine solution(ixI^L, ixO^L, x, time, sol)
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim), time
    double precision, intent(inout) :: sol(ixI^S)

    double precision :: D, t
    double precision :: G(ixI^S), distance2(ixI^S)

    t = time + 1.d-4
    D = m1_RSOL*m1_c/(3*rho*m1_kappa_e)

    distance2(ixI^S) = x(ixI^S,1)**2 + r0**2 - 2*x(ixI^S,1)*r0*cos(dpi/2-x(ixI^S,2))

    G(ixI^S) = 1/(4*dpi*D*t) * dexp(-distance2(ixI^S)/(4*D*t))

    sol(ixI^S) = Bg_E*(1+G(ixI^S))

  end subroutine solution

end module mod_usr
