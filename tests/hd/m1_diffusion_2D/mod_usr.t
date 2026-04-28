!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd
  use mod_m1

  implicit none

  double precision :: L, rho, Er, Bg_E
  integer :: i_Pxx,i_Pxy,i_Pyx,i_Pyy

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian")
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




    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

    i_Pxx = var_set_extravar('P_xx','P_xx')
    i_Pxy = var_set_extravar('P_xy','P_xy')
    i_Pyx = var_set_extravar('P_yx','P_yx')
    i_Pyy = var_set_extravar('P_yy','P_yy')

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

      print*, unit_density
      print*, unit_pressure

      w(ixI^S,rho_) = rho/unit_density
      w(ixI^S,mom(:)) = 0.d0
      w(ixI^S,Er_) = Bg_E*(1+Er*dexp(-sum(x(ixI^S,:)**2, dim=ndim+1)/0.01d0))/unit_pressure
      w(ixI^S,Fr(:)) = 0.d0

  end subroutine initial_conditions

  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    use mod_hd

    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: prad(ixI^S)

    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, prad)
    w(ixO^S,i_Pxx) = prad(ixO^S)/w(ixO^S,Er_)
    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 2, prad)
    w(ixO^S,i_Pxy) = prad(ixO^S)/w(ixO^S,Er_)
    call m1_get_prad(w, x, ixI^L, ixO^L, 2, 1, prad)
    w(ixO^S,i_Pyx) = prad(ixO^S)/w(ixO^S,Er_)
    call m1_get_prad(w, x, ixI^L, ixO^L, 2, 2, prad)
    w(ixO^S,i_Pyy) = prad(ixO^S)/w(ixO^S,Er_)

  end subroutine update_extravars
end module mod_usr
