!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd
  use mod_m1

  implicit none

  double precision :: Er_l, Fx_l, Fy_l, Er_r, Fx_r, Fy_r
  integer :: i_Pxx, i_Pxy, i_Pyy, i_chi, i_trap

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_1.5D")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    ! usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions

    !> Additional variables
    usr_process_grid => update_extravars

    ! Choose independent normalization units if using dimensionless variables.
    unit_length      = 1.d0 ! cm
    unit_velocity    = const_c ! K
    unit_density     = 1.d0 ! g cm^-3

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

    i_Pxx = var_set_extravar('Pxx','Pxx')
    i_Pxy = var_set_extravar('Pxy','Pxy')
    i_Pyy = var_set_extravar('Pyy','Pyy')
    i_chi = var_set_extravar('chi','chi')
    i_trap = var_set_extravar('theta','theta')

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /shock_list/ Er_l, Fx_l, Fy_l, Er_r, Fx_r, Fy_r

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, shock_list, end=111)
       111 close(unitpar)
    enddo

  end subroutine usr_params_read


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixI^S,rho_) = 1.d0
    w(ixI^S,mom(1)) = 0.d0

    where (x(ixI^S,1) .lt. 0.d0)
      w(ixI^S,Er_) = Er_l
      w(ixI^S,Fr(1)) = Fx_l
      w(ixI^S,Fr(2)) = Fy_l
    elsewhere
      w(ixI^S,Er_) = Er_r
      w(ixI^S,Fr(1)) = Fx_r
      w(ixI^S,Fr(2)) = Fy_r
    end where

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    use mod_hd

    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: prad(ixI^S), chi(ixI^S), theta(ixI^S)

    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, prad)
    w(ixO^S,i_Pxx) = prad(ixO^S)

    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 2, prad)
    w(ixO^S,i_Pxy) = prad(ixO^S)

    call m1_get_prad(w, x, ixI^L, ixO^L, 2, 2, prad)
    w(ixO^S,i_Pyy) = prad(ixO^S)

    call m1_get_chi(w, x, ixI^L, ixO^L, chi)
    w(ixO^S,i_chi) = chi(ixO^S)

    call m1_get_trapped_fraction(w, x, ixI^L, ixO^L, 1, theta)
    w(ixO^S,i_trap) = theta(ixO^S)

  end subroutine update_extravars

end module mod_usr
