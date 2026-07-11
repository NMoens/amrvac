!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd
  use mod_m1

  implicit none

  double precision :: rho1, v1, T1
  double precision :: p1, Er1

  integer :: i_Tgas, i_Trad

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system
    call set_coordinate_system("Cartesian")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions

    usr_process_grid => update_extravars

    ! Choose independent normalization units if using dimensionless variables.
    unit_velocity = v1
    unit_length = 7.d10
    unit_density = rho1

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

    i_Tgas = var_set_extravar('Tgas','Tgas')
    i_Trad = var_set_extravar('Trad','Trad')

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /shock_list/ rho1, v1, T1

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, shock_list, end=111)
       111 close(unitpar)
    enddo

  end subroutine usr_params_read

  subroutine initglobaldata_usr
    use mod_constants

    print*, 'unit_density', unit_density
    print*, 'unit_velocity', unit_velocity
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_pressure', unit_pressure
    print*, 'unit_time', unit_time
    print*, 'unit_length', unit_length

    p1 = const_kB*T1*rho1/(const_mp*hd_mu)
    p1 = p1/unit_pressure

    Er1 = const_rad_a*T1**4
    Er1 = Er1/unit_pressure

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixI^S,rho_) = 1.d0
    w(ixI^S,mom(1)) = -1.d0

    w(ixI^S,p_) = p1
    w(ixI^S,Er_) = Er1
    w(ixI^S,Fr(1)) = 0.d0



    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    use mod_hd

    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: pth(ixI^S)

    !Gas temperature
    call hd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,i_Tgas) = (pth(ixO^S)*unit_pressure*const_mp*hd_mu) &
                   /(w(ixO^S,rho_)*unit_density*const_kB)

    !Radiation temperature
    w(ixO^S,i_Trad) = (w(ixO^S,Er_)*unit_pressure/const_rad_a)**0.25d0

  end subroutine update_extravars

end module mod_usr
