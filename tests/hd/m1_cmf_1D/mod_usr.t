!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd
  use mod_m1

  implicit none

  double precision :: beta_max, F_red
  integer :: i_Eobs, i_Fobs, i_Pobs, i_Pcmf

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

    usr_internal_bc => Fix_hd
    usr_process_grid => update_extravars

    ! Choose independent normalization units if using dimensionless variables.
    unit_velocity    = const_c
    unit_time        = 1.d0/unit_velocity
    unit_density     = 1.d0

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

    i_Pcmf = var_set_extravar('Pr','Pr')

    i_Eobs = var_set_extravar('E_obs','E_obs')
    i_Fobs = var_set_extravar('F_obs','F_obs')
    i_Pobs = var_set_extravar('P_obs','P_obs')


  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /stream_list/ F_red, beta_max

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

      w(ixI^S,rho_) = 1.d0/(x(ixI^S,1)*beta_max*m1_c)
      w(ixI^S,mom(1)) = 1.d0
      w(ixI^S,Er_) = 1.d0
      w(ixI^S,Fr(1)) = F_red*m1_c*w(ixI^S,Er_)

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
    w(ixB^S,rho_) = 1.d0/(x(ixB^S,1)*beta_max*m1_c)
    w(ixB^S,mom(1)) = 1.d0
    w(ixB^S,Er_) = 1.d0
    w(ixB^S,Fr(1)) = F_red*m1_c*w(ixB^S,Er_)

  case default
    call mpistop('boundary not known')
  end select
end subroutine boundary_conditions

subroutine Fix_hd(level,qt,ixI^L,ixO^L,w,x)
  use mod_global_parameters
  use mod_physics, only: phys_get_pthermal

  integer, intent(in)             :: ixI^L,ixO^L,level
  double precision, intent(in)    :: qt
  double precision, intent(inout) :: w(ixI^S,1:nw)
  double precision, intent(in)    :: x(ixI^S,1:ndim)

  !> hd subsystem
  w(ixI^S,rho_) = 1.d0/(x(ixI^S,1)*beta_max*m1_c)
  w(ixI^S,mom(1)) = 1.d0

end subroutine Fix_hd

subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
  use mod_global_parameters
  use mod_hd

  integer, intent(in)             :: igrid,level,ixI^L,ixO^L
  double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
  double precision, intent(inout) :: w(ixI^S,1:nw)

  double precision :: v(ixI^S), prad(ixI^S)

  v(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
  call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, prad)

  w(ixI^S,i_Pcmf) = prad(ixI^S)

  w(ixI^S,i_Eobs) = w(ixI^S,Er_) + 2*v(ixI^S)/m1_c**2*w(ixI^S,Fr(1))
  w(ixI^S,i_Fobs) = v(ixI^S)*w(ixI^S,Er_) + w(ixI^S,Fr(1)) + v(ixI^S)*prad(ixI^S)
  w(ixI^S,i_Pobs) = prad(ixI^S) + 2*v(ixI^S)/m1_c**2*w(ixI^S,Fr(1))

end subroutine update_extravars

end module mod_usr
