!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: Er0 = 1.d0
  double precision :: rho0 = 0.025

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")


    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    usr_internal_bc => set_Erad

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: temp(ixI^S), pth(ixI^S)
    double precision :: kappa(ixO^S), fld_R(ixO^S), lambda(ixO^S)

    temp(ixI^S) = (Er0/const_rad_a)**0.25

    w(ixI^S,rho_) = rho0
    w(ixI^S,mom(:)) = 0.d0
    pth(ixI^S) = temp(ixI^S)*const_kB/(const_mp*fld_mu)*w(ixI^S,rho_)
    w(ixI^S,e_) = pth(ixI^S)/(rhd_gamma-1.d0)
    w(ixI^S,r_e) = Er0

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))

  end subroutine initial_conditions


  subroutine set_Erad(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    if (it>0) then
      where (x(ixI^S,1) .lt. 0.1d0)
        w(ixI^S,r_e) = 1.d0
      end where
    endif

  end subroutine set_Erad

end module mod_usr
