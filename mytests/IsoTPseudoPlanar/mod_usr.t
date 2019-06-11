!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd

  implicit none

  ! Custom variables can be defined here
  ! ...

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t

    ! Routine for setting special boundary conditions
    usr_special_bc => boundary_conditions

    ! Pseudo planar correction
    usr_source => PseudoPlanar

    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    ! Active the physics module
    call hd_activate()

  end subroutine usr_init

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixO^S, rho_) = one
    w(ixO^S, mom(1)) = one
    w(ixO^S, mom(2)) = half

  end subroutine initial_conditions

  ! Extra routines can be placed here
  ! ...

  subroutine boundary_conditions(qt,ixG^L,ixB^L,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixG^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    select case (iB)
    case(1)
      w(ixBmin1:ixBmax1,:,rho_) = one
      w(ixBmin1:ixBmax1,:,mom(1)) = one
      w(ixBmin1:ixBmax1,:,mom(2)) = half

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

  subroutine PseudoPlanar(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S)
    integer :: rdir, pdir

    rdir = 1
    pdir = 2

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    w(ixI^S,rho_) = w(ixI^S,rho_) - qdt*two*wCT(ixI^S,mom(rdir))/x(ixI^S,rdir)

    !> dm_r/dt = m_phi**2/r
    !> dm_phi/dt = - m_phi m_r/r
    w(ixI^S,mom(rdir)) = w(ixI^S,mom(rdir)) + qdt*(wCT(ixI^S,mom(pdir)))**two/(x(ixI^S,rdir)*wCT(ixI^S,rho_))
    w(ixI^S,mom(pdir)) = w(ixI^S,mom(pdir)) - qdt*wCT(ixI^S,mom(rdir))*wCT(ixI^S,mom(pdir))/(x(ixI^S,rdir)*wCT(ixI^S,rho_))
  end subroutine PseudoPlanar


end module mod_usr
