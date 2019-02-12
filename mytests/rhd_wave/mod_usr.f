!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

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

    ! Routine for setting special boundary conditions
    usr_special_bc => boundary_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! ...

    ! Choose independent normalization units if using dimensionless variables.
    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = 1.d0 ! cm
    unit_temperature   = 1.d0 ! K
    unit_numberdensity = 1.d0!/((1.d0+4.d0*He_abundance)*mp_cgs) ! cm-3,cm-3

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas
    use mod_fld

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: local_rad_e

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = two
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = zero
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = one
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = one

    call phys_get_tgas(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,Temp)

    !> Have to do this in 2 steps to define boundary e_rad
    local_rad_e = const_rad_a/unit_pressure*unit_temperature**4.d0 &
       *Temp(nghostcells+2,nghostcells+2)**4.d0

    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = local_rad_e

    call get_rad_extravars(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2)

  end subroutine initial_conditions

  ! Extra routines can be placed here
  ! ...
  subroutine boundary_conditions(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2,iB,w,x)
    use mod_global_parameters
    use mod_constants
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
        ixBmin1,ixBmin2,ixBmax1,ixBmax2, iB
    double precision, intent(in)    :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

    select case (iB)
    case(1)
        w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
           rho_) = two + dsin(global_time*two*3.1415)
    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

end module mod_usr
