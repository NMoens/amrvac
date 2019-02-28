!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    use mod_constants

    double precision :: rho_0 = one !1.d-7
    double precision :: t_0 = one !1.d-12
    double precision :: e_0 = one !1.d12

    call set_coordinate_system("Cartesian_2D")

    !Fix dimensionless stuff here
    unit_length        = dsqrt(e_0/rho_0)*t_0 !cm
    unit_numberdensity = rho_0/((1.d0+4.d0*He_abundance)*mp_cgs) !rho_0/(fld_mu*mp_cgs)                                      ! cm-3,cm-3
    unit_temperature   = e_0/(unit_numberdensity*(2.d0+&
       3.d0*He_abundance)*kB_cgs) !e_0/(unit_numberdensity*rhd_gamma*kB_cgs)                   ! K

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Keep the radiative energy constant with internal bound
    usr_internal_bc => constant_r_e

    ! Active the physics module
    call rhd_activate()

    print*, 'unit_time', unit_time
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_length', unit_length
    print*, 'unit_density', unit_density
    print*, 'unit_numberdensity', unit_numberdensity
    print*, 'unit_velocity', unit_velocity
    print*, 'unit_pressure', unit_pressure
    print*, '================================================================'

  end subroutine usr_init

!==========================================================================================

subroutine initglobaldata_usr
  use mod_global_parameters


end subroutine initglobaldata_usr

!==========================================================================================

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
     ixmax1,ixmax2, w, x)
    use mod_global_parameters
    use mod_constants
    use mod_rhd_phys, only: rhd_get_pthermal

    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
       ixmin2,ixmax1,ixmax2
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
        ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, nw)
    double precision :: e_eq

    ! Set initial values for w
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = 1.d-7
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(:)) = zero
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,r_e) = 1.d12

    e_eq = w(3,3,rho_)/(rhd_gamma-one) * (fld_speedofligt_0*w(3,3,&
       r_e)/(4.0d0*fld_sigma_0))**(1.d0/4.d0)

    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, e_) = 1.d-3*e_eq

    call get_rad_extravars(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
       ixmin2,ixmax1,ixmax2)

    print*, w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,e_)
    print*, w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,r_e)

    print*, "E", w(3,3,r_e)
    print*, "c", fld_speedofligt_0
    print*, "rho", w(3,3,rho_)
    print*, "kappa", fld_kappa0
    print*, "sigma", fld_sigma_0
    print*, "gamma", rhd_gamma
    print*, "e0", e_eq

  end subroutine initial_conditions

!==========================================================================================

  ! Extra routines can be placed here
  ! ...

!==========================================================================================

!> internal boundary, user defined
  !
  !> This subroutine can be used to artificially overwrite ALL conservative
  !> variables in a user-selected region of the mesh, and thereby act as
  !> an internal boundary region. It is called just before external (ghost cell)
  !> boundary regions will be set by the BC selection. Here, you could e.g.
  !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
  !> which can be used to identify the internal boundary region location.
  !> Its effect should always be local as it acts on the mesh.

  subroutine constant_r_e(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = 1.d-7
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = zero
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = 1d12

    print*, w(3,3,e_)

    ! if (it .eq. 0) open(1,file='energy')
    ! write(1,*) global_time*unit_time, w(3,3,r_e)*unit_pressure, w(3,3,e_)*unit_pressure
    ! if (global_time .ge. time_max - dt) close(1)

    print*, global_time*unit_time, w(3,3,r_e)*unit_pressure, w(3,3,e_)*unit_pressure

  end subroutine constant_r_e

end module mod_usr

!==========================================================================================
