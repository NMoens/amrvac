!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: e_eq

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

    fld_kappa0 = fld_kappa0/unit_opacity

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
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    use mod_constants
    use mod_rhd_phys, only: rhd_get_pthermal

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)

    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2), fld_R(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)


    ! Set initial values for w
    w(ixImin1:ixImax1,ixImin2:ixImax2, rho_) = 1.d-7
    w(ixImin1:ixImax1,ixImin2:ixImax2, mom(:)) = zero
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = 1.d12

    e_eq = (w(3,3,r_e)*unit_pressure/(const_rad_a))**(1.d0/4.d0) &
       *one/(rhd_gamma-one)*const_kB*w(3,3,&
       rho_)*unit_density /(fld_mu*const_mp)/unit_pressure

    w(ixImin1:ixImax1,ixImin2:ixImax2, e_) = 1.d2*e_eq

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, lambda, fld_R)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_diff_mg) = &
       (const_c/unit_velocity)*lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/(kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

    print*, w(ixImin1:ixImax1,ixImin2:ixImax2,e_)
    print*, w(ixImin1:ixImax1,ixImin2:ixImax2,r_e)

    print*, "E", w(3,3,r_e)
    print*, "c", (const_c/unit_velocity)
    print*, "rho", w(3,3,rho_)
    print*, "kappa", fld_kappa0
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

    ! w(ixI^S,rho_) = 1.d-7
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = zero
    !w(ixI^S,r_e) = 1d12

    ! print*, w(3,3,e_)

    if (it .eq. 0) open(1,file='Halley1_1.d2')
    write(1,*) global_time*unit_time, e_eq*unit_pressure, w(3,3,&
       r_e)*unit_pressure, w(3,3,e_)*unit_pressure
    if (global_time .ge. time_max - dt) close(1)

    ! print*, global_time*unit_time, w(3,3,r_e)*unit_pressure, w(3,3,e_)*unit_pressure

  end subroutine constant_r_e

end module mod_usr

!==========================================================================================
