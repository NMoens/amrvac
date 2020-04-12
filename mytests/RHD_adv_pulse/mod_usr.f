!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: rho0 = 1.2d0
  double precision :: v0 = 5.d7
  double precision :: T0 = 1.d7
  double precision :: T1 = 2.d7
  double precision :: wdth = 24.d0

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! Boundary conditions
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    unit_velocity = v0 !r_arr(nghostcells) ! cm
    unit_numberdensity = rho0/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = wdth

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+&
       3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    print*, unit_time, 's'

    rho0 = rho0/unit_density
    v0 = v0/unit_velocity
    T0 = T0/unit_temperature
    T1 = T1/unit_temperature
    wdth = wdth/unit_length
  end subroutine initglobaldata_usr


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: temp(ixImin1:ixImax1,ixImin2:ixImax2),&
        pth(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2), lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    v0 = 0.d0

    temp(ixImin1:ixImax1,ixImin2:ixImax2) = T0 + &
       (T1-T0)*dexp(-x(ixImin1:ixImax1,ixImin2:ixImax2,1)**2.d0/(2*wdth**2))
    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = rho0*T0/temp(ixImin1:ixImax1,&
       ixImin2:ixImax2) + const_rad_a*fld_mu*const_mp/(3.d0*const_kB) * &
       unit_temperature**3/unit_density * (T0**4.d0/temp(ixImin1:ixImax1,&
       ixImin2:ixImax2) - temp(ixImin1:ixImax1,ixImin2:ixImax2)**3.d0)
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_)*v0
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(2)) = 0.d0
    pth(ixImin1:ixImax1,ixImin2:ixImax2) = temp(ixImin1:ixImax1,&
       ixImin2:ixImax2)*w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = pth(ixImin1:ixImax1,&
       ixImin2:ixImax2)/(rhd_gamma-1.d0) + half*w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_)*v0**2
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = const_rad_a*(temp(ixImin1:ixImax1,&
       ixImin2:ixImax2)*unit_temperature)**4.d0/unit_pressure

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, lambda, fld_R)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_test) = lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_diff_mg) = &
       (const_c/unit_velocity)*lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/(kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2,iB,w,x)
    use mod_global_parameters
    use mod_fld


    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixBmin1,ixBmin2,ixBmax1,ixBmax2, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: temp(ixBmin1:ixBmax1,ixBmin2:ixBmax2),&
        pth(ixBmin1:ixBmax1,ixBmin2:ixBmax2)

    select case (iB)
    case(1,2)
      temp(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = T0 + &
         (T1-T0)*dexp(-x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,1)**2.d0/(2*wdth**2))
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_) = rho0*T0/temp(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2) + const_rad_a*fld_mu*const_mp/(3.d0*const_kB) * &
         unit_temperature**3/unit_density * (T0**4.d0/temp(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2) - temp(ixBmin1:ixBmax1,ixBmin2:ixBmax2)**3.d0)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1)) = w(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,rho_)*v0
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(2)) = 0.d0
      pth(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = temp(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)*w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,e_) = pth(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)/(rhd_gamma-1.d0) + half*w(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,rho_)*v0**2
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,r_e) = &
         const_rad_a*(temp(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)*unit_temperature)**4.d0/unit_pressure

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions


  subroutine mg_boundary_conditions(iB)

    use mod_global_parameters
    use mod_multigrid_coupling

    integer, intent(in)             :: iB

    select case (iB)
    case (1)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
        mg%bc(iB, mg_iphi)%bc_value = 0.d0
    case (2)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
        mg%bc(iB, mg_iphi)%bc_value = 0.d0

      case default
        print *, "Not a standard: ", trim(typeboundary(r_e, iB))
        error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions

end module mod_usr
