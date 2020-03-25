!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: rho1 = 1.d-2
  double precision :: rho2 = 0.0685847d0
  double precision :: v1 = 1.d9
  double precision :: v2 = 1.458d8
  double precision :: T1 = 1.d4
  double precision :: T2 = 4.29d7

  double precision :: p1,p2,eg1,eg2,Er1,Er2

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

    integer :: i

    p1 = const_kB*T1/(const_mp*fld_mu)*rho1
    p2 = const_kB*T2/(const_mp*fld_mu)*rho2

    eg1 = p1/(rhd_gamma) + half*v1**2*rho1
    eg2 = p2/(rhd_gamma) + half*v2**2*rho2

    Er1 = const_rad_a*T1**4
    Er2 = const_rad_a*T2**4

    unit_velocity = v1 !r_arr(nghostcells) ! cm
    unit_numberdensity = rho1/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = 1.d5

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+&
       3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    rho1 = rho1/unit_density
    rho2 = rho2/unit_density

    v1 = v1/unit_velocity
    v2 = v2/unit_velocity

    T1 = T1/unit_temperature
    T2 = T2/unit_temperature

    p1 = p1/unit_pressure
    p2 = p2/unit_pressure

    eg1 = eg1/unit_pressure
    eg2 = eg2/unit_pressure

    Er1 = Er1/unit_pressure
    Er2 = Er2/unit_pressure

  end subroutine initglobaldata_usr


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2), lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = rho1
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)) = rho1*v1
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(2)) = 0
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = eg1
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = Er1

    where (x(ixImin1:ixImax1,ixImin2:ixImax2,1) .gt. 0.d0)
      w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = rho2
      w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)) = rho2*v2
      w(ixImin1:ixImax1,ixImin2:ixImax2,mom(2)) = 0
      w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = eg2
      w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = Er2
    end where

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, lambda, fld_R)

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
    select case (iB)

    case(1)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_) = rho1
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1)) = rho1*v1
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(2)) = 0
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,e_) = eg1
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,r_e) = Er1

    case(2)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_) = rho2
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1)) = rho2*v2
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(2)) = 0
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,e_) = eg2
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,r_e) = Er2

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions


  subroutine mg_boundary_conditions(iB)

    use mod_global_parameters
    use mod_multigrid_coupling

    integer, intent(in)             :: iB

    integer :: ixOmax2

    select case (iB)
    case (1)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = Er1
    case (2)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = Er2

      case default
        print *, "Not a standard: ", trim(typeboundary(r_e, iB))
        error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions

end module mod_usr
