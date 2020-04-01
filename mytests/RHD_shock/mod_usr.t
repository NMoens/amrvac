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

    print*, 'M_1: ', v1/dsqrt(const_kB*T1/(const_mp*fld_mu))
    print*, 'M_2: ', v2/dsqrt(const_kB*T2/(const_mp*fld_mu))

    unit_velocity = v1 !r_arr(nghostcells) ! cm
    unit_numberdensity = rho1/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = 1.d5

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB)
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
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: kappa(ixO^S), fld_R(ixO^S), lambda(ixO^S)

    w(ixI^S,rho_) = rho1
    w(ixI^S,mom(1)) = rho1*v1
    w(ixI^S,mom(2)) = 0
    w(ixI^S,e_) = eg1
    w(ixI^S,r_e) = Er1

    where (x(ixI^S,1) .gt. 0.d0)
      w(ixI^S,rho_) = rho2
      w(ixI^S,mom(1)) = rho2*v2
      w(ixI^S,mom(2)) = 0
      w(ixI^S,e_) = eg2
      w(ixI^S,r_e) = Er2
    end where

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    use mod_fld


    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    select case (iB)

    case(1)
      w(ixB^S,rho_) = rho1
      w(ixB^S,mom(1)) = rho1*v1
      w(ixB^S,mom(2)) = 0
      w(ixB^S,e_) = eg1
      w(ixB^S,r_e) = Er1

    case(2)
      w(ixB^S,rho_) = rho2
      w(ixB^S,mom(1)) = rho2*v2
      w(ixB^S,mom(2)) = 0
      w(ixB^S,e_) = eg2
      w(ixB^S,r_e) = Er2

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
