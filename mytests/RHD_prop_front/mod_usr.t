!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: Er0 = 1.d-22
  double precision :: Er1 = 1.d0
  double precision :: rho0 = 0.025d0
  double precision :: l1 = 1.d-1
  double precision :: l2 = 1.d-5

  double precision :: p0, T0, p1, T1

  integer :: i_sol

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Boundary conditions
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions

    ! extra output
    usr_modify_output => output_routine

    ! Active the physics module
    call rhd_activate()

    i_sol = var_set_extravar("sol", "sol")

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    integer :: i

    T0 = (Er0/const_rad_a)**0.25
    p0 = const_kB*T0*rho0/(const_mp*fld_mu)

    T1 = (Er1/const_rad_a)**0.25
    p1 = const_kB*T1*rho0/(const_mp*fld_mu)

    unit_velocity = dsqrt(p1/rho0)
    unit_numberdensity = rho0/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = 1.d0

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    print*, unit_time, unit_velocity
    print*, '2.5d-11s in dimless is ', 2.5d-11/unit_time
    print*, 'speed of light dimless is', const_c/unit_velocity
    print*, 'dimless time for crossing is', xprobmax1/const_c*unit_velocity

    rho0 = rho0/unit_density
    p0 = p0/unit_pressure
    p1 = p1/unit_pressure
    T0 = T0/unit_temperature
    T1= T1/unit_temperature
    Er0 = Er0/unit_pressure
    Er1 = Er1/unit_pressure
    l1 = l1/unit_length
    l2 = l2/unit_length

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: kappa(ixO^S), fld_R(ixO^S), lambda(ixO^S), step(ixI^S)

    w(ixI^S,rho_) = rho0
    w(ixI^S,mom(:)) = 0.d0
    w(ixI^S,e_) = p0/(rhd_gamma-1.d0)

    step(ixI^S) = (  1.d0-erf(  (x(ixI^S,1)-l1)/l2    )  )/2.d0
    w(ixI^S,r_e) = Er0 + step(ixI^S)*Er1

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    w(ixO^S,i_test) = fld_R(ixO^S)
    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))
    w(ixI^S,i_sol) = Er0 + step(ixI^S)*Er1

  end subroutine initial_conditions


  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    use mod_fld


    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: i

    select case (iB)
    case(1)
      do i = ixBmax1,ixBmin1, -1
        w(i,:,rho_) =   w(i+1,:,rho_)
        w(i,:,mom(1)) = w(i+1,:,mom(1))
        w(i,:,mom(2)) = w(i+1,:,mom(2))
        w(i,:,e_) = w(i+1,:,e_)
        w(i,:,r_e) = Er1
      enddo

    case(2)
      do i = ixBmin1,ixBmax1
        w(i,:,rho_) = w(i-1,:,rho_)
        w(i,:,mom(1)) = w(i-1,:,mom(1))
        w(i,:,mom(2)) = w(i-1,:,mom(2))
        w(i,:,e_) = w(i-1,:,e_)
        w(i,:,r_e) = w(i-1,:,r_e)
      enddo

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
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = Er1
    case (2)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
    case default
      print *, "Not a standard: ", trim(typeboundary(r_e, iB))
      error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions

  subroutine output_routine(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: step(ixI^S)

    step(ixI^S) = (  1.d0-erf((x(ixI^S,1)-l1-global_time*const_c/unit_velocity)/l2    )  )/2.d0
    w(ixI^S,i_sol) = Er0 + step(ixI^S)*Er1

  end subroutine output_routine

end module mod_usr
