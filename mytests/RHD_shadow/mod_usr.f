!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: rho0 = 5.d-2 !1.d3
  double precision :: rho1 = 1.d0
  double precision :: x0 = 0.1d0
  double precision :: y0 = 0.06d0
  double precision :: xc = 0.5d0
  double precision :: yc = 0.d0
  double precision :: T0 = 290.d0
  double precision :: T1 = 6.d0*290.d0
  double precision :: kap0 = 1.d2

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

    ! Special Opacity
    usr_special_opacity => special_kramer

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    unit_temperature = T0
    unit_numberdensity = rho0/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = 1.d0

    ! !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=(2.d0+3.d0*He_abundance)&
       *unit_numberdensity*const_kB*unit_temperature
    unit_velocity=dsqrt(unit_pressure/unit_density)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)


    !> ademensionalise parameters
    x0 = x0/unit_length
    y0 = y0/unit_length
    xc = xc/unit_length
    yc = yc/unit_length

    rho0 = rho0/unit_density
    rho1 = rho1/unit_density

    T0 = T0/unit_density
    T1 = T1/unit_density

    kap0 = kap0/unit_opacity

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

    double precision :: rho(ixImin1:ixImax1,ixImin2:ixImax2),&
        Delta(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: kappa(ixImin1:ixImax1,ixImin2:ixImax2),&
        temp(ixImin1:ixImax1,ixImin2:ixImax2)

    Delta(ixImin1:ixImax1,ixImin2:ixImax2) = 10.d0*( ((x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1)-xc)/x0)**2 + ((x(ixImin1:ixImax1,ixImin2:ixImax2,&
       2)-yc)/y0)**2 - 1.d0 )
    rho(ixImin1:ixImax1,ixImin2:ixImax2) = rho0 + (rho1 - rho0)/(1.d0 + &
       dexp(Delta(ixImin1:ixImax1,ixImin2:ixImax2)))

    temp(ixImin1:ixImax1,ixImin2:ixImax2) = T1 + (T1 - T0)/(xprobmin1 - &
       xprobmax1)*x(ixImin1:ixImax1,ixImin2:ixImax2,1)

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = rho(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = 0.d0
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = rho(ixImin1:ixImax1,&
       ixImin2:ixImax2)*temp(ixImin1:ixImax1,ixImin2:ixImax2)/(rhd_gamma-1.d0)
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = const_rad_a*(temp(ixImin1:ixImax1,&
       ixImin2:ixImax2)*unit_temperature)**4/unit_pressure

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

    double precision :: rho(ixBmin1:ixBmax1,ixBmin2:ixBmax2),&
        Delta(ixBmin1:ixBmax1,ixBmin2:ixBmax2)

    select case (iB)
    case(1)
      Delta(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = 10.d0*( ((x(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,1)-xc)/x0)**2 + ((x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         2)-xc)/x0)**2 - 1.d0 )
      rho(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = rho0 + (rho1 - rho0)/(1.d0 + &
         dexp(Delta(ixBmin1:ixBmax1,ixBmin2:ixBmax2)))

      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_) = rho(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(:)) = w(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,mom(:))
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,e_) = rho(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)*T1/(rhd_gamma-1.d0)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,r_e) = &
         const_rad_a*(T1*unit_temperature)**4/unit_pressure

    case(2)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,r_e) = &
         const_rad_a*(T0*unit_temperature)**4/unit_pressure


    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions


  subroutine mg_boundary_conditions(iB)

    use mod_global_parameters
    use mod_multigrid_coupling

    integer, intent(in)             :: iB

    select case (iB)
    case(1)
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = const_rad_a*(T1*unit_temperature)**&
         4.d0/unit_pressure

    case(2)
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = const_rad_a*(T0*unit_temperature)**&
         4.d0/unit_pressure

    case default
      print *, "Not a standard: ", trim(typeboundary(r_e, iB))
      error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions

  subroutine special_kramer(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,kappa)
    use mod_physics, only: phys_get_trad
    use mod_global_parameters
    use mod_opacity
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    double precision :: temp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    temp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (rhd_gamma - &
       1.d0)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_) - half * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:))**2,&
        dim=ndim+1)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_))

    temp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = temp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = kap0*(temp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/T0)**(-3.5d0) *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)/rho0

  end subroutine special_kramer


end module mod_usr
