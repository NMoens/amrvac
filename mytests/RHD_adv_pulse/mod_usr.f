!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: rho1 = 1.2d0
  double precision :: v1 = 1.d6
  double precision :: T1 = 1.d7
  double precision :: T2 = 2.d7
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
    ! usr_special_bc => boundary_conditions
    ! usr_special_mg_bc => mg_boundary_conditions

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    unit_velocity = v1 !r_arr(nghostcells) ! cm
    unit_numberdensity = rho1/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = wdth

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+&
       3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

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

    temp(ixImin1:ixImax1,ixImin2:ixImax2) = T1 + &
       (T2-T1)*dexp(-(x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1)*unit_length)**2.d0/(2*wdth**2))
    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = rho1*T1/temp(ixImin1:ixImax1,&
       ixImin2:ixImax2) + const_rad_a*fld_mu*const_mp/(3.d0*const_kB) * &
       (T1**4.d0/temp(ixImin1:ixImax1,ixImin2:ixImax2) - temp(ixImin1:ixImax1,&
       ixImin2:ixImax2)**3.d0)
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_)*v1
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(2)) = 0.d0
    pth(ixImin1:ixImax1,ixImin2:ixImax2) = const_kB*temp(ixImin1:ixImax1,&
       ixImin2:ixImax2)*w(ixImin1:ixImax1,ixImin2:ixImax2,&
       rho_) /(const_mp*fld_mu)
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = pth(ixImin1:ixImax1,&
       ixImin2:ixImax2)/(rhd_gamma-1.d0) + half*w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_)*v1**2
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = const_rad_a*temp(ixImin1:ixImax1,&
       ixImin2:ixImax2)**4.d0

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_)/unit_density
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(:))/(unit_density*unit_velocity)
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
       e_)/unit_pressure
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
       r_e)/unit_pressure

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, lambda, fld_R)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_diff_mg) = &
       (const_c/unit_velocity)*lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/(kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

  end subroutine initial_conditions

  ! subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
  !   use mod_global_parameters
  !   use mod_fld
  !
  !
  !   integer, intent(in)             :: ixI^L, ixB^L, iB
  !   double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
  !   double precision, intent(inout) :: w(ixI^S,1:nw)
  !   select case (iB)
  !
  !   case(1)
  !     w(ixB^S,rho_) = rho1
  !     w(ixB^S,mom(1)) = rho1*v1
  !     w(ixB^S,mom(2)) = 0
  !     w(ixB^S,e_) = eg1
  !     w(ixB^S,r_e) = Er1
  !
  !   case(2)
  !     w(ixB^S,rho_) = rho2
  !     w(ixB^S,mom(1)) = rho2*v2
  !     w(ixB^S,mom(2)) = 0
  !     w(ixB^S,e_) = eg2
  !     w(ixB^S,r_e) = Er2
  !
  !   case default
  !     call mpistop('boundary not known')
  !   end select
  ! end subroutine boundary_conditions
  !
  !
  ! subroutine mg_boundary_conditions(iB)
  !
  !   use mod_global_parameters
  !   use mod_multigrid_coupling
  !
  !   integer, intent(in)             :: iB
  !
  !   integer :: ixOmax2
  !
  !   select case (iB)
  !   case (1)
  !       mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
  !       mg%bc(iB, mg_iphi)%bc_value = Er1
  !   case (2)
  !       mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
  !       mg%bc(iB, mg_iphi)%bc_value = Er2
  !
  !     case default
  !       print *, "Not a standard: ", trim(typeboundary(r_e, iB))
  !       error stop "Set special bound for this Boundary "
  !   end select
  ! end subroutine mg_boundary_conditions

end module mod_usr
