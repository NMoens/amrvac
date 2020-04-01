!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  ! Custom variables can be defined here
  ! ...
  double precision, parameter :: M_sun = 1.99d33
  double precision, parameter :: R_sun = 6.96d10
  double precision, parameter :: year = 365.25*24*60*60

  double precision :: R_star
  double precision :: rho_bound
  double precision :: M_dot
  double precision :: v_inf
  double precision :: E_0
  double precision :: T_0
  double precision :: L_star

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Pseudo planar correction
    usr_source => PseudoPlanar

    ! Routine for setting special boundary conditions
    ! usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions

    !> Keep pressure fixed
    ! usr_internal_bc => FixPressure

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters

    R_star = 20.d0*R_sun
    rho_bound = 1.d-8
    M_dot = 1.d-5
    v_inf = M_dot/(rho_bound*4*dpi*R_star**2)
    T_0 = 2.d4
    E_0 = const_rad_a*T_0**4
    L_star = 1.d8

    ! Choose independent normalization units if using dimensionless variables.
    !> Define units
    unit_numberdensity = rho_bound/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_temperature = T_0
    unit_length = R_star

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
    unit_pressure=(2.d0+3.d0*He_abundance)&
       *unit_numberdensity*kB_cgs*unit_temperature
    unit_velocity=dsqrt(unit_pressure/unit_density)
    unit_time=unit_length/unit_velocity
    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    R_star = R_star/unit_length
    rho_bound = rho_bound/unit_density
    M_dot = M_dot/(unit_density*unit_length**3/unit_time)
    v_inf = v_inf/unit_velocity
    T_0 = t_0/unit_temperature
    E_0 = E_0/unit_pressure
    L_star = L_star/(unit_radflux*unit_time**2)
    fld_kappa0 = fld_kappa0/unit_opacity


    print*, 'Unitless: ########################'
    print*, 'R_star', R_star
    print*, 'rho_bound', rho_bound
    print*, 'M_dot', M_dot
    print*, 'v_inf', v_inf
    print*, 'T_0', T_0
    print*, 'E_0', E_0
    print*, 'L_star', L_star

    print*, 'units: ###########################'
    print*, 'unit_numberdensity', unit_numberdensity
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_length', unit_length
    print*, 'unit_density', unit_density
    print*, 'unit_pressure', unit_pressure
    print*, 'unit_velocity', unit_velocity
    print*, 'unit_time', unit_time
    print*, 'unit_radflux', unit_radflux
    print*, 'unit_opacity', unit_opacity

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision :: pressure(ixImin1:ixImax1,ixImin2:ixImax2)

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = &
       M_dot/(v_inf*4.d0*dpi*x(ixImin1:ixImax1,ixImin2:ixImax2,2)**2)
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)) = zero
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(2)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_)*v_inf

    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = E_0 + &
       fld_kappa0*L_star*M_dot/((4*dpi)**2*(const_c/unit_velocity)*v_inf) * &
       (one/x(ixImin1:ixImax1,ixImin2:ixImax2,2)**3 - one)

    print*, '#######################################################'
    print*, fld_kappa0*L_star*M_dot/((4*dpi)**2*(const_c/unit_velocity)*v_inf)

    pressure(ixImin1:ixImax1,ixImin2:ixImax2) = const_kB*w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_) /(const_mp*fld_mu)*(w(ixImin1:ixImax1,&
       ixImin2:ixImax2,r_e)/const_rad_a)**0.25 &
       *unit_density*unit_temperature/unit_pressure

    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = pressure(ixImin1:ixImax1,&
       ixImin2:ixImax2)/(rhd_gamma - one) + half*(w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(1))**2 + w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(2))**2)/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)

    call get_rad_extravars(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2)

  end subroutine initial_conditions

  ! subroutine boundary_conditions(qt,ixG^L,ixB^L,iB,w,x)
  !   use mod_global_parameters
  !   integer, intent(in)             :: ixG^L, ixB^L, iB
  !   double precision, intent(in)    :: qt, x(ixG^S,1:ndim)
  !   double precision, intent(inout) :: w(ixG^S,1:nw)
  !
  !   select case (iB)
  !   case(3)
  !     w(:,ixBmin2:ixBmax2,rho_) = M_dot/(v_inf*4.d0*dpi*x(:,ixBmin2:ixBmax2,2)**2)
  !     w(:,ixBmin2:ixBmax2,mom(1)) = zero
  !     w(:,ixBmin2:ixBmax2,mom(2)) = w(:,ixBmin2:ixBmax2,rho_)*v_inf
  !     w(:,ixBmin2:ixBmax2,e_) = p_0/(rhd_gamma - one) &
  !     + half*M_dot*v_inf/(4.d0*dpi*x(:,ixBmin2:ixBmax2,2)**2)
  !     w(:,ixBmin2:ixBmax2,r_e) = const_rad_a*(p_0/w(:,ixBmin2:ixBmax2,rho_) &
  !     *fld_mu*const_mp/const_kB*unit_pressure/unit_density/unit_temperature)**4.d0
  !
  !   case default
  !     call mpistop('boundary not known')
  !   end select
  ! end subroutine boundary_conditions

  subroutine mg_boundary_conditions(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iB,w,x)

    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: grad4

    grad4 = sum((w(ixOmin1:ixOmax1,ixOmax2,r_e) - w(ixOmin1:ixOmax1,ixOmax2+1,&
        r_e)) / (x(ixOmin1:ixOmax1,ixOmax2,2) - x(ixOmin1:ixOmax1,ixOmax2+1,&
       2)))/(ixOmax1-ixOmin1)

    select case (iB)
      case (1)
         mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
      case (2)
         mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
      case (3)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = 7.3583819042386223 !7.6447315544263788
      case (4)
        ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        ! mg%bc(iB, mg_iphi)%bc_value = 0.77780683570039721
        ! mg%bc(iB, mg_iphi)%bc_value = grad4 !min(grad4,zero)
        ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous

        ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann

        if (sum(w(ixOmin1:ixOmax1,ixOmax2,r_e)) .le. sum(w(ixOmin1:ixOmax1,&
           ixOmax2+1, r_e))) then
           mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
           mg%bc(iB, mg_iphi)%bc_value = zero
        else
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
        endif

      case default
        print *, "Not a standard: ", trim(typeboundary(iw_r_e, iB))
        error stop "You have to set a user-defined boundary method"
    end select
  end subroutine mg_boundary_conditions

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine PseudoPlanar(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, iwmin,iwmax
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2)
    integer :: rdir, pdir

    rdir = 2
    pdir = 1

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_) - qdt*two*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(rdir))/x(ixImin1:ixImax1,ixImin2:ixImax2,rdir)

    !> dm_r/dt = m_phi**2/r
    !> dm_phi/dt = - m_phi m_r/r
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(rdir)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(rdir)) + qdt*(wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(pdir)))**two/(x(ixImin1:ixImax1,ixImin2:ixImax2,&
       rdir)*wCT(ixImin1:ixImax1,ixImin2:ixImax2,rho_))
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(pdir)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(pdir)) - qdt*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(rdir))*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(pdir))/(x(ixImin1:ixImax1,ixImin2:ixImax2,rdir)*wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_))

    !> de/dt = -2 (e+p)v_r/r
    call phys_get_pthermal(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,pth)
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
       e_) - qdt*two*(wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       e_)+pth(ixImin1:ixImax1,ixImin2:ixImax2))*wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(rdir))/(wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       rho_)*x(ixImin1:ixImax1,ixImin2:ixImax2,rdir))

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      call get_rad_extravars(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2)
      w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = w(ixImin1:ixImax1,&
         ixImin2:ixImax2,r_e) - qdt*two*(wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         e_)*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         mom(rdir))/wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         rho_) + w(ixImin1:ixImax1,ixImin2:ixImax2,&
         i_flux(rdir)))/x(ixImin1:ixImax1,ixImin2:ixImax2,rdir)
    else
      w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = w(ixImin1:ixImax1,&
         ixImin2:ixImax2,r_e) - qdt*two*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         e_)*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         mom(rdir))/(wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         rho_)*x(ixImin1:ixImax1,ixImin2:ixImax2,rdir))
    endif

  end subroutine PseudoPlanar

end module mod_usr
