!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision, parameter :: M_sun = 1.9891000d33
  double precision, parameter :: R_sun = 6.9599000d10
  double precision, parameter :: L_sun = 3.8268000d33
  double precision, parameter :: year = 365.25*24*60*60

  double precision :: StefBoltz

  double precision :: cak_Q, cak_a, cak_base, cak_x0, cak_x1
  integer :: it_start_cak
  double precision :: rho_bound, v_inf, Mdot
  double precision :: T_bound, R_star, M_star

  double precision :: kappa_e, L_bound, Gamma_e_bound, F_bound, gradE, E_out

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

    ! PseudoPlanar correction
    usr_source => PseudoPlanar

    ! Graviatational field
    usr_gravity => set_gravitation_field

    ! Special Opacity
    usr_special_opacity => OPAL_and_CAK

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    ! Timestep for PseudoPlanar
    ! usr_get_dt => pp_dt

    ! usr_refine_grid => refine_base

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_opacity, only: init_opal

    use mod_fld

    integer :: i

    !> Initialise Opal
    call init_opal(He_abundance)

    !> read usr par
    call params_read(par_files)

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = R_star
    unit_numberdensity = rho_bound/((1.d0+4.d0*He_abundance)*const_mp)
    unit_velocity = v_inf

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+&
       3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    R_star = R_star/unit_length
    M_star = M_star/unit_density/unit_length**3
    T_bound = T_bound/unit_temperature
    rho_bound = rho_bound/unit_density
    v_inf = v_inf/unit_velocity
    Mdot = Mdot/unit_density/unit_length**3*unit_time

    kappa_e = 0.34d0/unit_opacity
    F_bound = F_bound/unit_radflux
    L_bound = L_bound/(unit_radflux*unit_length**2)

    StefBoltz = const_rad_a*const_c/4.d0*(unit_temperature**&
       4.d0)/(unit_velocity*unit_pressure)

    !> Very bad initial guess for gradE using kappa_e
    gradE = -F_bound*3*kappa_e*rho_bound*unit_velocity/const_c

    ! print*, 'L_bound', L_bound*(unit_radflux*unit_length**2), log10(L_bound*(unit_radflux*unit_length**2)/L_sun)
    ! stop

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ cak_Q, cak_a, cak_base, cak_x0, cak_x1, rho_bound,&
        T_bound, R_star, M_star, v_inf, Mdot, Gamma_e_bound, it_start_cak

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wind_list, end=113)
113    close(unitpar)
    end do

    R_star = R_star*R_sun
    M_star = M_star*M_sun
    L_bound = Gamma_e_bound * 4.0 * dpi * const_G * M_star * const_c/0.34d0
    F_bound = L_bound/(4*dpi*R_star**2)
    Mdot = Mdot*M_sun/year

    ! print*, Gamma_e_bound, dpi, const_G, M_star, const_c

  end subroutine params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_constants

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2), fld_R(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    double precision :: vel(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: T_out, E_out, E_gauge
    double precision :: T_in, E_in, rr(ixImin1:ixImax1,ixImin2:ixImax2), bb

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)   = rho_bound
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = 0.d0

    where(x(ixImin1:ixImax1,ixImin2:ixImax2,1) .gt. 1.d0)
      vel(ixImin1:ixImax1,ixImin2:ixImax2) =  v_inf*(1 - 0.999d0*( &
         1.d0/x(ixImin1:ixImax1,ixImin2:ixImax2,1)))**0.5d0
      w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = &
         Mdot/(4.d0*dpi*vel(ixImin1:ixImax1,ixImin2:ixImax2)*x(ixImin1:ixImax1,&
         ixImin2:ixImax2,1)**2)
    endwhere

    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)) = vel(ixImin1:ixImax1,&
       ixImin2:ixImax2)*w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)

    !> Outer/Inner temperature
    T_out = 28445.836732569689/unit_temperature
    E_out = const_rad_a*(T_out*unit_temperature)**4/unit_pressure
    T_in = T_bound
    E_in = const_rad_a*(T_in*unit_temperature)**4/unit_pressure

    !>Very bad initial profile using constant gradE
    ! w(ixI^S,r_e) = E_out + gradE*(x(ixI^S,1)-xprobmax1)
    ! w(ixI^S,r_e) = E_in + gradE*(x(ixI^S,1)-xprobmin1)
    rr(ixImin1:ixImax1,ixImin2:ixImax2) = dsqrt(x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1)-xprobmin1)*(16*x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1)**2 + 8*x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1)+ 6) /(15*x(ixImin1:ixImax1,ixImin2:ixImax2,1)**(5.d0/2))
    bb = -kappa_e*L_bound*Mdot*unit_velocity*3/(16*dpi**2*const_c*v_inf)

    ! bb = bb*2

    ! rr(ixI^S) = 2*dsqrt(x(ixI^S,1)-xprobmin1)/dsqrt(x(ixI^S,1))
    ! bb = kappa_e*F_bound*Mdot*unit_velocity*3/(4*dpi*const_c*v_inf)
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = E_in + bb*rr(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    E_gauge = E_in + bb*dsqrt(xprobmax1-xprobmin1)*(16*xprobmax1**2 + &
       8*xprobmax1+ 6)/(15*xprobmax1**(5.d0/2))

    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
       r_e) - E_gauge + E_out

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
    use mod_physics, only: phys_get_trad
    use mod_global_parameters
    use mod_opacity
    use mod_fld

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixBmin1,ixBmin2,ixBmax1,ixBmax2, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: kappa(ixImin1:ixImax1,ixImin2:ixImax2),&
        Temp(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: Temp0, rho0, T_out, n
    double precision :: Local_gradE(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: Local_tauout(ixBmin1:ixBmax1,ixBmin2:ixBmax2)
    double precision :: Local_Tout(ixBmin1:ixBmax1,ixBmin2:ixBmax2)
    integer :: ix1,ix2

    select case (iB)

    case(1)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_) = rho_bound
      do ix1 = ixBmax1-1,ixBmin1,-1
        w(ix1,ixBmin2:ixBmax2,rho_) = dexp(2*dlog(w(ix1+1,ixBmin2:ixBmax2,&
           rho_)) - dlog(w(ix1+2,ixBmin2:ixBmax2,rho_)))
      enddo

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,ixBmin2:ixBmax2,mom(1)) = w(ix1+1,ixBmin2:ixBmax2,&
           mom(1)) *(x(ix1+1,ixBmin2:ixBmax2,1)/x(ix1,ixBmin2:ixBmax2,1))**2
        w(ix1,ixBmin2:ixBmax2,mom(2)) = w(ix1+1,ixBmin2:ixBmax2,mom(2))
      enddo

      where(w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1)) .lt. 0.d0)
        w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1)) = 0.d0
      endwhere

      where(w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1))/w(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,rho_) .gt. 0.5d0)
        w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1)) = 0.1d0*w(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2,rho_)
      endwhere

      ! w(ixB^S,r_e) = const_rad_a*(T_bound*unit_temperature)**4/unit_pressure

      call get_kappa_OPAL(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
         ixImax1,ixImax2,w,x,kappa)
      do ix1 = ixBmin1,ixBmax1
        kappa(ix1,ixBmin2:ixBmax2) = kappa(ixBmax1+1,ixBmin2:ixBmax2)
      enddo

      Local_gradE(ixImin1:ixImax1,ixImin2:ixImax2) = &
         -F_bound*3*kappa(ixImin1:ixImax1,ixImin2:ixImax2)*w(ixImin1:ixImax1,&
         ixImin2:ixImax2,rho_)*unit_velocity/const_c
      gradE = sum(Local_gradE(nghostcells,ixBmin2:ixBmax2))/(ixBmax2-ixBmin2)

      ! print*, gradE

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,ixBmin2:ixBmax2,r_e) = w(ix1+2,ixBmin2:ixBmax2,r_e) + (x(ix1,&
           ixBmin2:ixBmax2,1)-x(ix1+2,ixBmin2:ixBmax2,1))*Local_gradE(ix1+1,&
           ixBmin2:ixBmax2)
      enddo

      ! w(nghostcells,ixBmin2:ixBmax2,r_e) = dexp(half*(dlog(w(nghostcells-1,ixBmin2:ixBmax2,r_e))+dlog(w(nghostcells+1,ixBmin2:ixBmax2,r_e))))
      ! w(nghostcells,ixBmin2:ixBmax2,r_e) = half*(w(nghostcells-1,ixBmin2:ixBmax2,r_e)+w(nghostcells+1,ixBmin2:ixBmax2,r_e))

      ! print*, it, 'bottom------------------------------------'
      ! print*, w(1:5,5,r_e)
      ! print*, '********************', (w(3:6,5,r_e) - w(1:4,5,r_e))
      ! print*, '********************', (w(3:6,5,r_e) - w(1:4,5,r_e))/(w(2:5,5,rho_)*kappa(2:5,5))

      ! print*, F_bound, gradE*const_c/unit_velocity/(3*kappa(2,5)*w(2,5,rho_))
      print*, 4*dpi*unit_length**2*F_bound*unit_radflux/L_sun &
      , -4*dpi*unit_length**2*gradE*const_c/unit_velocity/(3*kappa(2,5)*w(2,5,&
         rho_))*unit_radflux/L_sun
      print*, -4*dpi*unit_length**2*(w(1:5,5,r_e)-w(3:7,5,r_e))/(x(1:5,5,1)-x(3:7,5,1))&
      *const_c/unit_velocity/(3*kappa(2:6,5)*w(2:6,5,rho_))*unit_radflux/L_sun
      print*,

    case(2)
      Local_tauout(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = &
         kappa_e*w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         rho_)*R_star**2/(3*x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,1))
      Local_Tout(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = &
         F_bound/StefBoltz*(3.d0/4.d0*Local_tauout(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2))**0.25d0

      !> one single NaN will kill the average and then we automatically take the floor temp
      where (Local_Tout(ixBmin1:ixBmax1,ixBmin2:ixBmax2) .ne. &
         Local_Tout(ixBmin1:ixBmax1,ixBmin2:ixBmax2))
        Local_Tout(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = 2.d4
      endwhere
      T_out = sum(Local_Tout(ixBmin2:ixBmax2,ixBmin1))/(ixBmax2-ixBmin2)

      T_out = max(1.d4/unit_temperature, T_out)
      E_out = const_rad_a*(T_out*unit_temperature)**4.d0/unit_pressure

      do ix1 = ixBmin1,ixBmax1
        w(ix1,:,r_e) = 2*w(ix1-1,:,r_e) - w(ix1-2,:,r_e)
      enddo

      ! w(ixB^S,r_e) = const_rad_a*(Local_Tout(ixB^S)*unit_temperature)**4.d0/unit_pressure

      ! print*, T_out*unit_temperature, E_out

      ! print*, it, 'top---------------------------------------'
      ! print*, Local_tauout(ixBmax1-5:ixBmax1,5)
      ! print*, Local_Tout(ixBmax1-5:ixBmax1,5)
      ! print*, w(ixBmax1-5:ixBmax1,5,r_e)

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
      mg%bc(iB, mg_iphi)%bc_value = gradE
      ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      ! mg%bc(iB, mg_iphi)%bc_value = const_rad_a*(T_bound*unit_temperature)**4/unit_pressure

    case (2)
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = E_out
      ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
      ! mg%bc(iB, mg_iphi)%bc_value = 0.d0

    case default
      print *, "Not a standard: ", trim(typeboundary(r_e, iB))
      error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions


  !> Calculate gravitational acceleration in each dimension
  subroutine set_gravitation_field(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,&
       ixImin2:ixImax2,ndim)

    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: mass

    radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,:) = zero
    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
       1) = -const_G*mass/radius(ixImin1:ixImax1,&
       ixImin2:ixImax2)**2*(unit_time**2/unit_length)

  end subroutine set_gravitation_field

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine PseudoPlanar(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, iwmin,iwmax
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision :: ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nw)

    call PseudoPlanarSource(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,wCT,x,ppsource)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_) + qdt*ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_) !> OK
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(1)) + qdt*ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(1)) !> OK
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2)) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(2)) + qdt*ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(2)) !> OK
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       r_e) + qdt*ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) !> TROUBLEMAKER

  end subroutine PseudoPlanar

  ! subroutine pp_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)
  !   use mod_global_parameters
  !   integer, intent(in)             :: ixI^L, ixO^L
  !   double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim)
  !   double precision, intent(in)    :: w(ixI^S,1:nw)
  !   double precision, intent(inout) :: dtnew
  !   double precision :: ppsource(ixO^S,1:nw)
  !
  !   call PseudoPlanarSource(ixI^L,ixO^L,w,x,ppsource)
  !
  !   dtnew = 1.d-1* minval(abs(w(ixO^S,r_e)/ppsource(ixO^S,r_e)))
  ! end subroutine pp_dt


  subroutine PseudoPlanarSource(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x,source)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out) :: source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1:nw)

    double precision :: rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndir)
    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),v(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)
    double precision :: radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
         pert(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    integer :: rdir, pdir

    source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nw) = zero

    rdir = 1
    pdir = 2

    v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rdir) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(rdir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,pdir) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(pdir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rdir) !+ half*block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rdir)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) = -two*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rdir)/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    call phys_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,pth)

    !> dm_r/dt = +(rho*v_p**2 + 2pth)/r -2 (rho*v_r**2 + pth)/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(rdir)) = - 2*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rdir)**two/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2) + w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       pdir)**two/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(pdir)) = - 3*v(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rdir)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       pdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)/radius(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      !> THIS BAD BOiii IS GIVING US SOME TROUBLE
      call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, rad_flux)
      source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = source(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,r_e) - two*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rdir)/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    endif

    if (rhd_radiation_advection) then
      source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = source(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,r_e) - two*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         r_e)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rdir)/radius(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
    endif

    ! Not sure about this one
    if (rhd_radiation_force) then
      source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = source(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,r_e) + two*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         r_e)/(3*radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    endif

  end subroutine PseudoPlanarSource


  subroutine OPAL_and_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
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

    double precision :: OPAL(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        CAK(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    !> Get OPAL opacities by reading from table
    call get_kappa_OPAL(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,x,OPAL)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    call get_kappa_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,w,x,CAK)

    !> Add OPAL and CAK for total opacity
    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = OPAL(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) + CAK(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    where(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .ne. kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))
      kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = kappa_e
    endwhere

  end subroutine OPAL_and_CAK


  subroutine get_kappa_OPAL(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
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

    integer :: ix1,ix2
    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: n, rho0, Temp0

    !> Get OPAL opacities by reading from table
    call phys_get_trad(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,Temp)
    do ix1=ixOmin1,ixOmax1
     do ix2=ixOmin2,ixOmax2
    
        rho0 = w(ix1,ix2,rho_)*unit_density
        Temp0 = Temp(ix1,ix2)*unit_temperature
        Temp0 = max(Temp0,1.d4)
        call set_opal_opacity(rho0,Temp0,n)
        kappa(ix1,ix2) = n/unit_opacity
    enddo
     enddo
    

    !> test without opal
    ! kappa(ixO^S) = kappa_e

  end subroutine get_kappa_OPAL

  subroutine get_kappa_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,kappa)
    use mod_global_parameters
    use mod_opacity
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    double precision :: vel(ixImin1:ixImax1,ixImin2:ixImax2),&
        gradv(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision :: xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixImin1:ixImax1,ixImin2:ixImax2) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(1))/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    call gradientO(vel,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,1,gradv)
    gradv(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs(gradv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))

    xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 1.d0-xprobmin1/x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)

    where (xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .le. cak_x0)
      alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cak_base
    elsewhere (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1) .le. cak_x1)
      alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cak_base + (cak_a - &
         cak_base)*(cak_x0 - xx(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))/(cak_x0 - cak_x1)
    elsewhere
      alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cak_a
    endwhere

    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
       kappa_e*cak_Q/(1-alpha(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)) *(gradv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*unit_velocity/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)*const_c*cak_Q*kappa_e))**alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if (it .le. it_start_cak) then
      kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = kappa(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*dexp(-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)*kappa_e)
    endif

    ! !> Limit cak for stability purposes
    ! where(kappa(ixO^S) .ge. 1.d3*kappa_e)
    !   kappa(ixO^S) = 1.d3*kappa_e
    ! endwhere

    !> test with no cak
    ! kappa(ixO^S) = 0.d0
  end subroutine get_kappa_CAK

  subroutine specialvar_output(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_physics
    use mod_fld

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: g_rad(ixImin1:ixImax1,&
       ixImin2:ixImax2), big_gamma(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                   :: g_grav(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                   :: Tgas(ixImin1:ixImax1,&
       ixImin2:ixImax2),Trad(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                   :: kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2), OPAL(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        CAK(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision                   :: vel(ixImin1:ixImax1,ixImin2:ixImax2),&
        gradv(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision                   :: rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:ndim), Lum(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision                   :: pp_rf(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    integer                            :: idim
    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: mass

    radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, kappa)
    call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, rad_flux)

    g_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1)/(const_c/unit_velocity)
    g_grav(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
       const_G*mass/radius(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2*(unit_time**2/unit_length)
    big_gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = g_rad(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/g_grav(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    call rhd_get_tgas(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, Tgas)
    call rhd_get_trad(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, Trad)

    call get_kappa_OPAL(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,x,OPAL)
    call get_kappa_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,w,x,CAK)

    vel(ixImin1:ixImax1,ixImin2:ixImax2) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(1))/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    call gradientO(vel,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,1,gradv)

    pp_rf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = two*rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)*dt

    Lum = 4*dpi*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)*(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)*unit_length)**2*unit_radflux/L_sun

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1) = kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/kappa_e
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2) = rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3) = Trad(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*unit_temperature
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+4) = big_gamma(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+5) = 4*dpi*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(1))*radius(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2 *unit_density*unit_velocity/M_sun*year
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+6) = OPAL(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/kappa_e
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+7) = CAK(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/kappa_e
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+8) = gradv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+9) = pp_rf(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+10) = Lum(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'kappa F1 Trad Gamma Mdot OPAL CAK gradv pp_rf L'
  end subroutine specialvarnames_output

end module mod_usr
