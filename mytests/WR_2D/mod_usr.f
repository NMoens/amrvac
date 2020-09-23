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
  logical :: fixed_lum, Cak_in_D

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
    ! usr_get_dt => get_dt_cak

    ! Refine mesh near base
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

    kappa_e = kappa_e/unit_opacity
    F_bound = F_bound/unit_radflux
    L_bound = L_bound/(unit_radflux*unit_length**2)

    StefBoltz = const_rad_a*const_c/4.d0*(unit_temperature**&
       4.d0)/(unit_velocity*unit_pressure)

    !> Very bad initial guess for gradE using kappa_e
    gradE = -F_bound*3*kappa_e*rho_bound*unit_velocity/const_c

    ! print*, 'L_bound (cgs)', L_bound*(unit_radflux*unit_length**2)
    ! print*, 'log10(L_bound)', log10(L_bound*(unit_radflux*unit_length**2)/L_sun)
    ! print*, 'L_bound', L_bound*(unit_radflux*unit_length**2)/L_sun
    ! ! stop
    ! print*, 'unit_density', unit_density
    ! print*, 'unit_time', unit_time
    ! print*, 'unit_pressure', unit_pressure

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ cak_Q, cak_a, cak_base, cak_x0, cak_x1, rho_bound,&
        kappa_e, T_bound, R_star, M_star, v_inf, Mdot, Gamma_e_bound,&
        it_start_cak, fixed_lum, Cak_in_D

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wind_list, end=113)
113    close(unitpar)
    end do

    R_star = R_star*R_sun
    M_star = M_star*M_sun
    L_bound = Gamma_e_bound * 4.0 * dpi * const_G * M_star * const_c/kappa_e
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
    double precision :: T_in, E_in, rr(ixImin1:ixImax1,ixImin2:ixImax2), bb,&
        temp(ixImin1:ixImax1,ixImin2:ixImax2)

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

    if (rhd_energy) then
      temp(ixImin1:ixImax1,ixImin2:ixImax2) = (w(ixImin1:ixImax1,&
         ixImin2:ixImax2,r_e)*unit_pressure/const_rad_a)**&
         0.25d0/unit_temperature
      w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = w(ixImin1:ixImax1,&
         ixImin2:ixImax2,rho_)*temp(ixImin1:ixImax1,&
         ixImin2:ixImax2)/(rhd_gamma-1.d0) + half*w(ixImin1:ixImax1,&
         ixImin2:ixImax2,mom(1))**2/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    endif

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, lambda, fld_R)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_diff_mg) = &
       (const_c/unit_velocity)*lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/(kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_diff_mg) = &
       (const_c/unit_velocity)/(3.d0*kappa(ixOmin1:ixOmax1,&
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
    double precision :: Local_gradE(ixImin1:ixImax1,ixImin2:ixImax2),&
        F_adv(ixBmin1:ixBmax1,ixBmin2:ixBmax2)
    double precision :: Local_tauout(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: Local_Tout(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision :: kappa_out

    integer :: ix1,ix2

    select case (iB)

    case(1)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_) = rho_bound
      do ix1 = ixBmax1-1,ixBmin1,-1
        w(ix1,ixBmin2:ixBmax2,rho_) = dexp(2*dlog(w(ix1+1,ixBmin2:ixBmax2,&
           rho_)) - dlog(w(ix1+2,ixBmin2:ixBmax2,rho_)))
      enddo

      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(2)) = 0.d0

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,ixBmin2:ixBmax2,mom(1)) = w(ix1+1,ixBmin2:ixBmax2,mom(1))
      enddo

      where(w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1)) .lt. 0.d0)
       w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1)) = 0.d0
      endwhere

      F_adv(ixBmax1,ixBmin2:ixBmax2) = 4.d0/3.d0*(w(ixBmax1,ixBmin2:ixBmax2,&
         mom(1))/w(ixBmax1,ixBmin2:ixBmax2,rho_))*w(ixBmax1,ixBmin2:ixBmax2,&
         r_e) * 4*dpi*xprobmin1**2

      where (F_adv(ixBmin1:ixBmax1,ixBmin2:ixBmax2) .ne. F_adv(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)) F_adv(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = 0.d0
      where (F_adv(ixBmin1:ixBmax1,ixBmin2:ixBmax2) .le. 0.d0) &
         F_adv(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = 0.d0

      do ix1 = ixImin1,ixImax1
        Local_gradE(ix1,ixBmin2:ixBmax2) = -(F_bound-F_adv(ixBmax1,&
           ixBmin2:ixBmax2))/w(nghostcells+1,ixBmin2:ixBmax2,i_diff_mg)
      enddo
      gradE = sum(Local_gradE(nghostcells,ixBmin2:ixBmax2))/(ixBmax2-ixBmin2)

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,ixBmin2:ixBmax2,r_e) = w(ix1+2,ixBmin2:ixBmax2,r_e) + (x(ix1,&
           ixBmin2:ixBmax2,1)-x(ix1+2,ixBmin2:ixBmax2,1))*Local_gradE(ix1+1,&
           ixBmin2:ixBmax2)
      enddo

      if (rhd_energy) then
        temp(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = (w(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2,r_e)*unit_pressure/const_rad_a)**&
           0.25d0/unit_temperature
        w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,e_) = w(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2,rho_)*temp(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2)/(rhd_gamma-1.d0) + half*w(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2,mom(1))**2/w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)
      endif

      ! print*, gradE

    case(2)

      !> Compute mean kappa in outer blocks
      call get_kappa_OPAL(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
         ixImax1,ixImax2,w,x,kappa)

      kappa_out = sum(kappa(ixImax1-nghostcells,&
         ixBmin2:ixBmax2))/((ixBmax2-ixBmin2))

      if (kappa_out .ne. kappa_out) kappa_out = kappa_e
      kappa_out = max(kappa_out,kappa_e)
      kappa_out = min(kappa_out,20*kappa_e)

      Local_tauout(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = &
         kappa_out*w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         rho_)*R_star**2/(3*x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,1))
      Local_Tout(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = &
         F_bound/StefBoltz*(3.d0/4.d0*Local_tauout(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2))**0.25d0

      T_out = sum(Local_Tout(ixBmin1,ixBmin2:ixBmax2))/(ixBmax2-ixBmin2)

      T_out = max(1.5d4/unit_temperature, T_out)
      E_out = const_rad_a*(T_out*unit_temperature)**4.d0/unit_pressure

      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,r_e) = &
         const_rad_a*(Local_Tout(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)*unit_temperature)**4.d0/unit_pressure

      ! print*, E_out
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

    case (2)
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = E_out

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
       ixImin2:ixImax2,1:ndim)

    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: mass

    radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,:) = 0.d0
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

    double precision :: k_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndim)

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
    if (rhd_energy) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,e_) + qdt*ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) !> OK
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       r_e) + qdt*ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) !> TROUBLEMAKER

    if (.not. Cak_in_D) then
      call get_kappa_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,wCT,x,k_cak)

      if (fixed_lum) then
        !> Fixed L = L_bound
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) = w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(1)) + qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           rho_)*L_bound/(4*dpi*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1)**2)/const_c*k_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*unit_velocity
        if (rhd_energy) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,e_) + qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(1))*L_bound/(4*dpi*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             1)**2)/const_c*k_cak(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)*unit_velocity
        endif
      else
        !> Local flux
        call fld_get_radflux(wCT, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, rad_flux)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) = w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(1)) + qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           rho_)*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1)/const_c*k_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*unit_velocity
        if (rhd_energy) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,e_) + qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(1))*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             1)/const_c*k_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*unit_velocity
        endif
      endif
    endif

    ! print*, k_cak(10,10)

  end subroutine PseudoPlanar

  subroutine get_dt_cak(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: dx1,dx2, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: mass

    double precision :: dt_cak
    double precision :: k_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndim)

    call get_kappa_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,w,x,k_cak)

    if (fixed_lum) then
      !> Fixed L = L_bound
      dt_cak = courantpar*minval(dsqrt(dxlevel(1)/(L_bound/(4*dpi*x(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)**2)/const_c*k_cak(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*unit_velocity -const_G*mass/radius(ixImin1:ixImax1,&
         ixImin2:ixImax2)**2*(unit_time**2/unit_length))))
    else
      !> Local flux
      call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, rad_flux)
      dt_cak = courantpar*minval(dsqrt(dxlevel(1)/abs(rad_flux(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1)/const_c*k_cak(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*unit_velocity -const_G*mass/radius(ixImin1:ixImax1,&
         ixImin2:ixImax2)**2*(unit_time**2/unit_length))))
    endif

    dtnew = min(dt_cak, dtnew)

  end subroutine get_dt_cak
  !

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
       ixOmin2:ixOmax2,1:ndim)
    double precision :: radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
         pert(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision :: edd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndim,1:ndim)
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
      call fld_get_eddington(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, edd)
      source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = source(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,r_e) + two*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e)*edd(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1,1)/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
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
    if (Cak_in_D) then
      call get_kappa_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x,CAK)
    else
      CAK(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.d0
    endif

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
    

    where(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .ne. kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))
      kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = kappa_e
    endwhere

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
        gradv(ixOmin1:ixOmax1,ixOmin2:ixOmax2), gradvI(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixImin1:ixImax1,ixImin2:ixImax2) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(1))/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    ! call gradientO(vel,x,ixI^L,ixO^L,1,gradv,1)

    call gradient(vel,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,1,gradvI)
    gradv(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = gradvI(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    !> Absolute value of gradient:
    gradv(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs(gradv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))

    xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 1.d0-xprobmin1/x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)

    alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cak_a

    where (xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .le. cak_x0)
      alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cak_base
    elsewhere (xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .le. cak_x1)
      alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cak_base + (cak_a - &
         cak_base)*(xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2) - cak_x0)/(cak_x1 - &
         cak_x0)
    endwhere

    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
       kappa_e*cak_Q/(1-alpha(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)) *(gradv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*unit_velocity/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)*const_c*cak_Q*kappa_e))**alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    ! if (it .le. it_start_cak) then
    !   kappa(ixO^S) = kappa(ixO^S)*dexp(-w(ixO^S,rho_)*kappa_e)
    ! endif

    !> test with no cak
    ! kappa(ixO^S) = 0.d0
  end subroutine get_kappa_CAK

  ! subroutine refine_base(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
  !   ! Enforce additional refinement or coarsening
  !   ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
  !   ! you must set consistent values for integers refine/coarsen:
  !   ! refine = -1 enforce to not refine
  !   ! refine =  0 doesn't enforce anything
  !   ! refine =  1 enforce refinement
  !   ! coarsen = -1 enforce to not coarsen
  !   ! coarsen =  0 doesn't enforce anything
  !   ! coarsen =  1 enforce coarsen
  !   use mod_global_parameters
  !
  !   integer, intent(in) :: igrid, level, ixG^L, ix^L
  !   double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
  !   integer, intent(inout) :: refine, coarsen
  !
  !   !> Refine close to base
  !   refine = -1
  !
  !   if (qt .gt. 1.d0) then
  !     if (any(x(ixG^S,1) < 1.d0)) refine=1
  !   endif
  !
  !   if (qt .gt. 2.d0) then
  !     if (any(x(ixG^S,1) < 1.d0)) refine=1
  !   endif
  !
  !   if (qt .gt. 4.d0) then
  !     if (any(x(ixG^S,1) < 1.d0)) refine=1
  !   endif
  !
  ! end subroutine refine_base


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
        gradv(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                   :: rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:ndim), Lum(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision                   :: pp_rf(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2), lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
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

    call rhd_get_tgas(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, Tgas)
    call rhd_get_trad(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, Trad)

    call get_kappa_OPAL(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,x,OPAL)
    call get_kappa_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,w,x,CAK)

    g_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (OPAL(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)+CAK(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1)/(const_c/unit_velocity)
    g_grav(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
       const_G*mass/radius(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2*(unit_time**2/unit_length)
    big_gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = g_rad(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/g_grav(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    vel(ixImin1:ixImax1,ixImin2:ixImax2) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(1))/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    call gradient(vel,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,1,gradv)

    pp_rf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = two*rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)*dt

    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, lambda, fld_R)

    Lum(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 4*dpi*rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1)*unit_length)**2*unit_radflux/L_sun

    ! w(ixO^S,nw+1) = Trad(ixO^S)*unit_temperature
    ! w(ixO^S,nw+2) = 4*dpi*w(ixO^S,mom(1))*radius(ixO^S)**2 &
    ! *unit_density*unit_velocity/M_sun*year
    ! w(ixO^S,nw+3) = OPAL(ixO^S)/kappa_e
    ! w(ixO^S,nw+4) = CAK(ixO^S)/kappa_e
    ! w(ixO^S,nw+5) = lambda(ixO^S)
    ! w(ixO^S,nw+6) = big_gamma(ixO^S)
    ! w(ixO^S,nw+7) = Lum(ixO^S)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1) = OPAL(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/kappa_e
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2) = CAK(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/kappa_e
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3) = big_gamma(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+4) = Lum(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    ! w(ixO^S,nw+3) = Trad(ixO^S)*unit_temperature
    ! w(ixO^S,nw+4) = Tgas(ixO^S)*unit_temperature
    ! w(ixO^S,nw+4) = lambda(ixO^S)



  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'OPAL CAK Gamma Lum'
    ! varnames = 'Gamma Lum Trad Tgas'
    ! varnames = 'Trad Mdot OPAL CAK lambda Gamma Lum'
  end subroutine specialvarnames_output

end module mod_usr
