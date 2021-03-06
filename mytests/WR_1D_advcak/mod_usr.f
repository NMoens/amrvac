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

  integer :: i_v1, i_v2, i_p
  integer :: i_Trad, i_Tgas, i_Mdot, i_Opal, i_CAK, i_lambda
  integer :: i_Gamma, i_Lum, i_F1, i_F2, i_CAK2

  double precision :: kappa_e, L_bound, Gamma_e_bound, F_bound, gradE, E_out
  logical :: fixed_lum, Cak_in_D, read_cak_table

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_1D")

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

    !> Additional variables
    usr_process_grid => update_extravars

    ! Output routines
    ! usr_aux_output    => specialvar_output
    ! usr_add_aux_names => specialvarnames_output

    ! Timestep for PseudoPlanar
    ! usr_get_dt => get_dt_cak

    ! Refine mesh near base
    usr_refine_grid => refine_base

    ! Active the physics module
    call rhd_activate()

    i_v1 = var_set_extravar("v1", "v1")
    i_p = var_set_extravar("p","p")
    i_Trad = var_set_extravar("Trad", "Trad")
    i_Tgas = var_set_extravar("Tgas", "Tgas")
    i_Mdot = var_set_extravar("Mdot", "Mdot")
    i_Opal = var_set_extravar("OPAL", "OPAL")
    i_CAK = var_set_extravar("CAK", "CAK")
    i_CAK2 = var_set_extravar("CAK2", "CAK2")
    i_lambda = var_set_extravar("lambda", "lambda")
    i_Gamma = var_set_extravar("Gamma", "Gamma")
    i_Lum = var_set_extravar("Lum", "Lum")
    i_F1 = var_set_extravar("F1", "F1")

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_opal_opacity, only: init_opal
    use mod_cak_opacity, only: init_cak

    use mod_fld

    integer :: i

    !> Initialise Opal
    call init_opal(He_abundance)

    !> Initialise CAK tables
    call init_cak

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

    print*, 'L_bound (cgs)', L_bound*(unit_radflux*unit_length**2)
    print*, 'log10(L_bound)', log10(L_bound*(unit_radflux*unit_length**2)/L_sun)
    print*, 'L_bound', L_bound*(unit_radflux*unit_length**2)/L_sun
    ! stop
    print*, 'unit_density', unit_density
    print*, 'unit_time', unit_time
    print*, 'unit_pressure', unit_pressure

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ cak_Q, cak_a, cak_base, cak_x0, cak_x1, rho_bound,&
        kappa_e, T_bound, R_star, M_star, v_inf, Mdot, Gamma_e_bound,&
        it_start_cak, fixed_lum, Cak_in_D, read_cak_table

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
  subroutine initial_conditions(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x)
    use mod_constants

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision :: kappa(ixOmin1:ixOmax1), lambda(ixOmin1:ixOmax1),&
        fld_R(ixOmin1:ixOmax1)
    double precision :: vel(ixImin1:ixImax1)
    double precision :: T_out, E_out, E_gauge
    double precision :: T_in, E_in, rr(ixImin1:ixImax1), bb,&
        temp(ixImin1:ixImax1)

    w(ixImin1:ixImax1,rho_)   = rho_bound
    w(ixImin1:ixImax1,mom(1)) = 0.d0

    where(x(ixImin1:ixImax1,1) .gt. 1.d0)
      vel(ixImin1:ixImax1) =  v_inf*(1 - 0.999d0*( 1.d0/x(ixImin1:ixImax1,&
         1)))**0.5d0
      w(ixImin1:ixImax1,rho_) = Mdot/(4.d0*dpi*vel(ixImin1:ixImax1)*x(&
         ixImin1:ixImax1,1)**2)
    endwhere

    w(ixImin1:ixImax1,mom(1)) = vel(ixImin1:ixImax1)*w(ixImin1:ixImax1,rho_)

    !> Outer/Inner temperature
    T_out = 28445.836732569689/unit_temperature
    E_out = const_rad_a*(T_out*unit_temperature)**4/unit_pressure
    T_in = T_bound
    E_in = const_rad_a*(T_in*unit_temperature)**4/unit_pressure

    !>Very bad initial profile using constant gradE
    ! w(ixI^S,r_e) = E_out + gradE*(x(ixI^S,1)-xprobmax1)
    ! w(ixI^S,r_e) = E_in + gradE*(x(ixI^S,1)-xprobmin1)
    rr(ixImin1:ixImax1) = dsqrt(x(ixImin1:ixImax1,&
       1)-xprobmin1)*(16*x(ixImin1:ixImax1,1)**2 + 8*x(ixImin1:ixImax1,&
       1)+ 6) /(15*x(ixImin1:ixImax1,1)**(5.d0/2))
    bb = -kappa_e*L_bound*Mdot*unit_velocity*3/(16*dpi**2*const_c*v_inf)

    ! bb = bb*2

    ! rr(ixI^S) = 2*dsqrt(x(ixI^S,1)-xprobmin1)/dsqrt(x(ixI^S,1))
    ! bb = kappa_e*F_bound*Mdot*unit_velocity*3/(4*dpi*const_c*v_inf)
    w(ixImin1:ixImax1,r_e) = E_in + bb*rr(ixImin1:ixImax1)

    E_gauge = E_in + bb*dsqrt(xprobmax1-xprobmin1)*(16*xprobmax1**2 + &
       8*xprobmax1+ 6)/(15*xprobmax1**(5.d0/2))

    w(ixImin1:ixImax1,r_e) = w(ixImin1:ixImax1,r_e) - E_gauge + E_out

    if (rhd_energy) then
      temp(ixImin1:ixImax1) = (w(ixImin1:ixImax1,&
         r_e)*unit_pressure/const_rad_a)**0.25d0/unit_temperature
      w(ixImin1:ixImax1,e_) = w(ixImin1:ixImax1,&
         rho_)*temp(ixImin1:ixImax1)/(rhd_gamma-1.d0) + half*w(ixImin1:ixImax1,&
         mom(1))**2/w(ixImin1:ixImax1,rho_)
    endif

    call fld_get_opacity(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, lambda,&
        fld_R)

    w(ixOmin1:ixOmax1,i_diff_mg) = (const_c/unit_velocity)*lambda(&
       ixOmin1:ixOmax1)/(kappa(ixOmin1:ixOmax1)*w(ixOmin1:ixOmax1,rho_))
    w(ixOmin1:ixOmax1,i_diff_mg) = (const_c/unit_velocity)/(3.d0*kappa(&
       ixOmin1:ixOmax1)*w(ixOmin1:ixOmax1,rho_))

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixImin1,ixImax1,ixBmin1,ixBmax1,iB,w,x)
    use mod_physics, only: phys_get_trad
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)             :: ixImin1,ixImax1, ixBmin1,ixBmax1, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision :: kappa(ixImin1:ixImax1), Temp(ixImin1:ixImax1)
    double precision :: Temp0, rho0, T_out, n
    double precision :: Local_gradE(ixImin1:ixImax1), F_adv
    double precision :: Local_tauout(ixBmin1:ixBmax1)
    double precision :: Local_Tout(ixBmin1:ixBmax1)

    double precision :: kappa_out

    integer :: ix1

    select case (iB)

    case(1)
      w(ixBmin1:ixBmax1,rho_) = rho_bound
      do ix1 = ixBmax1-1,ixBmin1,-1
        w(ix1,rho_) = dexp(2*dlog(w(ix1+1,rho_)) - dlog(w(ix1+2,rho_)))
      enddo

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,mom(1)) = w(ix1+1,mom(1))
      enddo

      where(w(ixBmin1:ixBmax1,mom(1)) .lt. 0.d0)
       w(ixBmin1:ixBmax1,mom(1)) = 0.d0
      endwhere

      call get_kappa_OPAL(ixImin1,ixImax1,ixImin1,ixImax1,w,x,kappa)
      do ix1 = ixBmin1,ixBmax1
        kappa(ix1) = kappa(ixBmax1+1)
      enddo

      F_adv = 4.d0/3.d0*(w(ixBmax1,mom(1))/w(ixBmax1,rho_))*w(ixBmax1,&
         r_e) * 4*dpi*xprobmin1**2
      if (F_adv .ne. F_adv) F_adv = 0.d0

      Local_gradE(ixImin1:ixImax1) = -(F_bound-F_adv)/w(nghostcells+1,&
         i_diff_mg)
      gradE = Local_gradE(nghostcells)

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,r_e) = w(ix1+2,r_e) + (x(ix1,1)-x(ix1+2,1))*Local_gradE(ix1+1)
      enddo

      if (rhd_energy) then
        temp(ixBmin1:ixBmax1) = (w(ixBmin1:ixBmax1,&
           r_e)*unit_pressure/const_rad_a)**0.25d0/unit_temperature
        w(ixBmin1:ixBmax1,e_) = w(ixBmin1:ixBmax1,&
           rho_)*temp(ixBmin1:ixBmax1)/(rhd_gamma-1.d0) + &
           half*w(ixBmin1:ixBmax1,mom(1))**2/w(ixBmin1:ixBmax1,rho_)
      endif

      ! print*, it

    case(2)

      !> Compute mean kappa in outer blocks
      call get_kappa_OPAL(ixImin1,ixImax1,ixImin1,ixImax1,w,x,kappa)

      kappa_out = kappa(ixImax1-nghostcells)
      ! kappa_out = kappa_e

      if (kappa_out .ne. kappa_out) kappa_out = kappa_e
      kappa_out = max(kappa_out,kappa_e)
      kappa_out = min(kappa_out,20*kappa_e)

      Local_tauout(ixBmin1:ixBmax1) = kappa_out*w(ixBmin1:ixBmax1,&
         rho_)*R_star**2/(3*x(ixBmin1:ixBmax1,1))
      Local_Tout(ixBmin1:ixBmax1) = F_bound/StefBoltz*(3.d0/4.d0*Local_tauout(&
         ixBmin1:ixBmax1))**0.25d0

      T_out = Local_Tout(ixBmin1)

      T_out = max(1.5d4/unit_temperature, T_out)
      ! T_out = max(3.5d4/unit_temperature, T_out)
      E_out = const_rad_a*(T_out*unit_temperature)**4.d0/unit_pressure

      do ix1 = ixBmin1,ixBmax1
        w(ix1,r_e) = 2*w(ix1-1,r_e) - w(ix1-2,r_e)
      enddo


      w(ixBmin1:ixBmax1,r_e) = const_rad_a*(Local_Tout(&
         ixBmin1:ixBmax1)*unit_temperature)**4.d0/unit_pressure

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

    case (2)
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = E_out

    case default
      print *, "Not a standard: ", trim(typeboundary(r_e, iB))
      error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions


  !> Calculate gravitational acceleration in each dimension
  subroutine set_gravitation_field(ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,x,&
     gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw)
    double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,ndim)

    double precision :: radius(ixImin1:ixImax1)
    double precision :: mass

    radius(ixOmin1:ixOmax1) = x(ixOmin1:ixOmax1,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    gravity_field(ixImin1:ixImax1,1) = -const_G*mass/radius(ixImin1:ixImax1)**&
       2*(unit_time**2/unit_length)

  end subroutine set_gravitation_field

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine PseudoPlanar(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,&
     wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,&
       iwmax
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision :: ppsource(ixOmin1:ixOmax1,1:nw)

    double precision :: k_cak(ixOmin1:ixOmax1), rad_flux(ixOmin1:ixOmax1,&
       1:ndim)

    call PseudoPlanarSource(ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,x,ppsource)
    w(ixOmin1:ixOmax1,rho_) = w(ixOmin1:ixOmax1,&
       rho_) + qdt*ppsource(ixOmin1:ixOmax1,rho_) !> OK
    w(ixOmin1:ixOmax1,mom(1)) = w(ixOmin1:ixOmax1,&
       mom(1)) + qdt*ppsource(ixOmin1:ixOmax1,mom(1)) !> OK
    if (rhd_energy) w(ixOmin1:ixOmax1,e_) = w(ixOmin1:ixOmax1,&
       e_) + qdt*ppsource(ixOmin1:ixOmax1,e_) !> OK
    w(ixOmin1:ixOmax1,r_e) = w(ixOmin1:ixOmax1,&
       r_e) + qdt*ppsource(ixOmin1:ixOmax1,r_e) !> TROUBLEMAKER

    if (.not. Cak_in_D) then
      if (read_cak_table) then
        call get_kappa_CAK2(ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,x,k_cak)
      else
        call get_kappa_CAK(ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,x,k_cak)
      endif

      if (fixed_lum) then
        !> Fixed L = L_bound
        w(ixOmin1:ixOmax1,mom(1)) = w(ixOmin1:ixOmax1,&
           mom(1)) + qdt*wCT(ixOmin1:ixOmax1,&
           rho_)*L_bound/(4*dpi*x(ixOmin1:ixOmax1,&
           1)**2)/const_c*k_cak(ixOmin1:ixOmax1)*unit_velocity
        if (rhd_energy) then
          w(ixOmin1:ixOmax1,e_) = w(ixOmin1:ixOmax1,&
             e_) + qdt*wCT(ixOmin1:ixOmax1,&
             mom(1))*L_bound/(4*dpi*x(ixOmin1:ixOmax1,&
             1)**2)/const_c*k_cak(ixOmin1:ixOmax1)*unit_velocity
        endif
      else
        !> Local flux
        call fld_get_radflux(wCT, x, ixImin1,ixImax1, ixOmin1,ixOmax1,&
            rad_flux)
        w(ixOmin1:ixOmax1,mom(1)) = w(ixOmin1:ixOmax1,&
           mom(1)) + qdt*wCT(ixOmin1:ixOmax1,rho_)*rad_flux(ixOmin1:ixOmax1,&
           1)/const_c*k_cak(ixOmin1:ixOmax1)*unit_velocity
        if (rhd_energy) then
          w(ixOmin1:ixOmax1,e_) = w(ixOmin1:ixOmax1,&
             e_) + qdt*wCT(ixOmin1:ixOmax1,mom(1))*rad_flux(ixOmin1:ixOmax1,&
             1)/const_c*k_cak(ixOmin1:ixOmax1)*unit_velocity
        endif
      endif
    endif

    ! print*, k_cak(10,10)

  end subroutine PseudoPlanar

  subroutine get_dt_cak(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: dx1, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: radius(ixImin1:ixImax1)
    double precision :: mass

    double precision :: dt_cak
    double precision :: k_cak(ixOmin1:ixOmax1), rad_flux(ixOmin1:ixOmax1,&
       1:ndim)

    if (read_cak_table) then
      call get_kappa_CAK2(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,k_cak)
    else
      call get_kappa_CAK(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,k_cak)
    endif

    if (fixed_lum) then
      !> Fixed L = L_bound
      dt_cak = courantpar*minval(dsqrt(dxlevel(1)/(L_bound/(4*dpi*x(&
         ixOmin1:ixOmax1,1)**2)/const_c*k_cak(ixOmin1:ixOmax1)*unit_velocity &
         -const_G*mass/radius(ixImin1:ixImax1)**2*(unit_time**&
         2/unit_length))))
    else
      !> Local flux
      call fld_get_radflux(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, rad_flux)
      dt_cak = courantpar*minval(dsqrt(dxlevel(1)/abs(rad_flux(ixOmin1:ixOmax1,&
         1)/const_c*k_cak(ixOmin1:ixOmax1)*unit_velocity &
         -const_G*mass/radius(ixImin1:ixImax1)**2*(unit_time**&
         2/unit_length))))
    endif

    dtnew = min(dt_cak, dtnew)

  end subroutine get_dt_cak
  !

  subroutine PseudoPlanarSource(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,source)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(out) :: source(ixOmin1:ixOmax1,1:nw)

    double precision :: rad_flux(ixOmin1:ixOmax1,1:ndir)
    double precision :: pth(ixImin1:ixImax1),v(ixOmin1:ixOmax1,1:ndim)
    double precision :: radius(ixOmin1:ixOmax1),  pert(ixOmin1:ixOmax1)

    double precision :: edd(ixOmin1:ixOmax1,1:ndim,1:ndim)

    integer :: rdir

    source(ixOmin1:ixOmax1,1:nw) = zero

    rdir = 1

    v(ixOmin1:ixOmax1,rdir) = w(ixOmin1:ixOmax1,mom(rdir))/w(ixOmin1:ixOmax1,&
       rho_)

    radius(ixOmin1:ixOmax1) = x(ixOmin1:ixOmax1,rdir) !+ half*block%dx(ixOmin1:ixOmax1,rdir)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    source(ixOmin1:ixOmax1,rho_) = -two*w(ixOmin1:ixOmax1,&
       rho_)*v(ixOmin1:ixOmax1,rdir)/radius(ixOmin1:ixOmax1)

    call phys_get_pthermal(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,pth)

    !> dm_r/dt = +(rho*v_p**2 + 2pth)/r -2 (rho*v_r**2 + pth)/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    source(ixOmin1:ixOmax1,mom(rdir)) = - 2*w(ixOmin1:ixOmax1,&
       rho_)*v(ixOmin1:ixOmax1,rdir)**two/radius(ixOmin1:ixOmax1)

    !> de/dt = -2 (e+p) v_r/r
    if (rhd_energy) source(ixOmin1:ixOmax1,e_) = -two*(w(ixOmin1:ixOmax1,&
       e_)+pth(ixOmin1:ixOmax1))*v(ixOmin1:ixOmax1,&
       rdir)/radius(ixOmin1:ixOmax1)

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      !> THIS BAD BOiii IS GIVING US SOME TROUBLE
      call fld_get_radflux(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, rad_flux)
      source(ixOmin1:ixOmax1,r_e) = source(ixOmin1:ixOmax1,&
         r_e) - two*rad_flux(ixOmin1:ixOmax1,rdir)/radius(ixOmin1:ixOmax1)
    endif

    if (rhd_radiation_advection) then
      source(ixOmin1:ixOmax1,r_e) = source(ixOmin1:ixOmax1,&
         r_e) - two*w(ixOmin1:ixOmax1,r_e)*v(ixOmin1:ixOmax1,&
         rdir)/radius(ixOmin1:ixOmax1)
    endif

    ! Not sure about this one
    if (rhd_radiation_force) then
      call fld_get_eddington(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, edd)
      source(ixOmin1:ixOmax1,r_e) = source(ixOmin1:ixOmax1,&
         r_e) + two*v(ixOmin1:ixOmax1,rdir)*w(ixOmin1:ixOmax1,&
         r_e)*edd(ixOmin1:ixOmax1,1,1)/radius(ixOmin1:ixOmax1)
    endif

  end subroutine PseudoPlanarSource


  subroutine OPAL_and_CAK(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,kappa)
    use mod_physics, only: phys_get_trad
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1,1:nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1)

    double precision :: OPAL(ixOmin1:ixOmax1), CAK(ixOmin1:ixOmax1)

    !> Get OPAL opacities by reading from table
    call get_kappa_OPAL(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,OPAL)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    if (Cak_in_D) then
      if (read_cak_table) then
        call get_kappa_CAK2(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,CAK)
      else
        call get_kappa_CAK(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,CAK)
      endif
    else
      CAK(ixOmin1:ixOmax1) = 0.d0
    endif

    !> Add OPAL and CAK for total opacity
    kappa(ixOmin1:ixOmax1) = OPAL(ixOmin1:ixOmax1) + CAK(ixOmin1:ixOmax1)

    where(kappa(ixOmin1:ixOmax1) .ne. kappa(ixOmin1:ixOmax1))
      kappa(ixOmin1:ixOmax1) = kappa_e
    endwhere

  end subroutine OPAL_and_CAK


  subroutine get_kappa_OPAL(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,kappa)
    use mod_physics, only: phys_get_trad, phys_get_tgas
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1,1:nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1)

    integer :: ix1
    double precision :: Temp(ixImin1:ixImax1)
    double precision :: n, rho0, Temp0

    !> Get OPAL opacities by reading from table
    if (rhd_energy) then
      call phys_get_tgas(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,Temp)
    else
      call phys_get_trad(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,Temp)
    endif

    do ix1=ixOmin1,ixOmax1
    
        rho0 = w(ix1,rho_)*unit_density
        Temp0 = Temp(ix1)*unit_temperature
        Temp0 = max(Temp0,1.d4)
        call set_opal_opacity(rho0,Temp0,n)
        kappa(ix1) = n/unit_opacity
    enddo
    

    where(kappa(ixOmin1:ixOmax1) .ne. kappa(ixOmin1:ixOmax1))
      kappa(ixOmin1:ixOmax1) = kappa_e
    endwhere

    !> test without opal
    ! kappa(ixO^S) = kappa_e

  end subroutine get_kappa_OPAL

  subroutine get_kappa_CAK(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,kappa)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1,1:nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1)

    double precision :: vel(ixImin1:ixImax1), gradv(ixOmin1:ixOmax1),&
        gradvI(ixImin1:ixImax1)
    double precision :: xx(ixOmin1:ixOmax1), alpha(ixOmin1:ixOmax1)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixImin1:ixImax1) = w(ixImin1:ixImax1,mom(1))/w(ixImin1:ixImax1,rho_)
    call gradientO(vel,x,ixImin1,ixImax1,ixOmin1,ixOmax1,1,gradv,nghostcells)

    ! call gradient(vel,ixI^L,ixO^L,1,gradvI)
    ! gradv(ixO^S) = gradvI(ixO^S)

    !> Absolute value of gradient:
    gradv(ixOmin1:ixOmax1) = abs(gradv(ixOmin1:ixOmax1))

    xx(ixOmin1:ixOmax1) = 1.d0-xprobmin1/x(ixOmin1:ixOmax1,1)

    alpha(ixOmin1:ixOmax1) = cak_a

    where (xx(ixOmin1:ixOmax1) .le. cak_x0)
      alpha(ixOmin1:ixOmax1) = cak_base
    elsewhere (xx(ixOmin1:ixOmax1) .le. cak_x1)
      alpha(ixOmin1:ixOmax1) = cak_base + (cak_a - &
         cak_base)*(xx(ixOmin1:ixOmax1) - cak_x0)/(cak_x1 - cak_x0)
    endwhere

    kappa(ixOmin1:ixOmax1) = kappa_e*cak_Q/(1-alpha(ixOmin1:ixOmax1)) &
       *(gradv(ixOmin1:ixOmax1)*unit_velocity/(w(ixOmin1:ixOmax1,&
       rho_)*const_c*cak_Q*kappa_e))**alpha(ixOmin1:ixOmax1)

    if (x(ixImax1,1) .ge. xprobmax1) then
      kappa(ixOmax1) = kappa(ixOmax1-1)
    endif

  end subroutine get_kappa_CAK

  subroutine get_kappa_CAK2(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,kappa)
    use mod_physics, only: phys_get_trad, phys_get_tgas
    use mod_global_parameters
    use mod_cak_opacity
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1,1:nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1)

    double precision :: Temp(ixImin1:ixImax1), rho0, temp0, gradv0, kap0
    integer :: ix1

    double precision :: alpha, Qbar, Q0, kappa_e_t
    double precision :: tau, M_t
    double precision :: vel(ixImin1:ixImax1), gradv(ixOmin1:ixOmax1),&
        gradvI(ixImin1:ixImax1)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixImin1:ixImax1) = w(ixImin1:ixImax1,mom(1))/w(ixImin1:ixImax1,rho_)
    call gradientO(vel,x,ixImin1,ixImax1,ixOmin1,ixOmax1,1,gradv,nghostcells)

    ! call gradient(vel,ixI^L,ixO^L,1,gradvI)
    ! gradv(ixO^S) = gradvI(ixO^S)

    !> Absolute value of gradient:
    gradv(ixOmin1:ixOmax1) = abs(gradv(ixOmin1:ixOmax1))

    !> Get CAK opacities by reading from table
    if (rhd_energy) then
      call phys_get_tgas(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,Temp)
    else
      call phys_get_trad(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,Temp)
    endif

    do ix1=ixOmin1,ixOmax1
    
        rho0 = w(ix1,rho_)*unit_density
        Temp0 = Temp(ix1)*unit_temperature
        Temp0 = max(Temp0,1.d4)
        gradv0 = gradv(ix1)*(unit_velocity/unit_length)
        call set_cak_opacity(rho0,Temp0,gradv0,alpha, Qbar, Q0, kappa_e_t)

        tau = (kappa_e*unit_opacity)*rho0*const_c/gradv0
        M_t = Qbar/(1-alpha)*((1+Q0*tau)**(1-alpha) - 1)/(Q0*tau)
        kap0 = (kappa_e*unit_opacity)*M_t

        kappa(ix1) = kap0/unit_opacity

        if (kappa(ix1) .ne. kappa(ix1)) kappa(ix1) = 0.d0

        kappa(ix1) = min(50*kappa_e,kappa(ix1))
    enddo
    


    if (x(ixImax1,1) .ge. xprobmax1) then
      kappa(ixOmax1) = kappa(ixOmax1-1)
    endif

  end subroutine get_kappa_CAK2


  subroutine ceil_diffcoef(w, wCT, x, ixImin1,ixImax1, ixOmin1,ixOmax1)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(inout) :: w(ixImin1:ixImax1, 1:nw)
    double precision, intent(in) :: wCT(ixImin1:ixImax1, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1, 1:ndim)


    where (w(ixOmin1:ixOmax1,i_diff_mg) .gt. 1.5d3) w(ixOmin1:ixOmax1,&
       i_diff_mg) = 1.5d3

  end subroutine ceil_diffcoef

  subroutine refine_base(igrid,level,ixGmin1,ixGmax1,ixmin1,ixmax1,qt,w,x,&
     refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixGmin1,ixGmax1, ixmin1,ixmax1
    double precision, intent(in) :: qt, w(ixGmin1:ixGmax1,1:nw),&
        x(ixGmin1:ixGmax1,1:ndim)
    integer, intent(inout) :: refine, coarsen

    double precision :: lim_1, lim_2, lim_3, lim_4

    lim_3 = 1.5d0
    lim_2 = 2.5d0
    lim_1 = 4.d0

    !refine= -1
    !coarsen= -1

    if (all(x(ixGmin1:ixGmax1,1) < lim_3)) then
      if (level > 4) coarsen=1
      if (level < 4) refine=1
    elseif (all(x(ixGmin1:ixGmax1,1) < lim_2)) then
      if (level > 3) coarsen=1
      if (level < 3) refine=1
    elseif (all(x(ixGmin1:ixGmax1,1) < lim_1)) then
      if (level > 2) coarsen=1
      if (level < 2) refine=1
    endif

  end subroutine refine_base



  subroutine update_extravars(igrid,level,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,&
     x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixImin1,ixImax1,ixOmin1,&
       ixOmax1
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision                   :: g_rad(ixImin1:ixImax1),&
        big_gamma(ixImin1:ixImax1)
    double precision                   :: g_grav(ixImin1:ixImax1)
    double precision                   :: Tgas(ixImin1:ixImax1),&
       Trad(ixImin1:ixImax1)
    double precision                   :: kappa(ixOmin1:ixOmax1),&
        OPAL(ixOmin1:ixOmax1), CAK(ixOmin1:ixOmax1), CAK2(ixOmin1:ixOmax1)
    double precision                   :: vel(ixImin1:ixImax1),&
        gradv(ixImin1:ixImax1)
    double precision                   :: rad_flux(ixOmin1:ixOmax1,1:ndim),&
        Lum(ixOmin1:ixOmax1)
    double precision                   :: pp_rf(ixOmin1:ixOmax1),&
        lambda(ixOmin1:ixOmax1), fld_R(ixOmin1:ixOmax1)
    integer                            :: idim
    double precision :: radius(ixImin1:ixImax1)
    double precision :: mass

    radius(ixOmin1:ixOmax1) = x(ixOmin1:ixOmax1,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    call fld_get_opacity(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, kappa)
    call fld_get_radflux(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, rad_flux)

    call rhd_get_tgas(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, Tgas)
    call rhd_get_trad(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, Trad)

    call get_kappa_OPAL(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,OPAL)
    call get_kappa_CAK(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,CAK)

    call get_kappa_CAK2(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,CAK2)

    g_rad(ixOmin1:ixOmax1) = (OPAL(ixOmin1:ixOmax1)+&
       CAK(ixOmin1:ixOmax1))*rad_flux(ixOmin1:ixOmax1,&
       1)/(const_c/unit_velocity)
    g_grav(ixOmin1:ixOmax1) = const_G*mass/radius(ixOmin1:ixOmax1)**&
       2*(unit_time**2/unit_length)
    big_gamma(ixOmin1:ixOmax1) = g_rad(ixOmin1:ixOmax1)/g_grav(&
       ixOmin1:ixOmax1)

    vel(ixImin1:ixImax1) = w(ixImin1:ixImax1,mom(1))/w(ixImin1:ixImax1,rho_)
    call gradient(vel,ixImin1,ixImax1,ixOmin1,ixOmax1,1,gradv)

    pp_rf(ixOmin1:ixOmax1) = two*rad_flux(ixOmin1:ixOmax1,1)/x(ixOmin1:ixOmax1,&
       1)*dt

    call fld_get_fluxlimiter(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, lambda,&
        fld_R)

    Lum(ixOmin1:ixOmax1) = 4*dpi*rad_flux(ixOmin1:ixOmax1,&
       1)*(x(ixOmin1:ixOmax1,1)*unit_length)**2*unit_radflux/L_sun

    w(ixOmin1:ixOmax1,i_v1) = w(ixOmin1:ixOmax1,mom(1))/w(ixOmin1:ixOmax1,&
       rho_)
    w(ixOmin1:ixOmax1,i_p) = (w(ixOmin1:ixOmax1,&
       e_) - 0.5d0 * sum(w(ixOmin1:ixOmax1, mom(:))**2,&
        dim=ndim+1) / w(ixOmin1:ixOmax1, rho_)) *(rhd_gamma - 1)

    w(ixOmin1:ixOmax1,i_Trad) = Trad(ixOmin1:ixOmax1)*unit_temperature
    w(ixOmin1:ixOmax1,i_Tgas) = Tgas(ixOmin1:ixOmax1)*unit_temperature
    w(ixOmin1:ixOmax1,i_Mdot) = 4*dpi*w(ixOmin1:ixOmax1,&
       mom(1))*radius(ixOmin1:ixOmax1)**2 &
       *unit_density*unit_velocity/M_sun*year
    w(ixOmin1:ixOmax1,i_Opal) = OPAL(ixOmin1:ixOmax1)/kappa_e
    w(ixOmin1:ixOmax1,i_CAK) = CAK(ixOmin1:ixOmax1)/kappa_e
    w(ixOmin1:ixOmax1,i_CAK2) = CAK2(ixOmin1:ixOmax1)/kappa_e
    w(ixOmin1:ixOmax1,i_lambda) = lambda(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,i_Gamma) = big_gamma(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,i_Lum) = Lum(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,i_F1) = rad_flux(ixOmin1:ixOmax1,1)/F_bound

  end subroutine update_extravars
end module mod_usr
