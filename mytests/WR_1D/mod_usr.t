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
  logical :: fixed_lum

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

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    ! Timestep for PseudoPlanar
    usr_get_dt => get_dt_cak

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
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB)
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

    StefBoltz = const_rad_a*const_c/4.d0*(unit_temperature**4.d0)/(unit_velocity*unit_pressure)

    !> Very bad initial guess for gradE using kappa_e
    gradE = -F_bound*3*kappa_e*rho_bound*unit_velocity/const_c

    print*, 'L_bound', L_bound*(unit_radflux*unit_length**2), log10(L_bound*(unit_radflux*unit_length**2)/L_sun)
    ! stop

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ cak_Q, cak_a, cak_base, cak_x0, cak_x1, rho_bound, kappa_e, &
    T_bound, R_star, M_star, v_inf, Mdot, Gamma_e_bound, it_start_cak, fixed_lum

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
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_constants

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)
    double precision :: vel(ixI^S)
    double precision :: T_out, E_out, E_gauge
    double precision :: T_in, E_in, rr(ixI^S), bb

    w(ixI^S,rho_)   = rho_bound
    w(ixI^S,mom(1)) = 0.d0

    where(x(ixI^S,1) .gt. 1.d0)
      vel(ixI^S) =  v_inf*(1 - 0.999d0*( 1.d0/x(ixI^S,1)))**0.5d0
      w(ixI^S,rho_) = Mdot/(4.d0*dpi*vel(ixI^S)*x(ixI^S,1)**2)
    endwhere

    w(ixI^S,mom(1)) = vel(ixI^S)*w(ixI^S,rho_)

    !> Outer/Inner temperature
    T_out = 28445.836732569689/unit_temperature
    E_out = const_rad_a*(T_out*unit_temperature)**4/unit_pressure
    T_in = T_bound
    E_in = const_rad_a*(T_in*unit_temperature)**4/unit_pressure

    !>Very bad initial profile using constant gradE
    ! w(ixI^S,r_e) = E_out + gradE*(x(ixI^S,1)-xprobmax1)
    ! w(ixI^S,r_e) = E_in + gradE*(x(ixI^S,1)-xprobmin1)
    rr(ixI^S) = dsqrt(x(ixI^S,1)-xprobmin1)*(16*x(ixI^S,1)**2 + 8*x(ixI^S,1)+ 6) &
                /(15*x(ixI^S,1)**(5.d0/2))
    bb = -kappa_e*L_bound*Mdot*unit_velocity*3/(16*dpi**2*const_c*v_inf)

    ! bb = bb*2

    ! rr(ixI^S) = 2*dsqrt(x(ixI^S,1)-xprobmin1)/dsqrt(x(ixI^S,1))
    ! bb = kappa_e*F_bound*Mdot*unit_velocity*3/(4*dpi*const_c*v_inf)
    w(ixI^S,r_e) = E_in + bb*rr(ixI^S)

    E_gauge = E_in + bb*dsqrt(xprobmax1-xprobmin1)&
              *(16*xprobmax1**2 + 8*xprobmax1+ 6)/(15*xprobmax1**(5.d0/2))

    w(ixI^S,r_e) = w(ixI^S,r_e) - E_gauge + E_out

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))
    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)/(3.d0*kappa(ixO^S)*w(ixO^S,rho_))

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_physics, only: phys_get_trad
    use mod_global_parameters
    use mod_opacity
    use mod_fld

    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: kappa(ixI^S), Temp(ixI^S)
    double precision :: Temp0, rho0, T_out, n
    double precision :: Local_gradE(ixI^S)
    double precision :: Local_tauout(ixB^S)
    double precision :: Local_Tout(ixB^S)

    double precision :: kappa_out

    integer :: ix^D

    select case (iB)

    case(1)
      w(ixB^S,rho_) = rho_bound
      do ix1 = ixBmax1-1,ixBmin1,-1
        w(ix1,rho_) = dexp(2*dlog(w(ix1+1,rho_)) - dlog(w(ix1+2,rho_)))
      enddo

      w(ixB^S,mom(2)) = 0.d0

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,mom(1)) = w(ix1+1,mom(1))
      enddo

      !where(w(ixB^S,mom(1)) .lt. 0.d0)
      !  w(ixB^S,mom(1)) = 0.d0
      !endwhere

      !where(w(ixB^S,mom(1))/w(ixB^S,rho_) .gt. 0.5d0)
      !  w(ixB^S,mom(1)) = 0.1d0*w(ixB^S,rho_)
      !endwhere

      ! w(ixB^S,r_e) = const_rad_a*(T_bound*unit_temperature)**4/unit_pressure

      call get_kappa_OPAL(ixI^L,ixI^L,w,x,kappa)
      do ix1 = ixBmin1,ixBmax1
        kappa(ix1) = kappa(ixBmax1+1)
      enddo

      Local_gradE(ixI^S) = -F_bound*3*kappa(ixI^S)*w(ixI^S,rho_)&
      *unit_velocity/const_c
      gradE = Local_gradE(nghostcells)

      ! print*, gradE

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,r_e) = w(ix1+2,r_e) &
        + (x(ix1,1)-x(ix1+2,1))*Local_gradE(ix1+1)
      enddo

    case(2)

      !> Compute mean kappa in outer blocks
      call OPAL_and_CAK(ixI^L,ixI^L,w,x,kappa)
      kappa_out = kappa(ixImax1-nghostcells)
      ! kappa_out = kappa_e

      if (kappa_out .ne. kappa_out) kappa_out = kappa_e
      kappa_out = max(kappa_out,kappa_e)
      kappa_out = min(kappa_out,20*kappa_e)

      Local_tauout(ixB^S) = kappa_out*w(ixB^S,rho_)*R_star**2/(3*x(ixB^S,1))
      Local_Tout(ixB^S) = F_bound/StefBoltz*(3.d0/4.d0*Local_tauout(ixB^S))**0.25d0

      T_out = Local_Tout(ixBmin1)

      T_out = max(1.5d4/unit_temperature, T_out)
      E_out = const_rad_a*(T_out*unit_temperature)**4.d0/unit_pressure

      do ix1 = ixBmin1,ixBmax1
        w(ix1,r_e) = 2*w(ix1-1,r_e) - w(ix1-2,r_e)
      enddo

      ! w(ixB^S,r_e) = const_rad_a*(Local_Tout(ixB^S)*unit_temperature)**4.d0/unit_pressure

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
  subroutine set_gravitation_field(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    double precision :: radius(ixI^S)
    double precision :: mass

    radius(ixO^S) = x(ixO^S,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    gravity_field(ixI^S,1) = -const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length)

  end subroutine set_gravitation_field

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine PseudoPlanar(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: ppsource(ixO^S,1:nw)

    double precision :: k_cak(ixO^S), rad_flux(ixO^S,1:ndim)

    call PseudoPlanarSource(ixI^L,ixO^L,wCT,x,ppsource)
    w(ixO^S,rho_) = w(ixO^S,rho_) + qdt*ppsource(ixO^S,rho_) !> OK
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt*ppsource(ixO^S,mom(1)) !> OK
    w(ixO^S,r_e) = w(ixO^S,r_e) + qdt*ppsource(ixO^S,r_e) !> TROUBLEMAKER

    call get_kappa_CAK(ixI^L,ixO^L,wCT,x,k_cak)

    if (fixed_lum) then
      !> Fixed L = L_bound
      w(ixO^S,mom(1)) = w(ixO^S,mom(1)) &
        + qdt*wCT(ixO^S,rho_)*L_bound/(4*dpi*x(ixO^S,1)**2*const_c)*k_cak(ixO^S)*unit_velocity
    else
      !> Local flux
      call fld_get_radflux(wCT, x, ixI^L, ixO^L, rad_flux)
      w(ixO^S,mom(1)) = w(ixO^S,mom(1)) &
        + qdt*wCT(ixO^S,rho_)*rad_flux(ixO^S,1)/const_c*k_cak(ixO^S)*unit_velocity
    endif

    ! print*, k_cak(10,10)

  end subroutine PseudoPlanar

  subroutine get_dt_cak(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: radius(ixI^S)
    double precision :: mass

    double precision :: dt_cak
    double precision :: k_cak(ixO^S), rad_flux(ixO^S,1:ndim)

    call get_kappa_CAK(ixI^L,ixO^L,w,x,k_cak)

    if (fixed_lum) then
      !> Fixed L = L_bound
      dt_cak = courantpar*minval(dsqrt(dxlevel(1)/(L_bound/(4*dpi*x(ixO^S,1)**2*const_c)*k_cak(ixO^S)*unit_velocity &
      -const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length))))
    else
      !> Local flux
      call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
      dt_cak = courantpar*minval(dsqrt(dxlevel(1)/abs(rad_flux(ixO^S,1)/const_c*k_cak(ixO^S)*unit_velocity &
      -const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length))))
    endif

    dtnew = min(dt_cak, dtnew)

  end subroutine get_dt_cak
  !

  subroutine PseudoPlanarSource(ixI^L,ixO^L,w,x,source)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: source(ixO^S,1:nw)

    double precision :: rad_flux(ixO^S,1:ndir)
    double precision :: pth(ixI^S),v(ixO^S,1:ndim)
    double precision :: radius(ixO^S),  pert(ixO^S)
    integer :: rdir

    source(ixO^S,1:nw) = zero

    rdir = 1

    v(ixO^S,rdir) = w(ixO^S,mom(rdir))/w(ixO^S,rho_)

    radius(ixO^S) = x(ixO^S,rdir) ! + half*block%dx(ixO^S,rdir)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    source(ixO^S,rho_) = -two*w(ixO^S,rho_)*v(ixO^S,rdir)/radius(ixO^S)

    call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)

    !> dm_r/dt = +(rho*v_p**2 + 2pth)/r -2 (rho*v_r**2 + pth)/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    source(ixO^S,mom(rdir)) = - 2*w(ixO^S,rho_)*v(ixO^S,rdir)**two/radius(ixO^S)

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      !> THIS BAD BOiii IS GIVING US SOME TROUBLE
      call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
      source(ixO^S,r_e) = source(ixO^S,r_e) - two*rad_flux(ixO^S,rdir)/radius(ixO^S)
    endif

    if (rhd_radiation_advection) then
      source(ixO^S,r_e) = source(ixO^S,r_e) - two*w(ixO^S,r_e)*v(ixO^S,rdir)/radius(ixO^S)
    endif

    ! Not sure about this one
    if (rhd_radiation_force) then
      source(ixO^S,r_e) = source(ixO^S,r_e) + two*v(ixO^S,rdir)*w(ixO^S,r_e)/(3*radius(ixO^S))
    endif

  end subroutine PseudoPlanarSource


  subroutine OPAL_and_CAK(ixI^L,ixO^L,w,x,kappa)
    use mod_physics, only: phys_get_trad
    use mod_global_parameters
    use mod_opacity
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    double precision :: OPAL(ixO^S), CAK(ixO^S)

    !> Get OPAL opacities by reading from table
    call get_kappa_OPAL(ixI^L,ixO^L,w,x,OPAL)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    ! call get_kappa_CAK(ixI^L,ixO^L,w,x,CAK)
    CAK(ixO^S) = 0.d0

    !> Add OPAL and CAK for total opacity
    kappa(ixO^S) = OPAL(ixO^S) + CAK(ixO^S)

    where(kappa(ixO^S) .ne. kappa(ixO^S))
      kappa(ixO^S) = kappa_e
    endwhere

  end subroutine OPAL_and_CAK


  subroutine get_kappa_OPAL(ixI^L,ixO^L,w,x,kappa)
    use mod_physics, only: phys_get_trad
    use mod_global_parameters
    use mod_opacity
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    integer :: ix^D
    double precision :: Temp(ixI^S)
    double precision :: n, rho0, Temp0

    !> Get OPAL opacities by reading from table
    call phys_get_trad(w,x,ixI^L,ixO^L,Temp)
    {do ix^D=ixOmin^D,ixOmax^D\ }
        rho0 = w(ix^D,rho_)*unit_density
        Temp0 = Temp(ix^D)*unit_temperature
        Temp0 = max(Temp0,1.d4)
        call set_opal_opacity(rho0,Temp0,n)
        kappa(ix^D) = n/unit_opacity
    {enddo\ }

    where(kappa(ixO^S) .ne. kappa(ixO^S))
      kappa(ixO^S) = kappa_e
    endwhere

    !> test without opal
    ! kappa(ixO^S) = kappa_e

  end subroutine get_kappa_OPAL

  subroutine get_kappa_CAK(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    use mod_opacity
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    double precision :: vel(ixI^S), gradv(ixO^S), gradvI(ixI^S)
    double precision :: xx(ixO^S), alpha(ixO^S)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    ! call gradientO(vel,x,ixI^L,ixO^L,1,gradv,1)

    call gradient(vel,ixI^L,ixO^L,1,gradvI)
    gradv(ixO^S) = gradvI(ixO^S)

    !> Absolute value of gradient:
    gradv(ixO^S) = abs(gradv(ixO^S))

    xx(ixO^S) = 1.d0-xprobmin1/x(ixO^S,1)

    alpha(ixO^S) = cak_a

    where (xx(ixO^S) .le. cak_x0)
      alpha(ixO^S) = cak_base
    elsewhere (x(ixO^S,1) .le. cak_x1)
      alpha(ixO^S) = cak_base + (cak_a - cak_base)&
      *(xx(ixO^S) - cak_x0)/(cak_x1 - cak_x0)
    endwhere

    kappa(ixO^S) = kappa_e*cak_Q/(1-alpha(ixO^S)) &
    *(gradv(ixO^S)*unit_velocity/(w(ixO^S,rho_)*const_c*cak_Q*kappa_e))**alpha(ixO^S)

    if (it .le. it_start_cak) then
      kappa(ixO^S) = kappa(ixO^S)*dexp(-w(ixO^S,rho_)*kappa_e)
    endif

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


  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_physics
    use mod_fld

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: g_rad(ixI^S), big_gamma(ixI^S)
    double precision                   :: g_grav(ixI^S)
    double precision                   :: Tgas(ixI^S),Trad(ixI^S)
    double precision                   :: kappa(ixO^S), OPAL(ixO^S), CAK(ixO^S)
    double precision                   :: vel(ixI^S), gradv(ixI^S)
    double precision                   :: rad_flux(ixO^S,1:ndim), Lum(ixO^S)
    double precision                   :: pp_rf(ixO^S), lambda(ixO^S), fld_R(ixO^S)
    integer                            :: idim
    double precision :: radius(ixI^S)
    double precision :: mass

    radius(ixO^S) = x(ixO^S,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)

    call rhd_get_tgas(w, x, ixI^L, ixO^L, Tgas)
    call rhd_get_trad(w, x, ixI^L, ixO^L, Trad)

    call get_kappa_OPAL(ixI^L,ixO^L,w,x,OPAL)
    call get_kappa_CAK(ixI^L,ixO^L,w,x,CAK)

    g_rad(ixO^S) = (OPAL(ixO^S)+CAK(ixO^S))*rad_flux(ixO^S,1)/(const_c/unit_velocity)
    g_grav(ixO^S) = const_G*mass/radius(ixO^S)**2*(unit_time**2/unit_length)
    big_gamma(ixO^S) = g_rad(ixO^S)/g_grav(ixO^S)

    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradient(vel,ixI^L,ixO^L,1,gradv)

    pp_rf(ixO^S) = two*rad_flux(ixO^S,1)/x(ixO^S,1)*dt

    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    Lum(ixO^S) = 4*dpi*rad_flux(ixO^S,1)*(x(ixO^S,1)*unit_length)**2*unit_radflux/L_sun

    ! w(ixO^S,nw+1) = Trad(ixO^S)*unit_temperature
    ! w(ixO^S,nw+2) = 4*dpi*w(ixO^S,mom(1))*radius(ixO^S)**2 &
    ! *unit_density*unit_velocity/M_sun*year
    ! w(ixO^S,nw+3) = OPAL(ixO^S)/kappa_e
    ! w(ixO^S,nw+4) = CAK(ixO^S)/kappa_e
    ! w(ixO^S,nw+5) = lambda(ixO^S)
    ! w(ixO^S,nw+6) = big_gamma(ixO^S)
    ! w(ixO^S,nw+7) = Lum(ixO^S)

    w(ixO^S,nw+1) = OPAL(ixO^S)/kappa_e
    w(ixO^S,nw+2) = CAK(ixO^S)/kappa_e
    w(ixO^S,nw+3) = big_gamma(ixO^S)
    w(ixO^S,nw+4) = Lum(ixO^S)


  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'OPAL CAK Gamma Lum'
    ! varnames = 'Trad Mdot OPAL CAK lambda Gamma Lum'
  end subroutine specialvarnames_output

end module mod_usr
