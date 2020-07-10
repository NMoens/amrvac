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

  double precision :: rho_bound, Gamma, T_bound, E_dot, L_bound, M_dot, R_star, M_star, v_inf, start_Edot
  double precision :: gradE, valE_out, F_bound, StefBoltz, E_out, T_out

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    {^IFONED call set_coordinate_system("Cartesian_1D")}
    {^IFTWOD call set_coordinate_system("Cartesian_2D")}

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

    ! Output routines
    ! usr_aux_output    => specialvar_output
    ! usr_add_aux_names => specialvarnames_output

    usr_refine_grid => refine_base

    ! Active the physics module
    call rhd_activate

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

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


    !> Make all parameters dimensionless
    M_star = M_star/(unit_density*unit_length**3.d0)
    R_star = R_star/unit_length
    M_dot = M_dot/(unit_density*unit_length**3.d0)*unit_time
    L_bound = L_bound/(unit_pressure*unit_length**3.d0)*unit_time
    F_bound = F_bound/unit_radflux
    rho_bound = rho_bound/unit_density
    v_inf = v_inf/unit_velocity
    T_bound = T_bound/unit_temperature

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ rho_bound, Gamma, T_bound, E_dot, L_bound, M_dot, R_star, M_star, v_inf, start_Edot

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wind_list, end=113)
113    close(unitpar)
    end do

    M_star = M_star*M_sun
    R_star = R_star*R_sun
    M_dot = M_dot*M_sun/year
    L_bound = L_bound*L_sun
    F_bound = L_bound/(4*dpi*R_star**2)


    start_Edot = start_Edot

  end subroutine params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: lambda(ixO^S), fld_R(ixO^S), vel(ixI^S)

    double precision :: rr(ixI^S), bb
    double precision :: E_in, E_gauge

    w(ixI^S,rho_)   = rho_bound
    w(ixI^S,mom(:)) = 0.d0

    where(x(ixI^S,1) .gt. 1.d0)
      vel(ixI^S) =  v_inf*(1.d0 - 0.99d0*(1.d0/x(ixI^S,1)))**0.5d0
      w(ixI^S,rho_) = M_dot/(4.d0*dpi*vel(ixI^S)*x(ixI^S,1)**2)
      w(ixI^S,mom(1)) = vel(ixI^S)*w(ixI^S,rho_)
    endwhere


    !> Outer/Inner temperature
    T_out = 28445.d0/unit_temperature
    E_out = const_rad_a*(T_out*unit_temperature)**4/unit_pressure
    E_in = const_rad_a*(T_bound*unit_temperature)**4/unit_pressure

    !>Very bad initial profile using constant gradE
    ! w(ixI^S,r_e) = E_out + gradE*(x(ixI^S,1)-xprobmax1)
    ! w(ixI^S,r_e) = E_in + gradE*(x(ixI^S,1)-xprobmin1)
    rr(ixI^S) = dsqrt(x(ixI^S,1)-xprobmin1)*(16*x(ixI^S,1)**2 + 8*x(ixI^S,1)+ 6) &
                /(15*x(ixI^S,1)**(5.d0/2))
    bb = -fld_kappa0*L_bound*M_dot*unit_velocity*3/(16*dpi**2*const_c*v_inf)

    ! bb = bb*2

    ! rr(ixI^S) = 2*dsqrt(x(ixI^S,1)-xprobmin1)/dsqrt(x(ixI^S,1))
    ! bb = kappa_e*F_bound*Mdot*unit_velocity*3/(4*dpi*const_c*v_inf)
    w(ixI^S,r_e) = E_in + bb*rr(ixI^S)

    E_gauge = E_in + bb*dsqrt(xprobmax1-xprobmin1)&
              *(16*xprobmax1**2 + 8*xprobmax1+ 6)/(15*xprobmax1**(5.d0/2))

    w(ixI^S,r_e) = w(ixI^S,r_e) - E_gauge + E_out

    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(fld_kappa0*w(ixO^S,rho_))

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_physics, only: phys_get_trad
    use mod_global_parameters
    use mod_opacity
    use mod_fld

    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Temp(ixI^S)
    double precision :: Local_gradE(ixI^S)
    double precision :: Local_tauout(ixB^S)
    double precision :: Local_Tout(ixB^S)

    integer :: ix^D

    select case (iB)

    case(1)
      w(ixB^S,rho_) = rho_bound
      do ix1 = ixBmax1-1,ixBmin1,-1
        w(ix1,rho_) = dexp(2*dlog(w(ix1+1,rho_)) - dlog(w(ix1+2,rho_)))
      enddo

      ! w(ixB^S,mom(:)) = 0.d0

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,mom(1)) = w(ix1+1,mom(1))
      enddo

      ! where(w(ixB^S,mom(1)) .lt. 0.d0)
      !   w(ixB^S,mom(1)) = 0.d0
      ! endwhere
      !
      ! where(w(ixB^S,mom(1))/w(ixB^S,rho_) .gt. 0.5d0)
      !   w(ixB^S,mom(1)) = 0.1d0*w(ixB^S,rho_)
      ! endwhere

      Local_gradE(ixI^S) = -F_bound*3*fld_kappa0*w(ixI^S,rho_)&
      *unit_velocity/const_c
      gradE = Local_gradE(ixBmin1)

      ! print*, gradE

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,r_e) = w(ix1+2,r_e) &
        + (x(ix1,1)-x(ix1+2,1))*Local_gradE(ix1+1)
      enddo

    case(2)
      Local_tauout(ixB^S) = fld_kappa0*w(ixB^S,rho_)*R_star**2/(3*x(ixB^S,1))
      Local_Tout(ixB^S) = F_bound/StefBoltz*(3.d0/4.d0*Local_tauout(ixB^S))**0.25d0

      !> one single NaN will kill the average and then we automatically take the floor temp
      where (Local_Tout(ixB^S) .ne. Local_Tout(ixB^S))
        Local_Tout(ixB^S) = 2.d4
      endwhere

      T_out = Local_Tout(ixBmin1)

      T_out = max(1.5d4/unit_temperature, T_out)
      E_out = const_rad_a*(T_out*unit_temperature)**4.d0/unit_pressure

      do ix1 = ixBmin1,ixBmax1
        w(ix1,r_e) = 2*w(ix1-1,r_e) - w(ix1-2,r_e)
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
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
      mg%bc(iB, mg_iphi)%bc_value = gradE

    case (2)
      ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
      ! mg%bc(iB, mg_iphi)%bc_value = gradE_out

      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = valE_out

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

    gravity_field(ixI^S,:) = zero
    gravity_field(ixI^S,1) = -const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length)

  end subroutine set_gravitation_field

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine PseudoPlanar(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: ppsource(ixO^S,1:nw)

    call PseudoPlanarSource(ixI^L,ixO^L,wCT,x,ppsource,.false.)
    w(ixO^S,rho_) = w(ixO^S,rho_) + qdt*ppsource(ixO^S,rho_)
    w(ixO^S,mom(:)) = w(ixO^S,mom(:)) + qdt*ppsource(ixO^S,mom(:))
    w(ixO^S,r_e) = w(ixO^S,r_e) + qdt*ppsource(ixO^S,r_e)

    if (qtC .gt. start_Edot) then
      where (x(ixI^S,1) .lt. 1.d0)
        w(ixO^S,r_e) = w(ixO^S,r_e) + qdt*E_dot
      end where
    endif

  end subroutine PseudoPlanar


  subroutine PseudoPlanarSource(ixI^L,ixO^L,w,x,source,boundary)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: source(ixO^S,1:nw)
    logical, intent(in) :: boundary

    double precision :: rad_flux(ixO^S,1:ndir)
    double precision :: pth(ixI^S),v(ixO^S,2)
    double precision :: radius(ixO^S),  pert(ixO^S)
    integer :: rdir, pdir

    source(ixO^S,1:nw) = zero

    rdir = 1
    pdir = 2

    v(ixO^S,rdir) = w(ixO^S,mom(rdir))/w(ixO^S,rho_)
    {^IFTWOD v(ixO^S,pdir) = w(ixO^S,mom(pdir))/w(ixO^S,rho_)}

    radius(ixO^S) = x(ixO^S,rdir) ! + half*block%dx(ixO^S,rdir)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    source(ixO^S,rho_) = -two*w(ixO^S,rho_)*v(ixO^S,rdir)/radius(ixO^S)

    ! pth(ixO^S) = (rhd_gamma-1) &
    ! *(w(ixO^S,e_) - half*sum(w(ixO^S, mom(:))**2, dim=ndim+1) / w(ixO^S, rho_))

    call phys_get_pthermal(w,x,ixI^L,ixO^L,pth)

    !> dm_r/dt = +(rho*v_p**2 + 2pth)/r -2 (rho*v_r**2 + pth)/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    source(ixO^S,mom(rdir)) = - 2*w(ixO^S,rho_)*v(ixO^S,rdir)**two/radius(ixO^S) {^IFTWOD + w(ixO^S,rho_)*v(ixO^S,pdir)**two/radius(ixO^S)}

    {^IFTWOD source(ixO^S,mom(pdir)) = - 3*v(ixO^S,rdir)*v(ixO^S,pdir)*w(ixO^S,rho_)/radius(ixO^S) }

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      if (boundary) then
        rad_flux(ixO^S,rdir) = L_bound/(4.d0*dpi*x(ixO^S,rdir)**2)
      else
        call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
      endif
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

  subroutine refine_base(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
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

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    !> Refine close to base
    refine = -1

    if (qt .gt. 1.d-1) then
      print*, 'refine 1'
      if (any(x(ixG^S,1) < 1.d0)) refine=1
    endif

    if (qt .gt. 2.d-1) then
      print*, 'refine 2'
      if (any(x(ixG^S,1) < 2.d0)) refine=1
    endif

    if (qt .gt. 4.d-1) then
      print*, 'refine 3'
      if (any(x(ixG^S,1) < 3.d0)) refine=1
    endif

  end subroutine refine_base

  ! subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  !   ! this subroutine can be used in convert, to add auxiliary variables to the
  !   ! converted output file, for further analysis using tecplot, paraview, ....
  !   ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !   !
  !   ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  !   ! corresponding normalization values (default value 1)
  !   use mod_global_parameters
  !   use mod_physics
  !   use mod_fld
  !
  !   integer, intent(in)                :: ixI^L,ixO^L
  !   double precision, intent(in)       :: x(ixI^S,1:ndim)
  !   double precision                   :: w(ixI^S,nw+nwauxio)
  !   double precision                   :: normconv(0:nw+nwauxio)
  !
  !   double precision                   :: g_rad(ixI^S), big_gamma(ixI^S)
  !   double precision                   :: g_grav(ixI^S)
  !   double precision                   :: Tgas(ixI^S),Trad(ixI^S)
  !   double precision                   :: kappa(ixO^S), OPAL(ixO^S), CAK(ixO^S)
  !   double precision                   :: vel(ixI^S), gradv(ixO^S)
  !   double precision                   :: rad_flux(ixO^S,1:ndim), Lum_cmf(ixO^S)
  !   double precision                   :: Lum_adv(ixO^S), Lum_obs(ixO^S)
  !   double precision                   :: Lum_wnd(ixO^S), Mdot(ixO^S)
  !   double precision                   :: pp_rf(ixO^S), lambda(ixO^S), fld_R(ixO^S)
  !   integer                            :: idim
  !   double precision :: radius(ixI^S)
  !   double precision :: mass
  !
  !   radius(ixO^S) = x(ixO^S,1)*unit_length
  !   mass = M_star*(unit_density*unit_length**3.d0)
  !
  !   call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
  !
  !   g_rad(ixO^S) = fld_kappa0*rad_flux(ixO^S,1)/(const_c/unit_velocity)
  !   g_grav(ixO^S) = const_G*mass/radius(ixO^S)**2*(unit_time**2/unit_length)
  !   big_gamma(ixO^S) = g_rad(ixO^S)/g_grav(ixO^S)
  !
  !   Mdot(ixO^S) =  4*dpi*w(ixO^S,mom(1))*radius(ixO^S)**2 &
  !   *unit_density*unit_velocity/M_sun*year
  !
  !   call rhd_get_trad(w, x, ixI^L, ixO^L, Trad)
  !
  !   vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
  !
  !   pp_rf(ixO^S) = two*rad_flux(ixO^S,1)/x(ixO^S,1)*dt
  !
  !   call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)
  !
  !   Lum_cmf(ixO^S) = 4*dpi*rad_flux(ixO^S,1)*(x(ixO^S,1)*unit_length)**2*unit_radflux/L_sun
  !   Lum_adv(ixO^S) = 4*dpi*4.d0/3.d0*vel(ixO^S)*w(ixO^S,r_e)*(x(ixO^S,1)*unit_length)**2*unit_radflux/L_sun
  !   Lum_obs(ixO^S) = Lum_cmf(ixO^S) + Lum_adv(ixO^S)
  !   Lum_wnd(ixO^S) =  4*dpi*w(ixO^S,mom(1))*radius(ixO^S)**2*unit_density*unit_velocity/L_sun&
  !                           *((vel(ixO^S)*unit_velocity)**2/2.d0 &
  !                           - const_G*mass/radius(ixO^S) &
  !                           + const_G*mass/unit_length)
  !
  !   w(ixO^S,nw+1) = kappa(ixO^S)*unit_opacity
  !   w(ixO^S,nw+2) = rad_flux(ixO^S,1)
  !   w(ixO^S,nw+3) = Trad(ixO^S)*unit_temperature
  !   w(ixO^S,nw+4) = big_gamma(ixO^S)
  !   w(ixO^S,nw+5) = Mdot(ixO^S)
  !   w(ixO^S,nw+6) = vel(ixO^S)*unit_velocity
  !   w(ixO^S,nw+7) = pp_rf(ixO^S)
  !   w(ixO^S,nw+8) = Lum_cmf(ixO^S)
  !   w(ixO^S,nw+9) = Lum_adv(ixO^S)
  !   w(ixO^S,nw+10) = Lum_obs(ixO^S)
  !   w(ixO^S,nw+11) = Lum_wnd(ixO^S)
  !
  ! end subroutine specialvar_output
  !
  ! subroutine specialvarnames_output(varnames)
  !   ! newly added variables need to be concatenated with the w_names/primnames string
  !   use mod_global_parameters
  !   character(len=*) :: varnames
  !
  !              !9     10 11   12    13   14  15    16     17    18    19    20
  !   varnames = 'kappa F1 Trad Gamma Mdot vel pp_rf L_cmf L_adv L_obs L_wnd'
  ! end subroutine specialvarnames_output

end module mod_usr
