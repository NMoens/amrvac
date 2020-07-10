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


  double precision :: dinflo,gradE,gradE_out, valE_out
  double precision :: error_b, kappa_0, kappa_b

  double precision :: rho_bound, T_b0, Gamma_0, Gamma_b, L_star, Mdot, v_inf, R_star, M_star

  integer :: int_r, int_v, int_e, int_re, int_dt


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
    usr_internal_bc => reset_egas

    ! Graviatational field
    usr_gravity => set_gravitation_field

    ! Special Opacity
    usr_special_opacity => Opacity_stepfunction

    ! Analysis
    usr_modify_output => time_average_values

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    usr_refine_grid => refine_base

    ! Active the physics module
    call rhd_activate()

    int_r = var_set_extravar("int_r", "int_r")
    int_v = var_set_extravar("int_v", "int_v")
    int_e = var_set_extravar("int_e", "int_e")
    int_re = var_set_extravar("int_re", "int_re")
    int_dt = var_set_extravar("int_dt", "int_dt")

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    integer :: i

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
    Mdot = Mdot/(unit_density*unit_length**3.d0)*unit_time
    L_star = L_star/(unit_pressure*unit_length**3.d0)*unit_time
    kappa_0 = kappa_0/unit_opacity
    kappa_b = kappa_b/unit_opacity
    rho_bound = rho_bound/unit_density
    v_inf = v_inf/unit_velocity
    T_b0 = T_b0/unit_temperature

    StefBoltz = const_rad_a*const_c/4.d0*(unit_temperature**4.d0)/(unit_velocity*unit_pressure)
  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ rho_bound, T_b0, Gamma_0, Gamma_b, L_star, Mdot, v_inf, R_star, M_star, error_b

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wind_list, end=113)
113    close(unitpar)
    end do

    R_star = R_star*R_sun
    M_star = M_star*M_sun
    L_star = L_star*L_sun
    Mdot = Mdot*M_sun/year

    kappa_0 = Gamma_0*4*dpi*const_G*M_star*const_c/L_star
    kappa_b = Gamma_b*4*dpi*const_G*M_star*const_c/L_star

  end subroutine params_read


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)
    double precision :: pert(ixO^S), vel(ixI^S)

    double precision :: bb, prim(ixI^S), Er_b0

    vel(ixI^S) = v_inf*(1-R_star/x(ixI^S,1))**0.5
    w(ixI^S,rho_) = Mdot/(4*dpi*x(ixI^S,1)**2*vel(ixI^S))
    w(ixI^S,mom(1)) =  w(ixI^S,rho_)*vel(ixI^S)

    Er_b0 = const_rad_a*(T_b0*unit_temperature)**4/unit_pressure
    prim(ixI^S) = dsqrt(x(ixI^S,1)-R_star) &
    *(16*x(ixI^S,1)**2+8*x(ixI^S,1)*R_star+6*R_star**2) &
    /(15*R_star**3*x(ixI^S,1)**(5.d0/2.d0))
    bb = 3*kappa_b*L_star*unit_velocity/(4*dpi*const_c) &
         *Mdot/(4*dpi*v_inf)

    w(ixI^S,r_e) = Er_b0 - bb*prim(ixI^S)

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

    double precision :: kappa(ixI^S), gradE_l(ixB^S), L_vE_l(ixB^S)
    double precision :: Temp(ixI^S), pth(ixI^S),pert(ixB^S)
    double precision :: tau_out, T_out, F_bound
    integer :: i,j

    select case (iB)

    case(1)

      {^IFONED w(nghostcells,rho_) = dinflo}
      {^IFTWOD w(nghostcells,:,rho_) = dinflo}

      do i = ixBmax1-1,ixBmin1,-1
        ! w(ixBmin1:ixBmax1,i,rho_) = 2*w(ixBmin1:ixBmax1,i+1,rho_) - w(ixBmin1:ixBmax1,i+2,rho_)
        {^IFONED w(i,rho_) = dexp(2*dlog(w(i+1,rho_)) - dlog(w(i+2,rho_))) }
        {^IFTWOD w(i,:,rho_) = dexp(2*dlog(w(i+1,:,rho_)) - dlog(w(i+2,:,rho_))) }
      enddo

      w(ixB^S,mom(:)) = 0.d0

      do i = ixBmax1,ixBmin1,-1
        {^IFONED w(i,mom(1)) = w(i+1,mom(1))*(x(i+1,1)/x(i,1))**2}
        {^IFTWOD w(i,:,mom(1)) = w(i+1,:,mom(1))*(x(i+1,:,1)/x(i,:,1))**2}
      enddo

      call fld_get_opacity(w, x, ixI^L, ixI^L, kappa)

      L_vE_l(ixB^S) = 4*dpi*x(ixB^S,1)**2*4.d0/3.d0*w(ixB^S,mom(1))/w(ixB^S,rho_)*w(ixB^S,r_e)
      gradE_l(ixB^S) = -w(ixB^S,rho_)*kappa(ixB^S)*(3.d0*unit_velocity/const_c)*(L_star-L_vE_l(ixB^S))/(4.d0*dpi*x(ixB^S,1)**2.d0)

      {^IFONED gradE = gradE_l(nghostcells)}
      {^IFTWOD gradE = sum(gradE_l(nghostcells,ixBmin2:ixBmax2))/(ixBmax2-ixBmin2)}

      do i = ixBmax1-1,ixBmin1,-1
        {^IFONED w(i,r_e) = w(i+2,r_e) + (x(i,1)-x(i+2,1))*gradE_l(i+1)}
        {^IFTWOD w(i,:,r_e) = w(i+2,:,r_e) + (x(i,:,1)-x(i+2,:,1))*gradE_l(i+1,:)}
      enddo

      {^IFONED w(nghostcells,r_e) = dexp(half*(dlog(w(nghostcells-1,r_e))+dlog(w(nghostcells+1,r_e))))}
      {^IFTWOD w(nghostcells,:,r_e) = dexp(half*(dlog(w(nghostcells-1,:,r_e))+dlog(w(nghostcells+1,:,r_e))))}

    case(2)
      do i = ixBmin1,ixBmax1
        !> Conserve gradE/rho
        {^IFONED w(i,r_e) = (x(i-1,1)**2*w(i-1,mom(1))/w(i-1,rho_))/(x(i,1)**2*w(i,mom(1))/w(i,rho_))*w(i-1,r_e)}
        {^IFTWOD w(i,:,r_e) = (x(i-1,:,1)**2*w(i-1,:,mom(1))/w(i-1,:,rho_))/(x(i,:,1)**2*w(i,:,mom(1))/w(i,:,rho_))*w(i-1,:,r_e)}
      enddo

      ! gradE_out = sum((w(ixBmin1:ixBmax1,ixBmin2-1,r_e)-w(ixBmin1:ixBmax1,ixBmin2-2,r_e))&
      ! /(x(ixBmin1:ixBmax1,ixBmin2-1,2) - x(ixBmin1:ixBmax1,ixBmin2-2,2)))/(ixBmax1-ixBmin1)

      {^IFONED
      gradE_out = (w(ixBmin1,r_e)-w(ixBmin1-1,r_e))/(x(ixBmin1,1) - x(ixBmin1-1,1))
      }
      {^IFTWOD
      gradE_out = sum((w(ixBmin1,ixBmin2:ixBmax2,r_e)-w(ixBmin1-1,ixBmin2:ixBmax2,r_e))&
      /(x(ixBmin2,ixBmin2:ixBmax2,1) - x(ixBmin1-1,ixBmin2:ixBmax2,1)))/(ixBmax2-ixBmin2)
      }

      F_bound = L_star/(4*dpi*R_star**2)
      tau_out = kappa_0*w(ixImax1-nghostcells,rho_)*R_star**2/(3*x(ixImax1-nghostcells,1))
      T_out = F_bound/StefBoltz*(3.d0/4.d0*tau_out)**0.25d0

      T_out = max(T_out,1.d4/unit_temperature)

      print*, w(ixImax1-nghostcells,rho_), tau_out, T_out*unit_temperature

      valE_out = const_rad_a*(T_out*unit_temperature)**4/unit_pressure

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
    w(ixO^S,e_) = w(ixO^S,e_) + qdt*ppsource(ixO^S,e_)
    w(ixO^S,r_e) = w(ixO^S,r_e) + qdt*ppsource(ixO^S,r_e)

  end subroutine PseudoPlanar

  !> internal boundary, user defined
  !> This subroutine can be used to artificially overwrite ALL conservative
  !> variables in a user-selected region of the mesh, and thereby act as
  !> an internal boundary region. It is called just before external (ghost cell)
  !> boundary regions will be set by the BC selection. Here, you could e.g.
  !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
  !> which can be used to identify the internal boundary region location.
  !> Its effect should always be local as it acts on the mesh.
  subroutine reset_egas(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision :: Temp(ixI^S),pth(ixI^S)

    if (rhd_energy_interact) call mpistop('Resetting e_gas but e_interact true')

    Temp(ixI^S) = (w(ixI^S,r_e)*unit_pressure/const_rad_a)**0.25d0/unit_temperature
    pth(ixI^S) = Temp(ixI^S)*w(ixI^S,rho_)
    w(ixI^S,e_) = pth(ixI^S)/(rhd_gamma-1.d0) &
    + half*sum(w(ixI^S, mom(:))**2, dim=ndim+1)/w(ixI^S, rho_)

  end subroutine reset_egas


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

    ! !> de/dt = -2 (e+p)v_r/r
    ! source(ixO^S,e_) = - two*(w(ixO^S,e_)+pth(ixO^S))*v(ixO^S,rdir)/radius(ixO^S)

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
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


  subroutine Opacity_stepfunction(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    kappa(ixO^S) = kappa_b + (1.d0+erf((x(ixO^S,1)-one)*error_b-error_b/2.d0))*(kappa_0-kappa_b)/2.d0

  end subroutine Opacity_stepfunction

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


  subroutine time_average_values(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    if (global_time .gt. 0.5d0) then
      w(ixI^S,int_r) = w(ixI^S,int_r) &
      + w(ixI^S,rho_)*dt
      w(ixI^S,int_v) =  w(ixI^S,int_v) &
      + w(ixI^S,mom(1))/w(ixI^S,rho_)*dt
      w(ixI^S,int_e) = w(ixI^S,int_e) &
      + w(ixI^S,e_)*dt
      w(ixI^S,int_re) = w(ixI^S,int_re) &
      + w(ixI^S,r_e)*dt

      w(ixI^S,int_dt) =  w(ixI^S,int_dt) + dt
    else
      w(ixI^S,int_r) = zero
      w(ixI^S,int_v) = zero
      w(ixI^S,int_e) = zero
      w(ixI^S,int_re) = zero
      w(ixI^S,int_dt) = zero
    endif

  end subroutine time_average_values

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
    double precision                   :: vel(ixI^S), gradv(ixO^S)
    double precision                   :: rad_flux(ixO^S,1:ndim), Lum_cmf(ixO^S)
    double precision                   :: Lum_adv(ixO^S), Lum_obs(ixO^S)
    double precision                   :: Lum_wnd(ixO^S), Mdot(ixO^S)
    double precision                   :: pp_rf(ixO^S), lambda(ixO^S), fld_R(ixO^S)
    integer                            :: idim
    double precision :: radius(ixI^S)
    double precision :: mass

    radius(ixO^S) = x(ixO^S,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)

    g_rad(ixO^S) = kappa(ixO^S)*rad_flux(ixO^S,1)/(const_c/unit_velocity)
    g_grav(ixO^S) = const_G*mass/radius(ixO^S)**2*(unit_time**2/unit_length)
    big_gamma(ixO^S) = g_rad(ixO^S)/g_grav(ixO^S)

    Mdot(ixO^S) =  4*dpi*w(ixO^S,mom(1))*radius(ixO^S)**2 &
    *unit_density*unit_velocity/M_sun*year

    call rhd_get_tgas(w, x, ixI^L, ixO^L, Tgas)
    call rhd_get_trad(w, x, ixI^L, ixO^L, Trad)

    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)

    pp_rf(ixO^S) = two*rad_flux(ixO^S,1)/x(ixO^S,1)*dt

    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    Lum_cmf(ixO^S) = 4*dpi*rad_flux(ixO^S,1)*(x(ixO^S,1)*unit_length)**2*unit_radflux/L_sun
    Lum_adv(ixO^S) = 4*dpi*4.d0/3.d0*vel(ixO^S)*w(ixO^S,r_e)*(x(ixO^S,1)*unit_length)**2*unit_radflux/L_sun
    Lum_obs(ixO^S) = Lum_cmf(ixO^S) + Lum_adv(ixO^S)
    Lum_wnd(ixO^S) =  4*dpi*w(ixO^S,mom(1))*radius(ixO^S)**2*unit_density*unit_velocity/L_sun&
                            *((vel(ixO^S)*unit_velocity)**2/2.d0 &
                            - const_G*mass/radius(ixO^S) &
                            + const_G*mass/unit_length)

    w(ixO^S,nw+1) = kappa(ixO^S)*unit_opacity
    w(ixO^S,nw+2) = rad_flux(ixO^S,1)
    w(ixO^S,nw+3) = Trad(ixO^S)*unit_temperature
    w(ixO^S,nw+4) = big_gamma(ixO^S)
    w(ixO^S,nw+5) = Mdot(ixO^S)
    w(ixO^S,nw+6) = vel(ixO^S)*unit_velocity
    w(ixO^S,nw+7) = pp_rf(ixO^S)
    w(ixO^S,nw+8) = Lum_cmf(ixO^S)
    w(ixO^S,nw+9) = Lum_adv(ixO^S)
    w(ixO^S,nw+10) = Lum_obs(ixO^S)
    w(ixO^S,nw+11) = Lum_wnd(ixO^S)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

               !9     10 11   12    13   14  15    16     17    18    19    20
    varnames = 'kappa F1 Trad Gamma Mdot vel pp_rf L_cmf L_adv L_obs L_wnd'
  end subroutine specialvarnames_output

end module mod_usr
