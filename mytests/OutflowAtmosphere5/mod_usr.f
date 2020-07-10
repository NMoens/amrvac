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

  double precision :: rho_bound, T_b0, Gamma_0, Gamma_b, L_star, Mdot, v_inf,&
      R_star, M_star

  integer :: int_r, int_v, int_e, int_re, int_dt


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
    unit_temperature=unit_pressure/((2.d0+&
       3.d0*He_abundance)*unit_numberdensity*const_kB)
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

    StefBoltz = const_rad_a*const_c/4.d0*(unit_temperature**&
       4.d0)/(unit_velocity*unit_pressure)
  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ rho_bound, T_b0, Gamma_0, Gamma_b, L_star, Mdot,&
        v_inf, R_star, M_star, error_b

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
  subroutine initial_conditions(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x)

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision :: kappa(ixOmin1:ixOmax1), lambda(ixOmin1:ixOmax1),&
        fld_R(ixOmin1:ixOmax1)
    double precision :: pert(ixOmin1:ixOmax1), vel(ixImin1:ixImax1)

    double precision :: bb, prim(ixImin1:ixImax1), Er_b0

    vel(ixImin1:ixImax1) = v_inf*(1-R_star/x(ixImin1:ixImax1,1))**0.5
    w(ixImin1:ixImax1,rho_) = Mdot/(4*dpi*x(ixImin1:ixImax1,&
       1)**2*vel(ixImin1:ixImax1))
    w(ixImin1:ixImax1,mom(1)) =  w(ixImin1:ixImax1,rho_)*vel(ixImin1:ixImax1)

    Er_b0 = const_rad_a*(T_b0*unit_temperature)**4/unit_pressure
    prim(ixImin1:ixImax1) = dsqrt(x(ixImin1:ixImax1,&
       1)-R_star) *(16*x(ixImin1:ixImax1,1)**2+8*x(ixImin1:ixImax1,&
       1)*R_star+6*R_star**2) /(15*R_star**3*x(ixImin1:ixImax1,&
       1)**(5.d0/2.d0))
    bb = 3*kappa_b*L_star*unit_velocity/(4*dpi*const_c) *Mdot/(4*dpi*v_inf)

    w(ixImin1:ixImax1,r_e) = Er_b0 - bb*prim(ixImin1:ixImax1)

    call fld_get_opacity(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, lambda,&
        fld_R)

    w(ixOmin1:ixOmax1,i_diff_mg) = (const_c/unit_velocity)*lambda(&
       ixOmin1:ixOmax1)/(kappa(ixOmin1:ixOmax1)*w(ixOmin1:ixOmax1,rho_))

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixImin1,ixImax1,ixBmin1,ixBmax1,iB,w,x)
    use mod_global_parameters
    use mod_fld


    integer, intent(in)             :: ixImin1,ixImax1, ixBmin1,ixBmax1, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision :: kappa(ixImin1:ixImax1), gradE_l(ixBmin1:ixBmax1),&
        L_vE_l(ixBmin1:ixBmax1)
    double precision :: Temp(ixImin1:ixImax1), pth(ixImin1:ixImax1),&
       pert(ixBmin1:ixBmax1)
    double precision :: tau_out, T_out, F_bound
    integer :: i,j

    select case (iB)

    case(1)

       w(nghostcells,rho_) = dinflo
      

      do i = ixBmax1-1,ixBmin1,-1
        ! w(ixBmin1:ixBmax1,i,rho_) = 2*w(ixBmin1:ixBmax1,i+1,rho_) - w(ixBmin1:ixBmax1,i+2,rho_)
         w(i,rho_) = dexp(2*dlog(w(i+1,rho_)) - dlog(w(i+2,rho_)))
        
      enddo

      w(ixBmin1:ixBmax1,mom(:)) = 0.d0

      do i = ixBmax1,ixBmin1,-1
         w(i,mom(1)) = w(i+1,mom(1))*(x(i+1,1)/x(i,1))**2
        
      enddo

      call fld_get_opacity(w, x, ixImin1,ixImax1, ixImin1,ixImax1, kappa)

      L_vE_l(ixBmin1:ixBmax1) = 4*dpi*x(ixBmin1:ixBmax1,&
         1)**2*4.d0/3.d0*w(ixBmin1:ixBmax1,mom(1))/w(ixBmin1:ixBmax1,&
         rho_)*w(ixBmin1:ixBmax1,r_e)
      gradE_l(ixBmin1:ixBmax1) = -w(ixBmin1:ixBmax1,&
         rho_)*kappa(ixBmin1:ixBmax1)*(3.d0*unit_velocity/const_c)*(L_star-&
         L_vE_l(ixBmin1:ixBmax1))/(4.d0*dpi*x(ixBmin1:ixBmax1,1)**2.d0)

       gradE = gradE_l(nghostcells)
      

      do i = ixBmax1-1,ixBmin1,-1
         w(i,r_e) = w(i+2,r_e) + (x(i,1)-x(i+2,1))*gradE_l(i+1)
        
      enddo

       w(nghostcells,r_e) = dexp(half*(dlog(w(nghostcells-1,&
          r_e))+dlog(w(nghostcells+1,r_e))))
      

    case(2)
      do i = ixBmin1,ixBmax1
        !> Conserve gradE/rho
         w(i,r_e) = (x(i-1,1)**2*w(i-1,mom(1))/w(i-1,rho_))/(x(i,1)**2*w(i,&
            mom(1))/w(i,rho_))*w(i-1,r_e)
        
      enddo

      ! gradE_out = sum((w(ixBmin1:ixBmax1,ixBmin2-1,r_e)-w(ixBmin1:ixBmax1,ixBmin2-2,r_e))&
      ! /(x(ixBmin1:ixBmax1,ixBmin2-1,2) - x(ixBmin1:ixBmax1,ixBmin2-2,2)))/(ixBmax1-ixBmin1)

      
      gradE_out = (w(ixBmin1,r_e)-w(ixBmin1-1,r_e))/(x(ixBmin1,&
         1) - x(ixBmin1-1,1))
     
      

      F_bound = L_star/(4*dpi*R_star**2)
      tau_out = kappa_0*w(ixImax1-nghostcells,&
         rho_)*R_star**2/(3*x(ixImax1-nghostcells,1))
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

    gravity_field(ixImin1:ixImax1,:) = zero
    gravity_field(ixImin1:ixImax1,1) = -const_G*mass/radius(ixImin1:ixImax1)**&
       2*(unit_time**2/unit_length)

  end subroutine set_gravitation_field

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine PseudoPlanar(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,&
     wCT,qt,w,x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,&
       iwmax
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision :: ppsource(ixOmin1:ixOmax1,1:nw)

    call PseudoPlanarSource(ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,x,ppsource,&
       .false.)
    w(ixOmin1:ixOmax1,rho_) = w(ixOmin1:ixOmax1,&
       rho_) + qdt*ppsource(ixOmin1:ixOmax1,rho_)
    w(ixOmin1:ixOmax1,mom(:)) = w(ixOmin1:ixOmax1,&
       mom(:)) + qdt*ppsource(ixOmin1:ixOmax1,mom(:))
    w(ixOmin1:ixOmax1,e_) = w(ixOmin1:ixOmax1,&
       e_) + qdt*ppsource(ixOmin1:ixOmax1,e_)
    w(ixOmin1:ixOmax1,r_e) = w(ixOmin1:ixOmax1,&
       r_e) + qdt*ppsource(ixOmin1:ixOmax1,r_e)

  end subroutine PseudoPlanar

  !> internal boundary, user defined
  !> This subroutine can be used to artificially overwrite ALL conservative
  !> variables in a user-selected region of the mesh, and thereby act as
  !> an internal boundary region. It is called just before external (ghost cell)
  !> boundary regions will be set by the BC selection. Here, you could e.g.
  !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
  !> which can be used to identify the internal boundary region location.
  !> Its effect should always be local as it acts on the mesh.
  subroutine reset_egas(level,qt,ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)

    double precision :: Temp(ixImin1:ixImax1),pth(ixImin1:ixImax1)

    if (rhd_energy_interact) call mpistop(&
       'Resetting e_gas but e_interact true')

    Temp(ixImin1:ixImax1) = (w(ixImin1:ixImax1,&
       r_e)*unit_pressure/const_rad_a)**0.25d0/unit_temperature
    pth(ixImin1:ixImax1) = Temp(ixImin1:ixImax1)*w(ixImin1:ixImax1,rho_)
    w(ixImin1:ixImax1,e_) = pth(ixImin1:ixImax1)/(rhd_gamma-1.d0) + &
       half*sum(w(ixImin1:ixImax1, mom(:))**2, dim=ndim+1)/w(ixImin1:ixImax1,&
        rho_)

  end subroutine reset_egas


  subroutine PseudoPlanarSource(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,source,&
     boundary)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(out) :: source(ixOmin1:ixOmax1,1:nw)
    logical, intent(in) :: boundary

    double precision :: rad_flux(ixOmin1:ixOmax1,1:ndir)
    double precision :: pth(ixImin1:ixImax1),v(ixOmin1:ixOmax1,2)
    double precision :: radius(ixOmin1:ixOmax1),  pert(ixOmin1:ixOmax1)
    integer :: rdir, pdir

    source(ixOmin1:ixOmax1,1:nw) = zero

    rdir = 1
    pdir = 2

    v(ixOmin1:ixOmax1,rdir) = w(ixOmin1:ixOmax1,mom(rdir))/w(ixOmin1:ixOmax1,&
       rho_)
    

    radius(ixOmin1:ixOmax1) = x(ixOmin1:ixOmax1,rdir) !+ half*block%dx(ixOmin1:ixOmax1,rdir)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    source(ixOmin1:ixOmax1,rho_) = -two*w(ixOmin1:ixOmax1,&
       rho_)*v(ixOmin1:ixOmax1,rdir)/radius(ixOmin1:ixOmax1)

    ! pth(ixO^S) = (rhd_gamma-1) &
    ! *(w(ixO^S,e_) - half*sum(w(ixO^S, mom(:))**2, dim=ndim+1) / w(ixO^S, rho_))

    call phys_get_pthermal(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,pth)

    !> dm_r/dt = +(rho*v_p**2 + 2pth)/r -2 (rho*v_r**2 + pth)/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    source(ixOmin1:ixOmax1,mom(rdir)) = - 2*w(ixOmin1:ixOmax1,&
       rho_)*v(ixOmin1:ixOmax1,rdir)**two/radius(ixOmin1:ixOmax1) 

    

    ! !> de/dt = -2 (e+p)v_r/r
    ! source(ixO^S,e_) = - two*(w(ixO^S,e_)+pth(ixO^S))*v(ixO^S,rdir)/radius(ixO^S)

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
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
      source(ixOmin1:ixOmax1,r_e) = source(ixOmin1:ixOmax1,&
         r_e) + two*v(ixOmin1:ixOmax1,rdir)*w(ixOmin1:ixOmax1,&
         r_e)/(3*radius(ixOmin1:ixOmax1))
    endif

  end subroutine PseudoPlanarSource


  subroutine Opacity_stepfunction(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,kappa)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1,1:nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1)

    kappa(ixOmin1:ixOmax1) = kappa_b + (1.d0+erf((x(ixOmin1:ixOmax1,&
       1)-one)*error_b-error_b/2.d0))*(kappa_0-kappa_b)/2.d0

  end subroutine Opacity_stepfunction

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

    !> Refine close to base
    refine = -1

    if (qt .gt. 1.d-1) then
      print*, 'refine 1'
      if (any(x(ixGmin1:ixGmax1,1) < 1.d0)) refine=1
    endif

    if (qt .gt. 2.d-1) then
      print*, 'refine 2'
      if (any(x(ixGmin1:ixGmax1,1) < 2.d0)) refine=1
    endif

    if (qt .gt. 4.d-1) then
      print*, 'refine 3'
      if (any(x(ixGmin1:ixGmax1,1) < 3.d0)) refine=1
    endif

  end subroutine refine_base


  subroutine time_average_values(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    if (global_time .gt. 0.5d0) then
      w(ixImin1:ixImax1,int_r) = w(ixImin1:ixImax1,int_r) + w(ixImin1:ixImax1,&
         rho_)*dt
      w(ixImin1:ixImax1,int_v) =  w(ixImin1:ixImax1,int_v) + w(ixImin1:ixImax1,&
         mom(1))/w(ixImin1:ixImax1,rho_)*dt
      w(ixImin1:ixImax1,int_e) = w(ixImin1:ixImax1,int_e) + w(ixImin1:ixImax1,&
         e_)*dt
      w(ixImin1:ixImax1,int_re) = w(ixImin1:ixImax1,&
         int_re) + w(ixImin1:ixImax1,r_e)*dt

      w(ixImin1:ixImax1,int_dt) =  w(ixImin1:ixImax1,int_dt) + dt
    else
      w(ixImin1:ixImax1,int_r) = zero
      w(ixImin1:ixImax1,int_v) = zero
      w(ixImin1:ixImax1,int_e) = zero
      w(ixImin1:ixImax1,int_re) = zero
      w(ixImin1:ixImax1,int_dt) = zero
    endif

  end subroutine time_average_values

  subroutine specialvar_output(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_physics
    use mod_fld

    integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
    double precision                   :: w(ixImin1:ixImax1,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: g_rad(ixImin1:ixImax1),&
        big_gamma(ixImin1:ixImax1)
    double precision                   :: g_grav(ixImin1:ixImax1)
    double precision                   :: Tgas(ixImin1:ixImax1),&
       Trad(ixImin1:ixImax1)
    double precision                   :: kappa(ixOmin1:ixOmax1),&
        OPAL(ixOmin1:ixOmax1), CAK(ixOmin1:ixOmax1)
    double precision                   :: vel(ixImin1:ixImax1),&
        gradv(ixOmin1:ixOmax1)
    double precision                   :: rad_flux(ixOmin1:ixOmax1,1:ndim),&
        Lum_cmf(ixOmin1:ixOmax1)
    double precision                   :: Lum_adv(ixOmin1:ixOmax1),&
        Lum_obs(ixOmin1:ixOmax1)
    double precision                   :: Lum_wnd(ixOmin1:ixOmax1),&
        Mdot(ixOmin1:ixOmax1)
    double precision                   :: pp_rf(ixOmin1:ixOmax1),&
        lambda(ixOmin1:ixOmax1), fld_R(ixOmin1:ixOmax1)
    integer                            :: idim
    double precision :: radius(ixImin1:ixImax1)
    double precision :: mass

    radius(ixOmin1:ixOmax1) = x(ixOmin1:ixOmax1,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    call fld_get_opacity(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, kappa)
    call fld_get_radflux(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, rad_flux)

    g_rad(ixOmin1:ixOmax1) = kappa(ixOmin1:ixOmax1)*rad_flux(ixOmin1:ixOmax1,&
       1)/(const_c/unit_velocity)
    g_grav(ixOmin1:ixOmax1) = const_G*mass/radius(ixOmin1:ixOmax1)**&
       2*(unit_time**2/unit_length)
    big_gamma(ixOmin1:ixOmax1) = g_rad(ixOmin1:ixOmax1)/g_grav(&
       ixOmin1:ixOmax1)

    Mdot(ixOmin1:ixOmax1) =  4*dpi*w(ixOmin1:ixOmax1,&
       mom(1))*radius(ixOmin1:ixOmax1)**2 &
       *unit_density*unit_velocity/M_sun*year

    call rhd_get_tgas(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, Tgas)
    call rhd_get_trad(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, Trad)

    vel(ixImin1:ixImax1) = w(ixImin1:ixImax1,mom(1))/w(ixImin1:ixImax1,rho_)

    pp_rf(ixOmin1:ixOmax1) = two*rad_flux(ixOmin1:ixOmax1,1)/x(ixOmin1:ixOmax1,&
       1)*dt

    call fld_get_fluxlimiter(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, lambda,&
        fld_R)

    Lum_cmf(ixOmin1:ixOmax1) = 4*dpi*rad_flux(ixOmin1:ixOmax1,&
       1)*(x(ixOmin1:ixOmax1,1)*unit_length)**2*unit_radflux/L_sun
    Lum_adv(ixOmin1:ixOmax1) = 4*dpi*4.d0/3.d0*vel(ixOmin1:ixOmax1)*w(&
       ixOmin1:ixOmax1,r_e)*(x(ixOmin1:ixOmax1,&
       1)*unit_length)**2*unit_radflux/L_sun
    Lum_obs(ixOmin1:ixOmax1) = Lum_cmf(ixOmin1:ixOmax1) + &
       Lum_adv(ixOmin1:ixOmax1)
    Lum_wnd(ixOmin1:ixOmax1) =  4*dpi*w(ixOmin1:ixOmax1,&
       mom(1))*radius(ixOmin1:ixOmax1)**2*unit_density*unit_velocity/L_sun*((&
       vel(ixOmin1:ixOmax1)*unit_velocity)**2/2.d0 - &
       const_G*mass/radius(ixOmin1:ixOmax1) + const_G*mass/unit_length)

    w(ixOmin1:ixOmax1,nw+1) = kappa(ixOmin1:ixOmax1)*unit_opacity
    w(ixOmin1:ixOmax1,nw+2) = rad_flux(ixOmin1:ixOmax1,1)
    w(ixOmin1:ixOmax1,nw+3) = Trad(ixOmin1:ixOmax1)*unit_temperature
    w(ixOmin1:ixOmax1,nw+4) = big_gamma(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+5) = Mdot(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+6) = vel(ixOmin1:ixOmax1)*unit_velocity
    w(ixOmin1:ixOmax1,nw+7) = pp_rf(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+8) = Lum_cmf(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+9) = Lum_adv(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+10) = Lum_obs(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+11) = Lum_wnd(ixOmin1:ixOmax1)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

               !9     10 11   12    13   14  15    16     17    18    19    20
    varnames = 'kappa F1 Trad Gamma Mdot vel pp_rf L_cmf L_adv L_obs L_wnd'
  end subroutine specialvarnames_output

end module mod_usr
