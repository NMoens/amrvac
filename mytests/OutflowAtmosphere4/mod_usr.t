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

  double precision, allocatable :: r_arr(:)
  double precision, allocatable :: rho_arr(:)
  double precision, allocatable :: v_arr(:)
  double precision, allocatable :: e_arr(:)
  double precision, allocatable :: Er_arr(:)
  double precision, allocatable :: T_arr(:)
  double precision, allocatable :: p_arr(:)

  double precision :: M_dot_ratio
  double precision :: Gamma_0, Gamma_b
  double precision :: kappa_0, kappa_b
  double precision :: L_0,L_vE
  double precision :: M_star
  double precision :: R_star
  double precision :: M_dot
  double precision :: rho_base
  double precision :: T_base

  double precision :: dinflo,gradE,gradE_out
  double precision :: error_b

  integer :: int_r, int_v, int_e, int_re, int_dt


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
    usr_special_opacity => Opacity_stepfunction

    ! Analysis
    usr_modify_output => time_average_values

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

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

    !> Set stellar mass and radius
    call ReadInParams(M_star,R_star,Gamma_0,M_dot_ratio,M_dot,L_0,rho_base,error_b)

    !> Gamma at the base is below one!
    Gamma_b = 0.9d0
    kappa_0 = Gamma_0*4*dpi*const_G*M_star*const_c/L_0
    kappa_b = Gamma_b*4*dpi*const_G*M_star*const_c/L_0

    allocate(r_arr(domain_nx2+2*nghostcells))
    allocate(rho_arr(domain_nx2+2*nghostcells))
    allocate(v_arr(domain_nx2+2*nghostcells))
    allocate(e_arr(domain_nx2+2*nghostcells))
    allocate(Er_arr(domain_nx2+2*nghostcells))
    allocate(T_arr(domain_nx2+2*nghostcells))
    allocate(p_arr(domain_nx2+2*nghostcells))

    ! call ReadInTable(r_arr,rho_arr,v_arr,e_arr,Er_arr,T_arr, p_arr)
    call ReadInTable(r_arr,rho_arr,v_arr,e_arr,Er_arr,T_arr, p_arr)

    rho_base = rho_arr(nghostcells+1)
    T_base = T_arr(nghostcells+1)

    ! if (mype .eq. 0) then
    !   print*, 'density at base', rho_base
    !   print*, 'Temperature base', T_base
    !   print*, 'cgs opacity', kappa_0
    ! endif

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = R_star !r_arr(nghostcells) ! cm
    unit_numberdensity = rho_base/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_temperature = T_base

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=(2.d0+3.d0*He_abundance)*unit_numberdensity*const_kb*unit_temperature
    unit_velocity=dsqrt(unit_pressure/unit_density)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    ! if (mype .eq. 0) then
    !   print*, 'M_star ', 'R_star ','M_dot_ratio ', 'M_dot ', 'L_0'
    !   print*, 'Gamma_0 ', 'kappa_0 ', 'kappa_b'
    !   print*, M_star, R_star,M_dot_ratio, M_dot, L_0
    !   print*, Gamma_0, kappa_0, kappa_b
    !
    !   print*, 'Flux at boundary: ', L_0/(4*dpi*R_star**2)
    !
    !   print*, 'unit_length', unit_length
    !   print*, 'unit_density', unit_density
    !   print*, 'unit_numberdensity', unit_numberdensity
    !   print*, 'unit_pressure', unit_pressure
    !   print*, 'unit_temperature', unit_temperature
    !   print*, 'unit_radflux', unit_radflux
    !   print*, 'unit_opacity', unit_opacity
    !   print*, 'unit_time', unit_time
    !   print*, 'unit_velocity', unit_velocity
    ! endif

    !> Make all parameters dimensionless
    M_star = M_star/(unit_density*unit_length**3.d0)
    R_star = R_star/unit_length
    M_dot = M_dot/(unit_density*unit_length**3.d0)*unit_time
    L_0 = L_0/(unit_pressure*unit_length**3.d0)*unit_time
    kappa_0 = kappa_0/unit_opacity
    kappa_b = kappa_b/unit_opacity
    rho_base = rho_base/unit_density
    T_base = T_base/unit_temperature

    !> Make initial profiles dimensionless
    r_arr = r_arr/unit_length
    rho_arr = rho_arr/unit_density
    v_arr = v_arr/unit_velocity
    e_arr = e_arr/unit_pressure
    Er_arr = Er_arr/unit_pressure
    T_arr = T_arr/unit_temperature
    p_arr = p_arr/unit_pressure

    ! if (mype .eq. 0) then
    !   print*, 'M_star ', 'R_star ','M_dot_ratio ', 'M_dot ', 'L_0'
    !   print*, 'Gamma_0 ', 'kappa_0 ', 'kappa_b'
    !   print*, M_star, R_star,M_dot_ratio, M_dot, L_0
    !   print*, Gamma_0, kappa_0, kappa_b
    !
    !   print*, 'Flux at boundary: ', L_0/(4*dpi*R_star**2)
    ! endif

    L_vE = 4*dpi*R_star**2*v_arr(nghostcells+1)*4.d0/3.d0*Er_arr(nghostcells+1)
    !>Set bottom density from massloss rate

    dinflo = 20.0d0*M_dot/(4*dpi*R_star**2)
    gradE = -dinflo*kappa_b*(L_0-L_vE)/(4*dpi*R_star**2*const_c/unit_velocity)

    print*, dinflo*unit_density
  end subroutine initglobaldata_usr

  subroutine ReadInParams(M_star,R_star,Gamma_0,M_dot_ratio,M_dot,L_0,rho_base,error_b)
    use mod_global_parameters
    double precision, intent(out) :: M_star,R_star,Gamma_0
    double precision, intent(out) :: M_dot_ratio,M_dot,L_0
    double precision, intent(out) :: rho_base, error_b
    character :: dum
    integer :: line

    OPEN(1,FILE='InputStan/params_G2_m0.2.txt')
    READ(1,*) dum, Gamma_0
    READ(1,*) dum, M_dot_ratio
    READ(1,*) dum, M_star
    READ(1,*) dum, L_0
    READ(1,*) dum, R_star
    READ(1,*)
    READ(1,*) dum, M_dot
    READ(1,*) dum, rho_base
    READ(1,*) dum, error_b
    CLOSE(1)

    M_star = M_star*M_sun
    L_0 = L_0*L_sun
    R_star = R_star*R_sun
    M_dot = M_dot*M_sun/year

  end subroutine ReadInParams

  subroutine ReadInTable(r_arr,rho_arr,v_arr,e_arr,Er_arr,T_arr,p_arr)
    use mod_global_parameters
    ! use mod_constants
    ! use mod_fld

    integer :: i
    double precision :: dr
    double precision, intent(out) :: r_arr(domain_nx2+2*nghostcells)
    double precision, intent(out) :: rho_arr(domain_nx2+2*nghostcells)
    double precision, intent(out) :: v_arr(domain_nx2+2*nghostcells)
    double precision, intent(out) :: e_arr(domain_nx2+2*nghostcells)
    double precision, intent(out) :: Er_arr(domain_nx2+2*nghostcells)

    double precision :: i_arr(domain_nx2+2*nghostcells)
    double precision, intent(out) :: T_arr(domain_nx2+2*nghostcells)
    double precision, intent(out) :: p_arr(domain_nx2+2*nghostcells)

    OPEN(1,FILE='InputStan/structure_amrvac_G2_m0.2.txt')
    do i = 1,domain_nx2+2*nghostcells
      READ(1,*) r_arr(i),v_arr(i),rho_arr(i),Er_arr(i)
    enddo
    CLOSE(1)

    T_arr = (Er_arr/const_rad_a)**0.25d0
    p_arr = kb_cgs/(fld_mu*mp_cgs)*T_arr*rho_arr
    e_arr = p_arr/(rhd_gamma - one) + half*rho_arr*v_arr**2

  end subroutine ReadInTable

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: i,j,b
    integer :: NumberOfBlocks
    double precision :: x_perc
    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)
    double precision :: pert(ixO^S)

    NumberOfBlocks = domain_nx2/block_nx2

    x_perc = (x(nghostcells,ixOmin2,2)-xprobmin2)/(xprobmax2-xprobmin2)
    b = floor(x_perc*NumberOfBlocks)

    do i = ixImin2,ixImax2
      j = i + b*block_nx2

      w(:,i,rho_) = rho_arr(j)
      w(:,i,mom(1)) = zero
      w(:,i,mom(2)) = rho_arr(j)*v_arr(j)
      w(:,i,e_) = e_arr(j)
      w(:,i,r_e) = Er_arr(j)
    enddo

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))

    pert(ixO^S) = dsin(2*dpi*x(ixO^S,1)/(xprobmax1-xprobmin1))*dsin(2*dpi*x(ixO^S,2))

    ! where ((x(ixO^S,2) .lt. 3.d0) .and. (x(ixO^S,2) .gt. 1.1d0))
      ! w(ixO^S,rho_) = w(ixO^S,rho_) * (1.d0 + 0.1d0*pert(ixO^S))
    ! end where

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    use mod_fld


    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: kappa(ixB^S), gradE_l(ixB^S), L_vE_l(ixB^S)
    double precision :: Temp(ixI^S), pth(ixI^S),pert(ixB^S)
    integer :: i,j

    double precision :: cool(ixI^S), heat(ixI^S)

    select case (iB)

    case(3)

      w(ixBmin1:ixBmax1,nghostcells,rho_) = dinflo

      do i = ixBmax2-1,ixBmin2,-1
        ! w(ixBmin1:ixBmax1,i,rho_) = 2*w(ixBmin1:ixBmax1,i+1,rho_) - w(ixBmin1:ixBmax1,i+2,rho_)
        w(ixBmin1:ixBmax1,i,rho_) = dexp(2*dlog(w(ixBmin1:ixBmax1,i+1,rho_)) - dlog(w(ixBmin1:ixBmax1,i+2,rho_)))
      enddo

      do i = ixBmax2,ixBmin2,-1
        w(ixBmin1:ixBmax1,i,mom(2)) = w(ixBmin1:ixBmax1,i+1,mom(2))*(x(ixBmin1:ixBmax1,i+1,2)/x(ixBmin1:ixBmax1,i,2))**2
        w(ixBmin1:ixBmax1,i,mom(1)) = w(ixBmin1:ixBmax1,i+1,mom(1))
      enddo

        ! pert(ixb^S) = dsin(2*dpi*x(ixB^S,1)/(xprobmax1-xprobmin1))*dsin(2*dpi*x(ixB^S,2)-2*dpi*global_time)
        ! w(ixB^S,rho_) = w(ixB^S,rho_) * (1.d0 + 0.1d0*pert(ixB^S))
        ! w(ixB^S,mom(2)) = w(ixB^S,mom(2)) * (1.d0 + 0.1d0*pert(ixB^S))


      where (w(ixB^S,mom(2)) .lt. zero)
         w(ixB^S,mom(2)) = zero ! abs(w(ixB^S,mom(2)))
      end where

      ! print*, w(5,1:5,mom(2))

      call fld_get_opacity(w, x, ixI^L, ixB^L, kappa)

      L_vE_l(ixB^S) = 4*dpi*x(ixB^S,2)**2*4.d0/3.d0*w(ixB^S,mom(2))/w(ixB^S,rho_)*w(ixB^S,r_e)
      gradE_l(ixB^S) = -w(ixB^S,rho_)*kappa(ixB^S)*(3.d0*unit_velocity/const_c)*(L_0-L_vE_l(ixB^S))/(4.d0*dpi*x(ixB^S,2)**2.d0)
      gradE = sum(gradE_l(ixBmin1:ixBmax1,nghostcells))/(ixBmax1-ixBmin1)

      do i = ixBmax2-1,ixBmin2,-1
        w(ixBmin1:ixBmax1,i,r_e) = w(ixBmin1:ixBmax1,i+2,r_e) &
        + (x(ixBmin1:ixBmax1,i,2)-x(ixBmin1:ixBmax1,i+2,2))*gradE_l(ixBmin1:ixBmax1,i+1)
      enddo
      w(ixBmin1:ixBmax1,nghostcells,r_e) = dexp(half*(dlog(w(ixBmin1:ixBmax1,nghostcells-1,r_e))+dlog(w(ixBmin1:ixBmax1,nghostcells+1,r_e))))

      ! Temp(ixB^S) = (w(ixB^S,r_e)*unit_pressure/const_rad_a)**0.25d0
      ! pth(ixB^S) = Temp(ixB^S)*w(ixB^S,rho_)*const_kb/(const_mp*fld_mu)*unit_density/unit_pressure
      ! w(ixB^S,e_) = pth(ixB^S)/(rhd_gamma-1) + half*(w(ixB^S,mom(1))**2+w(ixB^S,mom(2))**2)/w(ixB^S,rho_)

      do i = ixBmax2,ixBmin2,-1
        w(ixBmin1:ixBmax1,i,e_) = dexp(2*dlog(w(ixBmin1:ixBmax1,i+1,e_)) - dlog(w(ixBmin1:ixBmax1,i+2,e_)))
      enddo


      ! !>>>>>>>>>>>>>>
      ! if (x(1,1,2) .lt. 1) then
      !   pth(ixI^S) = (w(ixI^S,e_) - half/w(ixI^S,rho_)*w(ixI^S,mom(2))**2) &
      !   *(rhd_gamma -1)
      !   temp(ixI^S) = pth(ixI^S)/w(ixI^S,rho_)
      !   cool(ixI^S) = temp(ixI^S)
      !
      !   temp(ixI^S) = (w(ixI^S,r_e)*unit_pressure/const_rad_a)**(1.d0/4.d0)/unit_temperature
      !   heat(ixI^S) = temp(ixI^S)
      !
      !   print*, it,'-------------------------------------------------------------'
      !   print*, cool(5,1:6)
      !   print*, heat(5,1:6)
      !   print*, cool(5,1:6) - heat(5,1:6)
      ! endif

    case(4)
      do i = ixBmin2,ixBmax2
        !> Conserve gradE/rho
        w(ixBmin1:ixBmax1,i,r_e) = (x(ixBmin1:ixBmax1,i-1,2)**2*w(ixBmin1:ixBmax1,i-1,mom(2))/w(ixBmin1:ixBmax1,i-1,rho_))&
        /(x(ixBmin1:ixBmax1,i,2)**2*w(ixBmin1:ixBmax1,i,mom(2))/w(ixBmin1:ixBmax1,i,rho_))*w(ixBmin1:ixBmax1,i-1,r_e)
      enddo

      gradE_out = sum((w(ixBmin1:ixBmax1,ixBmin2-1,r_e)-w(ixBmin1:ixBmax1,ixBmin2-2,r_e))&
      /(x(ixBmin1:ixBmax1,ixBmin2-1,2) - x(ixBmin1:ixBmax1,ixBmin2-2,2)))/(ixBmax1-ixBmin1)

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

  subroutine mg_boundary_conditions(iB)
    use mod_global_parameters
    use mod_multigrid_coupling

    integer, intent(in)             :: iB

    select case (iB)
      case (3)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
        mg%bc(iB, mg_iphi)%bc_value = gradE

      case (4)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
        mg%bc(iB, mg_iphi)%bc_value = gradE_out

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

    radius(ixO^S) = x(ixO^S,2)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    gravity_field(ixI^S,:) = zero
    gravity_field(ixO^S,2) = -const_G*mass/radius(ixO^S)**2*(unit_time**2/unit_length)

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
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt*ppsource(ixO^S,mom(1))
    w(ixO^S,mom(2)) = w(ixO^S,mom(2)) + qdt*ppsource(ixO^S,mom(2))
    w(ixO^S,e_) = w(ixO^S,e_) + qdt*ppsource(ixO^S,e_)
    w(ixO^S,r_e) = w(ixO^S,r_e) + qdt*ppsource(ixO^S,r_e)

  end subroutine PseudoPlanar

  subroutine PseudoPlanarSource(ixI^L,ixO^L,w,x,source,boundary)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: source(ixO^S,1:nw)
    logical, intent(in) :: boundary

    double precision :: rad_flux(ixO^S,1:ndir)
    double precision :: pth(ixI^S),v(ixO^S,2)
    double precision :: radius(ixO^S),  pert(ixO^S)
    integer :: rdir, pdir

    source(ixO^S,1:nw) = zero

    rdir = 2
    pdir = 1

    v(ixO^S,rdir) = w(ixO^S,mom(rdir))/w(ixO^S,rho_)
    v(ixO^S,pdir) = w(ixO^S,mom(pdir))/w(ixO^S,rho_)

    radius(ixO^S) = x(ixO^S,rdir) ! + half*block%dx(ixO^S,rdir)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    source(ixO^S,rho_) = -two*w(ixO^S,rho_)*v(ixO^S,rdir)/radius(ixO^S)

    pth(ixO^S) = (rhd_gamma-1) &
    *(w(ixO^S,e_) - half*sum(w(ixO^S, mom(:))**2, dim=ndim+1) / w(ixO^S, rho_))

    !> dm_r/dt = +(rho*v_p**2 + 2pth)/r -2 (rho*v_r**2 + pth)/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    source(ixO^S,mom(rdir)) = - 2*w(ixO^S,rho_)*v(ixO^S,rdir)**two/radius(ixO^S) &
                              + w(ixO^S,rho_)*v(ixO^S,pdir)**two/radius(ixO^S)

    source(ixO^S,mom(pdir)) = - 3*v(ixO^S,rdir)*v(ixO^S,pdir)*w(ixO^S,rho_)/radius(ixO^S)

    !> de/dt = -2 (e+p)v_r/r
    source(ixO^S,e_) = - two*(w(ixO^S,e_)+pth(ixO^S))*v(ixO^S,rdir)/radius(ixO^S)

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      if (boundary) then
        rad_flux(ixO^S,rdir) = L_0/(4.d0*dpi*x(ixO^S,rdir)**2)
      else
        call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
      endif
      source(ixO^S,r_e) = source(ixO^S,r_e) - two*rad_flux(ixO^S,rdir)/radius(ixO^S)
    endif

    if (rhd_radiation_advection) then
      source(ixO^S,r_e) = source(ixO^S,r_e) - two*w(ixO^S,r_e)*v(ixO^S,rdir)/radius(ixO^S)
    endif

    ! Not sure about this one
    if (rhd_energy_interact) then
      source(ixO^S,r_e) = source(ixO^S,r_e) + two*v(ixO^S,rdir)*w(ixO^S,r_e)/(3*radius(ixO^S))
    endif

  end subroutine PseudoPlanarSource


  subroutine Opacity_stepfunction(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    kappa(ixO^S) = kappa_b + (1.d0+erf((x(ixO^S,2)-one)*error_b-error_b/2.d0))*(kappa_0-kappa_b)/2.d0

  end subroutine Opacity_stepfunction

  ! subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
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
  !   ! test with different levels of refinement enforced
  !   if (it .gt. 1) then
  !     if (any(x(ix^S,2) < 3.d0/4.d0 * xprobmax2)) refine=1
  !     if (any(x(ix^S,2) < 2.d0/4.d0 * xprobmax2)) refine=1
  !     if (any(x(ix^S,2) < 1.d0/4.d0 * xprobmax2)) refine=1
  !   endif
  !
  ! end subroutine specialrefine_grid


  subroutine time_average_values(ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    if (global_time .gt. 0.5d0) then
      w(ixI^S,int_r) = w(ixI^S,int_r) &
      + w(ixI^S,rho_)*dt
      w(ixI^S,int_v) =  w(ixI^S,int_v) &
      + w(ixI^S,mom(2))/w(ixI^S,rho_)*dt
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
    double precision                   :: kappa(ixO^S)
    double precision                   :: rad_flux(ixO^S,1:ndim)
    integer                            :: idim
    double precision :: radius(ixI^S)
    double precision :: mass

    radius(ixO^S) = x(ixO^S,2)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)

    g_rad(ixO^S) = kappa(ixO^S)*rad_flux(ixO^S,2)/(const_c/unit_velocity)
    g_grav(ixO^S) = const_G*mass/radius(ixO^S)**2*(unit_time**2/unit_length)
    big_gamma(ixO^S) = g_rad(ixO^S)/g_grav(ixO^S)

    call rhd_get_tgas(w, x, ixI^L, ixO^L, Tgas)
    call rhd_get_trad(w, x, ixI^L, ixO^L, Trad)

    w(ixO^S,nw+1) = Tgas(ixO^S)*unit_temperature
    w(ixO^S,nw+2) = Trad(ixO^S)*unit_temperature
    w(ixO^S,nw+3) = big_gamma(ixO^S)

    w(ixO^S,nw+4) = 4*dpi*w(ixO^S,mom(2))*radius(ixO^S)**2 &
    *unit_density*unit_velocity/M_sun*year

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'Tgas Trad Gamma Mdot'
  end subroutine specialvarnames_output

end module mod_usr
