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
  double precision :: Gamma_b, Gamma_0
  double precision :: kappa_b, kappa_0
  double precision :: L_0
  double precision :: M_star
  double precision :: R_star, R_0, R_b
  double precision :: M_dot
  double precision :: M_dot_max
  double precision :: rho_b
  double precision :: v_esc, v_inf, v_0, bb
  double precision :: lambda_0

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

    ! Reset gas energy to radiation temperature
    usr_internal_bc => set_gas_energy

    ! PseudoPlanar correction
    usr_source => PseudoPlanar

    ! Graviatational field
    usr_gravity => set_gravitation_field

    ! Special Opacity
    usr_special_opacity => Opacity_stepfunction

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters

    !> Set stellar mass and radius
    call ReadInParams(M_star,R_star,Gamma_0,M_dot_ratio,M_dot,L_0)

    R_b = xprobmin2*R_star
    R_0 = (xprobmin2+0.1d0)*R_star

    v_esc = dsqrt(2*const_G*M_star/R_star)
    v_inf = dsqrt(Gamma_0 - one)*v_esc
    v_0 = 1.d5

    bb = one - (v_0/v_inf)**two

    Gamma_b = 0.9d0

    rho_b = M_dot/(4*dpi*R_star**2*v_0)

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
    call GetArrays(r_arr,rho_arr,v_arr,e_arr,Er_arr,T_arr, p_arr)

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = R_star !r_arr(nghostcells) ! cm
    unit_velocity   = v_0
    unit_numberdensity = rho_b/((1.d0+4.d0*He_abundance)*mp_cgs)

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
    unit_pressure=(2.d0+3.d0*He_abundance)*unit_numberdensity*kB_cgs*unit_temperature
    unit_velocity=dsqrt(unit_pressure/unit_density)
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*kB_cgs)
    unit_time=unit_length/unit_velocity
    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    if (mype .eq. 0) then
      print*, 'M_star ', 'R_star ','M_dot_ratio ', 'M_dot ', 'L_0'
      print*, 'R_b ', 'R_0 ',  'Gamma_b ', 'Gamma_0 ', 'kappa_b ', 'kappa_0'
      print*, M_star, R_star,M_dot_ratio, M_dot, L_0
      print*, R_b, R_0,  Gamma_b, Gamma_0, kappa_b, kappa_0

      print*, 'Flux at boundary: ', L_0/(4*dpi*R_b**2)

      print*, 'unit_length', unit_length
      print*, 'unit_density', unit_density
      print*, 'unit_pressure', unit_pressure
      print*, 'unit_temperature', unit_temperature
      print*, 'unit_radflux', unit_radflux
      print*, 'unit_opacity', unit_opacity
      print*, 'unit_time', unit_time
      print*, 'unit_velocity', unit_velocity
    endif

    !> Make all parameters dimensionless
    M_star = M_star/(unit_density*unit_length**3.d0)
    R_star = R_star/unit_length
    M_dot = M_dot/(unit_density*unit_length**3.d0)*unit_time
    L_0 = L_0/(unit_pressure*unit_length**3.d0)*unit_time
    R_b = R_b/unit_length
    R_0 = R_0/unit_length
    kappa_0 = kappa_0/unit_opacity
    kappa_b = kappa_b/unit_opacity

    !> Make initial profiles dimensionless
    r_arr = r_arr/unit_length
    rho_arr = rho_arr/unit_density
    v_arr = v_arr/unit_velocity
    e_arr = e_arr/unit_pressure
    Er_arr = Er_arr/unit_pressure
    T_arr = T_arr/unit_temperature
    p_arr = p_arr/unit_pressure

    if (mype .eq. 0) then
      print*, 'M_star ', 'R_star ','M_dot_ratio ', 'M_dot ', 'L_0'
      print*, 'R_b ', 'R_0 ',  'Gamma_b ', 'Gamma_0 ', 'kappa_b ', 'kappa_0'
      print*, M_star, R_star,M_dot_ratio, M_dot, L_0
      print*, R_b, R_0,  Gamma_b, Gamma_0, kappa_b, kappa_0

      print*, 'Flux at boundary: ', L_0/(4*dpi*R_b**2)
    endif

  end subroutine initglobaldata_usr

  subroutine ReadInParams(M_star,R_star,Gamma_0,M_dot_ratio,M_dot,L_0)
    use mod_global_parameters
    double precision, intent(out) :: M_star,R_star,Gamma_0
    double precision, intent(out) :: M_dot_ratio,M_dot,L_0
    character :: dum
    integer :: line

    OPEN(1,FILE='InitialConditions/init_params_amrvac')
    READ(1,*) dum, Gamma_0
    READ(1,*) dum, M_dot_ratio
    READ(1,*) dum, M_star
    READ(1,*) dum, L_0
    READ(1,*) dum, R_star
    READ(1,*)
    READ(1,*) dum, M_dot
    CLOSE(1)

    M_star = M_star*M_sun
    L_0 = L_0*L_sun
    R_star = R_star*R_sun
    M_dot = M_dot*M_sun/year

  end subroutine ReadInParams

  subroutine GetArrays(r_arr,rho_arr,v_arr,e_arr,Er_arr,T_arr,p_arr)
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

    dr = (xprobmax2-xprobmin2)*R_star/domain_nx2
    do i = 1,domain_nx2+2*nghostcells
      r_arr(i) = R_b+(i-1-nghostcells)*dr
    enddo

    v_arr = v_inf*(one - bb*R_star/r_arr)**half
    rho_arr = M_dot/(4*dpi*r_arr**2*v_arr)

    Er_arr(domain_nx2+2*nghostcells) = M_dot*Gamma_0*const_G*M_star/(4*dpi*v_inf)&
    *(one/r_arr(domain_nx2+2*nghostcells))**3.d0

    do i = domain_nx2+2*nghostcells-1,1,-1
      Er_arr(i) = Er_arr(i+1)+(3*Gamma_0*const_G*M_star*(rho_arr(i)/r_arr(i)**2 + rho_arr(i+1)/r_arr(i+1)**2))*half*dr
    enddo

    T_arr = (Er_arr/const_rad_a)**0.25d0
    p_arr = const_kb/(fld_mu*const_mp)*T_arr*rho_arr
    e_arr = p_arr/(rhd_gamma - one) + half*rho_arr*v_arr**2.d0

    ! do i = domain_nx2+2*nghostcells-1,1,-1
    !   print*, i, r_arr(i)/R_star, T_arr(i), rho_arr(i)/rho_b, Er_arr(i), v_arr(i)/1.d5
    ! enddo
    !
    ! stop

  end subroutine GetArrays

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: i,j,b
    integer :: NumberOfBlocks
    double precision :: x_perc

    NumberOfBlocks = domain_nx2/block_nx2

    x_perc = (x(nghostcells,ixOmin2,2)-xprobmin2)/(xprobmax2-xprobmin2)
    b = floor(x_perc*NumberOfBlocks)

    ! print*, 'block number', b, 'x_perc', x_perc

    do i = ixImin2,ixImax2
      j = i + b*block_nx2
      ! if (b .ne. 0) j = j-nghostcells

      w(:,i,rho_) = rho_arr(j)
      w(:,i,mom(1)) = zero
      w(:,i,mom(2)) = rho_arr(j)*v_arr(j)
      w(:,i,e_) = e_arr(j)
      w(:,i,r_e) = Er_arr(j)
      ! print*, b,j,i,r_arr(j),v_arr(j),rho_arr(j)
    enddo

    call get_rad_extravars(w, x, ixI^L, ixO^L)

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: a(ixImin1:ixImax1),b(ixImin1:ixImax1),c(ixImin1:ixImax1),d(ixImin1:ixImax1)
    double precision :: Temp(ixI^S), Press(ixI^S)

    integer :: i,j

    select case (iB)

    case(3)

      if (mype .eq. 0) then
      print*,it, '#############', 'Before bound'
      do i = ixImin2,ixImax2
        print*, w(nghostcells+1,i,r_e),w(nghostcells+1,i,i_diff_mg),w(nghostcells+1,i,e_),w(nghostcells+1,i,rho_)
      enddo
      endif


      do i = ixBmax2,ixBmin2,-1
        w(ixImin1:ixImax1,i,rho_) = rho_arr(i)
        w(ixImin1:ixImax1,i,mom(1)) = w(ixImin1:ixImax1,i+1,mom(1))
        w(ixImin1:ixImax1,i,mom(2)) = rho_arr(i)*v_arr(i) !w(ixImin1:ixImax1,i+1,mom(2)) &
        ! *(x(ixImin1:ixImax1,i+1,2)/x(ixImin1:ixImax1,i,2))**2

        a(ixImin1:ixImax1) = L_0/(4.d0*dpi*x(ixImin1:ixImax1,i,2)**2)
        b(ixImin1:ixImax1) = w(ixImin1:ixImax1,i+1,r_e)*fld_speedofligt_0 &
        /(3.d0*(x(ixImin1:ixImax1,i+1,2)-x(ixImin1:ixImax1,i,2))*rho_arr(i)*kappa_b)
        c(ixImin1:ixImax1) =fld_speedofligt_0 &
        /(3.d0*(x(ixImin1:ixImax1,i+1,2)-x(ixImin1:ixImax1,i,2))*rho_arr(i)*kappa_b)
        d(ixImin1:ixImax1) = 4.d0/3.d0*abs(w(ixImin1:ixImax1,i,mom(2))/w(ixImin1:ixImax1,i,rho_))
        w(ixImin1:ixImax1,i,r_e) = (a(ixImin1:ixImax1) + b(ixImin1:ixImax1))/(c(ixImin1:ixImax1) + d(ixImin1:ixImax1))
        w(ixImin1:ixImax1,i,r_e) = max(w(ixImin1:ixImax1,i,r_e),half*Er_arr(i))
        ! w(ixImin1:ixImax1,i,r_e) = Er_arr(i)

        ! Temp(ixImin1:ixImax1,i) = (w(ixImin1:ixImax1,i,r_e)/const_rad_a)**0.25d0
        ! Press(ixImin1:ixImax1,i) = const_kb/(fld_mu*const_mp)*Temp(ixImin1:ixImax1,i)*w(ixImin1:ixImax1,i,rho_) &
        ! *unit_temperature*unit_density/unit_pressure
        ! w(ixImin1:ixImax1,i,e_) = Press(ixImin1:ixImax1,i)/(rhd_gamma - one) &
        ! + half*(w(ixImin1:ixImax1,i,mom(1))**2.d0+w(ixImin1:ixImax1,i,mom(2))**2.d0)/w(ixImin1:ixImax1,i,rho_)

        ! w(ixImin1:ixImax1,i,e_) = e_arr(i)

        w(ixImin1:ixImax1,i,e_) = w(ixImin1:ixImax1,i+1,e_)

      enddo

      if (mype .eq. 0) then
      print*,it, '#############', 'After bound'
      do i = ixImin2,ixImax2
        print*, w(nghostcells+1,i,r_e),w(nghostcells+1,i,i_diff_mg),w(nghostcells+1,i,e_),w(nghostcells+1,i,rho_)
      enddo
      endif

    case(4)
      do i = ixBmin2,ixBmax2
        !> Conserve gradE/rho
        w(ixImin1:ixImax1,i,r_e) = w(ixImin1:ixImax1,i-1,rho_)/w(ixImin1:ixImax1,i-2,rho_) &
        *(w(ixImin1:ixImax1,i-1,r_e) - w(ixImin1:ixImax1,i-2,r_e)) + w(ixImin1:ixImax1,i-1,r_e)
        do j = ixImin1,ixImax1
          w(j,i,r_e) = min(w(j,i,r_e),w(j,i-1,r_e))
        enddo
      enddo

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions


  subroutine mg_boundary_conditions(qt,ixI^L,ixO^L,iB,w,x)

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)

    double precision :: a(ixImin1:ixImax1),b(ixImin1:ixImax1),c(ixImin1:ixImax1),d(ixImin1:ixImax1)
    double precision :: mean_RE(ixImin1:ixImax1)
    integer :: i

    select case (iB)
      case (3)

        i = nghostcells

        ! a(ixImin1:ixImax1) = L_0/(4.d0*dpi*x(ixImin1:ixImax1,i,2)**2)
        ! b(ixImin1:ixImax1) = w(ixImin1:ixImax1,i+1,r_e)*fld_speedofligt_0 &
        ! /(3.d0*(x(ixImin1:ixImax1,i+1,2)-x(ixImin1:ixImax1,i,2))*w(ixImin1:ixImax1,i+1,rho_)*kappa_b)
        ! c(ixImin1:ixImax1) =fld_speedofligt_0 &
        ! /(3.d0*(x(ixImin1:ixImax1,i+1,2)-x(ixImin1:ixImax1,i,2))*w(ixImin1:ixImax1,i+1,rho_)*kappa_b)
        ! d(ixImin1:ixImax1) = 4.d0/3.d0*abs(w(ixImin1:ixImax1,i,mom(2))/w(ixImin1:ixImax1,i,rho_))
        ! mean_RE(ixImin1:ixImax1) = (a(ixImin1:ixImax1) + b(ixImin1:ixImax1))/(c(ixImin1:ixImax1) + d(ixImin1:ixImax1))

        mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
        ! mg%bc(iB, mg_iphi)%bc_value = sum(mean_RE(ixOmin1:ixOmax1))/(ixOmax1-ixOmin1)

      case (4)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous

      case default
        print *, "Not a standard: ", trim(typeboundary(iw_r_e, iB))
        error stop "You have to set a user-defined boundary method"
    end select
  end subroutine mg_boundary_conditions


  !> internal boundary, user defined
  !
  !> This subroutine can be used to artificially overwrite ALL conservative
  !> variables in a user-selected region of the mesh, and thereby act as
  !> an internal boundary region. It is called just before external (ghost cell)
  !> boundary regions will be set by the BC selection. Here, you could e.g.
  !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
  !> which can be used to identify the internal boundary region location.
  !> Its effect should always be local as it acts on the mesh.
  subroutine set_gas_energy(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_physics
    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision :: Trad(ixI^S),pgas(ixI^S)
    integer :: i

    if (mype .eq. 0) then
    print*,it, '#############', 'Before internal bound'
    do i = ixImin2,ixImax2
      print*, w(nghostcells+1,i,r_e),w(nghostcells+1,i,i_diff_mg),w(nghostcells+1,i,e_),w(nghostcells+1,i,rho_)
    enddo
    endif

    call phys_get_trad(w,x,ixI^L,ixO^L,Trad)
    pgas(ixI^S) = const_kB/(fld_mu*const_mp)*Trad(ixI^S)*w(ixI^S,rho_) &
    *unit_temperature*unit_density/unit_pressure

    w(ixI^S,e_) = pgas(ixI^S)/(rhd_gamma - 1) &
    + half*(w(ixI^S,mom(1))**2+w(ixI^S,mom(2))**2)/w(ixI^S,rho_)

    if (mype .eq. 0) then
    print*,it, '#############', 'After internal bound'
    do i = ixImin2,ixImax2
      print*, w(nghostcells+1,i,r_e),w(nghostcells+1,i,i_diff_mg),w(nghostcells+1,i,e_),w(nghostcells+1,i,rho_)
    enddo
    endif

  end subroutine set_gas_energy


  !> Calculate gravitational acceleration in each dimension
  subroutine set_gravitation_field(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    double precision :: radius(ixI^S)
    double precision :: mass

    radius(ixI^S) = x(ixI^S,2)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    gravity_field(ixI^S,1) = zero
    gravity_field(ixI^S,2) = -const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length)

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

    double precision :: wCCT(ixI^S,1:nw)
    double precision :: pth(ixI^S),v(ixI^S,2)
    double precision :: radius(ixI^S)
    integer :: rdir, pdir

    rdir = 2
    pdir = 1

    v(ixI^S,1) = wCT(ixI^S,mom(1))/wCT(ixI^S,rho_)
    v(ixI^S,2) = wCT(ixI^S,mom(2))/wCT(ixI^S,rho_)

    radius(ixO^S) = x(ixO^S,2)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    w(ixO^S,rho_) = w(ixO^S,rho_) - qdt*two*wCT(ixO^S,mom(rdir))/radius(ixO^S)

    !> dm_r/dt = +rho*v_p**2/r -rho*v_p**2/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    w(ixI^S,mom(rdir)) = w(ixI^S,mom(rdir)) + qdt*wCT(ixI^S,rho_)*v(ixI^S,pdir)**two/x(ixI^S,rdir) &
                                            - qdt*2*wCT(ixI^S,rho_)*v(ixI^S,rdir)**two/x(ixI^S,rdir)
    w(ixI^S,mom(pdir)) = w(ixI^S,mom(pdir)) - qdt*3*v(ixI^S,rdir)*v(ixI^S,pdir)*wCT(ixI^S,rho_)/x(ixI^S,rdir)


    !> de/dt = -2 (e+p)v_r/r
    call phys_get_pthermal(wCT,x,ixI^L,ixO^L,pth)
    w(ixO^S,e_) = w(ixO^S,e_) - qdt*two*(wCT(ixO^S,e_)+pth(ixO^S))*wCT(ixO^S,mom(rdir))/(wCT(ixO^S,rho_)*radius(ixO^S))

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      wCCT = wCT
      call get_rad_extravars(wCCT, x, ixI^L, ixO^L)
      w(ixO^S,r_e) = w(ixO^S,r_e) - qdt*two*wCCT(ixO^S,i_flux(rdir))/radius(ixO^S)
    endif

    if (rhd_radiation_advection) then
      w(ixO^S,r_e) = w(ixO^S,r_e) - qdt*two*wCT(ixO^S,r_e)*wCT(ixO^S,mom(rdir))/(wCT(ixO^S,rho_)*radius(ixO^S))
    endif

  end subroutine PseudoPlanar

  subroutine Opacity_stepfunction(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    double precision :: new_x(ixO^S)

    new_x(ixO^S) = (x(ixO^S,2) - R_0)*2d1
    kappa(ixO^S) = kappa_b + (kappa_0-kappa_b)*half*(one+erf(new_x(ixO^S)))

    ! kappa(ixO^S) = kappa_0
    ! where (x(ixO^S,2) .lt. R_0)
    !   kappa(ixO^S) = kappa_b
    ! endwhere

  end subroutine Opacity_stepfunction

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
    integer                            :: idim
    double precision :: radius(ixI^S)
    double precision :: mass

    radius(ixI^S) = x(ixI^S,2)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    call get_rad_extravars(w, x, ixI^L, ixO^L)

    g_rad(ixO^S) = w(ixO^S,i_op)*w(ixO^S,i_flux(2))/fld_speedofligt_0
    g_grav(ixI^S) = const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length)
    big_gamma(ixO^S) = g_rad(ixO^S)/g_grav(ixO^S)

    call rhd_get_tgas(w, x, ixI^L, ixO^L, Tgas)
    call rhd_get_trad(w, x, ixI^L, ixO^L, Trad)

    w(ixO^S,nw+1) = Tgas(ixO^S)*unit_temperature
    w(ixO^S,nw+2) = Trad(ixO^S)*unit_temperature
    w(ixO^S,nw+3) = big_gamma(ixO^S)
    w(ixO^S,nw+4) = 4*dpi*w(ixO^S,mom(2))*radius(ixI^S)**2 &
    *unit_density*unit_velocity/M_sun*year
    w(ixO^S,nw+5) = w(ixO^S,i_flux(2)) + w(ixO^S,r_e)*w(ixO^S,mom(2))/w(ixO^S,rho_)


  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'Tgas Trad Gamma Mdot Ftot'
  end subroutine specialvarnames_output

end module mod_usr
