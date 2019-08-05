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
  double precision :: Gamma_0
  double precision :: kappa_0
  double precision :: L_0
  double precision :: M_star
  double precision :: R_star
  double precision :: M_dot
  double precision :: sp_sos
  double precision :: sp_rho

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

    integer :: i

    !> Set stellar mass and radius
    call ReadInParams(M_star,R_star,Gamma_0,M_dot_ratio,M_dot,L_0,sp_sos)

    kappa_0 = Gamma_0*4*dpi*const_G*M_star*const_c/L_0

    allocate(r_arr(domain_nx2+2*nghostcells))
    allocate(rho_arr(domain_nx2+2*nghostcells))
    allocate(v_arr(domain_nx2+2*nghostcells))
    allocate(e_arr(domain_nx2+2*nghostcells))
    allocate(Er_arr(domain_nx2+2*nghostcells))
    allocate(T_arr(domain_nx2+2*nghostcells))
    allocate(p_arr(domain_nx2+2*nghostcells))

    ! call ReadInTable(r_arr,rho_arr,v_arr,e_arr,Er_arr,T_arr, p_arr)
    call ReadInTable(r_arr,rho_arr,v_arr,e_arr,Er_arr,T_arr, p_arr)

    sp_rho = M_dot/(4*dpi*R_star**2*sp_sos)

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = R_star !r_arr(nghostcells) ! cm
    unit_velocity   = sp_sos
    unit_numberdensity = sp_rho/((1.d0+4.d0*He_abundance)*mp_cgs)

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
    unit_pressure=(2.d0+3.d0*He_abundance)&
       *unit_numberdensity*kB_cgs*unit_temperature
    unit_velocity=dsqrt(unit_pressure/unit_density)
    unit_temperature=unit_pressure/((2.d0+&
       3.d0*He_abundance)*unit_numberdensity*kB_cgs)
    unit_time=unit_length/unit_velocity
    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    if (mype .eq. 0) then
      print*, 'M_star ', 'R_star ','M_dot_ratio ', 'M_dot ', 'L_0'
      print*, 'Gamma_0 ', 'kappa_0'
      print*, M_star, R_star,M_dot_ratio, M_dot, L_0
      print*, Gamma_0, kappa_0

      print*, 'Flux at boundary: ', L_0/(4*dpi*R_star**2)

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
    kappa_0 = kappa_0/unit_opacity

    !> Make initial profiles dimensionless
    r_arr = r_arr/unit_length
    rho_arr = rho_arr/unit_density
    v_arr = v_arr/unit_velocity
    e_arr = e_arr/unit_pressure
    Er_arr = Er_arr/unit_pressure
    T_arr = T_arr/unit_temperature
    p_arr = p_arr/unit_pressure


    if (mype .eq. 0) then
      do i = 1,domain_nx2+2*nghostcells
        print*, r_arr(i), rho_arr(i), v_arr(i), e_arr(i), Er_arr(i)
      enddo
    endif

    stop

    if (mype .eq. 0) then
      print*, 'M_star ', 'R_star ','M_dot_ratio ', 'M_dot ', 'L_0'
      print*, 'Gamma_0 ', 'kappa_0'
      print*, M_star, R_star,M_dot_ratio, M_dot, L_0
      print*, Gamma_0, kappa_0

      print*, 'Flux at boundary: ', L_0/(4*dpi*R_star**2)
    endif

  end subroutine initglobaldata_usr

  subroutine ReadInParams(M_star,R_star,Gamma_0,M_dot_ratio,M_dot,L_0,sp_sos)
    use mod_global_parameters
    double precision, intent(out) :: M_star,R_star,Gamma_0
    double precision, intent(out) :: M_dot_ratio,M_dot,L_0
    double precision, intent(out) :: sp_sos
    character :: dum
    integer :: line

    OPEN(1,FILE='InputLuka/params.txt')
    READ(1,*) dum, Gamma_0
    READ(1,*) dum, M_dot_ratio
    READ(1,*) dum, M_star
    READ(1,*) dum, L_0
    READ(1,*) dum, R_star
    READ(1,*)
    READ(1,*) dum, M_dot
    READ(1,*) dum, sp_sos
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

    OPEN(1,FILE='InputLuka/structure_amrvac.txt')
    do i = 1,domain_nx2+2*nghostcells
      READ(1,*) r_arr(i),v_arr(i),rho_arr(i),Er_arr(i)
    enddo
    CLOSE(1)

    T_arr = (Er_arr/const_rad_a)**0.25d0
    p_arr = const_kb/(fld_mu*const_mp)*T_arr*rho_arr
    e_arr = p_arr/(rhd_gamma - one) + half*rho_arr*v_arr**2.d0

  end subroutine ReadInTable

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
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

    call get_rad_extravars(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2)

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2,iB,w,x)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixBmin1,ixBmin2,ixBmax1,ixBmax2, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: a(ixImin1:ixImax1),b(ixImin1:ixImax1),&
       c(ixImin1:ixImax1),d(ixImin1:ixImax1)
    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2),&
        Press(ixImin1:ixImax1,ixImin2:ixImax2)

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
        ! *(x(ixImin1:ixImax1,i+1,2)/x(ixImin1:ixImax1,i,2))**2

        w(ixImin1:ixImax1,i,mom(2)) = rho_arr(i)*v_arr(i) !w(ixImin1:ixImax1,i+1,mom(2)) a(ixImin1:ixImax1) = L_0/(4.d0*dpi*x(ixImin1:ixImax1,i,2)**2)
        b(ixImin1:ixImax1) = w(ixImin1:ixImax1,i+1,&
           r_e)*fld_speedofligt_0 /(3.d0*(x(ixImin1:ixImax1,i+1,&
           2)-x(ixImin1:ixImax1,i,2))*rho_arr(i)*kappa_0)
        c(ixImin1:ixImax1) =fld_speedofligt_0 /(3.d0*(x(ixImin1:ixImax1,i+1,&
           2)-x(ixImin1:ixImax1,i,2))*rho_arr(i)*kappa_0)
        d(ixImin1:ixImax1) = 4.d0/3.d0*abs(w(ixImin1:ixImax1,i,&
           mom(2))/w(ixImin1:ixImax1,i,rho_))
        w(ixImin1:ixImax1,i,r_e) = (a(ixImin1:ixImax1) + &
           b(ixImin1:ixImax1))/(c(ixImin1:ixImax1) + d(ixImin1:ixImax1))
        w(ixImin1:ixImax1,i,r_e) = max(w(ixImin1:ixImax1,i,r_e),&
           half*Er_arr(i))
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
        w(ixImin1:ixImax1,i,r_e) = w(ixImin1:ixImax1,i-1,&
           rho_)/w(ixImin1:ixImax1,i-2,rho_) *(w(ixImin1:ixImax1,i-1,&
           r_e) - w(ixImin1:ixImax1,i-2,r_e)) + w(ixImin1:ixImax1,i-1,r_e)
        do j = ixImin1,ixImax1
          w(j,i,r_e) = min(w(j,i,r_e),w(j,i-1,r_e))
        enddo
      enddo

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions


  subroutine mg_boundary_conditions(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iB,w,x)

    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: a(ixImin1:ixImax1),b(ixImin1:ixImax1),&
       c(ixImin1:ixImax1),d(ixImin1:ixImax1)
    double precision :: mean_RE(ixImin1:ixImax1)
    integer :: i

    select case (iB)
      case (3)

        i = nghostcells

        ! a(ixImin1:ixImax1) = L_0/(4.d0*dpi*x(ixImin1:ixImax1,i,2)**2)
        ! b(ixImin1:ixImax1) = w(ixImin1:ixImax1,i+1,r_e)*fld_speedofligt_0 &
        ! /(3.d0*(x(ixImin1:ixImax1,i+1,2)-x(ixImin1:ixImax1,i,2))*w(ixImin1:ixImax1,i+1,rho_)*kappa_0)
        ! c(ixImin1:ixImax1) =fld_speedofligt_0 &
        ! /(3.d0*(x(ixImin1:ixImax1,i+1,2)-x(ixImin1:ixImax1,i,2))*w(ixImin1:ixImax1,i+1,rho_)*kappa_0)
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
  subroutine set_gas_energy(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    use mod_physics
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    double precision :: Trad(ixImin1:ixImax1,ixImin2:ixImax2),&
       pgas(ixImin1:ixImax1,ixImin2:ixImax2)
    integer :: i

    if (mype .eq. 0) then
    print*,it, '#############', 'Before internal bound'
    do i = ixImin2,ixImax2
      print*, w(nghostcells+1,i,r_e),w(nghostcells+1,i,i_diff_mg),w(nghostcells+1,i,e_),w(nghostcells+1,i,rho_)
    enddo
    endif

    call phys_get_trad(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,Trad)
    pgas(ixImin1:ixImax1,ixImin2:ixImax2) = &
       const_kB/(fld_mu*const_mp)*Trad(ixImin1:ixImax1,&
       ixImin2:ixImax2)*w(ixImin1:ixImax1,ixImin2:ixImax2,&
       rho_) *unit_temperature*unit_density/unit_pressure

    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = pgas(ixImin1:ixImax1,&
       ixImin2:ixImax2)/(rhd_gamma - 1) + half*(w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(1))**2+w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(2))**2)/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)

    if (mype .eq. 0) then
    print*,it, '#############', 'After internal bound'
    do i = ixImin2,ixImax2
      print*, w(nghostcells+1,i,r_e),w(nghostcells+1,i,i_diff_mg),w(nghostcells+1,i,e_),w(nghostcells+1,i,rho_)
    enddo
    endif

  end subroutine set_gas_energy


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

    radius(ixImin1:ixImax1,ixImin2:ixImax2) = x(ixImin1:ixImax1,&
       ixImin2:ixImax2,2)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,1) = zero
    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
       2) = -const_G*mass/radius(ixImin1:ixImax1,&
       ixImin2:ixImax2)**2*(unit_time**2/unit_length)

  end subroutine set_gravitation_field

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

    double precision :: wCCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),v(ixImin1:ixImax1,&
       ixImin2:ixImax2,2)
    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2)
    integer :: rdir, pdir

    rdir = 2
    pdir = 1

    v(ixImin1:ixImax1,ixImin2:ixImax2,1) = wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(1))/wCT(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    v(ixImin1:ixImax1,ixImin2:ixImax2,2) = wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(2))/wCT(ixImin1:ixImax1,ixImin2:ixImax2,rho_)

    radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_) - qdt*two*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(rdir))/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    !> dm_r/dt = +rho*v_p**2/r -rho*v_p**2/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(rdir)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(rdir)) + qdt*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       rho_)*v(ixImin1:ixImax1,ixImin2:ixImax2,pdir)**two/x(ixImin1:ixImax1,&
       ixImin2:ixImax2,rdir) - qdt*2*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       rho_)*v(ixImin1:ixImax1,ixImin2:ixImax2,rdir)**two/x(ixImin1:ixImax1,&
       ixImin2:ixImax2,rdir)
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(pdir)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(pdir)) - qdt*3*v(ixImin1:ixImax1,ixImin2:ixImax2,&
       rdir)*v(ixImin1:ixImax1,ixImin2:ixImax2,pdir)*wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_)/x(ixImin1:ixImax1,ixImin2:ixImax2,rdir)


    !> de/dt = -2 (e+p)v_r/r
    call phys_get_pthermal(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,pth)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_) - qdt*two*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_)+pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*wCT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(rdir))/(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)*radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2))

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      wCCT = wCT
      call get_rad_extravars(wCCT, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,r_e) - qdt*two*wCCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         i_flux(rdir))/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    endif

    if (rhd_radiation_advection) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,r_e) - qdt*two*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         r_e)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(rdir))/(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)*radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    endif

  end subroutine PseudoPlanar

  subroutine Opacity_stepfunction(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x,kappa)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = kappa_0

  end subroutine Opacity_stepfunction

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
    integer                            :: idim
    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: mass

    radius(ixImin1:ixImax1,ixImin2:ixImax2) = x(ixImin1:ixImax1,&
       ixImin2:ixImax2,2)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    call get_rad_extravars(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2)

    g_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       i_op)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_flux(2))/fld_speedofligt_0
    g_grav(ixImin1:ixImax1,ixImin2:ixImax2) = &
       const_G*mass/radius(ixImin1:ixImax1,&
       ixImin2:ixImax2)**2*(unit_time**2/unit_length)
    big_gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = g_rad(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/g_grav(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    call rhd_get_tgas(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, Tgas)
    call rhd_get_trad(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, Trad)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1) = Tgas(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*unit_temperature
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2) = Trad(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*unit_temperature
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3) = big_gamma(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+4) = 4*dpi*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(2))*radius(ixImin1:ixImax1,&
       ixImin2:ixImax2)**2 *unit_density*unit_velocity/M_sun*year
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+5) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,i_flux(2)) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       r_e)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2))/w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_)


  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'Tgas Trad Gamma Mdot Ftot'
  end subroutine specialvarnames_output

end module mod_usr
