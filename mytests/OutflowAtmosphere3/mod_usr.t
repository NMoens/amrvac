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
  double precision :: sp_T
  double precision :: sp_p

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

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    integer :: i

    !> Set stellar mass and radius
    call ReadInParams(M_star,R_star,Gamma_0,M_dot_ratio,M_dot,L_0,sp_rho,sp_sos)

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
    sp_T = fld_mu*mp_cgs/kB_cgs*sp_sos**2
    sp_p = sp_sos**2*sp_rho

    if (mype .eq. 0) then
      print*, 'soundspeed at sonic point', sp_sos
      print*, 'density at sonic point', sp_rho
      print*, 'pressure at sonic point', sp_p
      print*, 'Temperature at sonic point', sp_T
      print*, 'cgs opacity', kappa_0
    endif

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = R_star !r_arr(nghostcells) ! cm
    unit_numberdensity = sp_rho/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_temperature = sp_T

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
    unit_pressure=(2.d0+3.d0*He_abundance)*unit_numberdensity*kB_cgs*unit_temperature
    unit_velocity=dsqrt(unit_pressure/unit_density)
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
    sp_sos = sp_sos/unit_velocity
    sp_rho = sp_rho/unit_density
    sp_T = sp_T/unit_temperature
    sp_p = sp_p/unit_temperature

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
      print*, 'Gamma_0 ', 'kappa_0'
      print*, M_star, R_star,M_dot_ratio, M_dot, L_0
      print*, Gamma_0, kappa_0

      print*, 'Flux at boundary: ', L_0/(4*dpi*R_star**2)
    endif

  end subroutine initglobaldata_usr

  subroutine ReadInParams(M_star,R_star,Gamma_0,M_dot_ratio,M_dot,L_0,sp_rho,sp_sos)
    use mod_global_parameters
    double precision, intent(out) :: M_star,R_star,Gamma_0
    double precision, intent(out) :: M_dot_ratio,M_dot,L_0
    double precision, intent(out) :: sp_sos,sp_rho
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
    READ(1,*) dum, sp_rho
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

    call get_rad_extravars(w, x, ixI^L, ixO^L)

    ! do i = ixImin2, ixImax2
    !   print*, i, x(5,i,2), w(5,i,i_flux(2)), x(5,i,2)**2*w(5,i,i_flux(2))
    ! enddo
    ! stop

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    use mod_fld


    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: a(ixImin1:ixImax1),b(ixImin1:ixImax1),c(ixImin1:ixImax1),d(ixImin1:ixImax1)
    double precision :: Temp(ixI^S), Press(ixI^S), kbTmu(ixI^S)

    integer :: i,j

    select case (iB)

    case(3)

      !++++++++++++++++++++++++++++++++++++++++++++
      i = ixBmax2
        w(ixImin1:ixImax1,i,rho_) = x(ixImin1:ixImax1,nghostcells+1,2)**2/x(ixImin1:ixImax1,i,2)**2*w(ixImin1:ixImax1,nghostcells+1,mom(2))/sp_sos
        w(ixImin1:ixImax1,i,mom(1)) = w(ixImin1:ixImax1,i+1,mom(1))
        w(ixImin1:ixImax1,i,mom(2)) = x(ixImin1:ixImax1,i+1,2)**2/x(ixImin1:ixImax1,i,2)**2*w(ixImin1:ixImax1,i+1,mom(2))
        w(ixImin1:ixImax1,i,e_) = sp_sos**2*w(ixImin1:ixImax1,i,rho_)/(rhd_gamma-one) + half*(w(ixImin1:ixImax1,i,mom(1))**2 &
        + w(ixImin1:ixImax1,i,mom(2))**2)/w(ixImin1:ixImax1,i,rho_)
        w(ixImin1:ixImax1,i,r_e) = const_rad_a*sp_T**4*unit_temperature**4/unit_pressure
      !++++++++++++++++++++++++++++++++++++++++++++

      do i = ixBmax2 -1, ixBmin2, -1
        w(ixImin1:ixImax1,i,rho_) = w(ixImin1:ixImax1,i+1,rho_) - 2*w(ixImin1:ixImax1,i+2,rho_)
        w(ixImin1:ixImax1,i,mom(1)) = w(ixImin1:ixImax1,i+1,mom(1))
        w(ixImin1:ixImax1,i,mom(2)) = x(ixImin1:ixImax1,i+1,2)**2/x(ixImin1:ixImax1,i,2)**2*w(ixImin1:ixImax1,i+1,mom(2))
        do j = ixImin1,ixImax1
          w(j,i,r_e) = max(w(j,i+1,r_e) -2*w(j,i+2,r_e),w(j,i+1,r_e))
        enddo
        kbTmu(ixImin1:ixImax1,i) = (w(ixImin1:ixImax1,i,r_e)*unit_pressure/const_rad_a)**0.25*const_kB/(fld_mu*const_mp)
        kbTmu(ixImin1:ixImax1,i) = kbTmu(ixImin1:ixImax1,i)/unit_velocity**2
        w(ixImin1:ixImax1,i,e_) = kbTmu(ixImin1:ixImax1,i)*w(ixImin1:ixImax1,i,rho_)/(rhd_gamma-one) &
        + half*(w(ixImin1:ixImax1,i,mom(1))**2 + w(ixImin1:ixImax1,i,mom(2))**2)/w(ixImin1:ixImax1,i,rho_)
      enddo

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

    integer :: i

    select case (iB)
      case (3)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous !mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = sum(w(ixImin1:ixImax1,ixOmin2-1,r_e))/(ixImax1-ixImin1)!const_rad_a*sp_T**4*unit_temperature**4/unit_pressure
      case (4)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = sum(w(ixImin1:ixImax1,ixOmax2+1,r_e))/(ixImax1-ixImin1)
      case default
        print *, "Not a standard: ", trim(typeboundary(r_e, iB))
        error stop "You have to set a user-defined boundary method"
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

    kappa(ixO^S) = kappa_0

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

    g_rad(ixO^S) = w(ixO^S,i_op)*w(ixO^S,i_flux(2))/(const_c/unit_velocity)
    g_grav(ixI^S) = const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length)
    big_gamma(ixO^S) = g_rad(ixO^S)/g_grav(ixO^S)

    call rhd_get_tgas(w, x, ixI^L, ixO^L, Tgas)
    call rhd_get_trad(w, x, ixI^L, ixO^L, Trad)

    w(ixO^S,nw+1) = Tgas(ixO^S)*unit_temperature
    w(ixO^S,nw+2) = Trad(ixO^S)*unit_temperature
    w(ixO^S,nw+3) = big_gamma(ixO^S)
    w(ixO^S,nw+4) = 4*dpi*w(ixO^S,mom(2))*radius(ixO^S)**2 &
    *unit_density*unit_velocity/M_sun*year
    w(ixO^S,nw+5) = (w(ixO^S,i_flux(2)) + w(ixO^S,r_e)*w(ixO^S,mom(2))/w(ixO^S,rho_))*radius(ixO^S)**2

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'Tgas Trad Gamma Mdot FtotR2'
  end subroutine specialvarnames_output

end module mod_usr
