!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision, parameter :: M_sun = 1.99d33
  double precision, parameter :: R_sun = 6.96d10
  double precision, parameter :: L_sun = 3.9d33
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

    !> Set stellar mass and radius
    call ReadInParams(M_star,R_star,Gamma_0,M_dot_ratio,M_dot,L_0)

    R_b = R_star
    R_0 = (1.d0+2.d0)*R_star

    Gamma_b = 0.95d0

    kappa_0 = Gamma_0*4*dpi*const_G*M_star*const_c/L_0
    kappa_b = Gamma_b*4*dpi*const_G*M_star*const_c/L_0

    allocate(r_arr(domain_nx2+2*nghostcells))
    allocate(rho_arr(domain_nx2+2*nghostcells))
    allocate(v_arr(domain_nx2+2*nghostcells))
    allocate(e_arr(domain_nx2+2*nghostcells))
    allocate(Er_arr(domain_nx2+2*nghostcells))
    allocate(T_arr(domain_nx2+2*nghostcells))
    allocate(p_arr(domain_nx2+2*nghostcells))

    call ReadInTable(r_arr,rho_arr,v_arr,e_arr,Er_arr,T_arr, p_arr)

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = R_star !r_arr(nghostcells) ! cm
    unit_temperature   = T_arr(nghostcells)
    unit_numberdensity = rho_arr(nghostcells)/((1.d0+&
       4.d0*He_abundance)*mp_cgs)

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
    unit_pressure=(2.d0+3.d0*He_abundance)&
       *unit_numberdensity*kB_cgs*unit_temperature
    unit_velocity=dsqrt(unit_pressure/unit_density)
    unit_time=unit_length/unit_velocity
    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    if (mype .eq. 0) then
      print*, M_star, R_star, Gamma_0, M_dot_ratio, M_dot, L_0
      print*, R_b, R_0, kappa_b, kappa_0

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
    M_star = M_star/unit_density*unit_length**3.d0
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
    M_dot = M_dot*M_star/year

  end subroutine ReadInParams

  subroutine ReadInTable(r_arr,rho_arr,v_arr,e_arr,Er_arr,T_arr,p_arr)
    use mod_global_parameters
    ! use mod_constants
    ! use mod_fld

    integer :: line
    double precision, intent(out) :: r_arr(domain_nx2+2*nghostcells)
    double precision, intent(out) :: rho_arr(domain_nx2+2*nghostcells)
    double precision, intent(out) :: v_arr(domain_nx2+2*nghostcells)
    double precision, intent(out) :: e_arr(domain_nx2+2*nghostcells)
    double precision, intent(out) :: Er_arr(domain_nx2+2*nghostcells)

    double precision :: i_arr(domain_nx2+2*nghostcells)
    double precision, intent(out) :: T_arr(domain_nx2+2*nghostcells)
    double precision, intent(out) :: p_arr(domain_nx2+2*nghostcells)

    OPEN(2,FILE='InitialConditions/init_struc_amrvac')
    READ(2,*)
    READ(2,*)
    do line = 1,domain_nx2+2*nghostcells
      READ(2,*) i_arr(line), r_arr(line), v_arr(line), rho_arr(line),&
          Er_arr(line)
    enddo
    CLOSE(2)

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
    integer :: i

    do i = ixImin1,ixImax1
      w(i,:,rho_) = rho_arr(:)
      w(i,:,mom(1)) = zero
      w(i,:,mom(2)) = rho_arr(:)*v_arr(:)
      w(i,:,e_) = e_arr(:)
      w(i,:,r_e) = Er_arr(:)
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
    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2)

    integer :: i,j

    select case (iB)

    case(3)
      do i = ixBmax2,ixBmin2,-1
        w(ixImin1:ixImax1,i,rho_) = rho_arr(i)
        w(ixImin1:ixImax1,i,mom(:)) = w(ixImin1:ixImax1,i+1,mom(:))

        a(ixImin1:ixImax1) = L_0/(4.d0*dpi*x(ixImin1:ixImax1,i,2)**2)
        b(ixImin1:ixImax1) = w(ixImin1:ixImax1,i+1,&
           r_e)*fld_speedofligt_0 /(3.d0*(x(ixImin1:ixImax1,i+1,&
           2)-x(ixImin1:ixImax1,i,2))*rho_arr(i)*kappa_b)
        c(ixImin1:ixImax1) =fld_speedofligt_0 /(3.d0*(x(ixImin1:ixImax1,i+1,&
           2)-x(ixImin1:ixImax1,i,2))*rho_arr(i)*kappa_b)
        d(ixImin1:ixImax1) = 4.d0/3.d0*abs(w(ixImin1:ixImax1,i,&
           mom(2))/w(ixImin1:ixImax1,i,rho_))
        w(ixImin1:ixImax1,i,r_e) = (a(ixImin1:ixImax1) + &
           b(ixImin1:ixImax1))/(c(ixImin1:ixImax1) + d(ixImin1:ixImax1))

        Temp(ixImin1:ixImax1,i) = (w(ixImin1:ixImax1,i,&
           r_e)/const_rad_a)**0.25d0
        Temp(ixImin1:ixImax1,i) = const_kb/(fld_mu*const_mp)*Temp(&
           ixImin1:ixImax1,i)*w(ixImin1:ixImax1,i,rho_)
        w(ixImin1:ixImax1,i,e_) = Temp(ixImin1:ixImax1,&
           i)/(rhd_gamma - one)+ half*w(ixImin1:ixImax1,i,&
           mom(2))**2.d0/w(ixImin1:ixImax1,i,rho_)

         print*,it, a(10), b(10), c(10), w(10,i,r_e), w(10,nghostcells+1,r_e)
      enddo

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
        ! /(3.d0*(x(ixImin1:ixImax1,i+1,2)-x(ixImin1:ixImax1,i,2))*w(ixImin1:ixImax1,i,rho_)*kappa_b)
        ! c(ixImin1:ixImax1) =fld_speedofligt_0 &
        ! /(3.d0*(x(ixImin1:ixImax1,i+1,2)-x(ixImin1:ixImax1,i,2))*w(ixImin1:ixImax1,i,rho_)*kappa_b)
        ! d(ixImin1:ixImax1) = 4.d0/3.d0*abs(w(ixImin1:ixImax1,i,mom(2))/w(ixImin1:ixImax1,i,rho_))
        ! mean_RE(ixImin1:ixImax1) = (a(ixImin1:ixImax1) + b(ixImin1:ixImax1))/(c(ixImin1:ixImax1) + d(ixImin1:ixImax1))
        ! print*, sum(mean_RE(ixOmin1:ixOmax1))/(ixOmax1-ixOmin1)


        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = Er_arr(i) !sum(mean_RE(ixOmin1:ixOmax1))/(ixOmax1-ixOmin1)




      case (4)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
      case default
        print *, "Not a standard: ", trim(typeboundary(iw_r_e, iB))
        error stop "You have to set a user-defined boundary method"
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
       ixImin2:ixImax2,ndim)

    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: mass

    radius(ixImin1:ixImax1,ixImin2:ixImax2) = x(ixImin1:ixImax1,&
       ixImin2:ixImax2,2)*unit_length
    mass = M_star*(unit_density/unit_length**3.d0)

    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,1) = zero
    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
       2) = -const_G*mass/(radius(ixImin1:ixImax1,&
       ixImin2:ixImax2))**2*(unit_time**2/unit_length)
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
    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2)
    integer :: rdir, pdir

    rdir = 2
    pdir = 1

    radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_) - qdt*two*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(rdir))/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    !> dm_r/dt = m_phi**2/r
    !> dm_phi/dt = - m_phi m_r/r
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(rdir)) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(rdir)) + qdt*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(pdir)))**two/(radius(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(pdir)) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(pdir)) - qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(rdir))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(pdir))/(radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*wCT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_))

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
         ixOmin2:ixOmax2,r_e) - qdt*two*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         r_e)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(rdir))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_) + wCCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         i_flux(rdir)))/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    else
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

    double precision :: new_x(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    new_x(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       2) - R_0
    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = kappa_b + &
       (kappa_0-kappa_b)*half*(one+erf(new_x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)))

    ! kappa(ixO^S) = kappa_0
    ! where (x(ixO^S,2) .lt. R_0)
    !   kappa(ixO^S) = kappa_b
    ! endwhere

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
       ixImin2:ixImax2,1:ndim), big_gamma(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                   :: g_grav(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                   :: Tgas(ixImin1:ixImax1,&
       ixImin2:ixImax2),Trad(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                            :: idim
    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: mass

    radius(ixImin1:ixImax1,ixImin2:ixImax2) = x(ixImin1:ixImax1,&
       ixImin2:ixImax2,2)*unit_length
    mass = M_star*(unit_density/unit_length**3.d0)

    call get_rad_extravars(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2)

    g_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,i_op)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       i_flux(2))/fld_speedofligt_0
    g_grav(ixImin1:ixImax1,ixImin2:ixImax2) = &
       const_G*mass/(radius(ixImin1:ixImax1,&
       ixImin2:ixImax2))**2*(unit_time**2/unit_length)
    big_gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = g_rad(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)/g_grav(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    print*, w(10,1:7,r_e)
    print*, w(10,1:7,i_flux(2))

    call rhd_get_tgas(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, Tgas)
    call rhd_get_trad(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, Trad)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=Tgas(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*unit_temperature
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=Trad(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*unit_temperature
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3)=big_gamma(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'Tgas Trad Gamma'
  end subroutine specialvarnames_output

end module mod_usr
