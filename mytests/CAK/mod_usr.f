!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_HD

  implicit none

  ! Custom variables can be defined here
  double precision :: M_sun = 1.989d33
  double precision :: L_sun = 3.827d33
  double precision :: R_sun = 6.96d10

  double precision :: G_dp = 6.67d-8

  double precision :: M_star
  double precision :: L_star
  double precision :: T_star
  double precision :: R_star
  double precision :: B_star

  double precision :: rho_bound
  double precision :: c_sound
  double precision :: c_light

  double precision :: alpha
  double precision :: qbar
  double precision :: kappa_e
  double precision :: beta

  double precision :: Gamma_e
  double precision :: typical_speed
  double precision :: escape_speed
  double precision :: M_dot
  double precision :: eta_confine
  double precision :: v_inf

  double precision :: new_timestep

  integer :: i_g_cak, i_g_eff, i_f_fd

  integer :: dummy_it
  logical :: first_it, second_it, third_it

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("cylindrical")

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Special Boundary conditions
    usr_special_bc => special_bound

    ! CAK force
    usr_source => CAK_source

    ! adjusted timestep
    usr_get_dt => special_dt

    call usr_params_read(par_files)
    call initglobaldata_usr

    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = 1.d0 ! cm
    unit_temperature   = 1.d0 ! K
    unit_numberdensity = 1.d0 ! cm-3

    ! Active the physics module
    call HD_activate()

    i_g_cak = var_set_extravar("g_cak", "g_cak")
    i_g_eff = var_set_extravar("g_eff", "g_eff")
    i_f_fd = var_set_extravar("f_fd", "f_fd")

  end subroutine usr_init
  !==========================================================================================

  subroutine usr_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /star_list/ M_star, L_star, R_star, T_star, rho_bound, alpha,&
        qbar, kappa_e, beta

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       !open(unitpar, file=files, status="old")
       read(unitpar, star_list, end=111)
       111    close(unitpar)
    end do

    M_star = M_star*M_sun
    L_star = L_star*L_sun
    R_star = R_star*R_sun
    T_star = T_star

    print*, "In CGS: #####################################################"

    print*, 'M_star',M_star
    print*, 'L_star',L_star
    print*, 'R_star',R_star
    print*, 'T_star',T_star
    print*, 'rho_bound',rho_bound
    print*, 'alpha',alpha
    print*, 'qbar',qbar
    print*, 'kappa_e',kappa_e
    print*, 'beta',beta

  end subroutine usr_params_read
  !==========================================================================================

  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_constants

    c_sound = dsqrt(T_star*kB_cgs/(0.6*mp_cgs))

    Gamma_e = kappa_e*L_star/(4.d0*dpi*G_dp*M_star*const_c)
    escape_speed = (two*G_dp*M_star/R_star)**0.5
    M_dot = L_star/const_c**2.d0 * alpha/(1-alpha) &
       *(qbar*Gamma_e/(1.-Gamma_e))**((1.-alpha)/alpha)
    v_inf = (one-Gamma_e)**0.5d0*escape_speed*(alpha/(1-alpha))**0.5d0

    !###############################################

    unit_length        = R_star
    unit_numberdensity = rho_bound/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_temperature   = T_star

    unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
    unit_pressure=(2.d0+3.d0*He_abundance)&
       *unit_numberdensity*const_kB*unit_temperature
    unit_velocity=dsqrt(unit_pressure/unit_density)
    unit_time=unit_length/unit_velocity

    !###############################################

    rho_bound = rho_bound/unit_density
    c_sound = c_sound/unit_velocity
    escape_speed = escape_speed/unit_velocity
    M_dot = M_dot/(unit_density/unit_time*unit_length**3.d0)
    M_star = M_star/(unit_density*unit_length**3.d0)
    R_star = R_star/unit_length
    T_star = T_star/unit_temperature
    L_star = L_star/(unit_pressure/unit_time*unit_length**3.d0)
    kappa_e = kappa_e/(one/(unit_density*unit_length))
    v_inf = (one-Gamma_e)**0.5d0*escape_speed*(alpha/(1-alpha))**0.5d0

    c_light = const_c/unit_velocity
    G_dp = G_dp*(unit_density*unit_time**2.d0)

    new_timestep = dt

    print*, "UNITS: #####################################################"

    print*, 'unit_time', unit_time
    print*, 'unit_length', unit_length
    print*, 'unit_density', unit_density
    print*, 'unit_pressure', unit_pressure
    print*, 'unit_velocity', unit_velocity
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_numberdensity', unit_numberdensity

    print*, "Dimension-LESS: #####################################################"

    print*, 'M_star',M_star
    print*, 'L_star',L_star
    print*, 'R_star',R_star
    print*, 'T_star',T_star
    print*, 'rho_bound',rho_bound
    print*, 'alpha',alpha
    print*, 'qbar',qbar
    print*, 'kappa_e',kappa_e
    print*, 'Gamma_e',Gamma_e
    print*, 'beta',beta
    print*, 'c_sound',c_sound
    print*, 'escape_speed',escape_speed
    print*, 'v_inf',v_inf
    print*, 'eta_confine', eta_confine

    dummy_it = -1
    first_it = .true.
    second_it = .true.
    third_it = .true.

  end subroutine initglobaldata_usr


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision :: velocity_field(ixImin1:ixImax1),pth(ixImin1:ixImax1)
    double precision :: pert(ixImin1:ixImax1), amplitude, xi(ixImin1:ixImax1)
    integer :: i

    do i = ixImin1,ixImax1
      if (x(i,1) .ge. R_star) then
        !Set initial velocity acc to beta law
        velocity_field(i) = v_inf*(1-(R_star/x(i,1)))**beta

        !Set initial density
        w(i,rho_) = M_dot/(4*dpi*velocity_field(i)*x(i,1)**2)

        !set momentum
        w(i,mom(1)) = w(i,rho_)*velocity_field(i)
      else
        w(i,rho_) = rho_bound
        w(i,mom(1)) = zero
      endif
    end do

    w(ixImin1:ixImax1,i_g_cak) = zero
    w(ixImin1:ixImax1,i_g_eff) = zero
    w(ixImin1:ixImax1,i_f_fd) = zero

  end subroutine initial_conditions
  !==========================================================================================

  subroutine special_bound(qt,ixGmin1,ixGmax1,ixBmin1,ixBmax1,iB,w,x)

    use mod_global_parameters
    use mod_variables
    use mod_constants
    use mod_physics, only: phys_get_pthermal

    integer, intent(in) :: ixGmin1,ixGmax1, ixBmin1,ixBmax1, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)
    double precision :: pth(ixGmin1:ixGmax1)
    integer :: i,j

    select case (iB)

    case(1)
      w(ixBmin1:ixBmax1,rho_) = rho_bound

      do i = ixBmax1,ixBmin1,-1
        w(i,mom(1)) = w(i+1,mom(1))-(x(i+1,1)-x(i,1))*(w(i+2,mom(1)) - w(i+1,&
           mom(1)))/(x(i+2,1)-x(i+1,1))
      enddo

      do i = ixBmin1,ixBmax1
        w(i,mom(1)) = min(w(i,mom(1)),one)
        w(i,mom(1)) = max(w(i,mom(1)),-one)
      enddo


    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

!==========================================================================================

  subroutine CAK_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,&
     wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,&
       iwmax
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision :: grad_fwd(ixImin1:ixImax1),grad_bwd(ixImin1:ixImax1),&
       grad_centCT(ixImin1:ixImax1)
    double precision :: g_CAK(ixImin1:ixImax1), g_grav_sc(ixImin1:ixImax1),&
        ppp(300)
    double precision :: vel(ixImin1:ixImax1)
    double precision :: fin_disc(ixImin1:ixImax1), mu_star(ixImin1:ixImax1),&
        sigma_star(ixImin1:ixImax1)
    double precision :: beta_fd(ixImin1:ixImax1), F_fd(ixImin1:ixImax1)
    double precision :: dum(ixOmin1:ixOmax1)


    integer :: output_step
    character(len=11) :: filename
    character(len=4) :: it_nmb
    logical :: give_output



    integer :: i,j
    integer :: jxmin1,jxmax1, hxmin1,hxmax1

    !-----------------------------------------------------------------------------
    jxmin1=ixOmin1+kr(1,1);jxmax1=ixOmax1+kr(1,1);
    hxmin1=ixOmin1-kr(1,1);hxmax1=ixOmax1-kr(1,1);
    !--------------------------------------------------------------------------------------

    vel(ixImin1:ixImax1) = wCT(ixImin1:ixImax1,mom(1))/wCT(ixImin1:ixImax1,&
       rho_)
    !--------------------------------------------------------------------------------------
    ! calculate grad_v
    ! forward difference
    grad_fwd(ixOmin1:ixOmax1) = (vel(jxmin1:jxmax1)-&
       vel(ixOmin1:ixOmax1))/(x(jxmin1:jxmax1,1)-x(ixOmin1:ixOmax1,1))
    ! backward difference
    grad_bwd(ixOmin1:ixOmax1) = (vel(ixOmin1:ixOmax1)-&
       vel(hxmin1:hxmax1))/(x(ixOmin1:ixOmax1,1)-x(hxmin1:hxmax1,1))
    ! central difference
    grad_centCT(ixOmin1:ixOmax1) = half*abs(grad_fwd(ixOmin1:ixOmax1)+&
       grad_bwd(ixOmin1:ixOmax1))
    !--------------------------------------------------------------------------------------

    beta_fd(ixOmin1:ixOmax1) = (1 - (vel(ixOmin1:ixOmax1)/x(ixOmin1:ixOmax1,&
       1))/grad_centCT(ixOmin1:ixOmax1))*(R_star/x(ixOmin1:ixOmax1,1))**2
    do i=ixOmin1,ixOmax1
        if (beta_fd(i) .ge. 1.) then
          F_fd(i) = 1./(1+alpha)
        else if (beta_fd(i).lt.-1.d10) then
          F_fd(i) = abs(beta)**alpha/(1+alpha)
        else if (abs(beta_fd(i)).gt.1.d-3) then
          F_fd(i) = (1.-(1.-beta_fd(i))**(1+alpha))/(beta_fd(i)*(1+alpha))
        else
          F_fd(i) = 1.-0.5*alpha*beta_fd(i)*(1.+&
             0.333333*(1.-alpha)*beta_fd(i))
        end if
        if (F_fd(i) .lt. smalldouble) F_fd(i) = one
        if (F_fd(i) .gt. 5.d0) F_fd(i) = one
    enddo

    w(ixOmin1:ixOmax1,i_f_fd) = F_fd(ixOmin1:ixOmax1)

    !--------------------------------------------------------------------------------------

    ! calculate g_CAK
    g_CAK(ixImin1:ixImax1) = 1./(1.-alpha)*kappa_e*L_star*qbar/(4.*dpi*x(&
       ixImin1:ixImax1,1)**2*c_light) *(grad_centCT(ixImin1:ixImax1)/(wCT(&
       ixImin1:ixImax1,rho_)*c_light*qbar*kappa_e))**alpha

    w(ixImin1:ixImax1,i_g_cak) = g_CAK(ixImin1:ixImax1)

    ! finite disk correction
    g_CAK(ixImin1:ixImax1) = g_CAK(ixImin1:ixImax1) * F_fd(ixImin1:ixImax1)

    ! effective gravity
    g_grav_sc(ixImin1:ixImax1) = G_dp*M_star/(x(ixImin1:ixImax1,&
       1)**2)*(one-Gamma_e)

    w(ixImin1:ixImax1,i_g_eff) = g_grav_sc(ixImin1:ixImax1)

    !w = w + qdt*gsource
    w(ixOmin1:ixOmax1,mom(1)) = w(ixOmin1:ixOmax1,&
       mom(1)) + qdt * (g_CAK(ixOmin1:ixOmax1) - g_grav_sc(ixOmin1:ixOmax1)) * &
       wCT(ixOmin1:ixOmax1,rho_)

    dum = (x(jxmin1:jxmax1,1)-x(ixOmin1:ixOmax1,&
       1))/(g_CAK(ixOmin1:ixOmax1) - g_grav_sc(ixOmin1:ixOmax1))
    new_timestep = 0.3* minval( abs(dum))**half

  end subroutine CAK_source
  !==========================================================================================

  subroutine special_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)

    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: dx1, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(inout) :: dtnew

    if (it .ge. 1) dtnew = new_timestep

  end subroutine special_dt
!==========================================================================================


end module mod_usr
