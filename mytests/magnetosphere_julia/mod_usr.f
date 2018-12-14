!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_MHD

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

  double precision :: mag_mu

  double precision :: new_timestep

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("spherical_2D")

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Special Boundary conditions
    usr_special_bc => special_bound

    ! CAK force
    usr_source => CAK_source

    ! adjusted timestep
    usr_get_dt => special_dt

    ! get grid refined near equatorial
    usr_refine_grid => specialrefine_grid

    ! Extra Output
    usr_aux_output => specialvar_output
    usr_add_aux_names => specialvarnames_output

    call usr_params_read(par_files)
    call initglobaldata_usr


    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = 1.d0 ! cm
    unit_temperature   = 1.d0 ! K
    unit_numberdensity = 1.d0 ! cm-3,cm-3

    ! Active the physics module
    call MHD_activate()

  end subroutine usr_init
  !==========================================================================================

  subroutine usr_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /star_list/ M_star, L_star, R_star, T_star, B_star, rho_bound,&
        alpha, qbar, kappa_e, beta

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
    print*, 'B_star',B_star
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
    eta_confine = (B_star**2 / R_star**2)/(M_dot*v_inf)

    !###############################################

    unit_length        = R_star
    unit_numberdensity = rho_bound/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_temperature   = T_star

    unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
    unit_pressure=(2.d0+3.d0*He_abundance)&
       *unit_numberdensity*const_kB*unit_temperature
    unit_velocity=dsqrt(unit_pressure/unit_density)
    unit_time=unit_length/unit_velocity
    unit_magneticfield=sqrt(4.d0*dpi*unit_pressure)


    !###############################################

    rho_bound = rho_bound/unit_density
    c_sound = c_sound/unit_velocity
    escape_speed = escape_speed/unit_velocity
    M_dot = M_dot/(unit_density/unit_time*unit_length**3.0)
    M_star = M_star/(unit_density*unit_length**3.d0)
    R_star = R_star/unit_length
    T_star = T_star/unit_temperature
    L_star = L_star/(unit_pressure/unit_time*unit_length**3.0)
    kappa_e = kappa_e/(one/(unit_density*unit_length))
    B_star = B_star/unit_magneticfield
    v_inf = (one-Gamma_e)**0.5d0*escape_speed*(alpha/(1-alpha))**0.5d0
    eta_confine = (B_star**2 / R_star**2)/(M_dot*v_inf)

    c_light = const_c/unit_velocity
    G_dp = G_dp/((unit_density*unit_time**2))

    new_timestep = dt

    print*, "UNITS: #####################################################"

    print*, 'unit_time', unit_time
    print*, 'unit_length', unit_length
    print*, 'unit_density', unit_density
    print*, 'unit_pressure', unit_pressure
    print*, 'unit_velocity', unit_velocity
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_magneticfield', unit_magneticfield
    print*, 'unit_numberdensity', unit_numberdensity

    print*, "Dimension-LESS: #####################################################"

    print*, 'M_star',M_star
    print*, 'L_star',L_star
    print*, 'R_star',R_star
    print*, 'T_star',T_star
    print*, 'B_star',B_star
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

  end subroutine initglobaldata_usr


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision :: velocity_field(ixImin1:ixImax1),pth(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: pert(ixImin1:ixImax1,ixImin2:ixImax2), amplitude,&
        xi(ixImin1:ixImax1,ixImin2:ixImax2)
    integer :: i

    do i = ixImin1,ixImax1
      if (x(i,3,1) .ge. R_star) then
        !Set initial velocity acc to beta law
        velocity_field(i) = v_inf*(1-(R_star/x(i,3,1)))**beta

        !Set initial density
        w(i,:, rho_) = M_dot/(4*dpi*velocity_field(i)*x(i,:,1)**2)

        !set momentum
        w(i,:, mom(1)) = w(i,:,rho_)*velocity_field(i)
        w(i,:, mom(2)) = zero
      else
        w(i,:, rho_) = rho_bound
        w(i,:, mom(:)) = zero
      endif
    end do

    !Set initial magnetic field
    w(ixImin1:ixImax1,ixImin2:ixImax2,mag(1)) = B_star*dcos(x(ixImin1:ixImax1,&
       ixImin2:ixImax2,2))*(R_star/x(ixImin1:ixImax1,ixImin2:ixImax2,1))**3.0
    w(ixImin1:ixImax1,ixImin2:ixImax2,mag(2)) = &
       B_star*half*dsin(x(ixImin1:ixImax1,ixImin2:ixImax2,&
       2))*(R_star/x(ixImin1:ixImax1,ixImin2:ixImax2,1))**3.0

    if(mhd_glm) w(ixImin1:ixImax1,ixImin2:ixImax2,psi_)=0.d0

    ! !> Random perturbation on rho
    ! amplitude = 3.d-2
    ! !> perturb rho
    ! call RANDOM_NUMBER(pert)
    ! w(ixI^S, rho_) = w(ixI^S, rho_)*(one + amplitude*(pert(ixI^S)-half))

    !> Perturbation on Momentum
    xi(ixImin1:ixImax1,ixImin2:ixImax2) = half*(one-atan(abs(x(ixImin1:ixImax1,&
       ixImin2:ixImax2,2)-half*dpi/(dpi/16.d0))-one)*0.25d0)
    xi(ixImin1:ixImax1,ixImin2:ixImax2) = sin((x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1)-R_star)/R_star*two*dpi)*xi(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    w(ixImin1:ixImax1,ixImin2:ixImax2, mom(2)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2, mom(2)) + 0.001d0*xi(ixImin1:ixImax1,ixImin2:ixImax2)

  end subroutine initial_conditions
  !==========================================================================================

  subroutine special_bound(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,ixBmin2,&
     ixBmax1,ixBmax2,iB,w,x)

    use mod_global_parameters
    use mod_variables
    use mod_constants
    use mod_physics, only: phys_get_pthermal

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixBmin1,ixBmin2,&
       ixBmax1,ixBmax2, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    double precision :: pth(ixGmin1:ixGmax1,ixGmin2:ixGmax2)
    integer :: i,j

    select case (iB)

    case(1)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_) = rho_bound

      do i = ixBmax1,ixBmin1,-1
        w(i,:,mom(1)) = w(i+1,:,mom(1))-(x(i+1,:,1)-x(i,:,1))*(w(i+2,:,&
           mom(1)) - w(i+1,:,mom(1)))/(x(i+2,:,1)-x(i+1,:,1))
        ! w(i,:,mom(2)) = w(i+1,:,mom(2))-(x(i+1,:,1)-x(i,:,1))&
        ! *(w(i+2,:,mom(2)) - w(i+1,:,mom(2)))/(x(i+2,:,1)-x(i+1,:,1))
      enddo

      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(2)) = zero
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mag(1)) = &
         B_star*dcos(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         2))*(R_star/x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,1))**3.0

      !> Cont grad
      do i=ixBmax1,ixBmin1,-1
        w(i,:,mag(2)) = w(i+1,:,mag(2))-(x(i+1,:,1)-x(i,:,1))*(w(i+2,:,&
           mag(2)) - w(i+1,:,mag(2)))/(x(i+2,:,1)-x(i+1,:,1))
      enddo

      if(mhd_glm) w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,psi_)=0.d0

      do i = ixGmin2,ixGmax2
        do j = ixBmin1,ixBmax1
          w(j,i,mom(1)) = min(w(j,i,mom(1)),one)
          w(j,i,mom(2)) = min(w(j,i,mom(2)),one) !> obsolete
          w(j,i,mom(1)) = max(w(j,i,mom(1)),-one)
          w(j,i,mom(2)) = max(w(j,i,mom(2)),-one) !> obsolete

        enddo
      enddo

    ! case(2)
    !   w(ixB^S,mom(2)) = zero
    !

    ! case(3)
    !   do i=ixBmax2,ixBmin2, -1
    !     w(:,i,rho_) = w(:,i+1,rho_)
    !     w(:,i,mom(1)) = w(:,i+1,mom(1))
    !     w(:,i,mom(2)) = -w(:,i+1,mom(0))
    !     w(:,i,mag(1)) = half*B_star*dsin(x(:,i,2))
    !     w(:,i,mag(2)) = zero
    !   enddo
    !   if(mhd_glm) w(ixB^S,psi_)=0.d0
    !
    ! case(4)
    !   do i=ixBmin2,ixBmax2
    !     w(:,i,rho_) = w(:,i-1,rho_)
    !     w(:,i,mom(1)) = w(:,i-1,mom(1))
    !     w(:,i,mom(2)) = -w(:,i-1,mom(0))
    !     w(:,i,mag(1)) = half*B_star*dsin(x(:,i,2))
    !     w(:,i,mag(2)) = zero
    !   enddo
    !   if(mhd_glm) w(ixB^S,psi_)=0.d0

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

!==========================================================================================

  subroutine CAK_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, iwmin,iwmax
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: grad_fwd(ixImin1:ixImax1,ixImin2:ixImax2),&
       grad_bwd(ixImin1:ixImax1,ixImin2:ixImax2),grad_centCT(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: g_CAK(ixImin1:ixImax1,ixImin2:ixImax2),&
        g_grav_sc(ixImin1:ixImax1,ixImin2:ixImax2), ppp(300)
    double precision :: vel(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: fin_disc(ixImin1:ixImax1,ixImin2:ixImax2),&
        mu_star(ixImin1:ixImax1,ixImin2:ixImax2), sigma_star(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: beta_fd(ixImin1:ixImax1,ixImin2:ixImax2),&
        F_fd(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: dum(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    integer :: i,j
    integer :: jxmin1,jxmin2,jxmax1,jxmax2, hxmin1,hxmin2,hxmax1,hxmax2

    !-----------------------------------------------------------------------------
    jxmin1=ixOmin1+kr(1,1);jxmin2=ixOmin2+kr(1,2);jxmax1=ixOmax1+kr(1,1)
    jxmax2=ixOmax2+kr(1,2);
    hxmin1=ixOmin1-kr(1,1);hxmin2=ixOmin2-kr(1,2);hxmax1=ixOmax1-kr(1,1)
    hxmax2=ixOmax2-kr(1,2);
    !--------------------------------------------------------------------------------------

    vel(ixImin1:ixImax1,ixImin2:ixImax2) = wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(1))/wCT(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    !--------------------------------------------------------------------------------------
    ! calculate grad_v
    ! forward difference
    grad_fwd(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (vel(jxmin1:jxmax1,&
       jxmin2:jxmax2)-vel(ixOmin1:ixOmax1,ixOmin2:ixOmax2))/(x(jxmin1:jxmax1,&
       jxmin2:jxmax2,1)-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))
    ! backward difference
    grad_bwd(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (vel(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)-vel(hxmin1:hxmax1,hxmin2:hxmax2))/(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)-x(hxmin1:hxmax1,hxmin2:hxmax2,1))
    ! central difference
    grad_centCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
       half*abs(grad_fwd(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)+grad_bwd(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    !--------------------------------------------------------------------------------------

    beta_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (1 - (vel(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1))/grad_centCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*(R_star/x(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))**2
    do i=ixImin1,ixImax1
      do j=ixImin2,ixImax2
        if (beta_fd(i,j) .ge. 1.) then
          F_fd(i,j) = 1./(1+alpha)
        else if (beta_fd(i,j).lt.-1.d10) then
          F_fd(i,j) = abs(beta)**alpha/(1+alpha)
        else if (abs(beta_fd(i,j)).gt.1.d-3) then
          F_fd(i,j) = (1.-(1.-beta_fd(i,j))**(1+alpha))/(beta_fd(i,&
             j)*(1+alpha))
        else
          F_fd(i,j) = 1.-0.5*alpha*beta_fd(i,&
             j)*(1.+0.333333*(1.-alpha)*beta_fd(i,j))
        end if
        if (F_fd(i,j) .lt. smalldouble) F_fd(i,j) = one
        if (F_fd(i,j) .gt. 5.d0) F_fd(i,j) = one
      enddo
    enddo

    !--------------------------------------------------------------------------------------

    ! calculate g_CAK
    ! this is WHITOUT GEOMETRICAL CORRECTION!!!!!
    g_CAK(ixImin1:ixImax1,ixImin2:ixImax2) = &
       1./(1.-alpha)*kappa_e*L_star*qbar/(4.*dpi*x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1)**2*c_light) *(grad_centCT(ixImin1:ixImax1,&
       ixImin2:ixImax2)/(wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       rho_)*c_light*qbar*kappa_e))**alpha

    g_grav_sc(ixImin1:ixImax1,ixImin2:ixImax2) = &
       G_dp*M_star/(x(ixImin1:ixImax1,ixImin2:ixImax2,1)**2)*(one-Gamma_e)

    !FINITE DISC CORRECTION
    g_CAK(ixImin1:ixImax1,ixImin2:ixImax2) = g_CAK(ixImin1:ixImax1,&
       ixImin2:ixImax2) * F_fd(ixImin1:ixImax1,ixImin2:ixImax2)

    !w = w + qdt*gsource
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(1)) + qdt * (g_CAK(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) - g_grav_sc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    dum = (x(jxmin1:jxmax1,jxmin2:jxmax2,1)-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1))/(g_CAK(ixOmin1:ixOmax1,ixOmin2:ixOmax2) - g_grav_sc(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))
    new_timestep = 0.3* minval( abs(dum))**half

  end subroutine CAK_source
  !==========================================================================================

  subroutine special_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,dtnew,dx1,dx2,x)

    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: dx1,dx2, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: dtnew

    if (it .ge. 1) dtnew = new_timestep
    ! if (it .ge. 1504001+1) dtnew = new_timestep
    ! if (it .ge. 36838+1) dtnew = new_timestep

  end subroutine special_dt
  !==========================================================================================

  subroutine specialrefine_grid(igrid,level,ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
     ixmin1,ixmin2,ixmax1,ixmax2,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
  ! you must set consistent values for integers refine/coarsen:
  ! refine = -1 enforce to not refine
  ! refine =  0 doesn't enforce anything
  ! refine =  1 enforce refinement
  ! coarsen = -1 enforce to not coarsen
  ! coarsen =  0 doesn't enforce anything
  ! coarsen =  1 enforce coarsen
  integer, intent(in) :: igrid, level, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
     ixmin2,ixmax1,ixmax2
  double precision, intent(in) :: qt, w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw),&
      x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
  integer, intent(inout) :: refine, coarsen

  double precision :: ll, l, h, hh

  ll = (0.5d0-1.d0/16.d0)*dpi
  l = (0.5d0-1.d0/8.d0)*dpi
  h = (0.5d0+1.d0/8.d0)*dpi
  hh = (0.5d0+1.d0/16.d0)*dpi

  ! test with different levels of refinement enforced
  if (any(x(ixmin1:ixmax1,ixmin2:ixmax2,2) < h) .and. any(x(ixmin1:ixmax1,&
     ixmin2:ixmax2,2) > l)) refine=1
  if (any(x(ixmin1:ixmax1,ixmin2:ixmax2,2) < hh) .and. any(x(ixmin1:ixmax1,&
     ixmin2:ixmax2,2) > ll)) refine=1

end subroutine specialrefine_grid
!==========================================================================================

subroutine specialvar_output(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,w,x,normconv)
  use mod_global_parameters
  use mod_physics
  integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2
  double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim)
  double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)

  double precision :: B_vec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
  double precision :: Div_B(ixImin1:ixImax1,ixImin2:ixImax2)

  B_vec(ixImin1:ixImax1,ixImin2:ixImax2,:) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
     mag(:))
  call divvector(B_vec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,Div_B)

  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1) = Div_B(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)

end subroutine specialvar_output
!==========================================================================================

subroutine specialvarnames_output(varnames)
  use mod_global_parameters
  character(len=*) :: varnames

  varnames = 'Div_B'
end subroutine specialvarnames_output
end module mod_usr
