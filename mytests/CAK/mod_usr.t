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
    unit_numberdensity = 1.d0 ! cm^-3

    ! Active the physics module
    call HD_activate()

  end subroutine usr_init
  !==========================================================================================

  subroutine usr_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /star_list/ M_star, L_star, R_star, T_star, rho_bound, alpha, qbar, kappa_e, beta

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
    unit_pressure=(2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB*unit_temperature
    unit_velocity=dsqrt(unit_pressure/unit_density)
    unit_time=unit_length/unit_velocity


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

  end subroutine initglobaldata_usr


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: velocity_field(ixImin1:ixImax1),pth(ixI^S)
    double precision :: pert(ixI^S), amplitude, xi(ixI^S)
    integer :: i

    do i = ixImin1,ixImax1
      if (x(i,3,1) .ge. R_star) then
        !Set initial velocity acc to beta law
        velocity_field(i) = v_inf*(1-(R_star/x(i,3,1)))**beta

        !Set initial density
        w(i,:, rho_) = M_dot/(4*dpi*velocity_field(i)*x(i,:,1)**2)

        !set momentum
        w(i,:, mom(1)) = w(i,:,rho_)*velocity_field(i)
      else
        w(i,:, rho_) = rho_bound
        w(i,:, mom(:)) = zero
      endif
    end do

  end subroutine initial_conditions
  !==========================================================================================

  subroutine special_bound(qt,ixG^L,ixB^L,iB,w,x)

    use mod_global_parameters
    use mod_variables
    use mod_constants
    use mod_physics, only: phys_get_pthermal

    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision :: pth(ixG^S)
    integer :: i,j

    select case (iB)

    case(1)
      w(ixB^S,rho_) = rho_bound

      do i = ixBmax1,ixBmin1,-1
        w(i,:,mom(1)) = w(i+1,:,mom(1))-(x(i+1,:,1)-x(i,:,1))&
        *(w(i+2,:,mom(1)) - w(i+1,:,mom(1)))/(x(i+2,:,1)-x(i+1,:,1))
        ! w(i,:,mom(2)) = w(i+1,:,mom(2))-(x(i+1,:,1)-x(i,:,1))&
        ! *(w(i+2,:,mom(2)) - w(i+1,:,mom(2)))/(x(i+2,:,1)-x(i+1,:,1))
      enddo

      do i = ixGmin2,ixGmax2
        do j = ixBmin1,ixBmax1
          w(j,i,mom(1)) = min(w(j,i,mom(1)),one)
          w(j,i,mom(1)) = max(w(j,i,mom(1)),-one)
        enddo
      enddo

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

!==========================================================================================

  subroutine CAK_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: grad_fwd(ixI^S),grad_bwd(ixI^S),grad_centCT(ixI^S)
    double precision :: g_CAK(ixI^S), g_grav_sc(ixI^S), ppp(300)
    double precision :: vel(ixI^S)
    double precision :: fin_disc(ixI^S), mu_star(ixI^S), sigma_star(ixI^S)
    double precision :: beta_fd(ixI^S), F_fd(ixI^S)
    double precision :: dum(ixO^S)

    integer :: i,j
    integer :: jx^L, hx^L

    !-----------------------------------------------------------------------------
    jx^L=ixO^L+kr(1,^D);
    hx^L=ixO^L-kr(1,^D);
    !--------------------------------------------------------------------------------------

    vel(ixI^S) = wCT(ixI^S,mom(1))/wCT(ixI^S,rho_)
    !--------------------------------------------------------------------------------------
    ! calculate grad_v
    ! forward difference
    grad_fwd(ixO^S) = (vel(jx^S)-vel(ixO^S))/(x(jx^S,1)-x(ixO^S,1))
    ! backward difference
    grad_bwd(ixO^S) = (vel(ixO^S)-vel(hx^S))/(x(ixO^S,1)-x(hx^S,1))
    ! central difference
    grad_centCT(ixO^S) = half*max(grad_fwd(ixO^S)+grad_bwd(ixO^S),zero) !> FIX SYNTAX HERE
    !--------------------------------------------------------------------------------------

    beta_fd(ixO^S) = (1 - (vel(ixO^S)/x(ixO^S,1))/grad_centCT(ixO^S))*(R_star/x(ixO^S,1))**2
    do i=ixImin1,ixImax1
      do j=ixImin2,ixImax2
        if (beta_fd(i,j) .ge. 1.) then
          F_fd(i,j) = 1./(1+alpha)
        else if (beta_fd(i,j).lt.-1.d10) then
          F_fd(i,j) = abs(beta)**alpha/(1+alpha)
        else if (abs(beta_fd(i,j)).gt.1.d-3) then
          F_fd(i,j) = (1.-(1.-beta_fd(i,j))**(1+alpha))/(beta_fd(i,j)*(1+alpha))
        else
          F_fd(i,j) = 1.-0.5*alpha*beta_fd(i,j)*(1.+0.333333*(1.-alpha)*beta_fd(i,j))
        end if
        if (F_fd(i,j) .lt. smalldouble) F_fd(i,j) = one
        if (F_fd(i,j) .gt. 5.d0) F_fd(i,j) = one
      enddo
    enddo

    !--------------------------------------------------------------------------------------

    ! calculate g_CAK
    ! this is WHITOUT GEOMETRICAL CORRECTION!!!!!
    g_CAK(ixI^S) = 1./(1.-alpha)*kappa_e*L_star*qbar/(4.*dpi*x(ixI^S,1)**2*c_light) &
    *(grad_centCT(ixI^S)/(wCT(ixI^S,rho_)*c_light*qbar*kappa_e))**alpha

    g_grav_sc(ixI^S) = G_dp*M_star/(x(ixI^S,1)**2)*(one-Gamma_e)

    !FINITE DISC CORRECTION
    g_CAK(ixI^S) = g_CAK(ixI^S) * F_fd(ixI^S)

    !w = w + qdt*gsource
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt * (g_CAK(ixO^S) - g_grav_sc(ixO^S)) * wCT(ixO^S,rho_)

    dum = (x(jx^S,1)-x(ixO^S,1))/(g_CAK(ixO^S) - g_grav_sc(ixO^S))
    new_timestep = 0.3* minval( abs(dum))**half

  end subroutine CAK_source
  !==========================================================================================

  subroutine special_dt(w,ixI^L,ixO^L,dtnew,dx^D,x)

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    if (it .ge. 1) dtnew = new_timestep
    ! if (it .ge. 1504001+1) dtnew = new_timestep
    ! if (it .ge. 36838+1) dtnew = new_timestep

  end subroutine special_dt

!===========================================================================================
end module mod_usr
