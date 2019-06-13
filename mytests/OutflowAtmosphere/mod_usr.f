!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd

  implicit none

  double precision, parameter :: M_sun = 1.99d33
  double precision, parameter :: R_sun = 6.96d10
  double precision, parameter :: year = 365.25*24*60*60

  double precision, allocatable :: r_arr
  double precision, allocatable :: rho_arr
  double precision, allocatable :: v_arr
  double precision, allocatable :: e_arr
  double precision, allocatable :: Er_arr

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

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters

    !> Set stellar mass and radius
    call ReadInParams(M_star,R_star)

    M_star = 50.d0*M_sun
    R_star = 20.d0*R_sun

    R_b = R_star
    R_0 = 1.2d0*R_star

    !> Select mass loss parameters, this determines Luminosity
    M_dot_ratio = 0.5d0
    M_dot = 1d-6

    L_0 = M_dot*const_G*M_star/(M_dot_ratio*R_0)

    !> Set Gamma ratios for base and outer wind
    Gamma_0 = 3.d0
    Gamma_b = 0.8d0

    kappa_0 = Gamma_0*4*dpi*const_G*M_star*const_c/L_0
    kappa_b = Gamma_b*4*dpi*const_G*M_star*const_c/L_0


    call ReadInTable(r_arr,rho_arr,v_arr,e_arr,Er_arr)

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = R_star ! cm
    unit_velocity   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm-3,cm-3


    if (mype .eq. 0) then
      print*, 'unit_length', unit_length
      print*, 'unit_density', unit_density
      print*, 'unit_pressure', unit_pressure
      print*, 'unit_temperature', unit_temperature
      print*, 'unit_radflux', unit_radflux
      print*, 'unit_opacity', unit_opacity
      print*, 'unit_time', unit_time
      print*, 'unit_velocity', unit_velocity
    endif

  end subroutine initglobaldata_usr

  subroutine ReadInTable(rho_arr,v_arr,e_arr,Er_arr)
    use mod_global_parameters

  end subroutine ReadInTable

  subroutine ReadInTable(rho_arr,v_arr,e_arr,Er_arr)
    use mod_global_parameters

  end subroutine ReadInTable

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    ! Set initial values for w
    ! w(ixO^S, rho_) = ...

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
        ixBmin1,ixBmin2,ixBmax1,ixBmax2, iB
    double precision, intent(in)    :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

    double precision :: x_vac(ixGmin2:ixGmax2)
    double precision :: rho_vac(ixGmin2:ixGmax2)
    double precision :: v_vac(ixGmin2:ixGmax2)
    double precision :: pg_vac(ixGmin2:ixGmax2)
    double precision :: er_vac(ixGmin2:ixGmax2)
    double precision :: temp_vac(ixGmin2:ixGmax2)
    integer :: i,j

    select case (iB)

    case(3)

    case(4)
      do i = ixBmin2,ixBmax2
        !> Conserve gradE/rho
        w(ixGmin1:ixGmax1,i,r_e) = w(ixGmin1:ixGmax1,i-1,&
           rho_)/w(ixGmin1:ixGmax1,i-2,rho_) *(w(ixGmin1:ixGmax1,i-1,&
           r_e) - w(ixGmin1:ixGmax1,i-2,r_e)) + w(ixGmin1:ixGmax1,i-1,r_e)
        do j = ixGmin1,ixGmax1
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

    select case (iB)
      case (3)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = 7.3583819042386223 !7.6447315544263788
      case (4)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
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

    radius(ixImin1:ixImax1,ixImin2:ixImax2) = x(ixImin1:ixImax1,&
       ixImin2:ixImax2,2)*unit_length

    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,1) = zero
    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
       2) = -const_G*mstar/(radius(ixImin1:ixImax1,&
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
      call get_rad_extravars(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,r_e) - qdt*two*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         r_e)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(rdir))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
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

    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = kappa_0/unit_opacity

    where (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .lt. R_0)
      kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = kappa_b/unit_opacity
    endwhere

  end subroutine Opacity_stepfunction

end module mod_usr
