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
    unit_numberdensity = 1.d9 ! cm^-3


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
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixImin1:ixImax1,:,rho_) = rho_arr(:)
    w(ixImin1:ixImax1,:,mom(1)) = zero
    w(ixImin1:ixImax1,:,mom(2)) = rho_arr(:)*v_arr(:)
    w(ixImin1:ixImax1,:,e_) = e_arr(:)
    w(ixImin1:ixImax1,:,r_e) = Er_arr(:)

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: i,j

    select case (iB)

    case(3)
      do i = ixBmax2,ixBmin2,-1
        w(ixImin1:ixImax1,i,rho_) = rho_arr(i)
        w(ixImin1:ixImax1,i,mom(:)) = w(ixImin1:ixImax1,i+1,mom(:))
        w(ixImin1:ixImax1,i,e_) = e_arr(i)
        w(ixImin1:ixImax1,i,r_e) = Er_arr(i)
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
  subroutine set_gravitation_field(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    double precision :: radius(ixI^S)

    radius(ixI^S) = x(ixI^S,2)*unit_length

    gravity_field(ixI^S,1) = zero
    gravity_field(ixI^S,2) = -const_G*mstar/(radius(ixI^S))**2*(unit_time**2/unit_length)
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

    double precision :: pth(ixI^S)
    double precision :: radius(ixI^S)
    integer :: rdir, pdir

    rdir = 2
    pdir = 1

    radius(ixO^S) = x(ixO^S,2)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    w(ixO^S,rho_) = w(ixO^S,rho_) - qdt*two*wCT(ixO^S,mom(rdir))/radius(ixO^S)

    !> dm_r/dt = m_phi**2/r
    !> dm_phi/dt = - m_phi m_r/r
    w(ixO^S,mom(rdir)) = w(ixO^S,mom(rdir)) + qdt*(wCT(ixO^S,mom(pdir)))**two/(radius(ixO^S)*wCT(ixO^S,rho_))
    w(ixO^S,mom(pdir)) = w(ixO^S,mom(pdir)) - qdt*wCT(ixO^S,mom(rdir))*wCT(ixO^S,mom(pdir))/(radius(ixO^S)*wCT(ixO^S,rho_))

    !> de/dt = -2 (e+p)v_r/r
    call phys_get_pthermal(wCT,x,ixI^L,ixO^L,pth)
    w(ixO^S,e_) = w(ixO^S,e_) - qdt*two*(wCT(ixO^S,e_)+pth(ixO^S))*wCT(ixO^S,mom(rdir))/(wCT(ixO^S,rho_)*radius(ixO^S))

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      call get_rad_extravars(w, x, ixI^L, ixO^L)
      w(ixO^S,r_e) = w(ixO^S,r_e) - qdt*two*(wCT(ixO^S,r_e)*wCT(ixO^S,mom(rdir))/wCT(ixO^S,rho_) + w(ixO^S,i_flux(rdir)))/radius(ixO^S)
    else
      w(ixO^S,r_e) = w(ixO^S,r_e) - qdt*two*wCT(ixO^S,r_e)*wCT(ixO^S,mom(rdir))/(wCT(ixO^S,rho_)*radius(ixO^S))
    endif

  end subroutine PseudoPlanar

  subroutine Opacity_stepfunction(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    kappa(ixO^S) = kappa_0/unit_opacity

    where (x(ixO^S) .lt. R_0)
      kappa(ixO^S) = kappa_b/unit_opacity
    endwhere

  end subroutine Opacity_stepfunction

end module mod_usr
