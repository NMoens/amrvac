!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd

  implicit none

  ! Custom variables can be defined here
  double precision :: M_sun = 1.989d33
  double precision :: L_sun = 3.827d33
  double precision :: R_sun = 6.96d10

  double precision :: M_star
  double precision :: T_star
  double precision :: R_star

  double precision :: Omega
  double precision :: rho_bound

  double precision :: GM

  integer :: i_vt, i_vr

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("polar_2D")
    call usr_params_read(par_files)

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    usr_special_bc => boundary_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    usr_gravity => set_gravitation_field

    usr_modify_output => output_routine

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

    i_vt = var_set_extravar("v_t", "v_t")
    i_vr = var_set_extravar("v_r", "v_r")

  end subroutine usr_init

  subroutine usr_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /star_list/ M_star, R_star, T_star, rho_bound, Omega

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       !open(unitpar, file=files, status="old")
       read(unitpar, star_list, end=111)
       111    close(unitpar)
    end do

    print*, M_star
    print*, R_star
    print*, T_star
    print*, rho_bound
    print*, Omega

    M_star = M_star*M_sun
    R_star = R_star*R_sun

    GM = const_G*M_star

  end subroutine usr_params_read

  subroutine initglobaldata_usr
    use mod_global_parameters

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = R_star !r_arr(nghostcells) ! cm
    unit_numberdensity = rho_bound/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_temperature = T_star

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
    unit_pressure=(2.d0+3.d0*He_abundance)&
       *unit_numberdensity*kB_cgs*unit_temperature
    unit_velocity=dsqrt(unit_pressure/unit_density)
    unit_time=unit_length/unit_velocity

    M_star = M_star*unit_length**3/unit_density
    R_star = R_star/unit_length
    T_star = T_star/unit_temperature
    rho_bound = rho_bound/unit_density

    GM = GM*unit_time**2/unit_length**3

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: H_eff(ixImin1:ixImax1,ixImin2:ixImax2),&
        v_t(ixImin1:ixImax1,ixImin2:ixImax2)

    v_t(ixImin1:ixImax1,ixImin2:ixImax2) = Omega*dsqrt(GM/x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1))
    H_eff(ixImin1:ixImax1,ixImin2:ixImax2) = 1.d0/(GM/x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1)**2-v_t(ixImin1:ixImax1,&
       ixImin2:ixImax2)**2/x(ixImin1:ixImax1,ixImin2:ixImax2,1))

    ! cell-center density
    ! w(ixI^S, rho_) = rho_bound!*dexp(-x(ixI^S,1)/H_eff(ixI^S))
    w(ixImin1:ixImax1,ixImin2:ixImax2, rho_) = &
       rho_bound*dexp(-(one/x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1)-one)*GM*(Omega**2-one))
    ! cell-center momentum
    w(ixImin1:ixImax1,ixImin2:ixImax2, mom(1)) = 0.d0
    w(ixImin1:ixImax1,ixImin2:ixImax2, mom(2)) = v_t(ixImin1:ixImax1,&
       ixImin2:ixImax2)*w(ixImin1:ixImax1,ixImin2:ixImax2, rho_)

    ! print*, GM/x(1:6,5,1)**2
    ! print*, v_t(1:6,5)**2/x(1:6,5,1)
    ! print*, H_eff(1:6,5)
    ! print*, x(1:6,5,1)/H_eff(1:6,5)
    ! print*, x(1:6,5,1)/H_eff(1:6,5)
    ! print*, w(1:6,5,rho_)

    w(ixImin1:ixImax1,ixImin2:ixImax2,i_vt) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2, mom(2))/w(ixImin1:ixImax1,ixImin2:ixImax2, rho_)
    w(ixImin1:ixImax1,ixImin2:ixImax2,i_vr) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2, mom(1))/w(ixImin1:ixImax1,ixImin2:ixImax2, rho_)
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

    double precision :: H_eff(ixBmin1:ixBmax1,ixBmin2:ixBmax2),&
        v_t(ixBmin1:ixBmax1,ixBmin2:ixBmax2)

    select case (iB)

    case(1)

    v_t(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = Omega*dsqrt(GM/x(ixBmin1:ixBmax1,&
       ixBmin2:ixBmax2,1))
    H_eff(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = 1.d0/(GM/x(ixBmin1:ixBmax1,&
       ixBmin2:ixBmax2,1)**2 - v_t(ixBmin1:ixBmax1,&
       ixBmin2:ixBmax2)**2/x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,1))

    ! cell-center density
    ! w(ixB^S, rho_) = rho_bound*dexp(-x(ixB^S,1)/H_eff(ixB^S))
    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2, rho_) = &
       rho_bound*dexp(-(one/x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
       1)-one)*GM*(Omega**2-one))
    ! cell-center velocity
    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2, mom(1)) = 0.d0
    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2, mom(2)) = v_t(ixBmin1:ixBmax1,&
       ixBmin2:ixBmax2)*w(ixBmin1:ixBmax1,ixBmin2:ixBmax2, rho_)

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions


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

    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,:) = zero
    gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2) = -GM/x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)**2

  end subroutine set_gravitation_field
  ! Extra routines can be placed here
  ! ...

  subroutine output_routine(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    w(ixImin1:ixImax1,ixImin2:ixImax2,i_vt) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2, mom(2))/w(ixImin1:ixImax1,ixImin2:ixImax2, rho_)
    w(ixImin1:ixImax1,ixImin2:ixImax2,i_vr) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2, mom(1))/w(ixImin1:ixImax1,ixImin2:ixImax2, rho_)

  end subroutine output_routine

end module mod_usr
