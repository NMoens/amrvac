!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd

  implicit none

  double precision :: rho1, v1, T1
  double precision :: rho2, v2, T2
  double precision :: m1, m2, p1, p2, e1, e2, Er1, Er2

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions
    usr_special_bc => boundary_conditions

    ! Choose independent normalization units if using dimensionless variables.
    ! unit_length      = 1.d0 ! cm
    ! unit_temperature = 1.d0 ! K
    ! unit_density     = 1.d0 ! g cm^-3

    unit_velocity = v1
    unit_numberdensity = rho1/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = 1.d5

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity
    unit_magneticfield = sqrt(4.0*dpi * unit_pressure)

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /shock_list/ rho1, v1, T1, rho2, v2, T2

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, shock_list, end=111)
       111 close(unitpar)
    enddo

  end subroutine usr_params_read

  subroutine initglobaldata_usr
    use mod_constants

    m1 = rho1*v1
    m2 = rho2*v2
    p1 = const_kb*T1*rho1/(const_mp*0.5d0)
    p2 = const_kb*T2*rho2/(const_mp*0.5d0)
    Er1 = const_rad_a*T1**4
    Er2 = const_rad_a*T2**4

    print*, 'rho', rho1, '|', rho2
    print*, 'm', m1, '|', m2
    print*, 'e', p1/(hd_gamma-1), '|', p2/(hd_gamma-1)
    print*, 'Er', Er1, '|', Er2

    print*, 'unit_density', unit_density
    print*, 'unit_velocity', unit_velocity
    print*, 'unit_pressure', unit_pressure

    print*, 'mu', hd_mu

    rho1 = rho1/unit_density
    v1 = v1/unit_velocity
    p1 = p1/unit_pressure
    Er1 = Er1/unit_pressure

    rho2 = rho2/unit_density
    v2 = v2/unit_velocity
    p2 = p2/unit_pressure
    Er2 = Er2/unit_pressure

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    where (x(ixI^S,1) .lt. 0.d0)
      w(ixI^S,rho_) = rho1
      w(ixI^S,mom(1)) = v1
      w(ixI^S,p_) = p1
      w(ixI^S,Er_) = Er1
      w(ixI^S,Fr(1)) = 0.d0
    elsewhere
      w(ixI^S,rho_) = rho2
      w(ixI^S,mom(1)) = v2
      w(ixI^S,p_) = p2
      w(ixI^S,Er_) = Er2
      w(ixI^S,Fr(1)) = 0.d0
    end where

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
  use mod_global_parameters
  use mod_fld


  integer, intent(in)             :: ixI^L, ixB^L, iB
  double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
  double precision, intent(inout) :: w(ixI^S,1:nw)

  integer :: ii

  select case (iB)

  case(1)
    w(ixB^S,rho_) = rho1
    w(ixB^S,mom(1)) = v1
    w(ixB^S,e_) = p1
    w(ixB^S,Er_) = Er1
    w(ixB^S,Fr(1)) = 0.d0

    call hd_to_conserved(ixI^L,ixB^L,w,x)

  case(2)
    w(ixB^S,rho_) = rho2
    w(ixB^S,mom(1)) = v2
    w(ixB^S,e_) = p2
    w(ixB^S,Er_) = Er2
    w(ixB^S,Fr(1)) = 0.d0

    call hd_to_conserved(ixI^L,ixB^L,w,x)

  case default
    call mpistop('boundary not known')
  end select
end subroutine boundary_conditions

end module mod_usr
