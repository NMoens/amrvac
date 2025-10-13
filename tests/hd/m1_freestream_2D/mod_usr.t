!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd

  implicit none

  double precision :: L, rho, Er, F_red, Bg_E

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    ! usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions
    usr_special_bc => boundary_conditions

    ! Choose independent normalization units if using dimensionless variables.
    ! unit_length      = 1.d0 ! cm
    ! unit_temperature = 1.d0 ! K
    ! unit_density     = 1.d0 ! g cm^-3

    unit_length      = L    ! cm
    unit_temperature = (Er/const_rad_a)**0.25d0  ! K
    unit_density     = rho ! g cm^-3


    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /stream_list/ L, rho, Er, F_red, Bg_E

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, stream_list, end=111)
       111 close(unitpar)
    enddo

  end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

      w(ixI^S,rho_) = rho/unit_density
      w(ixI^S,mom(:)) = 0.d0
      w(ixI^S,Er_) = Bg_E*Er/unit_pressure
      w(ixI^S,Fr(:)) = 0.d0

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
  use mod_global_parameters
  use mod_hd

  integer, intent(in)             :: ixI^L, ixB^L, iB
  double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
  double precision, intent(inout) :: w(ixI^S,1:nw)

  integer :: ii

  select case (iB)

  case(1)

    select case (iprob)

    case(1)
      where (x(ixB^S,2) .gt. -0.25d0 .and. x(ixB^S,2) .lt. 0.25d0 )
        w(ixB^S,Er_) = Er/unit_pressure
        w(ixB^S,Fr(1)) = Er/unit_pressure*m1_c*F_red
      end where
    case(2)
      where (x(ixB^S,2) .lt. 0.25d0 )
        w(ixB^S,Er_) = Er/unit_pressure
        w(ixB^S,Fr(1)) = 1.d0/dsqrt(2.d0)*Er/unit_pressure*m1_c*F_red
        w(ixB^S,Fr(2)) = 1.d0/dsqrt(2.d0)*Er/unit_pressure*m1_c*F_red
      end where
    end select

  case(3)

    select case (iprob)

    case(1)
      where (x(ixB^S,1) .gt. -0.25d0 .and. x(ixB^S,1) .lt. 0.25d0 )
        w(ixB^S,Er_) = Er/unit_pressure
        w(ixB^S,Fr(2)) = Er/unit_pressure*m1_c*F_red
      end where
    case(2)
      where (x(ixB^S,1) .lt. 0.25d0 )
        w(ixB^S,Er_) = Er/unit_pressure
        w(ixB^S,Fr(1)) = 1.d0/dsqrt(2.d0)*Er/unit_pressure*m1_c*F_red
        w(ixB^S,Fr(2)) = 1.d0/dsqrt(2.d0)*Er/unit_pressure*m1_c*F_red
      end where
    end select

  case default
    call mpistop('boundary not known')
  end select
end subroutine boundary_conditions

end module mod_usr
