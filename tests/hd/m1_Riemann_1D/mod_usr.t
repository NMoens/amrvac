!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd
  use mod_m1

  implicit none

  double precision :: Er_l, Fx_l, Fy_l, Er_r, Fx_r, Fy_r

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_1.5D")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    ! usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions
    ! usr_special_bc => boundary_conditions

    ! Choose independent normalization units if using dimensionless variables.
    unit_length      = 1.d0 ! cm
    unit_velocity    = const_c ! K
    unit_density     = 1.d0 ! g cm^-3

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /shock_list/ Er_l, Fx_l, Fy_l, Er_r, Fx_r, Fy_r

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, shock_list, end=111)
       111 close(unitpar)
    enddo

  end subroutine usr_params_read


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixI^S,rho_) = 1.d0
    w(ixI^S,mom(1)) = 0.d0

    where (x(ixI^S,1) .lt. 0.d0)
      w(ixI^S,Er_) = Er_l
      w(ixI^S,Fr(1)) = Fx_l
      w(ixI^S,Fr(2)) = Fy_l
    elsewhere
      w(ixI^S,Er_) = Er_r
      w(ixI^S,Fr(1)) = Fx_r
      w(ixI^S,Fr(2)) = Fy_r
    end where

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

!   subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
!   use mod_global_parameters
!   use mod_fld
!
!
!   integer, intent(in)             :: ixI^L, ixB^L, iB
!   double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
!   double precision, intent(inout) :: w(ixI^S,1:nw)
!
!   integer :: ii
!
!   select case (iB)
!
!   case(1)
!     w(ixB^S,rho_) = rho1
!     w(ixB^S,mom(1)) = v1
!     w(ixB^S,e_) = p1
!     w(ixB^S,Er_) = Er1
!     w(ixB^S,Fr(1)) = 0.d0
!
!     call hd_to_conserved(ixI^L,ixB^L,w,x)
!
!   case(2)
!     w(ixB^S,rho_) = rho2
!     w(ixB^S,mom(1)) = v2
!     w(ixB^S,e_) = p2
!     w(ixB^S,Er_) = Er2
!     w(ixB^S,Fr(1)) = 0.d0
!
!     call hd_to_conserved(ixI^L,ixB^L,w,x)
!
!   case default
!     call mpistop('boundary not known')
!   end select
! end subroutine boundary_conditions

end module mod_usr
