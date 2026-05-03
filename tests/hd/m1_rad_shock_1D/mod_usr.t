!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd
  use mod_m1

  implicit none

  double precision :: rho1, rho2, v1, v2, p1, p2, Er1, Er2

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    ! usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions
    ! usr_special_bc => boundary_conditions

    ! Choose independent normalization units if using dimensionless variables.
    unit_length      = 1.d0 ! cm
    unit_velocity    = const_c !
    unit_density     = 1.d0 ! g cm^-3

    call hd_activate()

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /shock_list/ rho1, rho2, v1, v2, p1, p2, Er1, Er2

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

    where (x(ixI^S,1) .lt. 0.d0)
      w(ixI^S,rho_) = rho1
      w(ixI^S,mom(1)) = v1
      w(ixI^S,p_) = p1
      w(ixI^S,Er_) = Er1
      w(ixI^S,Fr(1)) = 1.d-2*Er1 - v1*Er1
    elsewhere
      w(ixI^S,rho_) = rho2
      w(ixI^S,mom(1)) = v2
      w(ixI^S,p_) = p2
      w(ixI^S,Er_) = Er2
      w(ixI^S,Fr(1)) = 1.d-2*Er2 - v2*Er2
    end where

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions
!
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
