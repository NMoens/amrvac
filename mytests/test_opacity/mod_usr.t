!> This is a template for a new user problem
module mod_usr

! Include a physics module
use mod_rhd
use mod_fld

implicit none

contains

!> This routine should set user methods, and activate the physics module
subroutine usr_init()
  use mod_global_parameters
  use mod_usr_methods
  use mod_constants

  call set_coordinate_system("Cartesian_2D")

  call initglobaldata_usr

  ! Initialize units
  usr_set_parameters => initglobaldata_usr

  ! A routine for initial conditions is always required
  usr_init_one_grid => initial_conditions

  ! Active the physics module
  call rhd_activate()

end subroutine usr_init

!==========================================================================================

subroutine initglobaldata_usr
use mod_global_parameters

  print*, 'unit_density', unit_density
  print*, 'unit_pressure', unit_pressure
  print*, 'unit_temperature', unit_temperature
  print*, 'unit_opacity', unit_opacity
  print*, 'unit_radflux', unit_radflux

end subroutine initglobaldata_usr

!==========================================================================================

!> A routine for specifying initial conditions
subroutine initial_conditions(ixG^L, ix^L, w, x)
  use mod_global_parameters
  use mod_constants
  use mod_variables
  use mod_physics, only: phys_get_pthermal
  
  integer, intent(in)             :: ixG^L, ix^L
  double precision, intent(in)    :: x(ixG^S, ndim)
  double precision, intent(inout) :: w(ixG^S, nw)

  integer :: i,j

  do i= ixGmin1, ixGmax1
  do i= ixGmin2, ixGmax2
    w(i,j, rho_) = 1d-8 * (1.d0 + x(i,j,1))
    w(i,j, mom(1)) = zero
    w(i,j, mom(2)) = zero
    w(i,j, e_) =
    w(i,j, r_e) =
  enddo
  enddo

end subroutine initial_conditions


!==========================================================================================

end module mod_usr

!==========================================================================================
