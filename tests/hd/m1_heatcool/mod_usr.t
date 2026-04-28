!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_m1

  implicit none

  double precision :: e_eq
  double precision :: rho0
  double precision :: t0
  double precision :: e0
  double precision :: E_r0

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods
    use mod_constants

    call set_coordinate_system("Cartesian")


    call usr_params_read(par_files)


    !Units THIS HAS BEEN A MESS AS AMRVAC DOESN'T ACCEPT "1" as a unit
    unit_numberdensity = rho0/((1d0 + 4d0 * He_abundance)*const_mp)
    unit_pressure = e_eq
    unit_length=sqrt(e_eq/rho0) !> hack to set unit_time=1

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Active the physics module
    call hd_activate()

  end subroutine usr_init

!==========================================================================================

subroutine usr_params_read(files)
  use mod_global_parameters, only: unitpar
  use mod_constants
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /test_list/ rho0, E_r0, t0

  do n = 1, size(files)
     open(unitpar, file=trim(files(n)), status="old")
     read(unitpar, test_list, end=111)
     111    close(unitpar)
  end do

  e_eq = (E_r0/(const_rad_a))**(1.d0/4.d0) &
        *one/(hd_gamma-one)*const_kB*rho0 &
        /(hd_mu*const_mp)

  E_r0 = E_r0/e_eq


end subroutine usr_params_read

!==========================================================================================

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_constants
    use mod_rhd_phys, only: rhd_get_pthermal

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S, ndim)
    double precision, intent(inout) :: w(ixI^S, nw)

    if (mype .eq. 0) then
      print*, 'unit_time', unit_time
      print*, 'unit_temperature', unit_temperature
      print*, 'unit_length', unit_length
      print*, 'unit_density', unit_density
      print*, 'unit_numberdensity', unit_numberdensity
      print*, 'unit_velocity', unit_velocity
      print*, 'unit_pressure', unit_pressure
      print*, '================================================================'
    end if

    ! Set initial values for w
    w(ixI^S, rho_) = 1.d0
    w(ixI^S, mom(:)) = zero
    w(ixI^S, e_) = t0
    w(ixI^S, Er_) = E_r0
    w(ixI^S, Fr(:)) = zero

  end subroutine initial_conditions

end module mod_usr

!==========================================================================================
