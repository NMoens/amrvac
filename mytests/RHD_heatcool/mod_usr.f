!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

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

    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Keep the radiative energy constant with internal bound
    usr_internal_bc => constant_r_e

    ! Write out energy levels and temperature
    usr_write_analysis => output_energy

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init

!==========================================================================================

subroutine initglobaldata_usr
  use mod_global_parameters

  !Units
  unit_numberdensity = 1.d0/((1.d0+4.d0*He_abundance)*const_mp)
  unit_pressure = 1.d0
  unit_time = 1.d0

  unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
  unit_velocity = dsqrt(unit_pressure/unit_density)
  unit_length = unit_time*unit_velocity
  unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)&
     *unit_numberdensity*const_kB)

  unit_radflux = unit_velocity*unit_pressure
  unit_opacity = one/(unit_density*unit_length)

  call usr_params_read(par_files)

  e_eq = (E_r0/(const_rad_a))**(1.d0/4.d0) *one/(rhd_gamma-one)*const_kB*rho0 &
     /(fld_mu*const_mp)

end subroutine initglobaldata_usr


subroutine usr_params_read(files)
  use mod_global_parameters, only: unitpar
  use mod_fld
  use mod_constants
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /test_list/ rho0, E_r0, t0

  do n = 1, size(files)
     open(unitpar, file=trim(files(n)), status="old")
     read(unitpar, test_list, end=111)
     111    close(unitpar)
  end do


end subroutine usr_params_read

!==========================================================================================

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    use mod_constants
    use mod_rhd_phys, only: rhd_get_pthermal

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)

    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2), fld_R(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)


    ! Set initial values for w
    w(ixImin1:ixImax1,ixImin2:ixImax2, rho_) = rho0
    w(ixImin1:ixImax1,ixImin2:ixImax2, mom(:)) = zero
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = E_r0
    w(ixImin1:ixImax1,ixImin2:ixImax2, e_) = t0*e_eq

    print*, 'unit_time', unit_time
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_length', unit_length
    print*, 'unit_density', unit_density
    print*, 'unit_numberdensity', unit_numberdensity
    print*, 'unit_velocity', unit_velocity
    print*, 'unit_pressure', unit_pressure
    print*, '================================================================'

  end subroutine initial_conditions

!==========================================================================================

  ! Extra routines can be placed here
  ! ...

!==========================================================================================

!> internal boundary, user defined
  !
  !> This subroutine can be used to artificially overwrite ALL conservative
  !> variables in a user-selected region of the mesh, and thereby act as
  !> an internal boundary region. It is called just before external (ghost cell)
  !> boundary regions will be set by the BC selection. Here, you could e.g.
  !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
  !> which can be used to identify the internal boundary region location.
  !> Its effect should always be local as it acts on the mesh.

  subroutine constant_r_e(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)


    !w(ixI^S,r_e) = 1d12

    ! if (it .eq. 0) open(1,file='Halley1_1.d2')
    ! write(1,*) global_time*unit_time, e_eq*unit_pressure, w(3,3,r_e)*unit_pressure, w(3,3,e_)*unit_pressure
    ! if (global_time .ge. time_max - dt) close(1)

    ! print*, global_time*unit_time, w(3,3,r_e)*unit_pressure, w(3,3,e_)*unit_pressure

  end subroutine constant_r_e

  subroutine output_energy()
    use mod_constants
    ! use mod_global_parameters

    double precision :: tmp_g, tmp_r, e_max, Er_max, rho_e_max
    double precision :: Tgas, Trad, tmp_T
    integer          :: iigrid, igrid, ierrmpi

    do iigrid = 1, igridstail
      igrid = igrids(iigrid)
      tmp_g = maxval(ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2, e_))
      tmp_r = maxval(ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2, r_e))
      tmp_T = maxval(ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
          e_)/ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2, rho_))
    end do

    call mpi_allreduce([tmp_g,tmp_r,tmp_T], [e_max,Er_max,rho_e_max], 3,&
        mpi_double_precision, mpi_max, icomm, ierrmpi)

    Tgas = tmp_T*(rhd_gamma - 1)*mp_cgs*fld_mu/kb_cgs*unit_temperature
    Trad = (tmp_r/const_rad_a)**0.25d0*unit_temperature

    ! print*, mype, global_time*unit_time, tmp_g*unit_pressure, tmp_r*unit_pressure, Tgas, Trad

    if (mype==0) then
      if (it .eq. 1) open(1,file = 'Instant1_1.d2',status = 'new')
      write(1,*) global_time*unit_time, tmp_g, tmp_r, e_eq, Tgas, Trad
      if (global_time .ge. time_max - dt) close(1)
    endif

  end subroutine output_energy

end module mod_usr

!==========================================================================================
