!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd

  implicit none

  double precision :: rho1, v1, T1
  double precision :: m1, p1, e1, Er1

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output

    ! Choose independent normalization units if using dimensionless variables.
    unit_velocity = v1
    unit_numberdensity = rho1/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = 7.d10

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

    namelist /shock_list/ rho1, v1, T1

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, shock_list, end=111)
       111 close(unitpar)
    enddo

  end subroutine usr_params_read

  subroutine initglobaldata_usr
    use mod_constants

    m1 = rho1*v1
    p1 = const_kb*T1*rho1/(const_mp*hd_mu)
    Er1 = const_rad_a*T1**4

    print*, 'rho', rho1
    print*, 'm', m1
    print*, 'e', p1/(hd_gamma-1)
    print*, 'Er', Er1

    print*, 'unit_density', unit_density
    print*, 'unit_velocity', unit_velocity
    print*, 'unit_pressure', unit_pressure
    print*, 'unit_time', unit_time

    rho1 = rho1/unit_density
    v1 = v1/unit_velocity
    p1 = p1/unit_pressure
    Er1 = Er1/unit_pressure

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixI^S,rho_) = rho1
    w(ixI^S,mom(1)) = -v1
    w(ixI^S,p_) = p1
    w(ixI^S,Er_) = Er1
    w(ixI^S,Fr(1)) = 0.d0

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: pth(ixI^S)

    !Gas temperature
    call hd_get_pthermal(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1) = (pth(ixO^S)*unit_pressure*const_mp*hd_mu) &
                   /(w(ixO^S,rho_)*unit_density*const_kB)

    !Radiation temperature
    w(ixO^S,nw+2) = (w(ixO^S,Er_)*unit_pressure/const_rad_a)**0.25d0

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames

    varnames='Tgas Trad'

  end subroutine specialvarnames_output

end module mod_usr
