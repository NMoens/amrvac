!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: g0
  double precision :: E0
  double precision :: kap0
  double precision :: a
  double precision :: b
  double precision :: m
  double precision :: kp
  double precision :: alpha

  double precision :: Xi0
  double precision :: Gam
  double precision :: Omega

  double precision :: delta_r


  integer :: ind_r
  integer :: ind_Tg, ind_Tr

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system
    ! Choose coordinate system as 2D Cartesian with three components for vectors
    {^IFONED call set_coordinate_system("Cartesian_1D")}
    {^IFTWOD call set_coordinate_system("Cartesian_2D")}
    {^IFTHREED call set_coordinate_system("Cartesian_3D")}

    usr_special_opacity => kramers_opacity

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Drive the wave using an internal boundary
    usr_internal_bc => Initialize_Wave

    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = 1.d0 ! cm
    unit_temperature   = 1.d0 ! K
    unit_numberdensity = 1.d0 ! cm^-3

    ! Active the physics module
    call rhd_activate()

    call usr_params_read(par_files)

    ind_r = var_set_extravar("r", "r")

  end subroutine usr_init


  subroutine usr_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_fld
    use mod_constants
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /sedov_list/ g0, E0, Xi0, a, b, m, delta_r

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, sedov_list, end=111)
       111    close(unitpar)
    end do

    m = -a
    b = ndim + 3
    kp = ((2*b-1)*ndim + 2)/(2*b - 2*a + 1)
    alpha = (2*b - 2*a + 1)/(2*b - (ndim+2)*a + ndim)

    Gam = 1 !const_kb/(const_mp*0.6d0)/(unit_pressure/(unit_time*unit_temperature))
    kap0 = 4.d0/3.d0*const_c*const_rad_a/Xi0/(unit_velocity*unit_pressure/unit_temperature**4)*(g0/unit_density)
    Omega = 2*Xi0/(Gam**(b+1) * g0**(1-a))*(E0/g0)**(b-1.d0/2.d0)

    if (mype ==0) then
      print*, 'g0', g0
      print*, 'E0', E0
      print*, 'Xi0', Xi0
      print*, 'kp', kp
      print*, 'alpha', alpha
      print*, 'Gamma', Gam
      print*, 'kappa0', kap0
      print*, 'Omega', Omega
    endif

  end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas,phys_get_trad
    use mod_fld

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Temp(ixI^S)
    double precision :: local_rad_e(ixI^S)

    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)
    double precision :: radius(ixI^S)

    radius(ixI^S) = dsqrt(sum(x(ixI^S,:)**2,dim=ndim+1))
    w(ixI^S,ind_r) = radius(ixI^S)

    where (radius(ixI^S) > delta_r/2)
      w(ixI^S,rho_) = g0*radius(ixI^S)**(-kp)
    else where (radius(ixI^S) <= delta_r/2)
      w(ixI^S,rho_) = g0*(delta_r/2)**(-kp)
    end where

    w(ixI^S,mom(:)) = zero
    w(ixI^S,e_) = 1d-3
    w(ixI^S,r_e) = w(ixI^S,e_)*w(ixI^S,rho_)/(rhd_gamma-1) &
    *const_rad_a/unit_pressure*unit_temperature**4

    w(ixI^S,e_) = w(ixI^S,e_) + E0/dsqrt(2*dpi*delta_r**2)*dexp(-radius(ixI^S)**2/(2*delta_r**2))

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)
    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))

  end subroutine initial_conditions

  subroutine Initialize_Wave(level,qt,ixI^L,ixO^L,w,x)
    use mod_global_parameters
    use mod_fld
    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision :: radius(ixI^S)

    where (radius(ixI^S) <= delta_r/2)
      w(ixI^S,e_) = E0 !w(ixI^S,e_) + E0/dsqrt(2*dpi*delta_r**2)*dexp(-radius(ixI^S)**2/(2*delta_r**2))
    end where

  end subroutine Initialize_Wave

  subroutine kramers_opacity(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    double precision :: Temp(ixI^S)

    call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)

    kappa(ixO^S) = kap0*(w(ixO^S,rho_)*unit_density)**m*(Temp(ixO^S)*unit_temperature)**(-ndim)*unit_opacity

  end subroutine kramers_opacity

end module mod_usr
