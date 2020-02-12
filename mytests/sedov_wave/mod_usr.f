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
    call set_coordinate_system("Cartesian_3D")

    usr_special_opacity => kramers_opacity

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = 1.d0 ! cm
    unit_temperature   = 1.d0 ! K
    unit_numberdensity = 1.d0 ! cm-3,cm-3,cm-3

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
    kap0 = 4.d0/3.d0*const_c*const_rad_a/Xi0/(&
       unit_velocity*unit_pressure/unit_temperature**4)*(g0/unit_density)
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
  subroutine initial_conditions(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas,phys_get_trad
    use mod_fld

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: local_rad_e(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)

    radius(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = dsqrt(sum(x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,:)**2,dim=ndim+1))
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       ind_r) = radius(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    where (radius(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) > delta_r/2)
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         rho_) = g0*radius(ixImin1:ixImax1,ixImin2:ixImax2,&
         ixImin3:ixImax3)**(-kp)
    else where (radius(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) <= delta_r/2)
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
         rho_) = g0*(delta_r/2)**(-kp)
    end where

    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mom(:)) = zero
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_) = 1d-3
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,r_e) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,e_)*w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,rho_)/(rhd_gamma-1) &
       *const_rad_a/unit_pressure*unit_temperature**4

    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,e_) + &
       E0/dsqrt(2*dpi*delta_r**2)*dexp(-radius(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)**2/(2*delta_r**2))

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, lambda,&
        fld_R)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_diff_mg) = (const_c/unit_velocity)*lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)/(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_))

  end subroutine initial_conditions

  subroutine kramers_opacity(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,kappa)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call phys_get_tgas(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Temp)

    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = kap0*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)*unit_density)**m*(Temp(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)*unit_temperature)**(-&
       ndim)*unit_opacity

  end subroutine kramers_opacity

end module mod_usr
