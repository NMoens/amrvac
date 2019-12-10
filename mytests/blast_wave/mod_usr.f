!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system
    call set_coordinate_system("Cartesian_2D")

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = 1.d0 ! cm
    unit_temperature   = 1.d0 ! K
    unit_numberdensity = 1.d0!/((1.d0+4.d0*He_abundance)*mp_cgs) ! cm-3,cm-3

    ! Active the physics module
    call rhd_activate()

    fld_kappa0 = fld_kappa0/unit_opacity

  end subroutine usr_init

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas
    use mod_fld

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: rbs,xc1,xc2
    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: local_rad_e

    if (mype .eq. 0) then
      print*, 'unit_length',unit_length
      print*, 'unit_velocity',unit_velocity
      print*, 'unit_density',unit_density
      print*, 'unit_pressure',unit_pressure
      print*, 'unit_flux',unit_radflux
    endif

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = one
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = zero
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = one

    call phys_get_tgas(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,Temp)

    !> Have to do this in 2 steps to define boundary e_rad
    local_rad_e = const_rad_a/unit_pressure*unit_temperature**4.d0 &
       *Temp(nghostcells+2,nghostcells+2)**4.d0

    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = local_rad_e

    xc1=(xprobmin1+xprobmax1)*0.5d0
    xc2=(xprobmin2+xprobmax2)*0.5d0
    rbs=0.2d0
    where((x(ixImin1:ixImax1,ixImin2:ixImax2,1)-xc1)**2+(x(ixImin1:ixImax1,&
       ixImin2:ixImax2,2)-xc2)**2<rbs**2)
      w(ixImin1:ixImax1,ixImin2:ixImax2,e_)=5.d0
    endwhere

    ! w(ixO^S,e_) = one + w(ixO^S,e_)*100.d0*dexp(-(x(ixO^S,1)**2 + x(ixO^S,2)**2)/rbs**2)

    call get_rad_extravars(w, w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2)

  end subroutine initial_conditions

end module mod_usr
