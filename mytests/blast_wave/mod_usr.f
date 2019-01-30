!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = 1.d0 ! cm
    unit_temperature   = 1.d0 ! K
    unit_numberdensity = 1.d0 ! cm-3,cm-3

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: rbs,xc1,xc2
    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2)

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = one
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = zero
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = one

    call phys_get_tgas(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,Temp)

    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = &
       const_rad_a/unit_pressure*unit_temperature**4.d0 *Temp(ixImin1:ixImax1,&
       ixImin2:ixImax2)**4.d0

    xc1=(xprobmin1+xprobmax1)*0.5d0
    xc2=(xprobmin2+xprobmax2)*0.5d0
    rbs=0.2d0
    where((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)-xc1)**2+(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)-xc2)**2<rbs**2)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=100.d0
    endwhere

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2)

    if (fld_diff_scheme .eq. 'mg') then
      call fld_get_diffcoef_central(w, x, ixImin1,ixImin2,ixImax1,ixImax2,&
          ixOmin1,ixOmin2,ixOmax1,ixOmax2)
      call set_mg_bounds()
    endif

  end subroutine initial_conditions

end module mod_usr
