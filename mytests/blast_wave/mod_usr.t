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
    unit_numberdensity = 1.d0 ! cm^-3

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas
    use mod_fld

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rbs,xc1,xc2
    double precision :: Temp(ixI^S)

    w(ixI^S,rho_) = one
    w(ixI^S,mom(:)) = zero
    w(ixI^S,e_) = one

    call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)

    w(ixI^S,r_e) = const_rad_a/unit_pressure*unit_temperature**4.d0 &
    *Temp(ixI^S)**4.d0

    xc1=(xprobmin1+xprobmax1)*0.5d0
    xc2=(xprobmin2+xprobmax2)*0.5d0
    rbs=0.2d0
    where((x(ixO^S,1)-xc1)**2+(x(ixO^S,2)-xc2)**2<rbs**2)
      w(ixO^S,e_)=10.d0
    endwhere

    ! w(ixO^S,r_e) = w(ixO^S,r_e)*100.d0*dexp(-(x(ixO^S,1)**2 + x(ixO^S,2)**2)/rbs**2)

    call fld_get_opacity(w, x, ixI^L, ixO^L)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L)
    call fld_get_radflux(w, x, ixI^L, ixO^L)

    if (fld_diff_scheme .eq. 'mg') then
      call fld_get_diffcoef_central(w, x, ixI^L, ixO^L)
      call set_mg_bounds()
    endif

  end subroutine initial_conditions

end module mod_usr
