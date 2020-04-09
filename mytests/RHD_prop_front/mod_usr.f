!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: Er0 = 1.d-1
  double precision :: rho0 = 0.025

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")


    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    usr_internal_bc => set_Erad

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: temp(ixImin1:ixImax1,ixImin2:ixImax2),&
        pth(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2), lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    temp(ixImin1:ixImax1,ixImin2:ixImax2) = (Er0/const_rad_a)**0.25

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = rho0
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = 0.d0
    pth(ixImin1:ixImax1,ixImin2:ixImax2) = temp(ixImin1:ixImax1,&
       ixImin2:ixImax2)*const_kB/(const_mp*fld_mu)*w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_)
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = pth(ixImin1:ixImax1,&
       ixImin2:ixImax2)/(rhd_gamma-1.d0)
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = Er0

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, lambda, fld_R)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_diff_mg) = &
       (const_c/unit_velocity)*lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/(kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

  end subroutine initial_conditions


  subroutine set_Erad(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    if (it>0) then
      where (x(ixImin1:ixImax1,ixImin2:ixImax2,1) .lt. 0.1d0)
        w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = 1.d0
      end where
    endif

  end subroutine set_Erad

end module mod_usr
