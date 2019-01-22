module mod_usr
  use mod_rho
  use mod_multigrid_coupling

  implicit none

  integer             :: i_sol
  integer             :: i_eps
  integer             :: i_err
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: diffusion_coeff = 0.2_dp
  real(dp), parameter :: solution_modes(2) = [1, 1]

contains

  subroutine usr_init()
    use mod_usr_methods
    use mod_variables
    use mod_physics

    usr_init_one_grid => initial_conditions

    use_multigrid = .true.
    phys_global_source => diffuse_density
    usr_process_grid => set_error
    mg_after_new_tree => set_epsilon

    mg%operator_type = mg_vhelmholtz
    mg%bc(:, mg_iphi)%bc_type = bc_neumann
    mg%bc(:, mg_iphi)%bc_value = 0.0d0

    call set_coordinate_system("Cartesian_2D")
    call rho_activate()

    i_sol = var_set_extravar("sol", "sol")
    i_err = var_set_extravar("err", "err")
    i_eps = var_set_extravar("eps", "eps")

  end subroutine usr_init

  subroutine initial_conditions(ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
     ixmax1,ixmax2, w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
       ixmin2,ixmax1,ixmax2
    double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
        1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
        1:nw)

    w(ixmin1:ixmax1,ixmin2:ixmax2, rho_) = solution(x(ixmin1:ixmax1,&
       ixmin2:ixmax2, 1), x(ixmin1:ixmax1,ixmin2:ixmax2, 2), 0.0d0)
    w(ixmin1:ixmax1,ixmin2:ixmax2, i_eps) = diffusion_coeff + 1.0 * &
       x(ixmin1:ixmax1,ixmin2:ixmax2, 1)

  end subroutine initial_conditions

  elemental function solution(x, y, t) result(sol)
    real(dp), intent(in) :: x, y, t
    real(dp)             :: sol, tmp(ndim)

    tmp = 2 * pi * solution_modes * [x, y]
    sol = 1 + product(cos(tmp)) * exp(-sum((2 * pi * solution_modes)**2) * &
       diffusion_coeff * t)
  end function solution

  subroutine diffuse_density(qdt, qt, active)
    use m_diffusion
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qt
    logical, intent(inout)       :: active
    double precision             :: max_res

    call mg_copy_to_tree(rho_, mg_iphi, .false., .false.)
    call diffusion_solve_vcoeff(mg, qdt, 1, 1d-4)
    call mg_copy_from_tree(mg_iphi, rho_)

    active = .true.
  end subroutine diffuse_density

  subroutine set_error(igrid,level,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,qt,w,x)
    integer, intent(in)             :: igrid,level,ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_sol) = solution(x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, 1), x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 2), qt)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_err) = abs(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_sol))
  end subroutine set_error

  subroutine set_epsilon()
    call mg_copy_to_tree(i_eps, mg_iveps, .true., .true.)
  end subroutine set_epsilon

end module mod_usr

