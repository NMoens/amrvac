!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: rho1
  double precision :: rho2
  double precision :: v1
  double precision :: v2
  double precision :: T1
  double precision :: T2

  double precision :: p1,p2,eg1,eg2,Er1,Er2

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
     call set_coordinate_system("Cartesian_1D")
    
    

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Manually refine grid near shock
    usr_refine_grid => refine_shock

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! Boundary conditions
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    double precision :: k1,k2, my_gamma
    integer :: i

    call params_read(par_files)

    p1 = const_kB*T1*rho1/(const_mp*fld_mu)
    p2 = const_kB*T2*rho2/(const_mp*fld_mu)

    eg1 = p1/(rhd_gamma-1.d0) + half*v1*v1*rho1
    eg2 = p2/(rhd_gamma-1.d0) + half*v2*v2*rho2

    Er1 = const_rad_a*T1**4.d0
    Er2 = const_rad_a*T2**4.d0

    print*, 'M_1: ', v1/dsqrt(rhd_gamma*p1/rho1)
    print*, 'M_2: ', v2/dsqrt(rhd_gamma*p2/rho2)

    print*, 'RHD-quantity: ', 'Left', ' | ', 'Right'
    print*, 'density', rho1, ' | ', rho2
    print*, 'velocity', v1, ' | ', v2
    print*, 'momentum', rho1*v1, ' | ', rho2*v2
    print*, 'gas pressure', p1, ' | ', p2
    print*, 'gas energy', eg1, ' | ', eg2
    print*, 'radiation energy', Er1, ' | ', Er2

    unit_velocity = v1 !r_arr(nghostcells) ! cm
    unit_numberdensity = rho1/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = 1.d5

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+&
       3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    print*, 'unit_numberdensity', unit_numberdensity
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_length', unit_length

    rho1 = rho1/unit_density
    rho2 = rho2/unit_density

    v1 = v1/unit_velocity
    v2 = v2/unit_velocity

    T1 = T1/unit_temperature
    T2 = T2/unit_temperature

    p1 = p1/unit_pressure
    p2 = p2/unit_pressure

    eg1 = eg1/unit_pressure
    eg2 = eg2/unit_pressure

    Er1 = Er1/unit_pressure
    Er2 = Er2/unit_pressure


    print*, 'RHD-fluxes: ', 'Left', ' | ', 'Right'
    print*, 'density', rho1*v1, ' | ', rho2*v2
    print*, 'momentum', rho1*v1*v1 + p1 + Er1/3, ' | ', rho2*v2*v2 + p2 + Er2/3
    print*, 'gas energy', p1*v1 + eg1*v1, ' | ', p2*v2 + eg2*v2
    print*, 'radiation energy', Er1*v1, ' | ', Er2*v2

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /shock_list/ rho1, rho2, v1, v2, T1, T2

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, shock_list, end=113)
113    close(unitpar)
    end do

  end subroutine params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x)

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision :: kappa(ixOmin1:ixOmax1), fld_R(ixOmin1:ixOmax1),&
        lambda(ixOmin1:ixOmax1)

    w(ixImin1:ixImax1,rho_) = rho1
    w(ixImin1:ixImax1,mom(:)) = 0.d0
    w(ixImin1:ixImax1,mom(1)) = rho1*v1
    w(ixImin1:ixImax1,e_) = eg1
    w(ixImin1:ixImax1,r_e) = Er1

    where (x(ixImin1:ixImax1,1) .gt. 0.d0)
      w(ixImin1:ixImax1,rho_) = rho2
      w(ixImin1:ixImax1,mom(1)) = rho2*v2
      w(ixImin1:ixImax1,e_) = eg2
      w(ixImin1:ixImax1,r_e) = Er2
    end where

    ! call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    ! call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)
    !
    ! w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixImin1,ixImax1,ixBmin1,ixBmax1,iB,w,x)
    use mod_global_parameters
    use mod_fld


    integer, intent(in)             :: ixImin1,ixImax1, ixBmin1,ixBmax1, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    select case (iB)

    case(1)
      w(ixBmin1:ixBmax1,rho_) = rho1
      w(ixBmin1:ixBmax1,mom(:)) = 0.d0
      w(ixBmin1:ixBmax1,mom(1)) = rho1*v1
      w(ixBmin1:ixBmax1,e_) = eg1
      w(ixBmin1:ixBmax1,r_e) = Er1

    case(2)
      w(ixBmin1:ixBmax1,rho_) = rho2
      w(ixBmin1:ixBmax1,mom(:)) = 0.d0
      w(ixBmin1:ixBmax1,mom(1)) = rho2*v2
      w(ixBmin1:ixBmax1,e_) = eg2
      w(ixBmin1:ixBmax1,r_e) = Er2

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions


  subroutine mg_boundary_conditions(iB)

    use mod_global_parameters
    use mod_multigrid_coupling

    integer, intent(in)             :: iB

    select case (iB)
    case (1)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = Er1
    case (2)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = Er2

      case default
        print *, "Not a standard: ", trim(typeboundary(r_e, iB))
        error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions

  subroutine refine_shock(igrid,level,ixGmin1,ixGmax1,ixmin1,ixmax1,qt,w,x,&
     refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixGmin1,ixGmax1, ixmin1,ixmax1
    double precision, intent(in) :: qt, w(ixGmin1:ixGmax1,1:nw),&
        x(ixGmin1:ixGmax1,1:ndim)
    integer, intent(inout) :: refine, coarsen

    !> Refine close to base
    coarsen = -1
    refine = -1

    if (it .gt. slowsteps) then
      if (any(x(ixGmin1:ixGmax1,1) < 2.d-1 .and. x(ixGmin1:ixGmax1,&
         1) > -2.d-1)) refine=1
    endif

    if (it .gt. slowsteps) then
      if (any(x(ixGmin1:ixGmax1,1) < 1.d-1 .and. x(ixGmin1:ixGmax1,&
         1) > -1.d-1)) refine=1
    endif

  end subroutine refine_shock

end module mod_usr
