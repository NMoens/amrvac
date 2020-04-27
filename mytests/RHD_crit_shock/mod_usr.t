!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: ri
  double precision :: ro
  double precision :: rho1
  double precision :: T1
  !subcritical:
  double precision :: v1
  ! !supercritical:
  ! double precision :: v1 = 16.d5/2.d0

  double precision :: Ti
  double precision :: To

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! Boundary conditions
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions
    ! usr_internal_bc => fix_v


    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters

    call params_read(par_files)

    unit_velocity = v1 !r_arr(nghostcells) ! cm
    unit_numberdensity = rho1/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = ri

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    rho1 = rho1/unit_density
    v1 = v1/unit_velocity
    T1 = T1/unit_temperature

    print*, v1, rho1, T1

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /shock_list/ ri, ro, rho1, v1, T1

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, shock_list, end=113)
113    close(unitpar)
    end do

  end subroutine params_read


  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: temp(ixI^S), pth(ixI^S)
    double precision :: kappa(ixO^S), fld_R(ixO^S), lambda(ixO^S)

    Ti = T1 + 75.d0*(xprobmin1-half*dxlevel(1))/(7.d10/unit_length)/unit_temperature
    To = T1 + 75.d0*(xprobmax1+half*dxlevel(1))/(7.d10/unit_length)/unit_temperature

    temp(ixI^S) = T1 + 75.d0*x(ixI^S,1)/(7.d10/unit_length)/unit_temperature

    w(ixI^S,rho_) = rho1
    w(ixI^S,mom(1)) = v1
    where (x(ixI^S,1) .gt. 1.d0)
      w(ixI^S,mom(1)) = 0.d0
    endwhere
    !> Smoothen initial conditions
    ! w(ixI^S,mom(1)) = w(ixI^S,mom(1))*(1.d0-dexp(-1.d3*(x(ixI^S,1) - 1d0)**2.d0))
    w(ixI^S,mom(2)) = 0.d0
    pth(ixI^S) = temp(ixI^S)*w(ixI^S,rho_)
    w(ixI^S,e_) = pth(ixI^S)/(rhd_gamma-1.d0) + half/rho1*w(ixI^S,mom(1))**2
    w(ixI^S,r_e) = const_rad_a*(temp(ixI^S)*unit_temperature)**4.d0/unit_pressure

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))
    w(ixI^S,i_test) = (1.d0-dexp(-1.d3*(x(ixI^S,1) - 1d0)**2.d0))

  end subroutine initial_conditions


  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: i
    double precision :: temp(ixB^S)

    select case (iB)
    case(1)
      ! Ti = 2*(w(ixImin1+nghostcells, nghostcells+1,r_e)*unit_pressure/const_rad_a)**0.25/unit_temperature

      do i = ixBmax1,ixBmin1,-1
        w(i,:,rho_) = w(i+1,:,rho_)
        w(i,:,mom(1)) = w(i+1,:,mom(1))
      enddo
      w(ixB^S,mom(2)) = 0.d0
      temp(ixB^S) = T1 + 75.d0*x(ixB^S,1)/(7.d10/unit_length)/unit_temperature
      w(ixB^S,e_) = temp(ixB^S)*rho1/(rhd_gamma-1) + half*w(ixB^S,mom(1))**2/w(ixB^S,rho_)
      w(ixB^S,r_e) = const_rad_a*(temp(ixB^S)*unit_temperature)**4.d0/unit_pressure

    case(2)
      ! To  = 2*(w(ixImax1-nghostcells, nghostcells+1,r_e)*unit_pressure/const_rad_a)**0.25/unit_temperature

      do i = ixBmin1,ixBmax1
        w(i,:,rho_) = w(i-1,:,rho_)
        w(i,:,mom(1)) = w(i-1,:,mom(1))
      enddo
      w(ixB^S,mom(2)) = 0.d0
      temp(ixB^S) = T1 + 75.d0*x(ixB^S,1)/(7.d10/unit_length)/unit_temperature
      w(ixB^S,e_) = temp(ixB^S)*rho1/(rhd_gamma-1) + half*w(ixB^S,mom(1))**2/w(ixB^S,rho_)
      w(ixB^S,r_e) = const_rad_a*(temp(ixB^S)*unit_temperature)**4.d0/unit_pressure


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
      mg%bc(iB, mg_iphi)%bc_value = const_rad_a*(Ti*unit_temperature)**4.d0/unit_pressure

      ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
      ! mg%bc(iB, mg_iphi)%bc_value = 0.d0
    case (2)

      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = const_rad_a*(To*unit_temperature)**4.d0/unit_pressure

      ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
      ! mg%bc(iB, mg_iphi)%bc_value = 0.d0

    case default
      print *, "Not a standard: ", trim(typeboundary(r_e, iB))
      error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions

  ! subroutine fix_v(level,qt,ixI^L,ixO^L,w,x)
  !   use mod_global_parameters
  !
  !   integer, intent(in)             :: ixI^L,ixO^L,level
  !   double precision, intent(in)    :: qt
  !   double precision, intent(inout) :: w(ixI^S,1:nw)
  !   double precision, intent(in)    :: x(ixI^S,1:ndim)
  !
  !   if (global_time .lt. 0.5d0) then
  !     w(ixI^S,mom(1)) = -rho1*v1
  !     where (x(ixI^S,1) .lt. 1d0)
  !       w(ixI^S,mom(1)) = rho1*v1
  !     endwhere
  !     !> Smoothen initial conditions
  !     w(ixI^S,mom(1)) = w(ixI^S,mom(1))*(1.d0-dexp(-1.d3*(x(ixI^S,1) - 1d0)**2.d0))
  !     w(ixI^S,mom(2)) = 0.d0
  !   endif
  !
  ! end subroutine fix_v


  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_physics
    use mod_fld

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: Tgas(ixI^S),Trad(ixI^S)
    double precision                   :: rad_flux(ixO^S,1:ndim)
    double precision                   :: fld_R(ixO^S), lambda(ixO^S)


    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
    call rhd_get_tgas(w, x, ixI^L, ixO^L, Tgas)
    call rhd_get_trad(w, x, ixI^L, ixO^L, Trad)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)

    w(ixO^S,nw+1) = Tgas(ixO^S)*unit_temperature
    w(ixO^S,nw+2) = Trad(ixO^S)*unit_temperature
    w(ixO^S,nw+3) = rad_flux(ixO^S,1)
    w(ixO^S,nw+4) = lambda(ixO^S)
    w(ixO^S,nw+5) = fld_R(ixO^S)
    w(ixO^S,nw+6) = w(ixO^S,mom(1))/w(ixO^S,rho_)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'Tgas Trad F1 lambda R v1'
  end subroutine specialvarnames_output

end module mod_usr
