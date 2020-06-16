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
    unit_temperature=unit_pressure/((2.d0+&
       3.d0*He_abundance)*unit_numberdensity*const_kB)
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

    Ti = T1 + 75.d0*(xprobmin1-half*dxlevel(1))/(&
       7.d10/unit_length)/unit_temperature
    To = T1 + 75.d0*(xprobmax1+half*dxlevel(1))/(&
       7.d10/unit_length)/unit_temperature

    temp(ixImin1:ixImax1,ixImin2:ixImax2) = T1 + 75.d0*x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1)/(7.d10/unit_length)/unit_temperature

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = rho1
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = 0.d0
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)) = v1
    where (x(ixImin1:ixImax1,ixImin2:ixImax2,1) .gt. 1.d0)
      w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)) = 0.d0
    endwhere
    !> Smoothen initial conditions
    ! w(ixI^S,mom(1)) = w(ixI^S,mom(1))*(1.d0-dexp(-1.d3*(x(ixI^S,1) - 1d0)**2.d0))
    pth(ixImin1:ixImax1,ixImin2:ixImax2) = temp(ixImin1:ixImax1,&
       ixImin2:ixImax2)*w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = pth(ixImin1:ixImax1,&
       ixImin2:ixImax2)/(rhd_gamma-1.d0) + half/rho1*w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(1))**2
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = const_rad_a*(temp(ixImin1:ixImax1,&
       ixImin2:ixImax2)*unit_temperature)**4.d0/unit_pressure

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, lambda, fld_R)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_diff_mg) = &
       (const_c/unit_velocity)*lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/(kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))
    w(ixImin1:ixImax1,ixImin2:ixImax2,i_test) = &
       (1.d0-dexp(-1.d3*(x(ixImin1:ixImax1,ixImin2:ixImax2,1) - 1d0)**2.d0))

  end subroutine initial_conditions


  subroutine boundary_conditions(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2,iB,w,x)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixBmin1,ixBmin2,ixBmax1,ixBmax2, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    integer :: i
    double precision :: temp(ixBmin1:ixBmax1,ixBmin2:ixBmax2)

    select case (iB)
    case(1)
      ! Ti = 2*(w(ixImin1+nghostcells, nghostcells+1,r_e)*unit_pressure/const_rad_a)**0.25/unit_temperature
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(:)) = 0.d0

      do i = ixBmax1,ixBmin1,-1
        w(i,:,rho_) = w(i+1,:,rho_)
        w(i,:,mom(1)) = w(i+1,:,mom(1))
      enddo
      temp(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = T1 + 75.d0*x(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,1)/(7.d10/unit_length)/unit_temperature
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,e_) = temp(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)*rho1/(rhd_gamma-1) + half*w(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,mom(1))**2/w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,r_e) = &
         const_rad_a*(temp(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)*unit_temperature)**4.d0/unit_pressure

    case(2)
      ! To  = 2*(w(ixImax1-nghostcells, nghostcells+1,r_e)*unit_pressure/const_rad_a)**0.25/unit_temperature
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(:)) = 0.d0

      do i = ixBmin1,ixBmax1
        w(i,:,rho_) = w(i-1,:,rho_)
        w(i,:,mom(1)) = w(i-1,:,mom(1))
      enddo
      temp(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = T1 + 75.d0*x(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,1)/(7.d10/unit_length)/unit_temperature
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,e_) = temp(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)*rho1/(rhd_gamma-1) + half*w(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,mom(1))**2/w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,r_e) = &
         const_rad_a*(temp(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)*unit_temperature)**4.d0/unit_pressure


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
      mg%bc(iB, mg_iphi)%bc_value = const_rad_a*(Ti*unit_temperature)**&
         4.d0/unit_pressure

      ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
      ! mg%bc(iB, mg_iphi)%bc_value = 0.d0
    case (2)

      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = const_rad_a*(To*unit_temperature)**&
         4.d0/unit_pressure

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


  subroutine specialvar_output(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_physics
    use mod_fld

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: Tgas(ixImin1:ixImax1,&
       ixImin2:ixImax2),Trad(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                   :: rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:ndim)
    double precision                   :: fld_R(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2), lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2)


    call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, rad_flux)
    call rhd_get_tgas(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, Tgas)
    call rhd_get_trad(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, Trad)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, lambda, fld_R)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1) = Tgas(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*unit_temperature
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2) = Trad(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*unit_temperature
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3) = rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+4) = lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+5) = fld_R(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+6) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(1))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'Tgas Trad F1 lambda R v1'
  end subroutine specialvarnames_output

end module mod_usr
