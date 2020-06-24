!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: Er0
  double precision :: Er1
  double precision :: rho0
  double precision :: l1
  double precision :: l2

  double precision :: p0, T0, p1, T1

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
     call set_coordinate_system("Cartesian_1D")
    
    

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Boundary conditions
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    integer :: i

    call params_read(par_files)

    T0 = (Er0/const_rad_a)**0.25
    p0 = const_kB*T0*rho0/(const_mp*fld_mu)

    T1 = (Er1/const_rad_a)**0.25
    p1 = const_kB*T1*rho0/(const_mp*fld_mu)

    unit_velocity = dsqrt(p1/rho0)
    unit_numberdensity = rho0/((1.d0+4.d0*He_abundance)*const_mp)
    unit_length = 1.d0

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


    print*, unit_time, unit_velocity
    print*, '2.5d-11s in dimless is ', 2.5d-11/unit_time
    print*, 'speed of light dimless is', const_c/unit_velocity
    print*, 'dimless time for crossing is', xprobmax1/const_c*unit_velocity

    rho0 = rho0/unit_density
    p0 = p0/unit_pressure
    p1 = p1/unit_pressure
    T0 = T0/unit_temperature
    T1= T1/unit_temperature
    Er0 = Er0/unit_pressure
    Er1 = Er1/unit_pressure
    l1 = l1/unit_length
    l2 = l2/unit_length

  end subroutine initglobaldata_usr


  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /front_list/ rho0, Er0, Er1, l1, l2

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, front_list, end=113)
113    close(unitpar)
    end do

  end subroutine params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision :: kappa(ixOmin1:ixOmax1), fld_R(ixOmin1:ixOmax1),&
        lambda(ixOmin1:ixOmax1)
    double precision :: rad_flux(ixOmin1:ixOmax1,1:ndim),&
        step(ixImin1:ixImax1)

    w(ixImin1:ixImax1,rho_) = rho0
    w(ixImin1:ixImax1,mom(:)) = 0.d0
    w(ixImin1:ixImax1,e_) = p0/(rhd_gamma-1.d0)

    w(ixImin1:ixImax1,r_e) = Er0
    step(ixImin1:ixImax1) = (  1.d0-erf(  (x(ixImin1:ixImax1,&
       1)-l1)/l2    )  )/2.d0
    w(ixImin1:ixImax1,r_e) = w(ixImin1:ixImax1,&
       r_e) + step(ixImin1:ixImax1)*Er1

    call fld_get_opacity(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, lambda,&
        fld_R)

    w(ixOmin1:ixOmax1,i_test) = fld_R(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,i_diff_mg) = (const_c/unit_velocity)*lambda(&
       ixOmin1:ixOmax1)/(kappa(ixOmin1:ixOmax1)*w(ixOmin1:ixOmax1,rho_))

  end subroutine initial_conditions


  subroutine boundary_conditions(qt,ixImin1,ixImax1,ixBmin1,ixBmax1,iB,w,x)
    use mod_global_parameters
    use mod_fld


    integer, intent(in)             :: ixImin1,ixImax1, ixBmin1,ixBmax1, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    integer :: i

    select case (iB)
    case(1)
        w(ixBmin1:ixBmax1,rho_) = rho0
        w(ixBmin1:ixBmax1,mom(:)) = 0.d0
        w(ixBmin1:ixBmax1,r_e) = Er1
        do i = ixBmax1,ixBmin1, -1
           w(i,e_) = w(i+1,e_)
           w(i,mom(1)) = w(i+1,mom(1))
          
          
        enddo

    case(2)
      do i = ixBmin1,ixBmax1
        
        w(i,rho_) = w(i-1,rho_)
        w(i,mom(1)) = w(i-1,mom(1))
        w(i,e_) = w(i-1,e_)
        w(i,r_e) = w(i-1,r_e)
       
        
      enddo

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
        ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous

        mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
        mg%bc(iB, mg_iphi)%bc_value = 0.d0
    case default
      print *, "Not a standard: ", trim(typeboundary(r_e, iB))
      error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions


  subroutine specialvar_output(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
    double precision                   :: w(ixImin1:ixImax1,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: step(ixImin1:ixImax1), rad_flux(ixOmin1:ixOmax1,&
       1:ndim)
    double precision :: lambda(ixOmin1:ixOmax1), fld_R(ixOmin1:ixOmax1),&
        kappa(ixOmin1:ixOmax1)

    double precision :: rad_e(ixImin1:ixImax1), normgrad2(ixOmin1:ixOmax1),&
        grad_r_e(ixImin1:ixImax1)
    double precision :: grE1(ixImin1:ixImax1), grE2(ixImin1:ixImax1)
    integer :: idir

    step(ixImin1:ixImax1) = (  1.d0-erf((x(ixImin1:ixImax1,&
       1)-l1-global_time*const_c/unit_velocity)/l2    )  )/2.d0
    w(ixImin1:ixImax1,nw+1) = step(ixImin1:ixImax1)

    call fld_get_radflux(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, rad_flux)
    w(ixOmin1:ixOmax1,nw+2) = rad_flux(ixOmin1:ixOmax1,1)

    w(ixOmin1:ixOmax1,nw+3) = const_c/unit_velocity*w(ixOmin1:ixOmax1,r_e)

    call fld_get_fluxlimiter(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, lambda,&
        fld_R)
    w(ixOmin1:ixOmax1,nw+4) = lambda(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+5) = fld_R(ixOmin1:ixOmax1)

    call fld_get_opacity(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, kappa)
    w(ixOmin1:ixOmax1,nw+6) = kappa(ixOmin1:ixOmax1)

    normgrad2(ixOmin1:ixOmax1) = 0.d0 !smalldouble

    rad_e(ixImin1:ixImax1) = w(ixImin1:ixImax1, r_e)
    do idir = 1,ndim
      call gradient(rad_e,ixImin1,ixImax1,ixOmin1,ixOmax1,idir,grad_r_e)
      normgrad2(ixOmin1:ixOmax1) = normgrad2(ixOmin1:ixOmax1) + &
         grad_r_e(ixOmin1:ixOmax1)**2
    end do

    w(ixOmin1:ixOmax1,nw+7) = normgrad2(ixOmin1:ixOmax1)

    call gradient(rad_e,ixImin1,ixImax1,ixOmin1,ixOmax1,1,grE1)

    w(ixOmin1:ixOmax1,nw+8) = grE1(ixOmin1:ixOmax1)

    ! if (x(1,1,1) .lt. xprobmin1) then
    !   print*, 'Er', w(1:5,5,r_e)
    !   print*, 'step', step(1:5,1)
    !   print*, 'lambda', lambda(1:5,5)
    !   print*, 'R', fld_R(1:5,5)
    !   print*, 'normgr', normgrad2(1:5,5)
    !   print*, 'gr1', grE1(1:5,5)
    !   print*, 'gr2', grE2(1:5,5)
    !   stop
    ! endif

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'step F1 cE lambda R kappa ngrd grE1'
  end subroutine specialvarnames_output

end module mod_usr
