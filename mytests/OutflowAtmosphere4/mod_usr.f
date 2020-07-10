!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision, parameter :: M_sun = 1.9891000d33
  double precision, parameter :: R_sun = 6.9599000d10
  double precision, parameter :: L_sun = 3.8268000d33
  double precision, parameter :: year = 365.25*24*60*60

  double precision, allocatable :: r_arr(:)
  double precision, allocatable :: rho_arr(:)
  double precision, allocatable :: v_arr(:)
  double precision, allocatable :: e_arr(:)
  double precision, allocatable :: Er_arr(:)
  double precision, allocatable :: T_arr(:)
  double precision, allocatable :: p_arr(:)

  double precision :: M_dot_ratio
  double precision :: Gamma_0, Gamma_b
  double precision :: kappa_0, kappa_b
  double precision :: L_0,L_vE
  double precision :: M_star
  double precision :: R_star
  double precision :: M_dot
  double precision :: rho_base
  double precision :: T_base

  double precision :: StefBoltz


  double precision :: dinflo,gradE,gradE_out, valE_out
  double precision :: error_b

  double precision :: inflow_density, inflow_gamma, pert_ampl
  character(len=std_len) :: init_struc

  integer :: int_r, int_v, int_e, int_re, int_dt


contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
     call set_coordinate_system("Cartesian_1D")
    

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! Boundary conditions
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions

    ! PseudoPlanar correction
    usr_source => PseudoPlanar
    usr_internal_bc => reset_egas

    ! Graviatational field
    usr_gravity => set_gravitation_field

    ! Special Opacity
    usr_special_opacity => Opacity_stepfunction

    ! Analysis
    usr_modify_output => time_average_values

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    usr_refine_grid => refine_base

    ! Active the physics module
    call rhd_activate()

    int_r = var_set_extravar("int_r", "int_r")
    int_v = var_set_extravar("int_v", "int_v")
    int_e = var_set_extravar("int_e", "int_e")
    int_re = var_set_extravar("int_re", "int_re")
    int_dt = var_set_extravar("int_dt", "int_dt")

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_fld

    integer :: i

    !> read usr par
    call params_read(par_files)
    !> Set stellar mass and radius
    call ReadInParams(M_star,R_star,Gamma_0,M_dot_ratio,M_dot,L_0,rho_base,&
       error_b)

    !> Gamma at the base is below one!
    Gamma_b = inflow_gamma
    kappa_0 = Gamma_0*4*dpi*const_G*M_star*const_c/L_0
    kappa_b = Gamma_b*4*dpi*const_G*M_star*const_c/L_0

    allocate(r_arr(domain_nx1+2*nghostcells))
    allocate(rho_arr(domain_nx1+2*nghostcells))
    allocate(v_arr(domain_nx1+2*nghostcells))
    allocate(e_arr(domain_nx1+2*nghostcells))
    allocate(Er_arr(domain_nx1+2*nghostcells))
    allocate(T_arr(domain_nx1+2*nghostcells))
    allocate(p_arr(domain_nx1+2*nghostcells))

    ! call ReadInTable(r_arr,rho_arr,v_arr,e_arr,Er_arr,T_arr, p_arr)
    call ReadInTable(r_arr,rho_arr,v_arr,e_arr,Er_arr,T_arr, p_arr)

    rho_base = rho_arr(nghostcells+1)
    T_base = T_arr(nghostcells+1)

    ! if (mype .eq. 0) then
    !   print*, 'density at base', rho_base
    !   print*, 'Temperature base', T_base
    !   print*, 'cgs opacity', kappa_0
    ! endif

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = 69599000000.0d0 !R_star
    unit_numberdensity = 2218242342924238.2d0 !rho_base/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_temperature = 296199.82122247218d0 !T_base

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=(2.d0+3.d0*He_abundance)&
       *unit_numberdensity*const_kb*unit_temperature
    unit_velocity=dsqrt(unit_pressure/unit_density)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    if (mype .eq. 0) then
      print*, 'M_star ', 'R_star ','M_dot_ratio ', 'M_dot ', 'L_0'
      print*, 'Gamma_0 ', 'kappa_0 ', 'kappa_b'
      print*, M_star, R_star,M_dot_ratio, M_dot, L_0
      print*, Gamma_0, kappa_0, kappa_b

      print*, 'Flux at boundary: ', L_0/(4*dpi*R_star**2)

      print*, 'unit_length', unit_length
      print*, 'unit_density', unit_density
      print*, 'unit_numberdensity', unit_numberdensity
      print*, 'unit_pressure', unit_pressure
      print*, 'unit_temperature', unit_temperature
      print*, 'unit_radflux', unit_radflux
      print*, 'unit_opacity', unit_opacity
      print*, 'unit_time', unit_time
      print*, 'unit_velocity', unit_velocity
    endif

    !> Make all parameters dimensionless
    M_star = M_star/(unit_density*unit_length**3.d0)
    R_star = R_star/unit_length
    M_dot = M_dot/(unit_density*unit_length**3.d0)*unit_time
    L_0 = L_0/(unit_pressure*unit_length**3.d0)*unit_time
    kappa_0 = kappa_0/unit_opacity
    kappa_b = kappa_b/unit_opacity
    rho_base = rho_base/unit_density
    T_base = T_base/unit_temperature

    !> Make initial profiles dimensionless
    r_arr = r_arr/unit_length
    rho_arr = rho_arr/unit_density
    v_arr = v_arr/unit_velocity
    e_arr = e_arr/unit_pressure
    Er_arr = Er_arr/unit_pressure
    T_arr = T_arr/unit_temperature
    p_arr = p_arr/unit_pressure

    ! if (mype .eq. 0) then
    !   print*, 'M_star ', 'R_star ','M_dot_ratio ', 'M_dot ', 'L_0'
    !   print*, 'Gamma_0 ', 'kappa_0 ', 'kappa_b'
    !   print*, M_star, R_star,M_dot_ratio, M_dot, L_0
    !   print*, Gamma_0, kappa_0, kappa_b
    !
    !   print*, 'Flux at boundary: ', L_0/(4*dpi*R_star**2)
    ! endif

    StefBoltz = const_rad_a*const_c/4.d0*(unit_temperature**&
       4.d0)/(unit_velocity*unit_pressure)

    L_vE = 4*dpi*R_star**2*v_arr(nghostcells+1)*4.d0/3.d0*Er_arr(nghostcells+&
       1)
    !>Set bottom density from massloss rate

    dinflo = inflow_density*M_dot/(4*dpi*R_star**2)
    gradE = -dinflo*kappa_b*(L_0-L_vE)/(4*dpi*R_star**2*const_c/unit_velocity)

    ! print*, dinflo*unit_density
  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ inflow_density, inflow_gamma, init_struc, pert_ampl

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wind_list, end=113)
113    close(unitpar)
    end do

  end subroutine params_read


  subroutine ReadInParams(M_star,R_star,Gamma_0,M_dot_ratio,M_dot,L_0,rho_base,&
     error_b)
    use mod_global_parameters
    double precision, intent(out) :: M_star,R_star,Gamma_0
    double precision, intent(out) :: M_dot_ratio,M_dot,L_0
    double precision, intent(out) :: rho_base, error_b
    character :: dum
    integer :: line, len

    len = len_trim(init_struc)
    print*, 'Structure from: ','InputStan/params_'//init_struc(1:len)//'.txt'


    OPEN(1,FILE='InputStan/params_'//init_struc(1:len)//'.txt')
    READ(1,*) dum, Gamma_0
    READ(1,*) dum, M_dot_ratio
    READ(1,*) dum, M_star
    READ(1,*) dum, L_0
    READ(1,*) dum, R_star
    READ(1,*)
    READ(1,*) dum, M_dot
    READ(1,*) dum, rho_base
    READ(1,*) dum, error_b
    CLOSE(1)

    M_star = M_star*M_sun
    L_0 = L_0*L_sun
    R_star = R_star*R_sun
    M_dot = M_dot*M_sun/year

  end subroutine ReadInParams

  subroutine ReadInTable(r_arr,rho_arr,v_arr,e_arr,Er_arr,T_arr,p_arr)
    use mod_global_parameters
    ! use mod_constants
    ! use mod_fld

    integer :: i
    double precision :: dr
    double precision, intent(out) :: r_arr(domain_nx1+2*nghostcells)
    double precision, intent(out) :: rho_arr(domain_nx1+2*nghostcells)
    double precision, intent(out) :: v_arr(domain_nx1+2*nghostcells)
    double precision, intent(out) :: e_arr(domain_nx1+2*nghostcells)
    double precision, intent(out) :: Er_arr(domain_nx1+2*nghostcells)

    double precision :: i_arr(domain_nx1+2*nghostcells)
    double precision, intent(out) :: T_arr(domain_nx1+2*nghostcells)
    double precision, intent(out) :: p_arr(domain_nx1+2*nghostcells)
    integer :: len

    len = len_trim(init_struc)

    print*, 'Parameters from: ', 'InputStan/structure_amrvac_'//init_struc(1:len)//'.txt'

    OPEN(1,FILE='InputStan/structure_amrvac_'//init_struc(1:len)//'.txt')
    do i = 1,domain_nx1+2*nghostcells
      READ(1,*) r_arr(i),v_arr(i),rho_arr(i),Er_arr(i)
    enddo
    CLOSE(1)

    T_arr = (Er_arr/const_rad_a)**0.25d0
    p_arr = kb_cgs/(fld_mu*mp_cgs)*T_arr*rho_arr
    e_arr = p_arr/(rhd_gamma - one) + half*rho_arr*v_arr**2

  end subroutine ReadInTable

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x)

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    integer :: i,j,b
    integer :: NumberOfBlocks
    double precision :: x_perc
    double precision :: kappa(ixOmin1:ixOmax1), lambda(ixOmin1:ixOmax1),&
        fld_R(ixOmin1:ixOmax1)
    double precision :: pert(ixOmin1:ixOmax1)

    NumberOfBlocks = domain_nx1/block_nx1

     x_perc = (x(ixOmin1,1)-xprobmin1)/(xprobmax1-xprobmin1)
    
    b = floor(x_perc*NumberOfBlocks)

    do i = ixImin1,ixImax1
      j = i + b*block_nx1

      
      w(i,rho_) = rho_arr(j)
      w(i,mom(:)) = zero
      w(i,mom(1)) = rho_arr(j)*v_arr(j)
      w(i,e_) = e_arr(j)
      w(i,r_e) = Er_arr(j)
     

      
    enddo

    call fld_get_opacity(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, lambda,&
        fld_R)

    w(ixOmin1:ixOmax1,i_diff_mg) = (const_c/unit_velocity)*lambda(&
       ixOmin1:ixOmax1)/(kappa(ixOmin1:ixOmax1)*w(ixOmin1:ixOmax1,rho_))

    ! pert(ixO^S) = dsin(2*dpi*x(ixO^S,1)/(xprobmax1-xprobmin1))*dsin(2*dpi*x(ixO^S,2))

    ! where ((x(ixO^S,2) .lt. 3.d0) .and. (x(ixO^S,2) .gt. 1.1d0))
      ! w(ixO^S,rho_) = w(ixO^S,rho_) * (1.d0 + 0.1d0*pert(ixO^S))
    ! end where

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixImin1,ixImax1,ixBmin1,ixBmax1,iB,w,x)
    use mod_global_parameters
    use mod_fld


    integer, intent(in)             :: ixImin1,ixImax1, ixBmin1,ixBmax1, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    double precision :: kappa(ixImin1:ixImax1), gradE_l(ixBmin1:ixBmax1),&
        L_vE_l(ixBmin1:ixBmax1)
    double precision :: Temp(ixImin1:ixImax1), pth(ixImin1:ixImax1),&
       pert(ixBmin1:ixBmax1)
    double precision :: tau_out, T_out, F_bound
    integer :: i,j

    select case (iB)

    case(1)

       w(nghostcells,rho_) = dinflo
      

      do i = ixBmax1-1,ixBmin1,-1
        ! w(ixBmin1:ixBmax1,i,rho_) = 2*w(ixBmin1:ixBmax1,i+1,rho_) - w(ixBmin1:ixBmax1,i+2,rho_)
         w(i,rho_) = dexp(2*dlog(w(i+1,rho_)) - dlog(w(i+2,rho_)))
        
      enddo

      w(ixBmin1:ixBmax1,mom(:)) = 0.d0

      do i = ixBmax1,ixBmin1,-1
         w(i,mom(1)) = w(i+1,mom(1))
        
      enddo

        ! pert(ixB^S) = dsin(3*2*dpi*x(ixB^S,2)/(xprobmax2-xprobmin2))*dsin(2*dpi*x(ixB^S,1)-3*2*dpi*global_time)
        ! w(ixB^S,rho_) = w(ixB^S,rho_) * (1.d0 + pert_ampl*pert(ixB^S))
        ! w(ixB^S,mom(1)) = w(ixB^S,mom(1)) * (1.d0 + pert_ampl*pert(ixB^S))

      !
      ! where (w(ixB^S,mom(1)) .lt. zero)
      !    w(ixB^S,mom(1)) = zero ! abs(w(ixB^S,mom(2)))
      ! end where

      ! print*, w(5,1:5,mom(2))

      call fld_get_opacity(w, x, ixImin1,ixImax1, ixImin1,ixImax1, kappa)

      L_vE_l(ixBmin1:ixBmax1) = 4*dpi*x(ixBmin1:ixBmax1,&
         1)**2*4.d0/3.d0*w(ixBmin1:ixBmax1,mom(1))/w(ixBmin1:ixBmax1,&
         rho_)*w(ixBmin1:ixBmax1,r_e)
      gradE_l(ixBmin1:ixBmax1) = -w(ixBmin1:ixBmax1,&
         rho_)*kappa(ixBmin1:ixBmax1)*(3.d0*unit_velocity/const_c)*(L_0-&
         L_vE_l(ixBmin1:ixBmax1))/(4.d0*dpi*x(ixBmin1:ixBmax1,1)**2.d0)

       gradE = gradE_l(nghostcells)
      

      do i = ixBmax1-1,ixBmin1,-1
         w(i,r_e) = w(i+2,r_e) + (x(i,1)-x(i+2,1))*gradE_l(i+1)
        
      enddo

       w(nghostcells,r_e) = dexp(half*(dlog(w(nghostcells-1,&
          r_e))+dlog(w(nghostcells+1,r_e))))
      

      Temp(ixBmin1:ixBmax1) = (w(ixBmin1:ixBmax1,&
         r_e)*unit_pressure/const_rad_a)**0.25d0
      pth(ixBmin1:ixBmax1) = Temp(ixBmin1:ixBmax1)*w(ixBmin1:ixBmax1,&
         rho_)*const_kb/(const_mp*fld_mu)*unit_density/unit_pressure
      w(ixBmin1:ixBmax1,e_) = pth(ixBmin1:ixBmax1)/(rhd_gamma-1) + &
         half*sum(w(ixBmin1:ixBmax1, mom(:))**2, dim=ndim+1)/w(ixBmin1:ixBmax1,&
         rho_)

    case(2)
      do i = ixBmin1,ixBmax1
        !> Conserve gradE/rho
         w(i,r_e) = (x(i-1,1)**2*w(i-1,mom(1))/w(i-1,rho_))/(x(i,1)**2*w(i,&
            mom(1))/w(i,rho_))*w(i-1,r_e)
        
      enddo

      ! gradE_out = sum((w(ixBmin1:ixBmax1,ixBmin2-1,r_e)-w(ixBmin1:ixBmax1,ixBmin2-2,r_e))&
      ! /(x(ixBmin1:ixBmax1,ixBmin2-1,2) - x(ixBmin1:ixBmax1,ixBmin2-2,2)))/(ixBmax1-ixBmin1)

      
      gradE_out = (w(ixBmin1,r_e)-w(ixBmin1-1,r_e))/(x(ixBmin1,&
         1) - x(ixBmin1-1,1))
     
      

      F_bound = L_0/(4*dpi*R_star**2)
      tau_out = kappa_0*w(ixImax1-nghostcells,&
         rho_)*R_star**2/(3*x(ixImax1-nghostcells,1))
      T_out = F_bound/StefBoltz*(3.d0/4.d0*tau_out)**0.25d0

      T_out = max(T_out,1.d3/unit_temperature)

      ! print*, w(ixImax1-nghostcells,rho_), tau_out, T_out*unit_temperature

      valE_out = const_rad_a*(T_out*unit_temperature)**4/unit_pressure

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
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
      mg%bc(iB, mg_iphi)%bc_value = gradE

    case (2)

      ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
      ! mg%bc(iB, mg_iphi)%bc_value = gradE_out

      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = valE_out

    case default
      print *, "Not a standard: ", trim(typeboundary(r_e, iB))
      error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions


  !> Calculate gravitational acceleration in each dimension
  subroutine set_gravitation_field(ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,x,&
     gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw)
    double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,ndim)

    double precision :: radius(ixImin1:ixImax1)
    double precision :: mass

    radius(ixOmin1:ixOmax1) = x(ixOmin1:ixOmax1,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    gravity_field(ixImin1:ixImax1,:) = zero
    gravity_field(ixImin1:ixImax1,1) = -const_G*mass/radius(ixImin1:ixImax1)**&
       2*(unit_time**2/unit_length)

  end subroutine set_gravitation_field

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine PseudoPlanar(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,&
     wCT,qt,w,x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,&
       iwmax
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision :: ppsource(ixOmin1:ixOmax1,1:nw)

    call PseudoPlanarSource(ixImin1,ixImax1,ixOmin1,ixOmax1,wCT,x,ppsource,&
       .false.)
    w(ixOmin1:ixOmax1,rho_) = w(ixOmin1:ixOmax1,&
       rho_) + qdt*ppsource(ixOmin1:ixOmax1,rho_)
    w(ixOmin1:ixOmax1,mom(:)) = w(ixOmin1:ixOmax1,&
       mom(:)) + qdt*ppsource(ixOmin1:ixOmax1,mom(:))
    w(ixOmin1:ixOmax1,e_) = w(ixOmin1:ixOmax1,&
       e_) + qdt*ppsource(ixOmin1:ixOmax1,e_)
    w(ixOmin1:ixOmax1,r_e) = w(ixOmin1:ixOmax1,&
       r_e) + qdt*ppsource(ixOmin1:ixOmax1,r_e)

  end subroutine PseudoPlanar

  !> internal boundary, user defined
  !> This subroutine can be used to artificially overwrite ALL conservative
  !> variables in a user-selected region of the mesh, and thereby act as
  !> an internal boundary region. It is called just before external (ghost cell)
  !> boundary regions will be set by the BC selection. Here, you could e.g.
  !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
  !> which can be used to identify the internal boundary region location.
  !> Its effect should always be local as it acts on the mesh.
  subroutine reset_egas(level,qt,ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)

    double precision :: Temp(ixImin1:ixImax1),pth(ixImin1:ixImax1)

    if (rhd_energy_interact) call mpistop(&
       'Resetting e_gas but e_interact true')

    Temp(ixImin1:ixImax1) = (w(ixImin1:ixImax1,&
       r_e)*unit_pressure/const_rad_a)**0.25d0/unit_temperature
    pth(ixImin1:ixImax1) = Temp(ixImin1:ixImax1)*w(ixImin1:ixImax1,rho_)
    w(ixImin1:ixImax1,e_) = pth(ixImin1:ixImax1)/(rhd_gamma-1.d0) + &
       half*sum(w(ixImin1:ixImax1, mom(:))**2, dim=ndim+1)/w(ixImin1:ixImax1,&
        rho_)

  end subroutine reset_egas


  subroutine PseudoPlanarSource(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,source,&
     boundary)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)           :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in)  :: w(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    double precision, intent(out) :: source(ixOmin1:ixOmax1,1:nw)
    logical, intent(in) :: boundary

    double precision :: rad_flux(ixOmin1:ixOmax1,1:ndir)
    double precision :: pth(ixImin1:ixImax1),v(ixOmin1:ixOmax1,2)
    double precision :: radius(ixOmin1:ixOmax1),  pert(ixOmin1:ixOmax1)
    integer :: rdir, pdir

    source(ixOmin1:ixOmax1,1:nw) = zero

    rdir = 1
    pdir = 2

    v(ixOmin1:ixOmax1,rdir) = w(ixOmin1:ixOmax1,mom(rdir))/w(ixOmin1:ixOmax1,&
       rho_)
    

    radius(ixOmin1:ixOmax1) = x(ixOmin1:ixOmax1,rdir) !+ half*block%dx(ixOmin1:ixOmax1,rdir)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    source(ixOmin1:ixOmax1,rho_) = -two*w(ixOmin1:ixOmax1,&
       rho_)*v(ixOmin1:ixOmax1,rdir)/radius(ixOmin1:ixOmax1)

    ! pth(ixO^S) = (rhd_gamma-1) &
    ! *(w(ixO^S,e_) - half*sum(w(ixO^S, mom(:))**2, dim=ndim+1) / w(ixO^S, rho_))

    call phys_get_pthermal(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,pth)

    !> dm_r/dt = +(rho*v_p**2 + 2pth)/r -2 (rho*v_r**2 + pth)/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    source(ixOmin1:ixOmax1,mom(rdir)) = - 2*w(ixOmin1:ixOmax1,&
       rho_)*v(ixOmin1:ixOmax1,rdir)**two/radius(ixOmin1:ixOmax1) 

    

    !> de/dt = -2 (e+p)v_r/r
    source(ixOmin1:ixOmax1,e_) = - two*(w(ixOmin1:ixOmax1,&
       e_)+pth(ixOmin1:ixOmax1))*v(ixOmin1:ixOmax1,&
       rdir)/radius(ixOmin1:ixOmax1)

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      if (boundary) then
        rad_flux(ixOmin1:ixOmax1,rdir) = L_0/(4.d0*dpi*x(ixOmin1:ixOmax1,&
           rdir)**2)
      else
        call fld_get_radflux(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, rad_flux)
      endif
      source(ixOmin1:ixOmax1,r_e) = source(ixOmin1:ixOmax1,&
         r_e) - two*rad_flux(ixOmin1:ixOmax1,rdir)/radius(ixOmin1:ixOmax1)
    endif

    if (rhd_radiation_advection) then
      source(ixOmin1:ixOmax1,r_e) = source(ixOmin1:ixOmax1,&
         r_e) - two*w(ixOmin1:ixOmax1,r_e)*v(ixOmin1:ixOmax1,&
         rdir)/radius(ixOmin1:ixOmax1)
    endif

    ! Not sure about this one
    if (rhd_radiation_force) then
      source(ixOmin1:ixOmax1,r_e) = source(ixOmin1:ixOmax1,&
         r_e) + two*v(ixOmin1:ixOmax1,rdir)*w(ixOmin1:ixOmax1,&
         r_e)/(3*radius(ixOmin1:ixOmax1))
    endif

  end subroutine PseudoPlanarSource


  subroutine Opacity_stepfunction(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,kappa)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: w(ixImin1:ixImax1,1:nw), x(ixImin1:ixImax1,&
       1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1)

    kappa(ixOmin1:ixOmax1) = kappa_b + (1.d0+erf((x(ixOmin1:ixOmax1,&
       1)-one)*error_b-error_b/2.d0))*(kappa_0-kappa_b)/2.d0

  end subroutine Opacity_stepfunction

  subroutine refine_base(igrid,level,ixGmin1,ixGmax1,ixmin1,ixmax1,qt,w,x,&
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
    refine = -1

    if (qt .gt. 1.d-1) then
      print*, 'refine 1'
      if (any(x(ixGmin1:ixGmax1,1) < 1.d0)) refine=1
    endif

    if (qt .gt. 2.d-1) then
      print*, 'refine 2'
      if (any(x(ixGmin1:ixGmax1,1) < 2.d0)) refine=1
    endif

    if (qt .gt. 4.d-1) then
      print*, 'refine 3'
      if (any(x(ixGmin1:ixGmax1,1) < 3.d0)) refine=1
    endif

  end subroutine refine_base


  subroutine time_average_values(ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)

    if (global_time .gt. 0.5d0) then
      w(ixImin1:ixImax1,int_r) = w(ixImin1:ixImax1,int_r) + w(ixImin1:ixImax1,&
         rho_)*dt
      w(ixImin1:ixImax1,int_v) =  w(ixImin1:ixImax1,int_v) + w(ixImin1:ixImax1,&
         mom(1))/w(ixImin1:ixImax1,rho_)*dt
      w(ixImin1:ixImax1,int_e) = w(ixImin1:ixImax1,int_e) + w(ixImin1:ixImax1,&
         e_)*dt
      w(ixImin1:ixImax1,int_re) = w(ixImin1:ixImax1,&
         int_re) + w(ixImin1:ixImax1,r_e)*dt

      w(ixImin1:ixImax1,int_dt) =  w(ixImin1:ixImax1,int_dt) + dt
    else
      w(ixImin1:ixImax1,int_r) = zero
      w(ixImin1:ixImax1,int_v) = zero
      w(ixImin1:ixImax1,int_e) = zero
      w(ixImin1:ixImax1,int_re) = zero
      w(ixImin1:ixImax1,int_dt) = zero
    endif

  end subroutine time_average_values

  subroutine specialvar_output(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,normconv)
    ! this subroutine can be used in convert, to add auxiliary variables to the
    ! converted output file, for further analysis using tecplot, paraview, ....
    ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
    !
    ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
    ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_physics
    use mod_fld

    integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
    double precision                   :: w(ixImin1:ixImax1,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision                   :: g_rad(ixImin1:ixImax1),&
        big_gamma(ixImin1:ixImax1)
    double precision                   :: g_grav(ixImin1:ixImax1)
    double precision                   :: Tgas(ixImin1:ixImax1),&
       Trad(ixImin1:ixImax1)
    double precision                   :: kappa(ixOmin1:ixOmax1),&
        OPAL(ixOmin1:ixOmax1), CAK(ixOmin1:ixOmax1)
    double precision                   :: vel(ixImin1:ixImax1),&
        gradv(ixOmin1:ixOmax1)
    double precision                   :: rad_flux(ixOmin1:ixOmax1,1:ndim),&
        Lum_cmf(ixOmin1:ixOmax1)
    double precision                   :: Lum_adv(ixOmin1:ixOmax1),&
        Lum_obs(ixOmin1:ixOmax1)
    double precision                   :: Lum_wnd(ixOmin1:ixOmax1),&
        Mdot(ixOmin1:ixOmax1)
    double precision                   :: pp_rf(ixOmin1:ixOmax1),&
        lambda(ixOmin1:ixOmax1), fld_R(ixOmin1:ixOmax1)
    integer                            :: idim
    double precision :: radius(ixImin1:ixImax1)
    double precision :: mass

    radius(ixOmin1:ixOmax1) = x(ixOmin1:ixOmax1,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    call fld_get_opacity(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, kappa)
    call fld_get_radflux(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, rad_flux)

    g_rad(ixOmin1:ixOmax1) = kappa(ixOmin1:ixOmax1)*rad_flux(ixOmin1:ixOmax1,&
       1)/(const_c/unit_velocity)
    g_grav(ixOmin1:ixOmax1) = const_G*mass/radius(ixOmin1:ixOmax1)**&
       2*(unit_time**2/unit_length)
    big_gamma(ixOmin1:ixOmax1) = g_rad(ixOmin1:ixOmax1)/g_grav(&
       ixOmin1:ixOmax1)

    Mdot(ixOmin1:ixOmax1) =  4*dpi*w(ixOmin1:ixOmax1,&
       mom(1))*radius(ixOmin1:ixOmax1)**2 &
       *unit_density*unit_velocity/M_sun*year

    call rhd_get_tgas(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, Tgas)
    call rhd_get_trad(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, Trad)

    vel(ixImin1:ixImax1) = w(ixImin1:ixImax1,mom(1))/w(ixImin1:ixImax1,rho_)

    pp_rf(ixOmin1:ixOmax1) = two*rad_flux(ixOmin1:ixOmax1,1)/x(ixOmin1:ixOmax1,&
       1)*dt

    call fld_get_fluxlimiter(w, x, ixImin1,ixImax1, ixOmin1,ixOmax1, lambda,&
        fld_R)

    Lum_cmf(ixOmin1:ixOmax1) = 4*dpi*rad_flux(ixOmin1:ixOmax1,&
       1)*(x(ixOmin1:ixOmax1,1)*unit_length)**2*unit_radflux/L_sun
    Lum_adv(ixOmin1:ixOmax1) = 4*dpi*4.d0/3.d0*vel(ixOmin1:ixOmax1)*w(&
       ixOmin1:ixOmax1,r_e)*(x(ixOmin1:ixOmax1,&
       1)*unit_length)**2*unit_radflux/L_sun
    Lum_obs(ixOmin1:ixOmax1) = Lum_cmf(ixOmin1:ixOmax1) + &
       Lum_adv(ixOmin1:ixOmax1)
    Lum_wnd(ixOmin1:ixOmax1) =  4*dpi*w(ixOmin1:ixOmax1,&
       mom(1))*radius(ixOmin1:ixOmax1)**2*unit_density*unit_velocity/L_sun*((&
       vel(ixOmin1:ixOmax1)*unit_velocity)**2/2.d0 - &
       const_G*mass/radius(ixOmin1:ixOmax1) + const_G*mass/unit_length)

    w(ixOmin1:ixOmax1,nw+1) = kappa(ixOmin1:ixOmax1)*unit_opacity
    w(ixOmin1:ixOmax1,nw+2) = rad_flux(ixOmin1:ixOmax1,1)
    w(ixOmin1:ixOmax1,nw+3) = Trad(ixOmin1:ixOmax1)*unit_temperature
    w(ixOmin1:ixOmax1,nw+4) = big_gamma(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+5) = Mdot(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+6) = vel(ixOmin1:ixOmax1)*unit_velocity
    w(ixOmin1:ixOmax1,nw+7) = pp_rf(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+8) = Lum_cmf(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+9) = Lum_adv(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+10) = Lum_obs(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+11) = Lum_wnd(ixOmin1:ixOmax1)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

               !9     10 11   12    13   14  15    16     17    18    19    20
    varnames = 'kappa F1 Trad Gamma Mdot vel pp_rf L_cmf L_adv L_obs L_wnd'
  end subroutine specialvarnames_output

end module mod_usr
