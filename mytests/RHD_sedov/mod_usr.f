!> This is a template for a new user problem
module mod_usr
  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision :: g0
  double precision :: E0
  double precision :: kap0
  double precision :: a
  double precision :: b
  double precision :: m
  double precision :: kp
  double precision :: alpha
  double precision :: Xi0
  double precision :: Gam
  double precision :: Omega
  double precision :: delta_r

  character(len=16) :: delta_type

  integer :: i_r, i_delta
  integer :: i_Tg, i_Tr
  integer :: i_kap, i_mom, i_v, i_chi, i_F, i_S

contains
  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    
    
     call set_coordinate_system("Cartesian_3D")

    !> Routine describing analytical opacity
    usr_special_opacity => kramers_opacity

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    !> Additional variables
    usr_process_grid => update_extravars

    !> Do not refine near boundaries
    usr_refine_threshold => coarse_bounds

    !> Set minimal value for diffusion coefficient
    usr_special_diffcoef => floor_diffcoef

    ! unit_velocity = 1.d0
    ! unit_numberdensity = 1.d0/2.3432130399999995E-024

    ! Active the physics module
    call rhd_activate()

    !> Set parameters
    call usr_params_read(par_files)

    i_r = var_set_extravar("r", "r")
    i_delta = var_set_extravar("delta","delta")
    i_Tg = var_set_extravar("Tgas", "Tgas")
    i_Tr = var_set_extravar("Trad", "Trad")
    i_kap = var_set_extravar("kappa","kappa")
    i_mom = var_set_extravar("m_r","m_r")
    i_v = var_set_extravar("v_r","v_r")
    i_chi = var_set_extravar("chi","chi")
    i_F = var_set_extravar("F_r","F_r")
    i_S = var_set_extravar("S_r","S_r")

  end subroutine usr_init

  subroutine usr_params_read(files)
    use mod_global_parameters, only: unitpar
    use mod_fld
    use mod_constants

    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /sedov_list/ g0, E0, Xi0, a, b, m, delta_r, delta_type

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, sedov_list, end=111)
       111    close(unitpar)
    end do

    !> dimensionless parameters
    m = -a
    b = ndim + 3
    kp = ((2*b-1)*ndim + 2)/(2*b - 2*a + 1)
    alpha = (2*b - 2*a + 1)/(2*b - (ndim+2)*a + ndim)

    ! g0 = g0/unit_density
    Gam = 1.d0

    !> Calculate kap0
    kap0 = 4.d0*const_c*const_rad_a/(3.d0*Xi0) &
       /(unit_velocity*unit_pressure/unit_temperature**4)

    !> shock strength (dimensionless)
    Omega = 2*Xi0/(Gam**(b+1) * g0**(1-a))*(E0/g0)**(b-1.d0/2.d0)

    if (mype ==0) then
      print*, 'g0', g0
      print*, 'E0', E0
      print*, 'Xi0', Xi0
      print*, 'kp', kp
      print*, 'alpha', alpha
      print*, 'Gamma', Gam
      print*, 'kappa0', kap0
      print*, 'Omega', Omega
    endif

    if (mype ==0) then
      print*, '-------------------------------'
      print*, 'unit_length', unit_length
      print*, 'unit_time', unit_time
      print*, 'unit_velocity', unit_velocity
      print*, 'unit_density', unit_density
      print*, 'unit_numberdensity', unit_numberdensity
      print*, 'unit_pressure', unit_pressure
      print*, 'unit_temperature', unit_temperature
      print*, 'unit_opacity', unit_opacity
    endif

  end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas,phys_get_trad
    use mod_fld

    integer, intent(in)             :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision :: local_rad_e(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)
    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), Tgas(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
        dirac(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    double precision :: rho_a, v_a, p_a, T_a, chi_a
    integer :: ix1,ix2,ix3

    radius(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3) = dsqrt(sum(x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,:)**2,dim=ndim+1))
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       i_r) = radius(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    dirac(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = 0.d0

    ! select case (delta_type)
    ! case('Gauss')
    ! !> Gaussian function to approximate dirac delta
    ! dirac(ixI^S) = 1.d0/dsqrt(2*dpi*delta_r**2)*dexp(-radius(ixI^S)**2/(2*delta_r**2))
    !
    ! case('Box')
    ! !> Box function  to approximate dirac delta
    ! where (radius(ixI^S) .lt. delta_r) &
    !   dirac(ixI^S) = 3.d0/(4*dpi*delta_r**3)
    !
    ! case('Wigner')
    ! !> Wigner function to approximate dirac delta
    ! where (radius(ixI^S) .lt. delta_r) &
    !   dirac(ixI^S) = 2.d0/(dpi*delta_r**2)*dsqrt(delta_r**2 - radius(ixI^S)**2)
    !
    ! case('Voxel')
    ! !> Fil only the eight voxels centered around 0,0,0 to approximate dirac delta
    ! where( radius(ixI^S) .lt. dsqrt(3* (dxlevel(1)/2)**2) + smalldouble) &
    !   dirac(ixI^S) = 1.d0/(2*dxlevel(1))**3
    !
    ! case('Triangle')
    ! !> Use a triangle function to approximate dirac delta
    ! where (radius(ixI^S) .lt. delta_r) &
    !   dirac(ixI^S) = (3.d0/(dpi*delta_r**3)) * (1- radius(ixI^S)/delta_r)
    !
    ! case default
    !   call mpistop("dirac delta not known")
    ! end select
    ! w(ixI^S,i_delta) = dirac(ixI^S)

    !> Set density
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       rho_) = g0*radius(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)**(-kp)
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,mom(:)) = zero
    ! w(ixI^S,e_) = small_pressure/(rhd_gamma-1.d0) + E0*dirac(ixI^S)
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3,rho_)*0.5d0/(rhd_gamma-1)
    call rhd_get_tgas(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, Tgas)
    w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
       r_e) = const_rad_a*(Tgas(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3)*unit_temperature)**4/unit_pressure

    do ix1=ixOmin1,ixOmax1
     do ix2=ixOmin2,ixOmax2
     do ix3=ixOmin3,ixOmax3
    
      if (radius(ix1,ix2,ix3) .lt. 0.27724) then
        rho_a = read_initial_conditions(radius(ix1,ix2,ix3),2)
        v_a = read_initial_conditions(radius(ix1,ix2,ix3),3)/unit_velocity
        p_a = read_initial_conditions(radius(ix1,ix2,ix3),4)
        T_a = read_initial_conditions(radius(ix1,ix2,ix3),5)
        ! chi_a = read_initial_conditions(radius(ix^D),6)

        w(ix1,ix2,ix3,rho_) = rho_a
        w(ix1,ix2,ix3,mom(:)) = rho_a*v_a*x(ix1,ix2,ix3,:)/radius(ix1,ix2,ix3)
        w(ix1,ix2,ix3,e_) = p_a/(rhd_gamma - 1) + rho_a*v_a**2/2
        w(ix1,ix2,ix3,r_e) = const_rad_a*T_a**4 * &
           (unit_temperature**4/unit_pressure)
        ! w(ix^D,i_chi) = chi_a
      endif
    enddo
     enddo
     enddo
    

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappa)
    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, lambda,&
        fld_R)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_diff_mg) = (const_c/unit_velocity)*lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)/(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_))

  end subroutine initial_conditions

  function read_initial_conditions(r_in,index) result(var)
    integer, intent(in) :: index
    double precision, intent(in) :: r_in
    double precision :: var

    double precision :: w(1:5), w_mo(1:5), w_po(1:5)
    integer :: ll

    w(:) = 0.d0

    open(unit=1, file='1D_sedov_in')
    ! read(1,*) !> header
    read(1,*) w !> first line of data
    do ll = 1,526
      w_mo = w
      read(1,*) w
        if (w(1) .gt. r_in) then
          w_po = w
          goto 8765
        endif
      w_po = w
    enddo

8765 CLOSE(1)
    var = w_mo(index) + (w_po(index) - w_mo(index))/(w_po(1) - w_mo(1))*(r_in &
       - w_mo(1))

  end function read_initial_conditions

  subroutine kramers_opacity(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,kappa)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)

    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

    call phys_get_tgas(w,x,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,Temp)

    ! where (Temp(ixO^S) .lt. 1.d0) &
    !   Temp(ixO^S) = 1.d0

    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = kap0 * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,rho_)**m * Temp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)**(-ndim)
    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_) !> KRUMHOLZ DEF OF KAPPA is muy kappa*rho

  end subroutine kramers_opacity

  subroutine floor_diffcoef(w, wCT, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
       ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3, 1:ndim)


    where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_diff_mg) .lt. 1.d-3) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,i_diff_mg) = 1.d-3

  end subroutine floor_diffcoef


  subroutine update_extravars(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixImin1,ixImin2,ixImin3,&
       ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    integer :: idir
    double precision                   :: Tgas(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3),Trad(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
    double precision                   :: kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    double precision                   :: rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndim), heat_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndim)
    double precision                   :: grad_T(ixImin1:ixImax1,&
       ixImin2:ixImax2,ixImin3:ixImax3)

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, kappa)
    call rhd_get_tgas(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, Tgas)
    call rhd_get_trad(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, Trad)
    call fld_get_radflux(w, x, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, rad_flux)


    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_kap) = kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_Tg) = Tgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*unit_temperature
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_Tr) = Trad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)*unit_temperature
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_chi) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)**a * Tgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)**b
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_mom) = dsqrt(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(1))**2 + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(2))**2 + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       mom(3))**2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,i_v) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,i_mom)/w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)*unit_velocity
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i_F) = dsqrt(rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)**2+ rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       2)**2+ rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       3)**2)*unit_radflux

  end subroutine update_extravars


  !> use different threshold in special regions for AMR to
  !> reduce/increase resolution there where nothing/something interesting happens.
  subroutine coarse_bounds(wlocal,xlocal,threshold,qt)
    use mod_global_parameters
    double precision, intent(in)    :: wlocal(1:nw),xlocal(1:ndim),qt
    double precision, intent(inout) :: threshold

    integer :: idim
    double precision :: rat, my_thr

    rat = 0.75d0
    my_thr = 1.d0

    if (xlocal(1) > rat* xprobmax1) threshold = my_thr
    if (xlocal(1) < rat* xprobmin1) threshold = my_thr

    if (xlocal(2) > rat* xprobmax2) threshold = my_thr
    if (xlocal(2) < rat* xprobmin2) threshold = my_thr

    if (xlocal(3) > rat* xprobmax3) threshold = my_thr
    if (xlocal(3) < rat* xprobmin3) threshold = my_thr

  end subroutine coarse_bounds


end module mod_usr
