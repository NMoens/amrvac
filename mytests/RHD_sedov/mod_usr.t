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
    {^IFONED call set_coordinate_system("Cartesian_1D")}
    {^IFTWOD call set_coordinate_system("Cartesian_2D")}
    {^IFTHREED call set_coordinate_system("Cartesian_3D")}

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
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas,phys_get_trad
    use mod_fld

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: Temp(ixI^S)
    double precision :: local_rad_e(ixI^S)
    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)
    double precision :: radius(ixI^S), Tgas(ixI^S), dirac(ixI^S)

    double precision :: rho_a, v_a, p_a, T_a, chi_a
    integer :: ix^D

    radius(ixI^S) = dsqrt(sum(x(ixI^S,:)**2,dim=ndim+1))
    w(ixI^S,i_r) = radius(ixI^S)
    dirac(ixI^S) = 0.d0

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
    w(ixI^S,rho_) = g0*radius(ixI^S)**(-kp)
    w(ixI^S,mom(:)) = zero
    ! w(ixI^S,e_) = small_pressure/(rhd_gamma-1.d0) + E0*dirac(ixI^S)
    w(ixI^S,e_) = w(ixI^S,rho_)*0.5d0/(rhd_gamma-1)
    call rhd_get_tgas(w, x, ixI^L, ixO^L, Tgas)
    w(ixI^S,r_e) = const_rad_a*(Tgas(ixI^S)*unit_temperature)**4/unit_pressure

    {do ix^D=ixOmin^D,ixOmax^D\ }
      if (radius(ix^D) .lt. 0.27724) then
        rho_a = read_initial_conditions(radius(ix^D),2)
        v_a = read_initial_conditions(radius(ix^D),3)
        p_a = read_initial_conditions(radius(ix^D),4)
        T_a = read_initial_conditions(radius(ix^D),5)
        ! chi_a = read_initial_conditions(radius(ix^D),6)

        w(ix^D,rho_) = rho_a
        w(ix^D,mom(:)) = rho_a*v_a*x(ix^D,:)/radius(ix^D)
        w(ix^D,e_) = p_a/(rhd_gamma - 1) + rho_a*v_a**2/2
        w(ix^D,r_e) = const_rad_a*T_a**4 * (unit_temperature**4/unit_pressure)
        ! w(ix^D,i_chi) = chi_a
      endif
    {enddo\ }

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)
    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))

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
    var = w_mo(index) + (w_po(index) - w_mo(index))/(w_po(1) - w_mo(1))*(r_in - w_mo(1))

  end function read_initial_conditions

  subroutine kramers_opacity(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    use mod_physics, only: phys_get_tgas
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    double precision :: Temp(ixI^S)

    call phys_get_tgas(w,x,ixI^L,ixO^L,Temp)

    ! where (Temp(ixO^S) .lt. 1.d0) &
    !   Temp(ixO^S) = 1.d0

    kappa(ixO^S) = kap0 * w(ixO^S,rho_)**m * Temp(ixO^S)**(-ndim)
    kappa(ixO^S) = kappa(ixO^S)/w(ixO^S,rho_) !> KRUMHOLZ DEF OF KAPPA is muy kappa*rho

  end subroutine kramers_opacity

  subroutine floor_diffcoef(w, wCT, x, ixI^L, ixO^L)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: wCT(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)


    where (w(ixO^S,i_diff_mg) .lt. 1.d-3) &
      w(ixO^S,i_diff_mg) = 1.d-3

  end subroutine floor_diffcoef


  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: idir
    double precision                   :: Tgas(ixI^S),Trad(ixI^S)
    double precision                   :: kappa(ixO^S)
    double precision                   :: rad_flux(ixO^S,1:ndim), heat_flux(ixO^S,1:ndim)
    double precision                   :: grad_T(ixI^S)

    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
    call rhd_get_tgas(w, x, ixI^L, ixO^L, Tgas)
    call rhd_get_trad(w, x, ixI^L, ixO^L, Trad)
    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)


    w(ixO^S,i_kap) = kappa(ixO^S)
    w(ixO^S,i_Tg) = Tgas(ixO^S)*unit_temperature
    w(ixO^S,i_Tr) = Trad(ixO^S)*unit_temperature
    w(ixO^S,i_chi) = w(ixO^S,rho_)**a * Tgas(ixO^S)**b
    w(ixO^S,i_mom) = dsqrt(w(ixO^S,mom(1))**2 + w(ixO^S,mom(2))**2 + w(ixO^S,mom(3))**2)
    w(ixO^S,i_v) = w(ixO^S,i_mom)/w(ixO^S,rho_)*unit_velocity
    w(ixO^S,i_F) = dsqrt(rad_flux(ixO^S,1)**2+ rad_flux(ixO^S,2)**2+ rad_flux(ixO^S,3)**2)*unit_radflux

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
