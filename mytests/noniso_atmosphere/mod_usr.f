!> This is a template for a new user problem
module mod_usr

! Include a physics module
use mod_rhd
use mod_fld
use mod_global_parameters

implicit none

  integer :: int_rho, int_r_e, int_p, int_m1, int_m2, int_t

  integer, parameter :: nyc = 136
  double precision, parameter :: M_sun = 1.99d33
  double precision, parameter :: R_sun = 6.96d10
  double precision, parameter :: year = 365.25*24*60*60


  integer :: i_is(1:nyc)
  double precision :: y_is(1:nyc)
  double precision :: rho_is(1:nyc)
  double precision :: tg_is(1:nyc)
  double precision :: pg_is(1:nyc)
  double precision :: er_is(1:nyc)
  double precision :: tau_is(1:nyc)
  double precision :: vel_is(1:nyc)
  double precision :: tr_is(1:nyc)

  double precision :: gamma0
  double precision :: mstar
  double precision :: rstar
  double precision :: yhe
  integer :: ny_vac
  double precision :: dy_vac
  double precision :: tau_max
  character :: kap_law

  double precision :: M_star, L_star, Gamma_e, R_core, R_up
  double precision :: log_g, T_gas_core, rho_core, tau_core
  double precision :: mdot, vinf, Hp

contains

!> This routine should set user methods, and activate the physics module
subroutine usr_init()
  use mod_global_parameters
  use mod_usr_methods
  use mod_constants

  call set_coordinate_system("Cartesian_2D")

  call initglobaldata_usr


  ! Initialize units
  usr_set_parameters => initglobaldata_usr

  ! A routine for initial conditions is always required
  usr_init_one_grid => initial_conditions

  ! Routine for setting special boundary conditions
  usr_special_bc => boundary_conditions

  ! Graviatational field
  usr_gravity => set_gravitation_field

  ! Get time integrated values
  usr_modify_output => time_average_values

  ! Output routines
  usr_aux_output    => specialvar_output
  usr_add_aux_names => specialvarnames_output

  ! Active the physics module
  call rhd_activate()

  int_rho = var_set_extravar('int_rho','int_rho')
  int_r_e = var_set_extravar('int_r_e','int_r_e')
  int_p = var_set_extravar('int_p','int_p')
  int_m1 = var_set_extravar('int_m1','int_m1')
  int_m2 = var_set_extravar('int_m2','int_m2')
  int_t = var_set_extravar('int_t','int_t')

end subroutine usr_init

!==========================================================================================

subroutine initglobaldata_usr
  use mod_global_parameters

  character(len=20) :: dum1
  double precision :: dum2,dum3,dum4
  integer :: i

  !> Reading in parameters
  OPEN(1,FILE='initial_conditions/init_indat')
   READ(1,*) gamma0,mstar,rstar,yhe,ny_vac,dy_vac,tau_max,kap_law
  CLOSE(1)

  mstar = mstar*M_sun
  rstar = rstar*R_sun

  !> Perform checks to see if init_indat and usr.par agree
  if (ny_vac .ne. nyc) call mpistop('Number of cells dont agree')

  OPEN(2,FILE='initial_conditions/init_struc_amrvac')
  do i = 1,domain_nx2+2*nghostcells
    READ(2,'(1i4,1e20.10,5e18.8)') i_is(i), y_is(i), rho_is(i), tg_is(i),&
        pg_is(i), er_is(i), tau_is(i)
  enddo
  CLOSE(2)

  OPEN(3,FILE='initial_conditions/init_params_amrvac')
    READ(3,*) dum1, M_star
    READ(3,*) dum1, L_star
    READ(3,*) dum1, Gamma_e
    READ(3,*) dum1, R_core
    READ(3,*) dum1, R_up
    READ(3,*)
    READ(3,*) dum1, log_g
    READ(3,*) dum1, T_gas_core
    READ(3,*) dum1, rho_core
    READ(3,*) dum1, Hp
    READ(3,*)
    READ(3,*) dum1, tau_core
  CLOSE(3)

  M_star = M_star*M_sun
  R_up = R_up*R_sun
  Hp = Hp*R_core

  !> Define units
  unit_numberdensity = rho_core/((1.d0+4.d0*He_abundance)*mp_cgs)
  unit_temperature = T_gas_core
  unit_length = Hp

  !> Remaining units
  unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
  unit_pressure=(2.d0+3.d0*He_abundance&
     )*unit_numberdensity*kB_cgs*unit_temperature
  unit_velocity=dsqrt(unit_pressure/unit_density)
  unit_time=unit_length/unit_velocity

  !> Make input dimensionless:
  y_is = y_is/unit_length
  rho_is = rho_is/unit_density
  tg_is = tg_is/unit_temperature
  pg_is = pg_is/unit_pressure
  er_is = er_is/unit_pressure

  y_is = y_is - y_is(1)
  tr_is  = (er_is/const_rad_a*unit_pressure/unit_temperature**4)**(1.0/4.0)

  if (mype .eq. 0) then
    print*, 'unit_length', unit_length
    print*, 'unit_density', unit_density
    print*, 'unit_pressure', unit_pressure
    print*, 'unit_temperature', unit_temperature
    print*, 'unit_radflux', unit_radflux
    print*, 'unit_opacity', unit_opacity
    print*, 'unit time', unit_time

    print*, minval(y_is), xprobmin2
    print*, maxval(y_is), xprobmax2

    if (minval(y_is) .gt. xprobmin2) call &
       mpistop("Simulation space not covered")
    if (maxval(y_is) .lt. xprobmax2) call &
       mpistop("Simulation space not covered")
  endif

end subroutine initglobaldata_usr

!==========================================================================================

!> A routine for specifying initial conditions
subroutine initial_conditions(ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
   ixmax1,ixmax2, w, x)
  use mod_global_parameters
  use mod_constants
  use mod_variables

  integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
     ixmin2,ixmax1,ixmax2
  double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2, ndim)
  double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, nw)

  double precision :: pert(ixGmin1:ixGmax1,ixGmin2:ixGmax2), amplitude,&
      y_res(1:nyc)
  integer :: i,j

  do i = ixGmin2,ixGmax2
    y_res(1:nyc) = y_is(1:nyc)-(x(1+nghostcells,i,2))
    j = minloc(abs(y_res), 1)

    w(ixGmin1: ixGmax1, i, rho_) = rho_is(j)
    w(ixGmin1: ixGmax1, i, mom(1)) = zero
    w(ixGmin1: ixGmax1, i, mom(2)) = zero
    w(ixGmin1: ixGmax1, i, r_e) = er_is(j)
    w(ixGmin1: ixGmax1, i, e_) = pg_is(j)/(rhd_gamma-1.0)
  enddo

  !> perturb rho
  amplitude = 0.0d-1
  call RANDOM_NUMBER(pert)
  do i = ixGmin2+10,ixGmax2
    w(ixGmin1:ixGmax1, i, rho_) = w(ixGmin1:ixGmax1, i,&
        rho_)*(one + amplitude*pert(ixGmin1:ixGmax1, i))
  enddo

  call get_rad_extravars(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
     ixmax1,ixmax2)
  call set_mg_bounds()

end subroutine initial_conditions

!==========================================================================================

subroutine boundary_conditions(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
   ixBmin2,ixBmax1,ixBmax2,iB,w,x)
  use mod_global_parameters
  integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2, iB
  double precision, intent(in)    :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
     1:ndim)
  double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
  double precision                :: w_rad(ixGmin1:ixGmax1,ixGmin2:ixGmax2)

  double precision :: temperature(ixGmin2:ixGmax2)
  double precision :: pressure(ixGmin2:ixGmax2)

  double precision :: y_res(1:nyc)
  integer :: i,j

  select case (iB)

  case(3)
    do i = ixBmin2,ixBmax2
      y_res(1:nyc) = y_is(1:nyc)-(x(ixBmin1+nghostcells,i,2))
      j = minloc(abs(y_res), 1)

      w(ixGmin1:ixGmax1,i,rho_) = rho_is(j)
      w(ixGmin1:ixGmax1,i,e_) = pg_is(j)/(rhd_gamma-1.0)
      w(ixGmin1:ixGmax1,i,r_e) = er_is(j)
    enddo

    ! do i = nghostcells,1,-1
    !   w(ixGmin1:ixGmax1,i,r_e) = w(ixGmin1:ixGmax1,i+1,rho_)/w(ixGmin1:ixGmax1,i+2,rho_) &
    !   *(w(ixGmin1:ixGmax1,i+1,r_e) - w(ixGmin1:ixGmax1,i+2,r_e)) + w(ixGmin1:ixGmax1,i+1,r_e)
    ! enddo

  case(4)
    do i = ixBmin2,ixBmax2
      !> Conserve gradE/rho
      w(ixGmin1:ixGmax1,i,r_e) = w(ixGmin1:ixGmax1,i-1,rho_)/w(ixGmin1:ixGmax1,&
         i-2,rho_) *(w(ixGmin1:ixGmax1,i-1,r_e) - w(ixGmin1:ixGmax1,i-2,&
         r_e)) + w(ixGmin1:ixGmax1,i-1,r_e)

      !> Conserve gradE
      ! w(ixGmin1:ixGmax1,i,r_e) = &
      ! (w(ixGmin1:ixGmax1,i-1,r_e) - w(ixGmin1:ixGmax1,i-2,r_e)) + w(ixGmin1:ixGmax1,i-1,r_e)
      ! do j = ixGmin1,ixGmax1
      !   w(j,i,r_e) = min(w(j,i-1,r_e), w(j,i,r_e))
      ! enddo

      !> Conserve vE + F
      ! w(ixGmin1:ixGmax1,i,r_e) = w(ixGmin1:ixGmax1,i-1,r_e) &
      ! + fld_speedofligt_0/(3*w(ixGmin1:ixGmax1,i,i_op)*w(ixGmin1:ixGmax1,i,rho_)) &
      ! *(w(ixGmin1:ixGmax1,i-1,mom(2))/w(ixGmin1:ixGmax1,i-1,rho_)*w(ixGmin1:ixGmax1,i-1,r_e) &
      ! - w(ixGmin1:ixGmax1,i,mom(2))/w(ixGmin1:ixGmax1,i,rho_)*w(ixGmin1:ixGmax1,i,r_e)) &
      ! + w(ixGmin1:ixGmax1,i-1,r_e) - w(ixGmin1:ixGmax1,i-2,r_e)

    enddo

  case default
    call mpistop('boundary not known')
  end select
end subroutine boundary_conditions

!==========================================================================================

!> Calculate gravitational acceleration in each dimension
subroutine set_gravitation_field(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,wCT,x,gravity_field)
  use mod_global_parameters
  integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2
  double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
  double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
  double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,&
     ixImin2:ixImax2,ndim)

  gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,1) = zero
    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
       2) = -const_G*mstar/rstar**2*(unit_time**2/unit_length)
end subroutine set_gravitation_field

!==========================================================================================

!> If defined, this routine is called before writing output, and it can
!> set/modify the variables in the w array.
subroutine time_average_values(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,qt,w,x)
  use mod_global_parameters
  use mod_physics, only: phys_get_pthermal

  integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2
  double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim)
  double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

  double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2)

  call phys_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,pth)

  if (global_time .eq. 1.d0) w(ixImin1:ixImax1,ixImin2:ixImax2,int_rho) = zero
  if (global_time .eq. 1.d0) w(ixImin1:ixImax1,ixImin2:ixImax2,int_r_e) = zero
  if (global_time .eq. 1.d0) w(ixImin1:ixImax1,ixImin2:ixImax2,int_p) = zero
  if (global_time .eq. 1.d0) w(ixImin1:ixImax1,ixImin2:ixImax2,int_m1) = zero
  if (global_time .eq. 1.d0) w(ixImin1:ixImax1,ixImin2:ixImax2,int_m2) = zero
  if (global_time .eq. 1.d0) w(ixImin1:ixImax1,ixImin2:ixImax2,int_t) = zero

  w(ixImin1:ixImax1,ixImin2:ixImax2,int_rho) = w(ixImin1:ixImax1,&
     ixImin2:ixImax2,int_rho) + w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)*dt
  w(ixImin1:ixImax1,ixImin2:ixImax2,int_r_e) = w(ixImin1:ixImax1,&
     ixImin2:ixImax2,int_r_e) + w(ixImin1:ixImax1,ixImin2:ixImax2,r_e)*dt
  w(ixImin1:ixImax1,ixImin2:ixImax2,int_p) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
     int_p) + pth(ixImin1:ixImax1,ixImin2:ixImax2)*dt
  w(ixImin1:ixImax1,ixImin2:ixImax2,int_m1) = w(ixImin1:ixImax1,&
     ixImin2:ixImax2,int_m1) + w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1))*dt
  w(ixImin1:ixImax1,ixImin2:ixImax2,int_m2) = w(ixImin1:ixImax1,&
     ixImin2:ixImax2,int_m2) + w(ixImin1:ixImax1,ixImin2:ixImax2,mom(2))*dt
  w(ixImin1:ixImax1,ixImin2:ixImax2,int_t) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
     int_t) + dt

end subroutine time_average_values




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
  double precision                   :: g_rad(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim), big_gamma(ixImin1:ixImax1,ixImin2:ixImax2)
  double precision                   :: Tgas(ixImin1:ixImax1,ixImin2:ixImax2),&
     Trad(ixImin1:ixImax1,ixImin2:ixImax2)
  integer                            :: idim

  do idim = 1,ndim
    g_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,i_op)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       i_flux(idim))/fld_speedofligt_0
  enddo
  big_gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = g_rad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,2)/(const_G*mstar/rstar**2*(unit_time**2/unit_length))

  call rhd_get_tgas(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, Tgas)
  call rhd_get_trad(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, Trad)

  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=Tgas(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)*unit_temperature
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=Trad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)*unit_temperature
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3)=big_gamma(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
end subroutine specialvar_output

subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
  use mod_global_parameters
  character(len=*) :: varnames

  varnames = 'Tgas Trad Gamma'
end subroutine specialvarnames_output

!==========================================================================================

end module mod_usr

!==========================================================================================
