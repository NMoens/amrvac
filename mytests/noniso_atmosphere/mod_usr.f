!> This is a template for a new user problem
module mod_usr

! Include a physics module
use mod_rhd
use mod_fld
use mod_global_parameters

implicit none

  integer, parameter :: nyc = 72
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

  ! Routine for setting radiation boundary conditions
  usr_radiation_bc => radiation_boundary

  ! Graviatational field
  usr_gravity => set_gravitation_field

  ! Attempt to use an innerbound to fix trouble with boundaries
  ! usr_internal_bc => fixleftright

  ! Output routines
  usr_aux_output    => specialvar_output
  usr_add_aux_names => specialvarnames_output

  ! Active the physics module
  call rhd_activate()

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

  mdot = 1.d-6

  mdot = mdot*year/M_sun
  vinf = vinf*1.d5

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

  vel_is = mdot/(4*dpi*R_core**2*rho_is*unit_density)/unit_velocity
  tr_is  = (er_is/const_rad_a*unit_pressure/unit_temperature**4)**(1.0/4.0)

  if (mype .eq. one) then
  print*, 'unit_length', unit_length
  print*, 'unit_density', unit_density
  print*, 'unit_pressure', unit_pressure
  print*, 'unit_temperature', unit_temperature
  endif

  if (minval(y_is) .gt. xprobmin2) call mpistop(&
     "Simulation space not covered")
  if (maxval(y_is) .lt. xprobmax2) call mpistop(&
     "Simulation space not covered")

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
    w(ixGmin1: ixGmax1, i, mom(1)) = vel_is(j)*rho_is(j)
    w(ixGmin1: ixGmax1, i, mom(2)) = zero
    w(ixGmin1: ixGmax1, i, e_) = pg_is(j)/(rhd_gamma-1.0)
    w(ixGmin1: ixGmax1, i, r_e) = er_is(j)
  enddo

  !> perturb rho
  amplitude = 0.5d-2
  call RANDOM_NUMBER(pert)

  where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
     2) .lt. 10d-3) .and. (x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,2) .gt. 5d-3))
  where((x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
     1) .lt. -5d-3) .and. (x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1) .gt. 5d-3))
    w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = w(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2, rho_)*(one + amplitude*pert(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2))
  endwhere
  endwhere

  call get_rad_extravars(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
     ixmax1,ixmax2)

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

  integer i,j

  select case (iB)
  case(1)
    call radiation_boundary(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,iB,w,w_rad,x)
    do i=ixGmin1,ixGmin1+nghostcells-1
      w(i,:,rho_) = rho_is(i)
      w(i,:,mom(:)) = zero
      w(i,:,e_) = pg_is(i)/(rhd_gamma-1.0)
      w(i,:,r_e) = w_rad(i,:)
    enddo

  case(2)
    call radiation_boundary(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,iB,w,w_rad,x)
    do i=ixGmax1-nghostcells+1,ixGmax1
      w(i,:,rho_) = rho_is(i)
      w(i,:,mom(:)) = zero
      w(i,:,e_) = pg_is(i)/(rhd_gamma-1.0)
      w(i,:,r_e) = w_rad(i,:)
    enddo

  case(3)
    do j = ixBmin2,ixBmax2
      w(ixGmin1:ixGmax1,j,rho_) = w(ixGmin1:ixGmax1,ixBmax2+1,rho_)
      w(ixGmin1:ixGmax1,j,mom(1)) = w(ixGmin1:ixGmax1,ixBmax2+1,mom(1))
      w(ixGmin1:ixGmax1,j,mom(2)) = w(ixGmin1:ixGmax1,ixBmax2+1,mom(2))
      w(ixGmin1:ixGmax1,j,e_) = w(ixGmin1:ixGmax1,ixBmax2+1,e_)
    enddo

    call radiation_boundary(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,iB,w,w_rad,x)
    do j = ixBmin2,ixBmax2
      w(ixGmin1:ixGmax1,j,r_e) = w_rad(ixGmin1:ixGmax1,j)
    enddo

  case(4)
    do j = ixBmin2,ixBmax2
      w(ixGmin1:ixGmax1,j,rho_) = w(ixGmin1:ixGmax1,ixBmin2-1,rho_)
      w(ixGmin1:ixGmax1,j,mom(1)) = w(ixGmin1:ixGmax1,ixBmin2-1,mom(1))
      w(ixGmin1:ixGmax1,j,mom(2)) = w(ixGmin1:ixGmax1,ixBmin2-1,mom(2))
      w(ixGmin1:ixGmax1,j,e_) = w(ixGmin1:ixGmax1,ixBmin2-1,e_)
    enddo

    call radiation_boundary(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,iB,w,w_rad,x)
    do j = ixBmin2,ixBmax2
      w(ixGmin1:ixGmax1,j,r_e) = w_rad(ixGmin1:ixGmax1,j)
    enddo

  case default
    call mpistop('boundary not known')
  end select
end subroutine boundary_conditions

!==========================================================================================

subroutine radiation_boundary(qt,ixImin1,ixImin2,ixImax1,ixImax2,iB,w,w_rad,x)
  use mod_global_parameters
  integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, iB
  double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim)
  double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
  double precision, intent(out)   :: w_rad(ixImin1:ixImax1,ixImin2:ixImax2)

  integer i,j

  select case (iB)
  case(1)
    do i=ixImin2,ixImin2+nghostcells-1
      w_rad(ixImin1:ixImin1+nghostcells-1,i) = w(ixImax2-nghostcells,i,r_e)
    enddo

  case(2)
    do i=ixImin2,ixImin2+nghostcells-1
      w_rad(ixImax1-nghostcells+1:ixImax1,i) = w(ixImax1+nghostcells,i,r_e)
    enddo


  case(3)
    do i=ixImin1,ixImax1
      ! w_rad(i,2) = 2.d0*w(i,3,r_e)-w(i,4,r_e)
      ! w_rad(i,1) = 2.d0*w(i,2,r_e)-w(i,3,r_e)
      ! print*, w_rad(i,1), w_rad(i,2)
      w_rad(i,ixImin2:ixImin2+nghostcells-1) = &
         er_is(ixImin2:ixImin2+nghostcells-1)
    enddo

  case(4)
    do i=ixImin1,ixImax1
      do j = ixImax2-nghostcells+1,ixImax2
        w_rad(i,j) = max(2.d0*w(i,j-1,r_e)-w(i,j-2,r_e), zero)
      enddo
    enddo

  case default
    call mpistop('boundary not known')
  end select
end subroutine radiation_boundary

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

!===============================================================================

subroutine fixleftright(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
   ixOmin2,ixOmax1,ixOmax2,w,x)

  use mod_global_parameters
  integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,level
  double precision, intent(in)    :: qt
  double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
  double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)

  integer :: i

  do i = ixOmin2, ixOmax2
    w(ixOmin1:ixOmin1+2,i,iw_r_e) =  er_is(i)
    w(ixOmax1-2:ixOmax1,i,iw_r_e) =  er_is(i)
  enddo

end subroutine fixleftright

!==========================================================================================

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

  integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2
  double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim)
  double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)
  double precision                   :: g_rad(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim), big_gamma(ixImin1:ixImax1,ixImin2:ixImax2)
  integer                            :: idim

  do idim = 1,ndim
    g_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,i_op)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       i_flux(idim))/fld_speedofligt_0
  enddo
  big_gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = g_rad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,2)/(6.67e-8*mstar/rstar**2*(unit_time**2/unit_length))

  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=g_rad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1)*unit_length/(unit_time**2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=g_rad(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,2)*unit_length/(unit_time**2)
  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3)=big_gamma(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2)
end subroutine specialvar_output

subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
  use mod_global_parameters
  character(len=*) :: varnames

  varnames = 'ar1 ar2 Gamma'

end subroutine specialvarnames_output

!==========================================================================================

end module mod_usr

!==========================================================================================
