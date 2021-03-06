!> This is a template for a new user problem
module mod_usr

! Include a physics module
use mod_rhd
use mod_fld

implicit none

double precision :: M_sun = 1.989d33
double precision :: L_sun = 3.827d33
double precision :: R_sun = 6.96d10

double precision :: M_star
double precision :: L_star
double precision :: T_star
double precision :: R_star

double precision :: Gamma, c_sound, Flux, g_eff, g_grav, H_eff, kappa, c_light

double precision :: Flux0, c_sound0, T_star0, kappa0, T_bound0
double precision :: c_light0, g0, geff0, heff0, Gamma_edd
double precision :: L_star0, R_star0, M_star0
double precision :: tau_bound,  P_bound, rho_bound, T_bound
double precision :: u_p, p_bound0, rho_bound0

double precision :: lower_bc_rho(2), lower_bc_e(2), lower_bc_re(2)

contains

!> This routine should set user methods, and activate the physics module
subroutine usr_init()
  use mod_global_parameters
  use mod_usr_methods
  use mod_constants

  call set_coordinate_system("Cartesian_2D")

  M_star = 80*M_sun !150*M_sun
  L_star = 10**6.2*L_sun !(M_star/M_sun)**3.d0*L_sun
  R_star = 60*R_sun !30*R_sun
  T_star = 9000.d0  !(L_star/(4d0*dpi*R_star**2*5.67051d-5))**0.25d0
  tau_bound = 100.d0

  call usr_params_read(par_files)
  call initglobaldata_usr

  ! Initialize units
  usr_set_parameters => initglobaldata_usr

  ! A routine for initial conditions is always required
  usr_init_one_grid => initial_conditions

  ! Special Boundary conditions
  usr_special_bc => special_bound

  ! Routine for setting radiation boundary conditions
  usr_radiation_bc => radiation_boundary

  ! Graviatational field
  usr_gravity => set_gravitation_field

  ! Output routines
  usr_aux_output    => specialvar_output
  usr_add_aux_names => specialvarnames_output

  ! Active the physics module
  call rhd_activate()

end subroutine usr_init

!===============================================================================

subroutine usr_params_read(files)
  use mod_global_parameters, only: unitpar
  use mod_constants
  character(len=*), intent(in) :: files(:)
  integer                      :: n

  namelist /star_list/ M_star, L_star, R_star, T_star, tau_bound

  do n = 1, size(files)
     open(unitpar, file=trim(files(n)), status="old")
     read(unitpar, star_list, end=111)
     111    close(unitpar)
  end do

  M_star = M_star*M_sun
  L_star = L_star*L_sun
  R_star = R_star*R_sun
  T_star = T_star
  tau_bound = tau_bound

  fld_kappa0 = fld_kappa0/unit_opacity

end subroutine usr_params_read

subroutine initglobaldata_usr
  use mod_global_parameters
  use mod_constants

  T_bound = (one +(3.d0/4.d0*tau_bound)**(1.d0/4.d0))*T_star
  c_sound =  dsqrt((const_kB*T_bound/(0.6*mp_cgs)))

  g_grav = const_G*M_star/R_star**two

  Flux =  L_star/(4.d0*dpi*R_star**2)
  kappa = 0.34d0
  Gamma_edd = (kappa*Flux)/(const_c*g_grav)
  g_eff = g_grav*(one - Gamma_edd)
  H_eff = c_sound**2.d0/g_eff

  p_bound = (one-Gamma_edd)/Gamma_edd*(4.d0*const_sigma)/(3.d0*const_c)*T_bound**4.d0
  rho_bound = p_bound/c_sound**two

  !###############################################

  unit_length        = H_eff
  unit_numberdensity = rho_bound/((1.d0+4.d0*He_abundance)*mp_cgs)
  unit_temperature   = T_bound

  unit_density=(1.d0+4.d0*He_abundance)*mp_cgs*unit_numberdensity
  unit_pressure=(2.d0+3.d0*He_abundance)*unit_numberdensity*const_kB*unit_temperature
  unit_velocity=dsqrt(unit_pressure/unit_density)
  unit_time=unit_length/unit_velocity

  !###############################################

  L_star0 = L_star*unit_time/(unit_pressure*unit_length**3.d0)
  R_star0 = R_star/unit_length
  M_star0 = M_star/(unit_density*unit_length**3.d0)

  Flux0 = L_star0/(4.d0*dpi*R_star0**2.d0)
  T_bound0 = T_bound/unit_temperature
  c_sound0 = dsqrt((const_kB*T_bound/(0.6d0*mp_cgs)))/unit_velocity

  kappa0 = fld_kappa0
  c_light0 = const_c/unit_velocity
  g0 =g_grav*(unit_density*unit_length/unit_pressure)
  geff0 = g0*(one - (kappa0*Flux0)/(c_light0*g0))
  heff0 = c_sound0**2.d0/geff0
  Gamma = (kappa0*Flux0)/(c_light0*g0)

  p_bound0 = p_bound/unit_pressure
  rho_bound0 = p_bound0/c_sound0**two

  print*, L_star0
  print*, R_star0

end subroutine initglobaldata_usr

!==========================================================================================

!> A routine for specifying initial conditions
subroutine initial_conditions(ixG^L, ix^L, w, x)
  use mod_global_parameters
  use mod_constants
  use mod_variables
  use mod_fld

  integer, intent(in)             :: ixG^L, ix^L
  double precision, intent(in)    :: x(ixG^S, ndim)
  double precision, intent(inout) :: w(ixG^S, nw)
  double precision :: density(ixG^S), pressure(ixG^S), pert(ixG^S), amplitude
  double precision :: rad_Flux(ix^S, 1:ndim)
  double precision :: opacity(ix^S), Gamma_dep(ix^S)

  double precision :: a,b,c

  double precision :: opt_depth(ixGmin2:ixGmax2)
  double precision :: temp_init(ixG^S)
  double precision :: temperature(ixG^S)
  double precision :: k1(ixG^S), k2(ixG^S)
  double precision :: k3(ixG^S), k4(ixG^S)

  integer :: i

  amplitude = 0.0d0 !5.d-1  !1.d-5 !3.d-2

  pressure(:,ixGmin2) = p_bound0

  ! Set pressure profile using RK
  a = Flux0*kappa0/(4.d0/3.d0*fld_sigma_0*geff0)
  b = T_bound0**4.d0 - a*p_bound0
  c = -geff0*mp_cgs *0.6d0/kB_cgs * unit_pressure/(unit_temperature*unit_density)

  pressure(:,ixGmin2) = p_bound0

  do i=ixGmin2,ixGmax2-1
    k1(:,i) = (x(:,1+1,2) - x(:,1,2))&
    *c*pressure(:,i)/((a*pressure(:,i)+b)**(1.d0/4.d0))
    k2(:,i) = (x(:,1+1,2) - x(:,1,2))&
    *c*(pressure(:,i)+half*k1(:,i))/((a*(pressure(:,i)+half*k1(:,i))+b)**(1.d0/4.d0))
    k3(:,i) = (x(:,1+1,2) - x(:,1,2))&
    *c*(pressure(:,i)+half*k2(:,i))/((a*(pressure(:,i)+half*k2(:,i))+b)**(1.d0/4.d0))
    k4(:,i) = (x(:,1+1,2) - x(:,1,2))&
    *c*(pressure(:,i)+k3(:,i))/((a*(pressure(:,i)+k3(:,i))+b)**(1.d0/4.d0))

    pressure(:,i+1) = pressure(:,i) + one/6.d0 * (k1(:,i) + two*k2(:,i) + two*k3(:,i) + k4(:,i))
  enddo

  w(ixG^S, e_) = pressure(ixG^S)/(rhd_gamma - one)
  density(ixG^S) = pressure(ixG^S)/((a*pressure(ixG^S)+b)**(1.d0/4.d0))&
  *mp_cgs*0.6d0/kB_cgs * unit_pressure/(unit_temperature*unit_density)
  temp_init(ixG^S) = (a*pressure(ixG^S) + b)**(1.d0/4.d0)

  do i = ixGmin2,ixGmax2
    opt_depth(i) =  tau_bound - fld_kappa0*sum(density(5,:i))*(x(5,2,2)-x(5,1,2))
  enddo

  ! Set initial values for w
  w(ixG^S, rho_) = density(ixG^S)
  w(ixG^S, mom(:)) = zero
  w(ixG^S, e_) = pressure(ixG^S)/(rhd_gamma - one)
  w(ixG^S,r_e) = 4.d0*fld_sigma_0/(const_c/unit_velocity)*(a*pressure(ixG^S) + b)


  lower_bc_rho(:) = w(5,1:2, rho_)
  lower_bc_e(:) = w(5,1:2, e_)
  lower_bc_re(:) = w(5,1:2, r_e)


  !> Write initial profile to to file
  open(1,file='initial_cond')
    write(1,*) 'x    ', 'rho   ', 'e_    ', 'r_e   ', 'T_i   ', 'tau   '
    do i=ixmin2,ixmax2
      write(1,222) x(5,i,2), w(5,i,rho_), w(5,i,e_), w(5,i,r_e), temp_init(5,i), opt_depth(i)
    enddo
    222 format(6e15.5E3)
  close(1)

  !> perturb rho
  call RANDOM_NUMBER(pert)
  w(ixG^S, rho_) = w(ixG^S, rho_)*(one + amplitude*(pert(ixG^S)-half))

  call fld_get_opacity(w,x,ixG^L,ix^L)


  print*, 'end of initial conditions', w(5,5:10,r_e)

end subroutine initial_conditions

!==========================================================================================

! Extra routines can be placed here
! ...

!==========================================================================================

subroutine special_bound(qt,ixG^L,ixB^L,iB,w,x)

  use mod_global_parameters
  use mod_variables
  use mod_constants
  use mod_fld
  use mod_physics, only: phys_get_pthermal

  integer, intent(in) :: ixG^L, ixB^L, iB
  double precision, intent(in) :: qt, x(ixG^S,1:ndim)
  double precision, intent(inout) :: w(ixG^S,1:nw)
  double precision :: w_rad(ixG^S)
  double precision :: velocity(ixG^S,1:ndir), pressure(ixG^S)
  double precision :: fld_R(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
  double precision :: Gamma_dep(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
  double precision :: fld_lambda(ixGmin1+2:ixGmax1-2,ixGmin2+2:ixGmax2-2)
  double precision :: kb0, mp0

  double precision :: num_flux, dens, op, dy
  double precision :: Flux_m
  double precision :: lb, ub, cb
  double precision :: e_ip, e_i, e_im

  integer :: i,j

  select case (iB)

  case(3)
    w(:,2, rho_) = lower_bc_rho(2)
    w(:,2,mom(2)) = zero !w(:,3,mom(2))

    w(:,2, r_e) = w(:,3, r_e) - 3.0*w(:,2, rho_)*w(:,2,i_op)*Flux0/(const_c/unit_velocity)*(x(:,2,2)-x(:,3,2))

    !w(:,2, r_e) = lower_bc_re(2)

    pressure(:,2) = const_kB/(0.6*mp_cgs)*w(:,2, rho_)/unit_pressure*(unit_temperature*unit_density)&
    *((const_c/unit_velocity)/(4.d0*fld_sigma_0)*w(:,2, r_e))**0.25d0
    w(:,2, e_) = pressure(:,2)/(rhd_gamma - one) + (w(:,2, mom(2))*w(:,2, mom(2))/(2*w(:,2, rho_)))

    !w(:,2, e_) = w(:,3, e_) + w(:,3, mom(2))*(half - one/(rhd_gamma - one))&
    !*(one/w(:,2,rho_) -one/w(:,3, rho_))


    w(:,1, rho_) = lower_bc_rho(1)
    w(:,1,mom(2)) = zero !w(:,3,mom(2))

    w(:,1, r_e) = w(:,2, r_e) - 3.0*w(:,1, rho_)*w(:,2,i_op)*Flux0/(const_c/unit_velocity)*(x(:,1,2)-x(:,2,2))

    !w(:,1, r_e) = lower_bc_re(1)

    pressure(:,1) = const_kB/(0.6*mp_cgs)*w(:,1, rho_)/unit_pressure*(unit_temperature*unit_density)&
    *((const_c/unit_velocity)/(4.d0*fld_sigma_0)*w(:,1, r_e))**0.25d0
    w(:,1, e_) = pressure(:,1)/(rhd_gamma - one) + (w(:,1, mom(2))*w(:,1, mom(2))/(2*w(:,1, rho_)))

    ! w(:,1, e_) = w(:,2, e_) + w(:,2, mom(2))*(half - one/(rhd_gamma - one))&
    ! *(one/w(:,1,rho_) -one/w(:,2, rho_))

    call radiation_boundary(qt,ixG^L,iB,w,w_rad,x)
    do j = ixBmin2,ixBmax2
      w(ixGmin1:ixGmax1,j,r_e) = w_rad(ixGmin1:ixGmax1,j)
    enddo


  case(4)
    do i = ixGmin1,ixGmax1
      !> Linear interpolation E_r/rho
      w(i, ixBmin2, r_e) = w(i, ixBmin2 - 1, r_e) + (w(i, ixBmin2, rho_)/w(i, ixBmin2-1, rho_))&
      *(w(i, ixBmin2-1, r_e) - w(i, ixBmin2-2, r_e))

      w(i, ixBmax2, r_e) = w(i, ixBmax2 - 1, r_e) + (w(i, ixBmax2, rho_)/w(i, ixBmax2-1, rho_))&
      *(w(i, ixBmax2-1, r_e) - w(i, ixBmax2-2, r_e))
    enddo

    call radiation_boundary(qt,ixG^L,iB,w,w_rad,x)
    do j = ixBmin2,ixBmax2
      w(ixGmin1:ixGmax1,j,r_e) = w_rad(ixGmin1:ixGmax1,j)
    enddo


    !> Corners?
    w(ixGmin1,ixBmin2:ixBmax2,r_e) = w(ixGmax1-3,ixBmin2:ixBmax2,r_e)
    w(ixGmin1+1,ixBmin2:ixBmax2,r_e) = w(ixGmax1-2,ixBmin2:ixBmax2,r_e)
    w(ixGmax1,ixBmin2:ixBmax2,r_e) = w(ixGmin1+3,ixBmin2:ixBmax2,r_e)
    w(ixGmax1-1,ixBmin2:ixBmax2,r_e) = w(ixGmin1+2,ixBmin2:ixBmax2,r_e)

  case default
    call mpistop("BC not specified")
  end select
end subroutine special_bound

!==========================================================================================

subroutine radiation_boundary(qt,ixI^L,iB,w,w_rad,x)
  use mod_global_parameters
  integer, intent(in)             :: ixI^L, iB
  double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
  double precision, intent(in)    :: w(ixI^S,1:nw)
  double precision, intent(out)   :: w_rad(ixI^S)

  integer i,j

  select case (iB)

  case(3)
    do i=ixImin1,ixImax1
      do j=ixImin2+nghostcells-1,ixImin2,-1
        w_rad(i,j) = 3.0*w(i,j, rho_)*w(i,j+1,i_op)*Flux0/(const_c/unit_velocity)*(x(i,j,2)-x(i,j+1,2))
      enddo
    enddo

  case(4)
    do i=ixImin1,ixImax1
      do j = ixImax2-nghostcells+1,ixImax2
        w_rad(i, j) = w(i, j - 1, r_e) + (w(i, j, rho_)/w(i, j-1, rho_))&
        *(w(i, j-1, r_e) - w(i, j-2, r_e))
      enddo
    enddo

  case default
    call mpistop('boundary not known')
  end select
end subroutine radiation_boundary

!==========================================================================================

!> Calculate gravitational acceleration in each dimension
subroutine set_gravitation_field(ixI^L,ixO^L,wCT,x,gravity_field)
  use mod_global_parameters
  integer, intent(in)             :: ixI^L, ixO^L
  double precision, intent(in)    :: x(ixI^S,1:ndim)
  double precision, intent(in)    :: wCT(ixI^S,1:nw)
  double precision, intent(out)   :: gravity_field(ixI^S,ndim)

  gravity_field(ixI^S,1) = zero
  gravity_field(ixI^S,2) = -const_G*M_star/(R_star**2)/unit_velocity*unit_time

  !print*, gravity_field(5,5,2), Flux0*fld_kappa0/(const_c/unit_velocity)

end subroutine set_gravitation_field


!==========================================================================================

subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
! corresponding normalization values (default value 1)
  use mod_global_parameters
  use mod_physics

  integer, intent(in)                :: ixI^L,ixO^L
  double precision, intent(in)       :: x(ixI^S,1:ndim)
  double precision                   :: w(ixI^S,nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)
  double precision                   :: rad_flux(ixO^S,1:ndim), fld_lambda(ixO^S), fld_R(ixO^S)
  double precision                   :: g_rad(ixI^S,1:ndim), big_gamma(ixI^S), D(ixI^S,1:ndim), Temp(ixI^S)
  integer                            :: idim

  call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
  call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)
  call fld_get_diffcoef(w, x, ixI^L, ixO^L, D)
  call phys_get_pthermal(w,x,ixI^L,ixO^L,Temp)

  do idim = 1,ndim
    g_rad(ixO^S,idim) = w(ixO^S,i_op)*rad_flux(ixO^S,idim)/c_light0
  enddo
  big_gamma(ixO^S) = g_rad(ixO^S,2)/(const_G*M_star/(R_star**2)/unit_velocity*unit_time)

  w(ixO^S,nw+1)=rad_flux(ixO^S,1)*(unit_pressure*unit_velocity)
  w(ixO^S,nw+2)=rad_flux(ixO^S,2)*(unit_pressure*unit_velocity)
  w(ixO^S,nw+3)=fld_lambda(ixO^S)
  w(ixO^S,nw+4)=fld_R(ixO^S)
  w(ixO^S,nw+5)=g_rad(ixO^S,1)*unit_length/(unit_time**2)
  w(ixO^S,nw+6)=g_rad(ixO^S,2)*unit_length/(unit_time**2)
  w(ixO^S,nw+7)=big_gamma(ixO^S)
  w(ixO^S,nw+8)=D(ixO^S,1)
  w(ixO^S,nw+9)=D(ixO^S,2)
  w(ixO^S,nw+10)=Temp(ixO^S)/w(ixO^S,rho_)
end subroutine specialvar_output

subroutine specialvarnames_output(varnames)
! newly added variables need to be concatenated with the w_names/primnames string
  use mod_global_parameters
  character(len=*) :: varnames

  varnames = 'F1 F2 lamnda fld_R ar1 ar2 Gamma D1 D2 Temp'

end subroutine specialvarnames_output

!==========================================================================================

end module mod_usr

!==========================================================================================
