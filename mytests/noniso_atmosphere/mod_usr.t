!> This is a template for a new user problem
module mod_usr

! Include a physics module
use mod_rhd
use mod_fld
use mod_global_parameters

implicit none

  integer, parameter :: nyc = 200

  integer :: i_is(1:nyc)
  double precision :: y_is(1:nyc)
  double precision :: rho_is(1:nyc)
  double precision :: tg_is(1:nyc)
  double precision :: pg_is(1:nyc)
  double precision :: er_is(1:nyc)
  double precision :: tau_is(1:nyc)

  double precision :: gamma0
  double precision :: mstar
  double precision :: rstar
  double precision :: yhe
  integer :: ny_vac
  double precision :: dy_vac
  double precision :: tau_max
  character :: kap_law

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

  ! Graviatational field
  usr_gravity => set_gravitation_field

  ! Output routines
  usr_aux_output    => specialvar_output
  usr_add_aux_names => specialvarnames_output

  ! Active the physics module
  call rhd_activate()

end subroutine usr_init

!==========================================================================================

subroutine initglobaldata_usr
  use mod_global_parameters

  integer :: i

  OPEN(1,FILE='init_struc_amrvac')
  do i = 1,domain_nx2
    READ(1,'(1i4,1e20.10,5e18.8)') i_is(i), y_is(i), rho_is(i), tg_is(i), pg_is(i), er_is(i), tau_is(i)
  enddo
  CLOSE(1)

  OPEN(2,FILE='init_indat')
   READ(2,*) gamma0,mstar,rstar,yhe,ny_vac,dy_vac,tau_max,kap_law
  CLOSE(2)

  !> Perform checks to see if init_indat and usr.par agree
  !if (ny_vac .ne. domain_nx2) call mpistop('Number of cells dont agree')

end subroutine initglobaldata_usr

!==========================================================================================

!> A routine for specifying initial conditions
subroutine initial_conditions(ixG^L, ix^L, w, x)
  use mod_global_parameters
  use mod_constants
  use mod_variables

  integer, intent(in)             :: ixG^L, ix^L
  double precision, intent(in)    :: x(ixG^S, ndim)
  double precision, intent(inout) :: w(ixG^S, nw)
  integer :: i

  do i = ixGmin2,ixGmax2
    w(ixGmin1: ixGmax1, i, rho_) = rho_is(i)
    w(ixGmin1: ixGmax1, i, mom(:)) = zero
    w(ixGmin1: ixGmax1, i, e_) = pg_is(i)/(rhd_gamma-1)
    w(ixGmin1: ixGmax1, i, r_e) = er_is(i)
  enddo

end subroutine initial_conditions

!==========================================================================================

! Extra routines can be placed here
! ...

!==========================================================================================

!> Calculate gravitational acceleration in each dimension
subroutine set_gravitation_field(ixI^L,ixO^L,wCT,x,gravity_field)
  use mod_global_parameters
  integer, intent(in)             :: ixI^L, ixO^L
  double precision, intent(in)    :: x(ixI^S,1:ndim)
  double precision, intent(in)    :: wCT(ixI^S,1:nw)
  double precision, intent(out)   :: gravity_field(ixI^S,ndim)

  gravity_field(ixI^S,1) = zero
  gravity_field(ixI^S,2) = -6.67e-8*mstar/rstar**2*(unit_time**2/unit_length)
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
  double precision                   :: rad_flux(ixO^S,1:ndim), rad_pressure(ixO^S), fld_lambda(ixO^S), fld_R(ixO^S), fld_kappa(ixO^S)
  double precision                   :: g_rad(ixI^S,1:ndim), big_gamma(ixI^S), D(ixI^S,1:ndim)
  integer                            :: idim

  call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux)
  call fld_get_radpress(w, x, ixI^L, ixO^L, rad_pressure)
  call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)
  call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)
  call fld_get_diffcoef(w, x, ixI^L, ixO^L, D)

  do idim = 1,ndim
    g_rad(ixO^S,idim) = fld_kappa(ixO^S)*rad_flux(ixO^S,idim)/const_c
  enddo
  big_gamma(ixO^S) = g_rad(ixO^S,2)/(6.67e-8*mstar/rstar**2*(unit_time**2/unit_length))

  w(ixO^S,nw+1)=rad_flux(ixO^S,1)*(unit_pressure*unit_velocity)
  w(ixO^S,nw+2)=rad_flux(ixO^S,2)*(unit_pressure*unit_velocity)
  w(ixO^S,nw+3)=rad_pressure(ixO^S)*unit_pressure
  w(ixO^S,nw+4)=fld_lambda(ixO^S)
  w(ixO^S,nw+5)=fld_R(ixO^S)
  w(ixO^S,nw+6)=g_rad(ixO^S,1)*unit_length/(unit_time**2)
  w(ixO^S,nw+7)=g_rad(ixO^S,2)*unit_length/(unit_time**2)
  w(ixO^S,nw+8)=big_gamma(ixO^S)
  w(ixO^S,nw+9)=D(ixO^S,1)
  w(ixO^S,nw+10)=D(ixO^S,2)
  w(ixO^S,nw+11)= fld_kappa(ixO^S)
end subroutine specialvar_output

subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
  use mod_global_parameters
  character(len=*) :: varnames

  varnames = 'F1 F2 P_rad lamnda fld_R ar1 ar2 Gamma D1 D2 kappa'

end subroutine specialvarnames_output

!==========================================================================================

end module mod_usr

!==========================================================================================
