!> Steady-state test of the M1 geometric and velocity source terms in spherical
!> coordinates (comoving-frame M1, zero opacity, frozen hydro).
!>
!> E, F and v are uniform *Cartesian* fields (or a rigid rotation), written in
!> spherical components. All Cartesian gradients vanish (or the rotation is
!> divergence- and shear-free), so the exact solution is that nothing changes:
!>   dE/dt = -div(vE + F) - P:grad v               = 0
!>   dF/dt = -div(vF + c^2 P) - (F.grad) v         = 0
!> In spherical components every term is non-zero, so every geometric correction
!> has to cancel the flux divergences. Any drift is a geometric-term error.
!>
!> iprob 1: v = V x_hat, F = f c E z_hat   (v perpendicular to F; 3D only)
!> iprob 2: v = V z_hat, F = f c E z_hat   (v parallel to F)
!> iprob 3: v = 0,       F = f c E z_hat   (static geometric terms only)
!> iprob 4: v_phi = Omega r sin(theta) (rigid rotation about z), F = f c E z_hat
!> iprob 5: v = 0,       F = f c E x_hat   (static, P_phi_theta /= 0; 3D only)
!> iprob 6: v = V z_hat, F = f c E x_hat   (v perpendicular to F, F_phi /= 0; 3D only)
module mod_usr

  use mod_hd
  use mod_m1

  implicit none

  double precision :: beta_max = 0.1d0, F_red = 0.5d0, E0 = 1.d0
  integer :: i_eE, i_eF, i_eFx, i_eFy, i_eFz, i_ev

contains

  subroutine usr_init()

    {^IFTWOD call set_coordinate_system("spherical_2.5D")}
    {^IFTHREED call set_coordinate_system("spherical")}
    call usr_params_read(par_files)

    usr_init_one_grid => initial_conditions
    usr_special_bc    => boundary_conditions
    usr_internal_bc   => Fix_hd
    usr_process_grid  => update_extravars

    !> c = 1 in code units
    unit_velocity = const_c
    unit_time     = 1.d0/unit_velocity
    unit_density  = 1.d0

    call hd_activate()

    !> Errors relative to the exact (initial) state
    i_eE  = var_set_extravar('err_E','err_E')     !< |E/E0 - 1|
    i_eF  = var_set_extravar('err_F','err_F')     !< |F - F_ex| / |F0|
    i_eFx = var_set_extravar('err_Fx','err_Fx')   !< |Fx - Fx_ex| / |F0|
    i_eFy = var_set_extravar('err_Fy','err_Fy')
    i_eFz = var_set_extravar('err_Fz','err_Fz')
    i_ev  = var_set_extravar('err_v','err_v')     !< frozen-hydro sanity check

  end subroutine usr_init

  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /cmf_list/ beta_max, F_red, E0

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, cmf_list, end=111)
       111 close(unitpar)
    enddo
  end subroutine usr_params_read

  !> Spherical components (ar, at, ap) of the uniform Cartesian vector (ax, ay, az)
  subroutine cart2sph(ixI^L, ixO^L, x, ax, ay, az, ar, at, ap)
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim), ax, ay, az
    double precision, intent(out) :: ar(ixI^S), at(ixI^S), ap(ixI^S)
    double precision              :: th(ixI^S), ph(ixI^S)

    th(ixO^S) = x(ixO^S,2)
    {^IFTWOD ph(ixO^S) = 0.d0}
    {^IFTHREED ph(ixO^S) = x(ixO^S,3)}

    ar(ixO^S) = ax*dsin(th(ixO^S))*dcos(ph(ixO^S)) + ay*dsin(th(ixO^S))*dsin(ph(ixO^S)) + az*dcos(th(ixO^S))
    at(ixO^S) = ax*dcos(th(ixO^S))*dcos(ph(ixO^S)) + ay*dcos(th(ixO^S))*dsin(ph(ixO^S)) - az*dsin(th(ixO^S))
    ap(ixO^S) = -ax*dsin(ph(ixO^S)) + ay*dcos(ph(ixO^S))
  end subroutine cart2sph

  !> Exact velocity and comoving flux, spherical components
  subroutine exact_fields(ixI^L, ixO^L, x, v, F)
    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S,1:3), F(ixI^S,1:3)
    double precision              :: V0, Omega

    V0 = beta_max*m1_c
    if (iprob == 5 .or. iprob == 6) then
      {^IFTWOD call mpistop('iprob 5/6 (F along x) are not axisymmetric: 3D only')}
      call cart2sph(ixI^L, ixO^L, x, F_red*m1_c*E0, 0.d0, 0.d0, F(ixI^S,1), F(ixI^S,2), F(ixI^S,3))
    else
      call cart2sph(ixI^L, ixO^L, x, 0.d0, 0.d0, F_red*m1_c*E0, F(ixI^S,1), F(ixI^S,2), F(ixI^S,3))
    end if

    select case (iprob)
    case (1)
      {^IFTWOD call mpistop('iprob 1 (v along x) is not axisymmetric: 3D only')}
      call cart2sph(ixI^L, ixO^L, x, V0, 0.d0, 0.d0, v(ixI^S,1), v(ixI^S,2), v(ixI^S,3))
    case (2, 6)
      call cart2sph(ixI^L, ixO^L, x, 0.d0, 0.d0, V0, v(ixI^S,1), v(ixI^S,2), v(ixI^S,3))
    case (3, 5)
      v(ixO^S,:) = 0.d0
    case (4)
      !> |v| <= beta_max c over the domain
      Omega = V0/xprobmax1
      v(ixO^S,1) = 0.d0
      v(ixO^S,2) = 0.d0
      v(ixO^S,3) = Omega*x(ixO^S,1)*dsin(x(ixO^S,2))
    case default
      call mpistop('iprob not known')
    end select
  end subroutine exact_fields

  subroutine set_exact(ixI^L, ixO^L, w, x, hydro_only)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    logical, intent(in)             :: hydro_only
    double precision                :: v(ixI^S,1:3), F(ixI^S,1:3)
    integer                         :: idir

    call exact_fields(ixI^L, ixO^L, x, v, F)
    w(ixO^S,rho_) = 1.d0
    do idir = 1, ndir
      w(ixO^S,mom(idir)) = w(ixO^S,rho_)*v(ixO^S,idir)
    end do
    if (hydro_only) return
    w(ixO^S,Er_) = E0
    do idir = 1, ndir
      w(ixO^S,Fr(idir)) = F(ixO^S,idir)
    end do
  end subroutine set_exact

  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    call set_exact(ixI^L, ixI^L, w, x, .false.)
  end subroutine initial_conditions

  !> All boundaries: exact state
  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    call set_exact(ixI^L, ixB^L, w, x, .false.)
  end subroutine boundary_conditions

  !> Frozen hydro
  subroutine Fix_hd(level,qt,ixI^L,ixO^L,w,x)
    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    call set_exact(ixI^L, ixI^L, w, x, .true.)
  end subroutine Fix_hd

  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: v(ixI^S,1:3), F(ixI^S,1:3), dF(ixI^S,1:3)
    double precision :: th(ixI^S), ph(ixI^S), dx_(ixI^S), dy_(ixI^S), dz_(ixI^S), F0
    integer          :: idir

    F0 = F_red*m1_c*E0
    call exact_fields(ixI^L, ixO^L, x, v, F)
    do idir = 1, 3
      dF(ixO^S,idir) = w(ixO^S,Fr(idir)) - F(ixO^S,idir)
    end do

    !> Back to Cartesian components
    th(ixO^S) = x(ixO^S,2)
    {^IFTWOD ph(ixO^S) = 0.d0}
    {^IFTHREED ph(ixO^S) = x(ixO^S,3)}
    dx_(ixO^S) = dF(ixO^S,1)*dsin(th(ixO^S))*dcos(ph(ixO^S)) + dF(ixO^S,2)*dcos(th(ixO^S))*dcos(ph(ixO^S)) &
               - dF(ixO^S,3)*dsin(ph(ixO^S))
    dy_(ixO^S) = dF(ixO^S,1)*dsin(th(ixO^S))*dsin(ph(ixO^S)) + dF(ixO^S,2)*dcos(th(ixO^S))*dsin(ph(ixO^S)) &
               + dF(ixO^S,3)*dcos(ph(ixO^S))
    dz_(ixO^S) = dF(ixO^S,1)*dcos(th(ixO^S)) - dF(ixO^S,2)*dsin(th(ixO^S))

    w(ixO^S,i_eE)  = dabs(w(ixO^S,Er_)/E0 - 1.d0)
    w(ixO^S,i_eFx) = dabs(dx_(ixO^S))/F0
    w(ixO^S,i_eFy) = dabs(dy_(ixO^S))/F0
    w(ixO^S,i_eFz) = dabs(dz_(ixO^S))/F0
    w(ixO^S,i_eF)  = dsqrt(dx_(ixO^S)**2 + dy_(ixO^S)**2 + dz_(ixO^S)**2)/F0
    w(ixO^S,i_ev)  = 0.d0
    do idir = 1, ndir
      w(ixO^S,i_ev) = w(ixO^S,i_ev) + dabs(w(ixO^S,mom(idir))/w(ixO^S,rho_) - v(ixO^S,idir))
    end do

  end subroutine update_extravars

end module mod_usr
