!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd
  use mod_m1

  implicit none

  double precision :: L, rho, Er, F_red, Bg_E
  integer :: i_Pxx,i_Pxy,i_Pyx,i_Pyy
  integer :: i_r, i_ff, i_chi, i_tau

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("spherical")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    ! usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions
    usr_special_bc => boundary_conditions

    usr_refine_grid => refine_beam
    usr_process_grid => update_extravars


    ! Choose independent normalization units if using dimensionless variables.
    ! unit_length      = 1.d0 ! cm
    ! unit_temperature = 1.d0 ! K
    ! unit_density     = 1.d0 ! g cm^-3

    unit_length      = L    ! cm
    unit_temperature = (Er/const_rad_a)**0.25d0  ! K
    unit_density     = rho ! g cm^-3

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

    i_Pxx = var_set_extravar('P_rr','P_rr')
    i_Pxy = var_set_extravar('P_rt','P_rt')
    i_Pyx = var_set_extravar('P_tr','P_tr')
    i_Pyy = var_set_extravar('P_tt','P_tt')
    i_r = var_set_extravar('r','r')
    i_ff = var_set_extravar('ff','ff')
    i_chi = var_set_extravar('chi','chi')
    i_tau = var_set_extravar('tau_c','tau_c')

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /stream_list/ L, rho, Er, F_red, Bg_E

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, stream_list, end=111)
       111 close(unitpar)
    enddo

  end subroutine usr_params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    select case (iprob)

    case(1,2)

      w(ixI^S,rho_) = rho/unit_density
      w(ixI^S,mom(:)) = 0.d0
      w(ixI^S,Er_) = Bg_E/unit_pressure
      w(ixI^S,Fr(:)) = 0.d0

    case(3)

      w(ixI^S,rho_) = rho/unit_density
      w(ixI^S,mom(:)) = 0.d0
      w(ixI^S,Er_) = Er/unit_pressure
      w(ixI^S,Fr(:)) = 0.d0
      w(ixI^S,Fr(1)) = Er/unit_pressure*m1_c*F_red/x(ixI^S,1)**2

    end select

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixI^L,ixB^L,iB,w,x)
  use mod_global_parameters
  use mod_hd

  integer, intent(in)             :: ixI^L, ixB^L, iB
  double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
  double precision, intent(inout) :: w(ixI^S,1:nw)

  integer :: ii

  select case (iprob)

  case(1) !> F Parallel to r

    select case (iB)

    case(1)
      w(ixB^S,Er_) = Bg_E/unit_pressure
      w(ixB^S,Fr(:)) = 0.d0
      where (x(ixB^S,2) .gt. 0.05d0*2*dpi .and. x(ixB^S,2) .lt. 0.2d0*2*dpi )
        w(ixB^S,Er_) = Er/unit_pressure
        w(ixB^S,Fr(1)) = Er/unit_pressure*m1_c*F_red
      end where

    case default
      call mpistop('boundary not known')
    end select

  case(2) !> F Parallel to z

    select case (iB)

    case(4)
      w(ixB^S,Er_) = Bg_E/unit_pressure
      w(ixB^S,Fr(:)) = 0.d0
      where (x(ixB^S,1) .gt.2 .and. x(ixB^S,1) .lt. 3 )
        w(ixB^S,Er_) = Er/unit_pressure
        w(ixB^S,Fr(2)) = -Er/unit_pressure*m1_c*F_red
      end where

    case default
      call mpistop('boundary not known')
    end select


  end select
end subroutine boundary_conditions

!> Enforce additional refinement or coarsening
!> One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
!> you must set consistent values for integers refine/coarsen:
!> refine = -1 enforce to not refine
!> refine =  0 doesn't enforce anything
!> refine =  1 enforce refinement
!> coarsen = -1 enforce to not coarsen
!> coarsen =  0 doesn't enforce anything
!> coarsen =  1 enforce coarsen
!> e.g. refine for negative first coordinate x < 0 as
!> if (any(x(ix^S,1) < zero)) refine=1
subroutine refine_beam(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
  use mod_global_parameters
  integer, intent(in)          :: igrid, level, ixI^L, ixO^L
  double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
  integer, intent(inout)       :: refine, coarsen

  double precision :: xx(ixI^S), yy(ixI^S)

  xx(ixI^S) = x(ixI^S,1)*sin(dpi-x(ixI^S,2))
  yy(ixI^S) = x(ixI^S,1)*cos(dpi-x(ixI^S,2))

  if (any(xx .gt. 2 .and. xx .lt. 3 .and. yy .lt. qt)) then
    refine=1
  else
    coarsen=1
  end if

end subroutine refine_beam

subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
  use mod_global_parameters
  use mod_hd

  integer, intent(in)             :: igrid,level,ixI^L,ixO^L
  double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
  double precision, intent(inout) :: w(ixI^S,1:nw)

  double precision :: prad(ixI^S), r0

  r0 = 3.d0
  w(ixI^S,i_r) = dsqrt(x(ixI^S,1)**2 + r0**2 - 2*x(ixI^S,1)*r0*cos(dpi/2-x(ixI^S,2)))

  w(ixI^S,i_ff) = sum(w(ixI^S, Fr(:))**2, dim=ndim+1)/(m1_c*w(ixI^S,Er_))**2
  w(ixI^S,i_ff) = min(w(ixI^S,i_ff),1.d0)
  w(ixI^S,i_ff) = max(w(ixI^S,i_ff),smalldouble)

  w(ixI^S,i_chi) = (3.d0+4*w(ixI^S,i_ff))/(5.d0+2*dsqrt(4.d0-3*w(ixI^S,i_ff)))

  call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, prad)
  w(ixI^S,i_Pxx) = prad(ixI^S)/w(ixI^S,Er_)
  call m1_get_prad(w, x, ixI^L, ixO^L, 1, 2, prad)
  w(ixI^S,i_Pxy) = prad(ixI^S)/w(ixI^S,Er_)
  call m1_get_prad(w, x, ixI^L, ixO^L, 2, 1, prad)
  w(ixI^S,i_Pyx) = prad(ixI^S)/w(ixI^S,Er_)
  call m1_get_prad(w, x, ixI^L, ixO^L, 2, 2, prad)
  w(ixI^S,i_Pyy) = prad(ixI^S)/w(ixI^S,Er_)

  w(ixI^S,i_tau) = block%dx(ixI^S,1)*w(ixI^S,kF_)*w(ixI^S,rho_)

end subroutine update_extravars

end module mod_usr
