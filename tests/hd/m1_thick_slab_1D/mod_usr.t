!> This is a template for a new user problem of mhd
module mod_usr

  ! Include a physics module: mod_rho, mod_hd, mod_mhd ...
  use mod_hd
  use mod_m1

  implicit none

  ! Custom global variables can be defined here
  ! User input parameters
  double precision :: rho, kappa_step, Fbound, Ebound

  integer :: i_prad, i_pht, i_chi, i_ff, i_FLD, i_theta

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian")
    call usr_params_read(par_files)

    ! A routine for initial conditions is always required
    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see src/mod_usr_methods.t
    usr_special_bc => special_bound

    ! Set opacity
    usr_special_opacity_m1 => step_opacity

    !> Additional variables
    usr_process_grid => update_extravars

    ! Choose independent normalization units if using dimensionless variables.
    unit_length      = 1.d0    ! cm
    unit_velocity    = const_c
    unit_density     = 1.d0/m1_kappa_e ! g cm^-3

    ! Active the physics module: rho_activate(), hd_activate(), mhd_activate()
    call hd_activate()

    !> Extra variables
    i_Prad = var_set_extravar('P_rr','P_rr')
    i_ff = var_set_extravar('ff','ff')
    i_chi = var_set_extravar('chi','chi')
    i_fld = var_set_extravar('FLD','FLD')
    i_theta = var_set_extravar('theta','theta')

  end subroutine usr_init

  subroutine usr_params_read(files)

    character(len=*), intent(in) :: files(:)
    integer :: n

    namelist /star_list/ rho, Fbound, kappa_step

    do n = 1,size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, star_list, end=111)
       111 close(unitpar)
    enddo


  end subroutine usr_params_read


  subroutine initglobaldata_usr
    use mod_global_parameters

    print*, 'unit_density', unit_density
    print*, 'unit_velocity', unit_velocity
    print*, 'unit_length', unit_length
    print*, 'unit_time', unit_time

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixI^S, rho_) = rho
    w(ixI^S, mom(:)) = 0.d0
    Ebound = 10*max(kappa_step,1.d0)*(xprobmax1-xprobmin1)*m1_kappa_e*rho*Fbound/m1_c
    w(ixI^S, Er_) = 1.d-3 * Ebound !Ebound-(x(ixI^S,1)-xprobmin1)*3*w(ixI^S, rho_)*m1_kappa_e*Fbound/m1_c
    w(ixI^S, Fr(1)) = 0.d0 !Fbound

  end subroutine initial_conditions


  subroutine special_bound(qt,ixI^L,ixB^L,iB,w,x)

    ! Subroutine arguments
    integer, intent(in)    :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    ! Local variable
    integer :: ir, ix1, ixOmin1
    double precision :: prad(ixI^S)
    double precision :: gradE
    double precision :: ff(ixI^S), chi(ixI^S)

    select case (iB)
    case(1)
      w(ixB^S,kF_) = m1_kappa_e
      w(ixB^S,kE_) = m1_kappa_e
      w(ixB^S,kB_) = m1_kappa_e

      w(ixB^S, Er_) = Ebound !-(x(ixB^S,1)-xprobmin1)*3*w(ixB^S, rho_)*w(ixB^S,kF_)*Fbound/m1_c
      ! w(ixB^S, Fr(1)) = Fbound

      do ix1 = ixBmax1, ixBmin1, -1
        w(ix1, Fr(1)) = w(ix1+1, Fr(1))
      end do

    case(2)

      w(ixB^S,kF_) = kappa_step *m1_kappa_e
      w(ixB^S,kE_) = kappa_step *m1_kappa_e
      w(ixB^S,kB_) = kappa_step *m1_kappa_e

      do ix1 = ixBmin1, ixBmax1
        w(ix1, Er_) = 2*w(ix1-1, Er_) - w(ix1-2, Er_)
        w(ix1, Fr(1)) = w(ix1-1, Fr(1))
      end do
      ! w(ixB^S, Er_) = Ebound-(x(ixB^S,1)-xprobmin1)*3*w(ixB^S, rho_)*w(ixB^S,kF_)*Fbound/m1_c
      ! w(ixB^S, Fr(1)) = Fbound

    case default
      call mpistop("BC not specified")
    end select

  end subroutine special_bound

  subroutine step_opacity(ixI^L, ixO^L, w, x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout)    :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    where (x(ixI^S,1) .lt. (xprobmax1+xprobmin1)/2)
      w(ixI^S,kF_) = m1_kappa_e
      w(ixI^S,kE_) = m1_kappa_e
      w(ixI^S,kB_) = m1_kappa_e
    elsewhere
      w(ixI^S,kF_) = kappa_step * m1_kappa_e
      w(ixI^S,kE_) = kappa_step * m1_kappa_e
      w(ixI^S,kB_) = kappa_step * m1_kappa_e
    end where

  end subroutine step_opacity


  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    use mod_hd

    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: grav(ixI^S), grad(ixI^S), prad(ixI^S), v(ixI^S), gradv(ixI^S), pth(ixI^S)
    double precision :: gradE(ixI^S), theta(ixI^S)

    !> Photon tiring + Radflux CMF correction
    v(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call m1_get_prad(w, x, ixI^L, ixO^L, 1, 1, prad, str=.true.)
    w(ixI^S,i_Prad) = prad(ixI^S)


    w(ixI^S,i_ff) = sum(w(ixI^S, Fr(:))**2, dim=ndim+1)/(m1_c*w(ixI^S,Er_))**2
    w(ixI^S,i_chi) = (3.d0+4*w(ixI^S,i_ff))/(5.d0+2*dsqrt(4.d0-3*w(ixI^S,i_ff)))

    call gradientS(w(ixI^S,Er_), ixI^L, ixO^L, 1, gradE)
    w(ixI^S,i_FLD) = -m1_c/(3*w(ixI^S,rho_)*w(ixI^S,kF_)) * gradE(ixI^S)

    call m1_get_trapped_fraction(w, x, ixI^L, ixO^L, 1, theta)
    w(ixI^S,i_theta) = theta(ixI^S)

  end subroutine update_extravars

end module mod_usr
