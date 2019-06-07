!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  ! Custom variables can be defined here
  ! ...

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Pseudo planar correction
    usr_source => PseudoPlanar

    ! Routine for setting special boundary conditions
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions

    ! User defined opacities
    usr_special_opacity => MyOpacity

    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    w(ixI^S,rho_) = one !/x(ixI^S,2)**2
    w(ixI^S,mom(1)) = zero
    w(ixI^S,mom(2)) = one !/x(ixI^S,2)**2
    w(ixI^S,e_) = one !/x(ixI^S,2)**2
    w(ixI^S,r_e) = one !/x(ixI^S,2)**2

    call get_rad_extravars(w, x, ixI^L, ixO^L)

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixG^L,ixB^L,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixG^L, ixB^L, iB
    double precision, intent(in)    :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    select case (iB)
    case(3)
      w(:,ixBmin2:ixBmax2,rho_) = one
      w(:,ixBmin2:ixBmax2,mom(1)) = zero
      w(:,ixBmin2:ixBmax2,mom(2)) = one
      w(:,ixBmin2:ixBmax2,e_) = one
      w(:,ixBmin2:ixBmax2,r_e) = one

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

  subroutine mg_boundary_conditions(qt,ixI^L,ixO^L,iB,w,x)

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iB
    double precision, intent(in)    :: qt, x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)

    double precision :: grad4

    grad4 = sum((w(ixOmin1:ixOmax1,ixOmax2,r_e) - w(ixOmin1:ixOmax1,ixOmax2+1, r_e)) &
    / (x(ixOmin1:ixOmax1,ixOmax2,2) - x(ixOmin1:ixOmax1,ixOmax2+1,2)))/(ixOmax1-ixOmin1)

    select case (iB)
      case (1)
         mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
      case (2)
         mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
      case (3)
        mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        mg%bc(iB, mg_iphi)%bc_value = 7.3583819042386223 !7.6447315544263788
      case (4)
        ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
        ! mg%bc(iB, mg_iphi)%bc_value = 0.77780683570039721
        ! mg%bc(iB, mg_iphi)%bc_value = grad4 !min(grad4,zero)
        ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous

        ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann

        if (sum(w(ixOmin1:ixOmax1,ixOmax2,r_e)) .le. &
           sum(w(ixOmin1:ixOmax1,ixOmax2+1, r_e))) then
           mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
           mg%bc(iB, mg_iphi)%bc_value = zero
        else
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
        endif

      case default
        print *, "Not a standard: ", trim(typeboundary(iw_r_e, iB))
        error stop "You have to set a user-defined boundary method"
    end select
  end subroutine mg_boundary_conditions

  subroutine MyOpacity(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    kappa(ixO^S) = 0.34d0/unit_opacity

  end subroutine MyOpacity


  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine PseudoPlanar(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: pth(ixI^S)
    integer :: rdir, pdir

    rdir = 2
    pdir = 1

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    w(ixI^S,rho_) = w(ixI^S,rho_) - qdt*two*wCT(ixI^S,mom(rdir))/x(ixI^S,rdir)

    !> dm_r/dt = m_phi**2/r
    !> dm_phi/dt = - m_phi m_r/r
    w(ixI^S,mom(rdir)) = w(ixI^S,mom(rdir)) + qdt*(wCT(ixI^S,mom(pdir)))**two/x(ixI^S,rdir)
    w(ixI^S,mom(pdir)) = w(ixI^S,mom(pdir)) - qdt*wCT(ixI^S,mom(rdir))*wCT(ixI^S,mom(pdir))/x(ixI^S,rdir)

    !> de/dt = -2 (e+p)v_r/r
    call phys_get_pthermal(wCT,x,ixI^L,ixO^L,pth)
    w(ixI^S,e_) = w(ixI^S,e_) + qdt*two*(wCT(ixI^S,e_)+pth(ixI^S))*wCT(ixI^S,mom(rdir))/(wCT(ixI^S,rho_)*x(ixI^S,rdir))

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      call get_rad_extravars(w, x, ixI^L, ixO^L)
      w(ixI^S,r_e) = w(ixI^S,r_e) - qdt*two*(wCT(ixI^S,e_)*wCT(ixI^S,mom(rdir))/wCT(ixI^S,rho_) + w(ixI^S,i_flux(rdir)))/x(ixI^S,rdir)
    else
      w(ixI^S,r_e) = w(ixI^S,r_e) - qdt*two*wCT(ixI^S,e_)*wCT(ixI^S,mom(rdir))/(wCT(ixI^S,rho_)*x(ixI^S,rdir))
    endif
  end subroutine PseudoPlanar

end module mod_usr
