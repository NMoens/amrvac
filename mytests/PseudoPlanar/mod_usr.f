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
    unit_numberdensity = 1.d9 ! cm-3,cm-3

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr
    use mod_global_parameters

  end subroutine initglobaldata_usr

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = one !/x(ixImin1:ixImax1,ixImin2:ixImax2,2)**2
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)) = zero
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(2)) = one !/x(ixImin1:ixImax1,ixImin2:ixImax2,2)**2
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = one !/x(ixImin1:ixImax1,ixImin2:ixImax2,2)**2
    w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = one !/x(ixImin1:ixImax1,ixImin2:ixImax2,2)**2

    call get_rad_extravars(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2)

  end subroutine initial_conditions

  subroutine boundary_conditions(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
        ixBmin1,ixBmin2,ixBmax1,ixBmax2, iB
    double precision, intent(in)    :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

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

  subroutine mg_boundary_conditions(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iB,w,x)

    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: grad4

    grad4 = sum((w(ixOmin1:ixOmax1,ixOmax2,r_e) - w(ixOmin1:ixOmax1,ixOmax2+1,&
        r_e)) / (x(ixOmin1:ixOmax1,ixOmax2,2) - x(ixOmin1:ixOmax1,ixOmax2+1,&
       2)))/(ixOmax1-ixOmin1)

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

        if (sum(w(ixOmin1:ixOmax1,ixOmax2,r_e)) .le. sum(w(ixOmin1:ixOmax1,&
           ixOmax2+1, r_e))) then
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

  subroutine MyOpacity(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,w,x,kappa)
    use mod_global_parameters
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.34d0/unit_opacity

  end subroutine MyOpacity


  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine PseudoPlanar(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, iwmin,iwmax
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2)
    integer :: rdir, pdir

    rdir = 2
    pdir = 1

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_) - qdt*two*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(rdir))/x(ixImin1:ixImax1,ixImin2:ixImax2,rdir)

    !> dm_r/dt = m_phi**2/r
    !> dm_phi/dt = - m_phi m_r/r
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(rdir)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(rdir)) + qdt*(wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(pdir)))**two/x(ixImin1:ixImax1,ixImin2:ixImax2,rdir)
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(pdir)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(pdir)) - qdt*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(rdir))*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(pdir))/x(ixImin1:ixImax1,ixImin2:ixImax2,rdir)

    !> de/dt = -2 (e+p)v_r/r
    call phys_get_pthermal(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,pth)
    w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
       e_) + qdt*two*(wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       e_)+pth(ixImin1:ixImax1,ixImin2:ixImax2))*wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(rdir))/(wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       rho_)*x(ixImin1:ixImax1,ixImin2:ixImax2,rdir))

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      call get_rad_extravars(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2)
      w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = w(ixImin1:ixImax1,&
         ixImin2:ixImax2,r_e) - qdt*two*(wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         e_)*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         mom(rdir))/wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         rho_) + w(ixImin1:ixImax1,ixImin2:ixImax2,&
         i_flux(rdir)))/x(ixImin1:ixImax1,ixImin2:ixImax2,rdir)
    else
      w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = w(ixImin1:ixImax1,&
         ixImin2:ixImax2,r_e) - qdt*two*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         e_)*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         mom(rdir))/(wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
         rho_)*x(ixImin1:ixImax1,ixImin2:ixImax2,rdir))
    endif
  end subroutine PseudoPlanar

end module mod_usr
