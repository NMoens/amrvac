!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd

  implicit none

  ! Custom variables can be defined here
  ! ...

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t

    ! Routine for setting special boundary conditions
    usr_special_bc => boundary_conditions

    ! Pseudo planar correction
    usr_source => PseudoPlanar

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm-3,cm-3

    ! Active the physics module
    call hd_activate()

  end subroutine usr_init

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_) = one
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(1)) = one
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(2)) = half

  end subroutine initial_conditions

  ! Extra routines can be placed here
  ! ...

  subroutine boundary_conditions(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2,iB,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
        ixBmin1,ixBmin2,ixBmax1,ixBmax2, iB
    double precision, intent(in)    :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

    select case (iB)
    case(1)
      w(ixBmin1:ixBmax1,:,rho_) = one
      w(ixBmin1:ixBmax1,:,mom(1)) = one
      w(ixBmin1:ixBmax1,:,mom(2)) = half

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

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

    rdir = 1
    pdir = 2

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_) - qdt*two*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(rdir))/x(ixImin1:ixImax1,ixImin2:ixImax2,rdir)

    !> dm_r/dt = m_phi**2/r
    !> dm_phi/dt = - m_phi m_r/r
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(rdir)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(rdir)) + qdt*(wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(pdir)))**two/(x(ixImin1:ixImax1,ixImin2:ixImax2,&
       rdir)*wCT(ixImin1:ixImax1,ixImin2:ixImax2,rho_))
    w(ixImin1:ixImax1,ixImin2:ixImax2,mom(pdir)) = w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(pdir)) - qdt*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(rdir))*wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(pdir))/(x(ixImin1:ixImax1,ixImin2:ixImax2,rdir)*wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,rho_))
  end subroutine PseudoPlanar

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

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'radius'
  end subroutine specialvarnames_output

end module mod_usr
