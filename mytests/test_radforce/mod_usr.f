!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd

  implicit none

  ! Custom variables can be defined here
  ! ...
  double precision, parameter :: gradE = -1.d0

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    use mod_constants

    call set_coordinate_system("Cartesian_2D")

    unit_velocity = one
    unit_numberdensity = one/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_length = one

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Set boundary conditions
    usr_special_bc => boundary_conditions

    ! Keep the radiative energy constant with internal bound
    usr_internal_bc => constant_var

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init

  !==========================================================================================

  subroutine initglobaldata_usr
    use mod_global_parameters


  end subroutine initglobaldata_usr

  !==========================================================================================

    !> A routine for specifying initial conditions
    subroutine initial_conditions(ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
       ixmin2,ixmax1,ixmax2, w, x)
      use mod_global_parameters
      use mod_constants
      use mod_fld

      integer, intent(in)             :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
          ixmin1,ixmin2,ixmax1,ixmax2
      double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          ndim)
      double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          nw)

      ! Set initial values for w
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, rho_) = one
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(1)) = zero
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(2)) = zero
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, e_) = one
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,r_e) = gradE*x(ixGmin1:ixGmax1,&
         ixGmin2:ixGmax2,1)


      call fld_get_opacity(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
         ixmin2,ixmax1,ixmax2)
      call fld_get_fluxlimiter(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
         ixmin2,ixmax1,ixmax2)
      call fld_get_radflux(w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,&
         ixmin2,ixmax1,ixmax2)

    end subroutine initial_conditions

  !==========================================================================================

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
      ! Set initial values for w
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2, rho_) = one
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2, mom(1)) = zero
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2, mom(2)) = zero
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2, e_) = one
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,r_e) = gradE*x(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,1)
    case(2)
      ! Set initial values for w
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2, rho_) = one
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2, mom(1)) = zero
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2, mom(2)) = zero
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2, e_) = one
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,r_e) = gradE*x(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,1)


    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

  !==========================================================================================

  !> internal boundary, user defined
    !
    !> This subroutine can be used to artificially overwrite ALL conservative
    !> variables in a user-selected region of the mesh, and thereby act as
    !> an internal boundary region. It is called just before external (ghost cell)
    !> boundary regions will be set by the BC selection. Here, you could e.g.
    !> want to introduce an extra variable (nwextra, to be distinguished from nwaux)
    !> which can be used to identify the internal boundary region location.
    !> Its effect should always be local as it acts on the mesh.

    subroutine constant_var(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,w,x)
      use mod_global_parameters
      use mod_fld
      integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2,level
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nw)
      double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:ndim)

      w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = one
      w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = one
      w(ixImin1:ixImax1,ixImin2:ixImax2,r_e) = gradE*x(ixImin1:ixImax1,&
         ixImin2:ixImax2,1)

    end subroutine constant_var

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
    use mod_fld

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: theoretical(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                   :: residual(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    theoretical(ixImin1:ixImax1,ixImin2:ixImax2) = &
       -gradE*1.d0/3.d0*global_time
    residual(ixImin1:ixImax1,ixImin2:ixImax2) = &
       abs(theoretical(ixImin1:ixImax1,ixImin2:ixImax2) - w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mom(1)))/theoretical(ixImin1:ixImax1,ixImin2:ixImax2)

    if (it .eq. 0) open(1,file='f1_1')
    write(1,*) global_time, w(100,100,mom(1)), -gradE*1.d0/3.d0*global_time,&
        w(100,100,i_flux(1))
    if (global_time .ge. time_max - dt) close(1)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1) = theoretical(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2) = residual(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
    use mod_global_parameters
    character(len=*) :: varnames

    varnames = 'theoretical residual'

  end subroutine specialvarnames_output

  !==========================================================================================

  end module mod_usr

  !==========================================================================================
