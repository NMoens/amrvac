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
    subroutine initial_conditions(ixG^L, ix^L, w, x)
      use mod_global_parameters
      use mod_constants
      use mod_fld

      integer, intent(in)             :: ixG^L, ix^L
      double precision, intent(in)    :: x(ixG^S, ndim)
      double precision, intent(inout) :: w(ixG^S, nw)

      ! Set initial values for w
      w(ixG^S, rho_) = one
      w(ixG^S, mom(1)) = one
      w(ixG^S, mom(2)) = zero
      w(ixG^S, e_) = one
      w(ixG^S,r_e) = gradE*x(ixG^S,1)


      call fld_get_opacity(w, x, ixG^L, ix^L)
      call fld_get_fluxlimiter(w, x, ixG^L, ix^L)
      call fld_get_radflux(w, x, ixG^L, ix^L)

    end subroutine initial_conditions

  !==========================================================================================

    ! Extra routines can be placed here
    ! ...

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

    subroutine constant_var(level,qt,ixI^L,ixO^L,w,x)
      use mod_global_parameters
      use mod_fld
      integer, intent(in)             :: ixI^L,ixO^L,level
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixI^S,1:nw)
      double precision, intent(in)    :: x(ixI^S,1:ndim)

      w(ixI^S,rho_) = one
      w(ixI^S,e_) = one
      w(ixI^S,r_e) = gradE*x(ixI^S,1)

    end subroutine constant_var

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
    use mod_fld

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: theoretical(ixI^S)
    double precision                   :: residual(ixI^S)

    theoretical(ixI^S) = -gradE*1.d0/3.d0*global_time
    residual(ixI^S) = abs(theoretical(ixI^S) - w(ixI^S,mom(1)))/theoretical(ixI^S)

    if (it .eq. 0) open(1,file='f1_1')
    write(1,*) global_time, w(100,100,mom(1)), -gradE*1.d0/3.d0*global_time, w(100,100,i_flux(1))
    if (global_time .ge. time_max - dt) close(1)
    w(ixO^S,nw+1) = theoretical(ixO^S)
    w(ixO^S,nw+2) = residual(ixO^S)

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
