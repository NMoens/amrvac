!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd

  implicit none

  ! Custom variables can be defined here
  ! ...

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    use mod_constants

    ! Choose coordinate system as n-D Cartesian
    {^IFONED call set_coordinate_system("Cartesian_1D")}
    {^IFTWOD call set_coordinate_system("Cartesian_2D")}
    {^IFTHREED call set_coordinate_system("Cartesian_3D")}

    unit_velocity = 1.d0
    unit_numberdensity = 1.d0/((1.d0+4.d0*He_abundance)*mp_cgs)
    unit_length = 1.d0

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Keep the radiative energy constant with internal bound
    usr_internal_bc => constant_var

    ! Output routines
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    ! Refine mesh for testing purposes
    usr_refine_grid => my_refinement

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
      w(ixG^S, mom(:)) = zero
      w(ixG^S, e_) = one
      w(ixG^S,r_e) =  spotpattern(x,ixG^L,0.d0)

      ! call fld_get_opacity(w, x, ixG^L, ix^L)
      ! call fld_get_fluxlimiter(w, x, ixG^L, ix^L)
      ! call fld_get_radflux(w, x, ixG^L, ix^L)
      ! call fld_get_eddington(w, x, ixG^L, ix^L)

      if (fld_diff_scheme .eq. 'mg') then
        call fld_get_diffcoef_central(w, w, x, ixG^L, ix^L)
        ! call set_mg_bounds(w, x, ixG^L, ix^L)
      endif

      print*, mype, npe

    end subroutine initial_conditions

    function spotpattern(x,ixG^L,t1) result(e0)
      use mod_global_parameters

      integer, intent(in) :: ixG^L
      double precision, intent(in) :: x(ixG^S, ndim), t1
      double precision :: e0(ixG^S)
      integer ix^D

      {^IFONED e0(ixG^S) =  2 + dexp(-8.d0 *dpi**two*t1*unit_time)*sin(two*dpi*x(ixG^S,1))}
      {^IFTWOD e0(ixG^S) =  2 + dexp(-8.d0 *dpi**two*t1*unit_time)*sin(two*dpi*x(ixG^S,1))*sin(two*dpi*x(ixG^S,2))}

    end function spotpattern

  !==========================================================================================

  subroutine my_refinement(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    !> Refine close to base
    coarsen = -1

    if (any(x(ix^S,1) < 0.5d0)) refine=1
    ! if (any(x(ix^S,1) < 1.0d0)) refine=1
    ! if (any(x(ix^S,1) < 1.5d0)) refine=1

  end subroutine my_refinement

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
      integer, intent(in)             :: ixI^L,ixO^L,level
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixI^S,1:nw)
      double precision, intent(in)    :: x(ixI^S,1:ndim)

      integer :: i
      w(ixI^S,rho_) = one
      w(ixI^S,mom(:)) = zero
      w(ixI^S,e_) = one

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

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: theoretical(ixI^S)
    double precision                   :: residual(ixI^S)

    theoretical(ixI^S) = spotpattern(x,ixI^L,global_time)
    residual(ixI^S) = abs(theoretical(ixI^S) - w(ixI^S,r_e))/theoretical(ixI^S)

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
