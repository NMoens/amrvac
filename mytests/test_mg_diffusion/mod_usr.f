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
    
     call set_coordinate_system("Cartesian_2D")
    

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
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, mom(:)) = zero
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2, e_) = one
      w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,r_e) =  spotpattern(x,ixGmin1,ixGmin2,&
         ixGmax1,ixGmax2,0.d0)

      ! call fld_get_opacity(w, x, ixG^L, ix^L)
      ! call fld_get_fluxlimiter(w, x, ixG^L, ix^L)
      ! call fld_get_radflux(w, x, ixG^L, ix^L)
      ! call fld_get_eddington(w, x, ixG^L, ix^L)

      if (fld_diff_scheme .eq. 'mg') then
        call fld_get_diffcoef_central(w, w, x, ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
            ixmin1,ixmin2,ixmax1,ixmax2)
        ! call set_mg_bounds(w, x, ixG^L, ix^L)
      endif

      print*, mype, npe

    end subroutine initial_conditions

    function spotpattern(x,ixGmin1,ixGmin2,ixGmax1,ixGmax2,t1) result(e0)
      use mod_global_parameters

      integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2
      double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2, ndim),&
          t1
      double precision :: e0(ixGmin1:ixGmax1,ixGmin2:ixGmax2)
      integer ix1,ix2

      
       e0(ixGmin1:ixGmax1,ixGmin2:ixGmax2) =  2 + dexp(-8.d0 &
          *dpi**two*t1*unit_time)*sin(two*dpi*x(ixGmin1:ixGmax1,&
          ixGmin2:ixGmax2,1))*sin(two*dpi*x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
          2))

    end function spotpattern

  !==========================================================================================

  subroutine my_refinement(igrid,level,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,&
     ixmin2,ixmax1,ixmax2,qt,w,x,refine,coarsen)
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

    integer, intent(in) :: igrid, level, ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
        ixmin1,ixmin2,ixmax1,ixmax2
    double precision, intent(in) :: qt, w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:nw), x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
    integer, intent(inout) :: refine, coarsen

    !> Refine close to base
    coarsen = -1

    if (any(x(ixmin1:ixmax1,ixmin2:ixmax2,1) < 0.5d0)) refine=1
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

    subroutine constant_var(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2,level
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nw)
      double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:ndim)

      integer :: i
      w(ixImin1:ixImax1,ixImin2:ixImax2,rho_) = one
      w(ixImin1:ixImax1,ixImin2:ixImax2,mom(:)) = zero
      w(ixImin1:ixImax1,ixImin2:ixImax2,e_) = one

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

    theoretical(ixImin1:ixImax1,ixImin2:ixImax2) = spotpattern(x,ixImin1,&
       ixImin2,ixImax1,ixImax2,global_time)
    residual(ixImin1:ixImax1,ixImin2:ixImax2) = &
       abs(theoretical(ixImin1:ixImax1,ixImin2:ixImax2) - w(ixImin1:ixImax1,&
       ixImin2:ixImax2,r_e))/theoretical(ixImin1:ixImax1,ixImin2:ixImax2)

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
