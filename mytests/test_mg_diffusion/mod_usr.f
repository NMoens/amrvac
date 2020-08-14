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
     call set_coordinate_system("Cartesian_1D")
    
    

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

    ! Active the physics module
    call rhd_activate()

  end subroutine usr_init

  !==========================================================================================

  subroutine initglobaldata_usr
    use mod_global_parameters


  end subroutine initglobaldata_usr

  !==========================================================================================

    !> A routine for specifying initial conditions
    subroutine initial_conditions(ixGmin1,ixGmax1, ixmin1,ixmax1, w, x)
      use mod_global_parameters
      use mod_constants
      use mod_fld

      integer, intent(in)             :: ixGmin1,ixGmax1, ixmin1,ixmax1
      double precision, intent(in)    :: x(ixGmin1:ixGmax1, ndim)
      double precision, intent(inout) :: w(ixGmin1:ixGmax1, nw)

      ! Set initial values for w
      w(ixGmin1:ixGmax1, rho_) = one
      w(ixGmin1:ixGmax1, mom(:)) = zero
      w(ixGmin1:ixGmax1, e_) = one
      w(ixGmin1:ixGmax1,r_e) =  spotpattern(x,ixGmin1,ixGmax1,0.d0)

      ! call fld_get_opacity(w, x, ixG^L, ix^L)
      ! call fld_get_fluxlimiter(w, x, ixG^L, ix^L)
      ! call fld_get_radflux(w, x, ixG^L, ix^L)
      ! call fld_get_eddington(w, x, ixG^L, ix^L)

      if (fld_diff_scheme .eq. 'mg') then
        call fld_get_diffcoef_central(w, w, x, ixGmin1,ixGmax1, ixmin1,ixmax1)
        ! call set_mg_bounds(w, x, ixG^L, ix^L)
      endif

      print*, mype, npe

    end subroutine initial_conditions

    function spotpattern(x,ixGmin1,ixGmax1,t1) result(e0)
      use mod_global_parameters

      integer, intent(in) :: ixGmin1,ixGmax1
      double precision, intent(in) :: x(ixGmin1:ixGmax1, ndim), t1
      double precision :: e0(ixGmin1:ixGmax1)
      integer ix1

       e0(ixGmin1:ixGmax1) =  2 + dexp(-8.d0 &
          *dpi**two*t1*unit_time)*sin(two*dpi*x(ixGmin1:ixGmax1,1))
      

    end function spotpattern

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

    subroutine constant_var(level,qt,ixImin1,ixImax1,ixOmin1,ixOmax1,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixImin1,ixImax1,ixOmin1,ixOmax1,level
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
      double precision, intent(in)    :: x(ixImin1:ixImax1,1:ndim)

      integer :: i
      w(ixImin1:ixImax1,rho_) = one
      w(ixImin1:ixImax1,mom(:)) = zero
      w(ixImin1:ixImax1,e_) = one

    end subroutine constant_var

  !==========================================================================================

  subroutine specialvar_output(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  !
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_global_parameters
    use mod_physics

    integer, intent(in)                :: ixImin1,ixImax1,ixOmin1,ixOmax1
    double precision, intent(in)       :: x(ixImin1:ixImax1,1:ndim)
    double precision                   :: w(ixImin1:ixImax1,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: theoretical(ixImin1:ixImax1)
    double precision                   :: residual(ixImin1:ixImax1)

    theoretical(ixImin1:ixImax1) = spotpattern(x,ixImin1,ixImax1,global_time)
    residual(ixImin1:ixImax1) = abs(theoretical(ixImin1:ixImax1) - &
       w(ixImin1:ixImax1,r_e))/theoretical(ixImin1:ixImax1)

    w(ixOmin1:ixOmax1,nw+1) = theoretical(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,nw+2) = residual(ixOmin1:ixOmax1)

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
