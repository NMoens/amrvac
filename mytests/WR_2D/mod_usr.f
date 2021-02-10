!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_rhd
  use mod_fld

  implicit none

  double precision, parameter :: M_sun = 1.9891000d33
  double precision, parameter :: R_sun = 6.9599000d10
  double precision, parameter :: L_sun = 3.8268000d33
  double precision, parameter :: year = 365.25*24*60*60

  double precision :: StefBoltz

  double precision :: cak_Q, cak_a, cak_base, cak_x0, cak_x1
  integer :: it_start_cak
  double precision :: rho_bound, v_inf, Mdot
  double precision :: T_bound, R_star, M_star

  double precision :: kappa_e, L_bound, Gamma_e_bound, F_bound, gradE, E_out
  logical :: fixed_lum, Cak_in_D
  logical :: read_cak_table = .false.

  integer :: i_v1, i_v2, i_p
  integer :: i_Trad, i_Tgas, i_Mdot, i_Opal, i_CAK, i_CAK2, i_lambda
  integer :: i_Gamma, i_Lum, i_F1, i_F2

  double precision :: sum_time
  double precision, allocatable :: vr_sumt(:), rho_sumt(:), rho2_sumt(:),&
      sumt(:)

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_2D")

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! Boundary conditions
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions

    !> lasy fix for inexplicable pressure
    usr_internal_bc => Fix_pressure

    ! PseudoPlanar correction
    usr_source => PseudoPlanar

    ! Graviatational field
    usr_gravity => set_gravitation_field

    ! Special Opacity
    usr_special_opacity => OPAL_and_CAK

    ! Write out energy levels and temperature
    usr_write_analysis => collapse_to_1D

    !> Additional variables
    usr_process_grid => update_extravars

    ! Output routines
    ! usr_aux_output    => specialvar_output
    ! usr_add_aux_names => specialvarnames_output

    ! Timestep for PseudoPlanar
    ! usr_get_dt => get_dt_cak

    ! Refine mesh near base
    ! usr_refine_grid => refine_base

    ! Active the physics module
    call rhd_activate()

    i_v1 = var_set_extravar("v1", "v1")
    i_v2 = var_set_extravar("v2", "v2")
    if (rhd_energy) i_p = var_set_extravar("p","p")
    i_Trad = var_set_extravar("Trad", "Trad")
    if (rhd_energy) i_Tgas = var_set_extravar("Tgas", "Tgas")
    i_Mdot = var_set_extravar("Mdot", "Mdot")
    i_Opal = var_set_extravar("OPAL", "OPAL")
    i_CAK = var_set_extravar("CAK", "CAK")
    i_CAK2 = var_set_extravar("CAK2", "CAK2")
    i_lambda = var_set_extravar("lambda", "lambda")
    i_Gamma = var_set_extravar("Gamma", "Gamma")
    i_Lum = var_set_extravar("Lum", "Lum")
    i_F1 = var_set_extravar("F1", "F1")
    i_F2 = var_set_extravar("F2", "F2")

  end subroutine usr_init


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_opal_opacity, only: init_opal
    use mod_cak_opacity, only: init_cak

    use mod_fld

    integer :: i

    !> Initialise Opal
    call init_opal(He_abundance)

    call init_cak()

    !> read usr par
    call params_read(par_files)

    ! Choose independent normalization units if using dimensionless variables.
    unit_length  = R_star
    unit_numberdensity = rho_bound/((1.d0+4.d0*He_abundance)*const_mp)
    unit_velocity = v_inf

    !> Remaining units
    unit_density=(1.d0+4.d0*He_abundance)*const_mp*unit_numberdensity
    unit_pressure=unit_density*unit_velocity**2
    unit_temperature=unit_pressure/((2.d0+&
       3.d0*He_abundance)*unit_numberdensity*const_kB)
    unit_time=unit_length/unit_velocity

    unit_radflux = unit_velocity*unit_pressure
    unit_opacity = one/(unit_density*unit_length)

    R_star = R_star/unit_length
    M_star = M_star/unit_density/unit_length**3
    T_bound = T_bound/unit_temperature
    rho_bound = rho_bound/unit_density
    v_inf = v_inf/unit_velocity
    Mdot = Mdot/unit_density/unit_length**3*unit_time

    kappa_e = kappa_e/unit_opacity
    F_bound = F_bound/unit_radflux
    L_bound = L_bound/(unit_radflux*unit_length**2)

    StefBoltz = const_rad_a*const_c/4.d0*(unit_temperature**&
       4.d0)/(unit_velocity*unit_pressure)

    !> Very bad initial guess for gradE using kappa_e
    gradE = -F_bound*3*kappa_e*rho_bound*unit_velocity/const_c

    ! print*, 'L_bound (cgs)', L_bound*(unit_radflux*unit_length**2)
    ! print*, 'log10(L_bound)', log10(L_bound*(unit_radflux*unit_length**2)/L_sun)
    ! print*, 'L_bound', L_bound*(unit_radflux*unit_length**2)/L_sun
    ! ! stop
    ! print*, 'unit_density', unit_density
    ! print*, 'unit_time', unit_time
    ! print*, 'unit_pressure', unit_pressure

    sum_time = 0.d0

    allocate(vr_sumt(domain_nx1))
    vr_sumt = 0.d0
    allocate(rho_sumt(domain_nx1))
    rho_sumt = 0.d0
    allocate(rho2_sumt(domain_nx1))
    rho2_sumt = 0.d0

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ cak_Q, cak_a, cak_base, cak_x0, cak_x1, rho_bound,&
        kappa_e, T_bound, R_star, M_star, v_inf, Mdot, Gamma_e_bound,&
        it_start_cak, fixed_lum, Cak_in_D, read_cak_table

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wind_list, end=113)
113    close(unitpar)
    end do

    R_star = R_star*R_sun
    M_star = M_star*M_sun
    L_bound = Gamma_e_bound * 4.0 * dpi * const_G * M_star * const_c/kappa_e
    F_bound = L_bound/(4*dpi*R_star**2)
    Mdot = Mdot*M_sun/year

    ! print*, Gamma_e_bound, dpi, const_G, M_star, const_c

  end subroutine params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, x)
    use mod_constants
    use mod_physics

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2), fld_R(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    double precision :: vel(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: T_out, E_out, E_gauge
    double precision :: T_in, E_in, rr(ixImin1:ixImax1,ixImin2:ixImax2), bb,&
        temp(ixImin1:ixImax1,ixImin2:ixImax2)

    integer :: ii

!    w(ixI^S,rho_)   = rho_bound
!    w(ixI^S,mom(:)) = 0.d0
!
!    where(x(ixI^S,1) .gt. 1.d0)
!      vel(ixI^S) =  v_inf*(1 - 0.999d0*( 1.d0/x(ixI^S,1)))**0.5d0
!      w(ixI^S,rho_) = Mdot/(4.d0*dpi*vel(ixI^S)*x(ixI^S,1)**2)
!    endwhere
!
!    w(ixI^S,mom(1)) = vel(ixI^S)*w(ixI^S,rho_)
!
!    !> Outer/Inner temperature
!    T_out = 28445.836732569689/unit_temperature
!    E_out = const_rad_a*(T_out*unit_temperature)**4/unit_pressure
!    T_in = T_bound
!    E_in = const_rad_a*(T_in*unit_temperature)**4/unit_pressure
!
!    !>Very bad initial profile using constant gradE
!    ! w(ixI^S,r_e) = E_out + gradE*(x(ixI^S,1)-xprobmax1)
!    ! w(ixI^S,r_e) = E_in + gradE*(x(ixI^S,1)-xprobmin1)
!    rr(ixI^S) = dsqrt(x(ixI^S,1)-xprobmin1)*(16*x(ixI^S,1)**2 + 8*x(ixI^S,1)+ 6) &
!                /(15*x(ixI^S,1)**(5.d0/2))
!    bb = -kappa_e*L_bound*Mdot*unit_velocity*3/(16*dpi**2*const_c*v_inf)
!
!    ! bb = bb*2
!
!    ! rr(ixI^S) = 2*dsqrt(x(ixI^S,1)-xprobmin1)/dsqrt(x(ixI^S,1))
!    ! bb = kappa_e*F_bound*Mdot*unit_velocity*3/(4*dpi*const_c*v_inf)
!    w(ixI^S,r_e) = E_in + bb*rr(ixI^S)
!
!    E_gauge = E_in + bb*dsqrt(xprobmax1-xprobmin1)&
!              *(16*xprobmax1**2 + 8*xprobmax1+ 6)/(15*xprobmax1**(5.d0/2))
!
!    w(ixI^S,r_e) = w(ixI^S,r_e) - E_gauge + E_out
!
!    if (rhd_energy) then
!      temp(ixI^S) = (w(ixI^S,r_e)*unit_pressure/const_rad_a)**0.25d0/unit_temperature
!      w(ixI^S,e_) = w(ixI^S,rho_)*temp(ixI^S)/(rhd_gamma-1.d0) + half*w(ixI^S,mom(1))**2/w(ixI^S,rho_)
!    endif
!
!
!    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
!    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R)
!
!    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))
!    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)/(3.d0*kappa(ixO^S)*w(ixO^S,rho_))
!

    do ii = ixOmin1,ixOmax1
      w(ii,:,rho_) = read_initial_conditions(x(ii,3,1),2)
      w(ii,:,mom(1)) = read_initial_conditions(x(ii,3,1),3)
      w(ii,:,mom(2)) = 0.d0 !1.d-2*dsin(x(ii,:,1)-1.d0)
      if (rhd_energy) w(ii,:,e_) = read_initial_conditions(x(ii,3,1),4)
      w(ii,:,r_e) = read_initial_conditions(x(ii,3,1),5)
      w(ii,:,i_diff_mg) = read_initial_conditions(x(ii,3,1),6)
    enddo

    call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,x)

  end subroutine initial_conditions



  function read_initial_conditions(r_in,index) result(var)
    integer, intent(in) :: index
    double precision, intent(in) :: r_in
    double precision :: var

    double precision :: w(1:6), w_mo(1:6), w_po(1:6)
    integer :: ll

    w(:) = 0.d0

    open(unit=1, file='1D_stable.blk')
    read(1,*) !> header
    read(1,*) !>header
    read(1,*) !>header
    read(1,*) w !> first line of data
    do ll = 1,1024
      w_mo = w
      read(1,*) w
        if (w(1) .gt. min(9.d0,r_in)) then
          w_po = w
          goto 8765
        endif
      w_po = w
    enddo

8765 CLOSE(1)
    var = w_mo(index) + (w_po(index) - w_mo(index))/(w_po(1) - w_mo(1))*(r_in &
       - w_mo(1))

  end function read_initial_conditions


  subroutine boundary_conditions(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2,iB,w,x)
    use mod_physics, only: phys_get_trad
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixBmin1,ixBmin2,ixBmax1,ixBmax2, iB
    double precision, intent(in)    :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: kappa(ixImin1:ixImax1,ixImin2:ixImax2),&
        Temp(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: Temp0, rho0, T_out, n
    double precision :: Local_gradE(ixImin1:ixImax1,ixImin2:ixImax2),&
        F_adv(ixBmin1:ixBmax1,ixBmin2:ixBmax2)
    double precision :: Local_tauout(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: Local_Tout(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision :: kappa_out

    integer :: ix1,ix2

    select case (iB)

    case(1)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_) = rho_bound
      do ix1 = ixBmax1-1,ixBmin1,-1
        w(ix1,ixBmin2:ixBmax2,rho_) = dexp(2*dlog(w(ix1+1,ixBmin2:ixBmax2,&
           rho_)) - dlog(w(ix1+2,ixBmin2:ixBmax2,rho_)))
      enddo

      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(2)) = 0.d0

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,ixBmin2:ixBmax2,mom(1)) = w(ix1+1,ixBmin2:ixBmax2,mom(1))
       ! w(ix1,ixBmin2:ixBmax2,mom(2)) = w(ix1+1,ixBmin2:ixBmax2,mom(2))
      enddo

      where(w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1)) .lt. 0.d0)
       w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1)) = 0.d0
      endwhere

      F_adv(ixBmax1,ixBmin2:ixBmax2) = 4.d0/3.d0*(w(ixBmax1,ixBmin2:ixBmax2,&
         mom(1))/w(ixBmax1,ixBmin2:ixBmax2,rho_))*w(ixBmax1,ixBmin2:ixBmax2,&
         r_e) * 4*dpi*xprobmin1**2

      where (F_adv(ixBmin1:ixBmax1,ixBmin2:ixBmax2) .ne. F_adv(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)) F_adv(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = 0.d0
      where (F_adv(ixBmin1:ixBmax1,ixBmin2:ixBmax2) .le. 0.d0) &
         F_adv(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = 0.d0

      do ix1 = ixImin1,ixImax1
        Local_gradE(ix1,ixBmin2:ixBmax2) = -(F_bound-F_adv(ixBmax1,&
           ixBmin2:ixBmax2))/w(nghostcells+1,ixBmin2:ixBmax2,i_diff_mg)
      enddo
      gradE = sum(Local_gradE(nghostcells,ixBmin2:ixBmax2))/(ixBmax2-ixBmin2)

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,ixBmin2:ixBmax2,r_e) = w(ix1+2,ixBmin2:ixBmax2,r_e) + (x(ix1,&
           ixBmin2:ixBmax2,1)-x(ix1+2,ixBmin2:ixBmax2,1))*Local_gradE(ix1+1,&
           ixBmin2:ixBmax2)
      enddo

      if (rhd_energy) then
        temp(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = (w(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2,r_e)*unit_pressure/const_rad_a)**&
           0.25d0/unit_temperature
        w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,e_) = w(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2,rho_)*temp(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2)/(rhd_gamma-1.d0) + half*w(ixBmin1:ixBmax1,&
           ixBmin2:ixBmax2,mom(1))**2/w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)
      endif


   case(2)

      !> Compute mean kappa in outer blocks
      call get_kappa_OPAL(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
         ixImax1,ixImax2,w,x,kappa)

      kappa_out = sum(kappa(ixImax1-nghostcells,&
         ixBmin2:ixBmax2))/((ixBmax2-ixBmin2))

      if (kappa_out .ne. kappa_out) kappa_out = kappa_e
      kappa_out = max(kappa_out,kappa_e)
      kappa_out = min(kappa_out,2*kappa_e)

      Local_tauout(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = &
         kappa_out*w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         rho_)*R_star**2/(3*x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,1))
      Local_Tout(ixBmin1:ixBmax1,ixBmin2:ixBmax2) = &
         F_bound/StefBoltz*(3.d0/4.d0*Local_tauout(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2))**0.25d0

      T_out = sum(Local_Tout(ixBmin1,ixBmin2:ixBmax2))/(ixBmax2-ixBmin2)

!      T_out = max(1.5d4/unit_temperature, T_out)
!      T_out = max(4.2d4/unit_temperature, T_out)


      E_out = const_rad_a*(T_out*unit_temperature)**4.d0/unit_pressure

      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,r_e) = &
         const_rad_a*(Local_Tout(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2)*unit_temperature)**4.d0/unit_pressure

      ! print*, 'E_out', E_out
    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

  subroutine mg_boundary_conditions(iB)
    use mod_global_parameters
    use mod_multigrid_coupling

    integer, intent(in)             :: iB
    integer          :: iigrid, igrid
    double precision :: snd_g, rcv_g, snd_e, rcv_e

    select case (iB)
    case (1)
      !> communicate boundary
      ! do iigrid=1,igridstail; igrid=igrids(iigrid);
      !   if (node(pig1_,igrid) > 1) &
      !     gradE = 0.d0
      ! enddo
      !
      ! snd_g = gradE
      ! call mpi_allreduce(snd_g,rcv_g,1,&
      ! MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)
      ! gradE = rcv_g

      mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
      mg%bc(iB, mg_iphi)%bc_value = gradE


    case (2)
      !> communicate boundary
      ! do iigrid=1,igridstail; igrid=igrids(iigrid);
      !   if (node(pig1_,igrid) < domain_nx1/block_nx1*node(plevel_,igrid)) &
      !     E_out = bigdouble
      ! enddo
      !
      ! snd_e = E_out
      ! call mpi_allreduce(snd_e,rcv_e,1,&
      ! MPI_DOUBLE_PRECISION,MPI_MIN,icomm,ierrmpi)
      ! E_out = rcv_e


      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
      mg%bc(iB, mg_iphi)%bc_value = E_out

    case default
      print *, "Not a standard: ", trim(typeboundary(r_e, iB))
      error stop "Set special bound for this Boundary "
    end select
  end subroutine mg_boundary_conditions


  subroutine Fix_pressure(level,qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,level
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2), mean_p

    !> fix density as well
    where ((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_) .lt. 1.d-5) .and. (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1) .lt. 4.d0))
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) = 1.d-5
    endwhere

    if (.not. rhd_energy) call mpistop("no energy equation, no pressure fix")

    call phys_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,pth)

    !if (any(press(ixO^S) .lt. 0.d0)) then
      mean_p = sum(pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2))/(block_nx1*block_nx2)
      where (pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .le. small_pressure + &
         smalldouble)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) = mean_p/(rhd_gamma - 1) +  0.5d0 * sum(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, mom(:))**2, dim=ndim+1)/w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,rho_)
      end where
    !endif

    !> fix density as well
!    where (pth(ixO^S)/w(ixO^S,rho_) .gt. 3.d-5)
!      w(ixO^S,rho_) = pth(ixO^S)/
!    endwhere



  end subroutine Fix_pressure


  !> Calculate gravitational acceleration in each dimension
  subroutine set_gravitation_field(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)

    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: mass

    radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,:) = 0.d0
    gravity_field(ixImin1:ixImax1,ixImin2:ixImax2,&
       1) = -const_G*mass/radius(ixImin1:ixImax1,&
       ixImin2:ixImax2)**2*(unit_time**2/unit_length)

  end subroutine set_gravitation_field

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine PseudoPlanar(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, iwmin,iwmax
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision :: ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nw)

    double precision :: k_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndim)

    call PseudoPlanarSource(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,wCT,x,ppsource)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_) + qdt*ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_) !> OK
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(1)) + qdt*ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(1)) !> OK
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2)) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(2)) + qdt*ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(2)) !> OK
    if (rhd_energy) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,e_) + qdt*ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) !> OK
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       r_e) + qdt*ppsource(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) !> TROUBLEMAKER

    if (.not. Cak_in_D) then
      if (read_cak_table) then
        call get_kappa_CAK2(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,wCT,x,k_cak)
      else
        call get_kappa_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,wCT,x,k_cak)
      endif

      if (fixed_lum) then
        !> Fixed L = L_bound
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) = w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(1)) + qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           rho_)*L_bound/(4*dpi*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1)**2)/const_c*k_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*unit_velocity
        if (rhd_energy) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,e_) + qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(1))*L_bound/(4*dpi*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             1)**2)/const_c*k_cak(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)*unit_velocity
        endif
      else
        !> Local flux
        call fld_get_radflux(wCT, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, rad_flux)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) = w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(1)) + qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           rho_)*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1)/const_c*k_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*unit_velocity
        if (rhd_energy) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,e_) + qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(1))*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             1)/const_c*k_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*unit_velocity
        endif
      endif

    endif

  end subroutine PseudoPlanar

  subroutine get_dt_cak(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: dx1,dx2, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: mass

    double precision :: dt_cak
    double precision :: k_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndim)

    if (read_cak_table) then
      call get_kappa_CAK2(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x,k_cak)
    else
      call get_kappa_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x,k_cak)
    endif

    if (fixed_lum) then
      !> Fixed L = L_bound
      dt_cak = courantpar*minval(dsqrt(dxlevel(1)/(L_bound/(4*dpi*x(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)**2)/const_c*k_cak(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*unit_velocity -const_G*mass/radius(ixImin1:ixImax1,&
         ixImin2:ixImax2)**2*(unit_time**2/unit_length))))
    else
      !> Local flux
      call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, rad_flux)
      dt_cak = courantpar*minval(dsqrt(dxlevel(1)/abs(rad_flux(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1)/const_c*k_cak(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*unit_velocity -const_G*mass/radius(ixImin1:ixImax1,&
         ixImin2:ixImax2)**2*(unit_time**2/unit_length))))
    endif

    dtnew = min(dt_cak, dtnew)

  end subroutine get_dt_cak
  !

  subroutine PseudoPlanarSource(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,x,source)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out) :: source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1:nw)

    double precision :: rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndir)
    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),v(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:ndim)
    double precision :: radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
         pert(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision :: edd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndim,1:ndim)
    integer :: rdir, pdir

    source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nw) = zero

    rdir = 1
    pdir = 2

    v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rdir) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(rdir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,pdir) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(pdir))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rdir) !+ half*block%dx(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rdir)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) = -two*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rdir)/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    call phys_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,pth)

    !> dm_r/dt = +(rho*v_p**2 + 2pth)/r -2 (rho*v_r**2 + pth)/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(rdir)) = - 2*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rdir)**two/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2) + w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rho_)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       pdir)**two/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(pdir)) = - 3*v(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,rdir)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       pdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)/radius(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    !> de/dt = -2 (e + p)v_r/r
    if (rhd_energy) source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_) = -two*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)+pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rdir)/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    !> dEr/dt = -2 (E v_r + F_r)/r
    if (rhd_radiation_diffusion) then
      call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, rad_flux)
      source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = source(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,r_e) - two*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rdir)/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    endif

    if (rhd_radiation_advection) then
      source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = source(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,r_e) - two*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         r_e)*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rdir)/radius(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
    endif

    ! Not sure about this one
    if (rhd_radiation_force) then
      call fld_get_eddington(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, edd)
      source(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = source(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,r_e) + two*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rdir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e)*edd(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1,1)/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    endif

  end subroutine PseudoPlanarSource


  subroutine OPAL_and_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,kappa)
    use mod_physics, only: phys_get_trad
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    double precision :: OPAL(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        CAK(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    !> Get OPAL opacities by reading from table
    call get_kappa_OPAL(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,x,OPAL)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    if (Cak_in_D) then
      if (read_cak_table) then
        call get_kappa_CAK2(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,w,x,CAK)
      else
        call get_kappa_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,w,x,CAK)
      endif
    else
      CAK(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.d0
    endif

    !> Add OPAL and CAK for total opacity
    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = OPAL(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) + CAK(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    where(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .ne. kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))
      kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = kappa_e
    endwhere

  end subroutine OPAL_and_CAK


  subroutine get_kappa_OPAL(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,kappa)
    use mod_physics, only: phys_get_trad, phys_get_tgas
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    integer :: ix1,ix2
    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: n, rho0, Temp0

    !> Get OPAL opacities by reading from table
    if (rhd_energy) then
      call phys_get_tgas(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,Temp)
    else
      call phys_get_trad(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,Temp)
    endif

    do ix1=ixOmin1,ixOmax1
     do ix2=ixOmin2,ixOmax2
    
        rho0 = w(ix1,ix2,rho_)*unit_density
        Temp0 = Temp(ix1,ix2)*unit_temperature
        Temp0 = max(Temp0,1.d4)
        call set_opal_opacity(rho0,Temp0,n)
        kappa(ix1,ix2) = n/unit_opacity
    enddo
     enddo
    

    where(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .ne. kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))
      kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = kappa_e
    endwhere

    !> Lower limit is electron scattering:
    where(kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .lt. kappa_e)
      kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = kappa_e
    endwhere


    !> test without opal
    ! kappa(ixO^S) = kappa_e

  end subroutine get_kappa_OPAL

  subroutine get_kappa_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,kappa)
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    integer :: ix1,ix2
    double precision :: vel(ixImin1:ixImax1,ixImin2:ixImax2),&
        gradv(ixOmin1:ixOmax1,ixOmin2:ixOmax2), gradvI(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixImin1:ixImax1,ixImin2:ixImax2) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(1))/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    call gradientO(vel,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,1,gradv,nghostcells)

    ! call gradient(vel,ixI^L,ixO^L,1,gradvI)
    ! gradv(ixO^S) = gradvI(ixO^S)

    !> Absolute value of gradient:
    gradv(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs(gradv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))

    xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 1.d0-xprobmin1/x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)

    alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cak_a

    where (xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .le. cak_x0)
      alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cak_base
    elsewhere (xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .le. cak_x1)
      alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cak_base + (cak_a - &
         cak_base)*(xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2) - cak_x0)/(cak_x1 - &
         cak_x0)
    endwhere

    kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
       kappa_e*cak_Q/(1-alpha(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)) *(gradv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*unit_velocity/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)*const_c*cak_Q*kappa_e))**alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if (x(ixImax1,nghostcells,1) .ge. xprobmax1) then
      kappa(ixOmax1,ixOmin2:ixOmax2) = kappa(ixOmax1-1,ixOmin2:ixOmax2)
    endif

    ! if (it .le. it_start_cak) then
    !   kappa(ixO^S) = kappa(ixO^S)*dexp(-w(ixO^S,rho_)*kappa_e)
    ! endif

    do ix1=ixOmin1,ixOmax1
     do ix2=ixOmin2,ixOmax2
    
        if (xx(ix1,ix2) .lt. cak_x0) then
          kappa(ix1,ix2) = min(2*kappa_e,kappa(ix1,ix2))
        else if (xx(ix1,ix2) .lt. cak_x1) then
          kappa(ix1,ix2) = min( (2 + (10.d0 - 2.d0)/(cak_x1 - cak_x0)*(xx(ix1,&
             ix2)-cak_x0) )*kappa_e, kappa(ix1,ix2))
        else
          kappa(ix1,ix2) = min(10*kappa_e,kappa(ix1,ix2))
        endif

      ! if (kappa(ix^D) .gt. 3*kappa_e) &
      !   kappa(ix^D) = 3*kappa_e + half*(kappa(ix^D)-3*kappa_e)
    enddo
     enddo
    

    !> test with no cak
    ! kappa(ixO^S) = 0.d0
  end subroutine get_kappa_CAK

  subroutine get_kappa_CAK2(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,kappa)
    use mod_physics, only: phys_get_trad, phys_get_tgas
    use mod_global_parameters
    use mod_cak_opacity
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: kappa(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    double precision :: Temp(ixImin1:ixImax1,ixImin2:ixImax2), rho0, temp0,&
        gradv0, kap0
    integer :: ix1,ix2

    double precision :: vel(ixImin1:ixImax1,ixImin2:ixImax2),&
        gradv(ixOmin1:ixOmax1,ixOmin2:ixOmax2), gradvI(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: xx(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixImin1:ixImax1,ixImin2:ixImax2) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(1))/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    call gradientO(vel,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,1,gradv,nghostcells)

    ! call gradient(vel,ixI^L,ixO^L,1,gradvI)
    ! gradv(ixO^S) = gradvI(ixO^S)

    !> Absolute value of gradient:
    gradv(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs(gradv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))

    !> Get CAK opacities by reading from table
    if (rhd_energy) then
      call phys_get_tgas(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,Temp)
    else
      call phys_get_trad(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,Temp)
    endif

    do ix1=ixOmin1,ixOmax1
     do ix2=ixOmin2,ixOmax2
    
        rho0 = w(ix1,ix2,rho_)*unit_density
        Temp0 = Temp(ix1,ix2)*unit_temperature
        Temp0 = max(Temp0,1.d4)
        gradv0 = gradv(ix1,ix2)*(unit_velocity/unit_length)
        call set_cak_opacity(rho0,Temp0,gradv0,kap0)
        kappa(ix1,ix2) = kap0/unit_opacity
    enddo
     enddo
    

    if (x(ixImax1,nghostcells,1) .ge. xprobmax1) then
      kappa(ixOmax1,:) = kappa(ixOmax1-1,:)
    endif


  end subroutine get_kappa_CAK2

  subroutine collapse_to_1D()
    use mod_global_parameters

    integer          :: iigrid, igrid, jj_blk,ii, jj, nbx, i, il, ih, dn_mdot
    integer          :: lvl, ibx
    integer          :: np_mdot, nc
    double precision :: ratio, sf_mdot, dx_l1

    double precision, allocatable :: rp_mdot(:), p_mdot(:)
    double precision, allocatable :: mdot_S(:), mdot_R(:)
    integer, allocatable :: jp_mdot(:)

    double precision :: rr(1:domain_nx1), rr_S(1:domain_nx1),&
        rr_R(1:domain_nx1)
    double precision :: vr(1:domain_nx1), vr_S(1:domain_nx1),&
        vr_R(1:domain_nx1)
    double precision :: mdot(1:domain_nx1)
    double precision :: rho(1:domain_nx1), rho_S(1:domain_nx1),&
        rho_R(1:domain_nx1)
    double precision :: rho2(1:domain_nx1), rho2_S(1:domain_nx1),&
        rho2_R(1:domain_nx1)

    integer :: lvl_h(1:domain_nx1), lvl_h_S(1:domain_nx1),&
        lvl_h_R(1:domain_nx1)
    integer :: lvl_l(1:domain_nx1), lvl_l_S(1:domain_nx1),&
        lvl_l_R(1:domain_nx1)

      ! if (refine_max_level .ne. 1) &
      ! call mpistop("collapse_to_1D doesnt work YET with mpi")

    !> #R_star -1 in simulation
    np_mdot = floor((xprobmax1-xprobmin1)/R_star)

    allocate(rp_mdot(1:np_mdot))
    allocate(p_mdot(1:np_mdot))
    allocate(jp_mdot(1:np_mdot))
    allocate(mdot_S(1:np_mdot))
    allocate(mdot_R(1:np_mdot))

    rr = 0.d0
    rp_mdot = 0.d0
    vr = 0.0d0
    rho = 0.d0
    rho2 = 0.d0
    mdot = 0.d0
    p_mdot = 0.d0
    lvl_h = 0.d0
    lvl_l = 2* refine_max_level

    !> Reconstruct radius at level 1
    do jj = 1,domain_nx1
      rr(jj) = xprobmin1 + (jj-0.5d0)*(xprobmax1-xprobmin1)/domain_nx1
    enddo

    !> Choose radii at which to save mdot
    !> This is done at every stellar radii
    do jj = 1,np_mdot
      ! rp_mdot(jj) = xprobmin1 + (jj-1)* R_star
      dx_l1 = (xprobmax1-xprobmin1)/domain_nx1
      dn_mdot = floor(1.d0*domain_nx1/np_mdot+smalldouble)
      rp_mdot(jj) = xprobmin1 + dx_l1/2 + (jj-1)*dn_mdot*dx_l1
    enddo

    !> Find cells that are closest to radii where we want to track mdot
    do jj = 1,domain_nx1
      do ii = 1,np_mdot
        if (rr(jj) - rp_mdot(ii) .le. smalldouble) jp_mdot(ii) = jj
          ! rp_mdot(ii) = rr(jj)
      enddo
    enddo

    !> loop over all grids for this proc
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)

      !> widht in cells per cell of lvl 1: On lvl 1: 1, lvl 2: 2, lvl 3: 4
      nc = 2**(node(plevel_,igrid)-1)

      !> Map the block index to the global index
      !> nbx = nr of blocks between i=0 and current block
      jj_blk = (node(pig1_,igrid)-1)*block_nx1/nc

      if (nc > block_nx1) call mpistop&
         ("collapse_to_1D doesnt work, reduce amr")

      !> For all cells in the current block, average velocity over lateral direction.
      !> Take into account possible amr
      do ii = 1,block_nx1/nc
        !> relate global index jj to local index ii on grid
        il = nghostcells + 1 + (ii-1)*nc
        ih = nghostcells + 1 + (ii)*nc -1
        jj = jj_blk + ii

        !> Velocity
        vr(jj) = vr(jj) + sum(block%w(il:ih,ixMlo2:ixMhi2,&
           mom(1))/block%w(il:ih,ixMlo2:ixMhi2,rho_))/domain_nx2/nc**2 !> average value on lvl 1

        !> Density
        rho(jj) = rho(jj) + sum(block%w(il:ih,ixMlo2:ixMhi2,&
           rho_))/domain_nx2/nc**2 !> average value on lvl 1

        !> Density squared
        rho2(jj) = rho2(jj) + sum(block%w(il:ih,ixMlo2:ixMhi2,&
           rho_)**2)/domain_nx2/nc**2 !> average value on lvl 1

        !> Mass loss rate
        mdot(jj) = mdot(jj) + 4*dpi*sum(block%x(il:ih,ixMlo2:ixMhi2,&
           1)**2*block%w(il:ih,ixMlo2:ixMhi2,mom(1)))/domain_nx2/nc**2 !> average value on lvl 1

        !> Highest/Lowest amr level
        lvl_h(jj) = node(plevel_,igrid)
        lvl_l(jj) = node(plevel_,igrid)
      enddo

    enddo

    !> communicate velocity
    vr_S=vr
    call mpi_allreduce(vr_S,vr_R,domain_nx1,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,&
       ierrmpi)
    vr=vr_R

    !> communicate density
    rho_S=rho
    call mpi_allreduce(rho_S,rho_R,domain_nx1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       icomm,ierrmpi)
    rho=rho_R

    !> communicate density squared
    rho2_S=rho2
    call mpi_allreduce(rho2_S,rho2_R,domain_nx1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       icomm,ierrmpi)
    rho2=rho2_R

    !> Only keep mdot in interested radii
    do ii = 1,np_mdot
      jj = jp_mdot(ii)
      p_mdot(ii) = mdot(jj)
    enddo

    !> communicate mdot array
    mdot_S=p_mdot
    call mpi_reduce(mdot_S,mdot_R,np_mdot,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,&
       ierrmpi)
    p_mdot=mdot_R*unit_density*unit_velocity*unit_length**2/M_sun*year

    !> communicate highest amr level
    lvl_h_S=lvl_h
    call mpi_allreduce(lvl_h_S,lvl_h_R,domain_nx1,MPI_INTEGER,MPI_MAX,icomm,&
       ierrmpi)
    lvl_h = lvl_h_R
    lvl_l_S=lvl_l
    call mpi_allreduce(lvl_l_S,lvl_l_R,domain_nx1,MPI_INTEGER,MPI_MIN,icomm,&
       ierrmpi)
    lvl_l = lvl_l_R

    !> integrate over time
    sum_time = sum_time + dt
    vr_sumt(1:domain_nx1) = vr_sumt(1:domain_nx1) + vr(1:domain_nx1)*dt
    rho_sumt(1:domain_nx1) = rho_sumt(1:domain_nx1) + rho(1:domain_nx1)*dt
    rho2_sumt(1:domain_nx1) = rho2_sumt(1:domain_nx1) + rho2(1:domain_nx1)*dt

    !> Write out average velocity profile
    if (mype==0) then
      !> Always update file to give the last mean snapshot
      open(unit=unitanalysis,file=trim(base_filename)//'_vr',status='replace')
      write(unitanalysis,*)&
          'r | vr_int_theta | vr_int_theta_int_t | lvl_h lvl_l'
      do i=1,domain_nx1
        write(unitanalysis,'(3f11.7,3i4)') rr(i), vr(i), vr_sumt(i)/sum_time,&
            lvl_h(i), lvl_l(i)
      enddo
    close(unitanalysis)
    endif

    !> Write out average density profile
    if (mype==0) then
      !> Always update file to give the last mean snapshot
      open(unit=unitanalysis,file=trim(base_filename)//'_rho',&
         status='replace')
      write(unitanalysis,*) 'r | rho_int_theta | rho_int_theta_int_t'
      do i=1,domain_nx1
        write(unitanalysis,'(6f11.7)') rr(i), rho(i), rho_sumt(i)/sum_time,&
            rho2(i), rho2_sumt(i)/sum_time,&
            sum_time*rho2_sumt(i)/rho_sumt(i)**2
      enddo
    close(unitanalysis)
    endif

    !> Write out average mdot at different radii
    if (mype==0) then
       if (global_time<smalldouble) then ! if very 1st iteration
         open(unit=unitanalysis,file=trim(base_filename)//'_mdot',&
            status='replace')
         write(unitanalysis,'(a16,*(f7.2))') 'Mdot at t | r= ', rp_mdot
         write(unitanalysis,'(f10.3,*(E15.6))') global_time, p_mdot
       else
         open(unit=unitanalysis,file=trim(base_filename)//'_mdot',&
            access='append')
         write(unitanalysis,'(f10.3,*(E15.6))') global_time, p_mdot
       endif
       close(unitanalysis)
     endif

  end subroutine collapse_to_1D


  subroutine update_extravars(igrid,level,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision                   :: g_rad(ixImin1:ixImax1,&
       ixImin2:ixImax2), big_gamma(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                   :: g_grav(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                   :: Tgas(ixImin1:ixImax1,&
       ixImin2:ixImax2),Trad(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                   :: kappa(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2), OPAL(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        CAK(ixOmin1:ixOmax1,ixOmin2:ixOmax2), CAK2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    double precision                   :: vel(ixImin1:ixImax1,ixImin2:ixImax2),&
        gradv(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision                   :: rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:ndim), Lum(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision                   :: pp_rf(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2), lambda(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        fld_R(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    integer                            :: idim
    double precision :: radius(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: mass

    radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = x(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    call fld_get_opacity(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, kappa)
    call fld_get_radflux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, rad_flux)

    if (rhd_energy)    call rhd_get_tgas(w, x, ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, Tgas)
    call rhd_get_trad(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, Trad)

    call get_kappa_OPAL(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,x,OPAL)
    call get_kappa_CAK(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,w,x,CAK)
    call get_kappa_CAK2(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,x,CAK2)

    g_rad(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (OPAL(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)+CAK(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))*rad_flux(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1)/(const_c/unit_velocity)
    g_grav(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
       const_G*mass/radius(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2*(unit_time**2/unit_length)
    big_gamma(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = g_rad(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/g_grav(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    vel(ixImin1:ixImax1,ixImin2:ixImax2) = w(ixImin1:ixImax1,ixImin2:ixImax2,&
       mom(1))/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    call gradient(vel,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,1,gradv)

    pp_rf(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = two*rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)*dt

    call fld_get_fluxlimiter(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, lambda, fld_R)

    Lum(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 4*dpi*rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1)*unit_length)**2*unit_radflux/L_sun

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_v1) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(1))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_v2) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(2))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    if (rhd_energy) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       i_p) = (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_) - 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:))**2,&
        dim=ndim+1) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        rho_)) *(rhd_gamma - 1)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_Trad) = Trad(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*unit_temperature
    if (rhd_energy) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       i_Tgas) = Tgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*unit_temperature
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_Mdot) = 4*dpi*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(1))*radius(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2 *unit_density*unit_velocity/M_sun*year
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_Opal) = OPAL(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/kappa_e
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_CAK) = CAK(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/kappa_e
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_CAK2) = CAK2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/kappa_e
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_lambda) = lambda(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_Gamma) = big_gamma(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_Lum) = Lum(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_F1) = rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1)/F_bound
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i_F2) = rad_flux(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,2)/F_bound

  end subroutine update_extravars

end module mod_usr
