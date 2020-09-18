!> This is a template for a new user problem
MODULE mod_usr

  ! Include a physics module
  USE mod_HD

  IMPLICIT NONE

  ! Custom dimentional variables
  ! Mass of the sun
  DOUBLE PRECISION :: M_sun_CGS = 1.988550d33
  ! Luminocity of the sun
  DOUBLE PRECISION :: L_sun_CGS = 3.828000d33
  ! Photospheric radiuse of the sun
  DOUBLE PRECISION :: R_sun_CGS = 6.957000d10
  ! Stefan-Boltzman constant
  DOUBLE PRECISION :: S_SB_CGS  = 5.670367d-5
  ! Newtons' gravitation constant
  DOUBLE PRECISION :: G_dp_CGS  = 6.671910d-8

  ! Stellar parameters in CGS units, Mass, Luminosity, Effective temperature
  ! Radiuse
  DOUBLE PRECISION :: M_star_CGS
  DOUBLE PRECISION :: L_star_CGS
  DOUBLE PRECISION :: T_star_CGS
  DOUBLE PRECISION :: R_star_CGS

  ! Thompson electron scattering opacity
  DOUBLE PRECISION :: kappa_e_CGS
  ! Computational domain lower boundary sound speed
  ! to limit the lower boundary momentum in the special_bound subroutine
  DOUBLE PRECISION :: Lower_bound_Cs_CGS
  ! Boundary dencity in CGS units, = input boundary dencity * unit_density
  DOUBLE PRECISION :: rho_bound_CGS

  ! Usefule normalisation units
  DOUBLE PRECISION :: unit_luminosity_my
  DOUBLE PRECISION :: unit_mass_my
  DOUBLE PRECISION :: unit_opacity_my


  ! Custom dimention-LESS variables
  ! see the comments above
  DOUBLE PRECISION :: M_sun_DML
  DOUBLE PRECISION :: L_sun_DML
  DOUBLE PRECISION :: R_sun_DML
  DOUBLE PRECISION :: S_SB_DML
  DOUBLE PRECISION :: G_dp_DML

  DOUBLE PRECISION :: M_star_DML
  DOUBLE PRECISION :: L_star_DML
  DOUBLE PRECISION :: T_star_DML
  DOUBLE PRECISION :: R_star_DML

  DOUBLE PRECISION :: rho_bound_DML
  ! The speed of light in unites of normalisation
  DOUBLE PRECISION :: const_c_DML
  DOUBLE PRECISION :: kappa_e_DML

  DOUBLE PRECISION :: Lower_bound_Cs_DML
  ! Mass-loss rate in unites of normalisation - used to compurute initial dencity
  DOUBLE PRECISION :: M_dot_DML
  ! terminal velocity = escape velocity - used for initial velocity
  DOUBLE PRECISION :: v_inf_DML

  ! intut variables form  *.par file
  DOUBLE PRECISION :: M_star
  DOUBLE PRECISION :: R_star
  DOUBLE PRECISION :: rho_bound
  DOUBLE PRECISION :: kappa_e
  DOUBLE PRECISION :: mu_molec
  DOUBLE PRECISION :: Gamma_e
  DOUBLE PRECISION :: Q_cak
  DOUBLE PRECISION :: a_cak_min
  DOUBLE PRECISION :: a_cak_max
  DOUBLE PRECISION :: a_cak_x0
  DOUBLE PRECISION :: a_cak_x1
  DOUBLE PRECISION :: G_init
  DOUBLE PRECISION :: ecut_time
  DOUBLE PRECISION :: avst_time !include this one in the .par file and assign value
  CHARACTER(len=std_len) :: TABLE_name

  DOUBLE PRECISION :: new_timestep

  DOUBLE PRECISION :: Opal_K_list(19,70)
  DOUBLE PRECISION :: Opal_T_list(70)
  DOUBLE PRECISION :: Opal_R_list(19)

  ! Arrays for statistical computations, make global to update properly
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: tmp1, tmp2, tmp3, tmp4 !temporary
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xc, dvav, vav, Tav, Gav

  DOUBLE PRECISION :: avst_time_DML !dimensionless form of stat time to be used in code
  DOUBLE PRECISION :: tnorm !
  DOUBLE PRECISION , DIMENSION(:), ALLOCATABLE :: Temperature_test_arrey
  DOUBLE PRECISION , DIMENSION(:), ALLOCATABLE :: Gcak_test_arrey
  ! ! Identifiers for the output vectors
  INTEGER :: i_g_cak, i_temper, i_kap, i_tau_sp, i_a_cak, i_P_thrm

CONTAINS

  !> This routine should set user methods, and activate the physics module
  SUBROUTINE usr_init()
    USE mod_global_parameters
    USE mod_usr_methods

    ! Choose coordinate system as 1D spherical with three components for vectors
    CALL set_coordinate_system("spherical")

    !afdasfas
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Special Boundary conditions
    usr_special_bc => special_bound

    ! CAK force
    usr_source => tot_Force

    ! adjusted timestep
    usr_get_dt => special_dt

    !menual ovverride of the p thermal
    usr_set_pthermal => p_thermal

    !menual ovverride of the thermal speed
    ! usr_set_csound2 => manual_csound2

    ! Invoke every 'ditsave_custom' a routine to compute statistics (not store!)
    usr_write_analysis => Stat_analysis
    ! Store the statistical quantities every .dat file save time
    usr_add_aux_names  => myextravarnames_output
    usr_aux_output     => myextravar_output

    CALL usr_params_read(par_files)

    ! Choose independent normalization units if using dimensionless variables.
    unit_length        = R_sun_CGS ! cm
    unit_temperature   = 1.0d5 ! K
    unit_numberdensity = rho_bound_CGS / ((1.d0+4.d0*He_abundance)*mp_cgs ) !cm-3

    ! Active the physics module
    CALL HD_activate()

    ! output variables
    ! g_cak CAK acceleration
    i_g_cak  = var_set_extravar("g_cak" , "g_cak" )
    ! Temperature
    i_temper = var_set_extravar("Temp"  , "Temp"  )
    ! Effective opacity
    i_kap    = var_set_extravar("Kap"   , "Kap"   )
    ! sphericaly modified Effective optical depth
    i_tau_sp = var_set_extravar("tau_sp", "tau_sp")
    ! Radial CAK power index
    i_a_cak  = var_set_extravar("a_cak" , "a_cak" )
    ! Pressure
    i_P_thrm = var_set_extravar("P_trm" , "P_trm" )
  END SUBROUTINE usr_init
  !=============================================================================









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Read the used defined parameters form the *.par file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE usr_params_read(files)
    USE mod_global_parameters, ONLY: unitpar
    USE mod_constants
    CHARACTER(len=*), INTENT(in) :: files(:)
    INTEGER                      :: n

    NAMELIST /star_list/ M_star, R_star, rho_bound, Gamma_e, kappa_e, mu_molec,&
       G_init, Q_cak, a_cak_min, a_cak_max, a_cak_x0, a_cak_x1, ecut_time,&
        avst_time, TABLE_name

    DO n = 1, SIZE(files)
       OPEN(unitpar, file=TRIM(files(n)), status="old")
       !open(unitpar, file=files, status="old")
       READ(unitpar, star_list, END=111)
111    CLOSE(unitpar)
    END DO

    M_star_CGS    = M_star * M_sun_CGS
    R_star_CGS    = R_star * R_sun_CGS
    kappa_e_CGS   = kappa_e
    rho_bound_CGS = rho_bound
    L_star_CGS    = Gamma_e * 4.0 * dpi * G_dp_CGS * M_star_CGS * &
       const_c/kappa_e_CGS
    T_star_CGS    = ( L_star_CGS/(4.0 * dpi * R_star_CGS**2 * S_SB_CGS) &
       )**0.25

    PRINT*,&
        "Input in CGS: #####################################################"

    PRINT*, 'M_star   ', M_star_CGS
    PRINT*, 'R_star   ', R_star_CGS
    PRINT*, 'L_star   ', L_star_CGS
    PRINT*, 'T_star   ', T_star_CGS
    PRINT*, 'rho_bound', rho_bound_CGS
    PRINT*, 'Gamma_e  ', Gamma_e
    PRINT*, 'kappa_e  ', kappa_e
    PRINT*, 'mu_molec ', mu_molec
    PRINT*, 'Gamma_i  ', G_init
    PRINT*, 'Q_cak    ', Q_cak
    PRINT*, 'a_cak_max', a_cak_min
    PRINT*, 'a_cak_min', a_cak_max
    PRINT*, 'a_cak_x0 ', a_cak_x0
    PRINT*, 'a_cak_x1 ', a_cak_x1
  END SUBROUTINE usr_params_read
  !=============================================================================








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Create dimentionall less quantities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE initglobaldata_usr
    USE mod_global_parameters
    USE mod_constants

    ! define usefule dimentions
    unit_mass_my          = unit_density * unit_length**3
    unit_opacity_my       = one / ( unit_density * unit_length )
    unit_luminosity_my    = unit_mass_my * unit_velocity**2 / unit_time

    ! create dimentionless quantities
    rho_bound_DML  = rho_bound_CGS /unit_density
    M_star_DML     = M_star_CGS    /unit_mass_my
    R_star_DML     = R_star_CGS    /unit_length
    T_star_DML     = T_star_CGS    /unit_temperature
    L_star_DML     = L_star_CGS    /unit_luminosity_my
    kappa_e_DML    = kappa_e_CGS   / unit_opacity_my
    const_c_DML    = const_c       / unit_velocity
    G_dp_DML       = G_dp_CGS * unit_time**2.0d0 * unit_mass_my / &
       unit_length**3.0d0
    v_inf_DML      = SQRT(2*M_star_CGS*G_dp_CGS*(G_init - 1)/R_star_CGS)/ &
       unit_velocity
    M_dot_DML      = 4.0d0 * dpi * rho_bound_DML * &
       R_star_DML**2.0d0*1.0d-3*v_inf_DML
    ! v_inf_DML      = 2.0d0 * dsqrt(M_star_DML*G_dp_DML/(2 * R_star_DML) ) !2.0d8/unit_velocity

    ! time mark for tharting the computation of avarege values
    avst_time_DML= avst_time/unit_time
    new_timestep = dt

    PRINT*, "UNITS: #####################################################"
    PRINT*, 'mp_cgs   ', mp_cgs
    PRINT*, 'const_kB ', const_kB
    PRINT*, 'H_abound ', He_abundance
    PRINT*, "UNITS: #####################################################"

    PRINT*, 'unit_time        ', unit_time
    PRINT*, 'unit_length      ', unit_length
    PRINT*, 'unit_density     ', unit_density
    PRINT*, 'unit_pressure    ', unit_pressure
    PRINT*, 'unit_velocity    ', unit_velocity
    PRINT*, 'unit_temperature ', unit_temperature

    ! Allocate space for statistical arrays
    ALLOCATE (tmp1(domain_nx1))
    ALLOCATE (tmp2(domain_nx1))
    ALLOCATE (tmp3(domain_nx1))
    ALLOCATE (tmp4(domain_nx1))

    ALLOCATE (xc(domain_nx1))
    ALLOCATE (dvav(domain_nx1))
    ALLOCATE (vav(domain_nx1))
    ALLOCATE (Tav(domain_nx1))
    ALLOCATE (Gav(domain_nx1))
    ALLOCATE (Temperature_test_arrey(domain_nx1))
    ALLOCATE (Gcak_test_arrey(domain_nx1))

  END SUBROUTINE initglobaldata_usr
  !===============================================================================





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Assign initial conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE initial_conditions(ixImin1,ixImax1, ixOmin1,ixOmax1, w, x)
    USE mod_global_parameters
    USE mod_physics, ONLY: phys_get_pthermal

    INTEGER         , INTENT(in)    :: ixImin1,ixImax1, ixOmin1,ixOmax1
    DOUBLE PRECISION, INTENT(in)    :: x(ixImin1:ixImax1,1:ndim)
    DOUBLE PRECISION, INTENT(inout) :: w(ixImin1:ixImax1,1:nw)

    DOUBLE PRECISION :: velocity_field(ixImin1:ixImax1),pth(ixImin1:ixImax1)
    DOUBLE PRECISION :: pert(ixImin1:ixImax1), amplitude, xi(ixImin1:ixImax1)

    DOUBLE PRECISION :: dtau(ixImin1:ixImax1), tau(ixImin1:ixImax1), tau_0
    DOUBLE PRECISION :: d_CGS(ixImin1:ixImax1), r_CGS(ixImin1:ixImax1),&
        T_CGS(ixImin1:ixImax1)
    DOUBLE PRECISION :: scale_hieght

    INTEGER :: i
    CHARACTER(100) :: line

    ! len = len_trim(TABLE_name)
    OPEN(123, file=TABLE_name(1:len_trim(TABLE_name)))

    READ(123,'(a)') line
    READ(123,'(a)') line
    READ(123,'(a)') line
    READ(123,'(a)') line
    ! READ(123,'(a)') line

    READ(123,*) Opal_R_list

    READ(123,'(a)') line

    DO i = 1,70
       READ(123,*) Opal_T_list(i), Opal_K_list(:,i)
    ENDDO

    CLOSE(123)

    ! asigne the beta velocity and corresponding dencity
    DO i = ixImin1,ixImax1
       IF (x(i,1) >= R_star_DML) THEN
          ! Set initial velocity acc to beta law
          velocity_field(i) =  v_inf_DML * (1 - 0.999d0*( R_star_DML/x(i,&
             1) + (1 - R_star_DML/x(1,1) )) )**0.5
           ! PRINT*,i, x(i,1), R_star_DML/x(i,1) + (1 - R_star_DML/x(1,1) ), velocity_field(i)*unit_velocity/1.0d5
          ! Set initial density  using the boundary dencuty
          w(i,rho_) = M_dot_DML/(4 * dpi * velocity_field(i) * x(i,1)**2 )

          !set momentum
          w(i,mom(1)) = w(i,rho_)*velocity_field(i)
       ELSE
          w(i,rho_)   = rho_bound_DML
          w(i,mom(1)) = zero
       ENDIF
    END DO
     ! STOP
    ! get the cgs unit dencity and radiuse to compute the optical depth
    r_CGS(ixImin1:ixImax1) = x(ixImin1:ixImax1,1)*unit_length
    d_CGS(ixImin1:ixImax1) = w(ixImin1:ixImax1,rho_)*unit_density

    ! outer boundary optical depth tau_0
    tau_0 = kappa_e_CGS * d_CGS(ixOmax1) * R_star_CGS**2 /(3*r_CGS(ixOmax1))
    ! the differentioal of the optical depth
    dtau(1:ixOmax1) = -d_CGS(1:ixOmax1) * G_init*kappa_e_CGS/Gamma_e  * &
       (R_star_CGS/r_CGS(1:ixOmax1))**2

    ! Do trapez integration to integrate the optical depth
    tau(ixOmax1) = tau_0
    DO i = ixOmax1-1,ixImin1, -1
       tau(i) = half*(dtau(i) + dtau(i+1))*(r_CGS(i) - r_CGS(i+1)) + tau(i+1)
    END DO

    ! call routine to compute the temperature
    CALL Compute_Lucy_temperature(T_CGS,tau,x,ixImin1,ixImax1,ixOmin1,ixOmax1)
    Lower_bound_Cs_DML = dsqrt(const_kB/(mu_molec* &
       mp_cgs)*T_CGS(ixOmin1))/unit_velocity

    ! save variables
    w(ixOmin1:ixOmax1,i_temper) = T_CGS(ixOmin1:ixOmax1)/unit_temperature
    w(ixOmin1:ixOmax1,i_tau_sp) = tau(ixOmin1:ixOmax1)
    w(ixOmin1:ixOmax1,i_kap)    = G_init*kappa_e_CGS/Gamma_e
    w(ixOmin1:ixOmax1,i_g_cak)  = zero
    tmp1(ixOmin1:ixOmax1)   = 0.0d0
    tmp2(ixOmin1:ixOmax1)   = 0.0d0
    tmp3(ixOmin1:ixOmax1)   = 0.0d0
    tmp4(ixOmin1:ixOmax1)   = 0.0d0

    dvav(ixOmin1:ixOmax1)   = 0.0d0
    vav(ixOmin1:ixOmax1)    = 0.0d0
    Tav(ixOmin1:ixOmax1)    = 0.0d0
    Gav(ixOmin1:ixOmax1)    = 0.0d0
    tnorm         = 1.0d0
  END SUBROUTINE initial_conditions
  !===============================================================================






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Subrutine for the bpundaries
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE special_bound(qt,ixGmin1,ixGmax1,ixBmin1,ixBmax1,iB,w,x)

    USE mod_global_parameters
    USE mod_variables
    USE mod_constants
    USE mod_physics, ONLY: phys_get_pthermal

    INTEGER         , INTENT(in)    :: ixGmin1,ixGmax1, ixBmin1,ixBmax1, iB
    DOUBLE PRECISION, INTENT(in)    :: qt, x(ixGmin1:ixGmax1,1:ndim)
    DOUBLE PRECISION, INTENT(inout) :: w(ixGmin1:ixGmax1,1:nw)

    DOUBLE PRECISION :: pth(ixGmin1:ixGmax1)
    DOUBLE PRECISION :: scale_hieght

    INTEGER :: i,j

    ! compute the scale hight of the system under for the lower boundary
    scale_hieght = M_star_DML * G_dp_DML *(1 - Gamma_e) / &
       Lower_bound_Cs_DML**2

    SELECT CASE (iB)

    CASE(1)

       w(ixBmin1:ixBmax1,rho_) = rho_bound_DML !*EXP(scale_hieght/R_star * (1 - R_star/x(i,1)))
       !PRINT*, 'Ghost dencity', w(ixB^S,rho_)

       DO i = ixBmax1,ixBmin1,-1
          w(i,mom(1)) = w(i+1,mom(1))*x(i+1,1)/x(i,1)
       ENDDO

       DO i = ixBmin1,ixBmax1
          w(i,mom(1)) = MIN(w(i,mom(1)),rho_bound_DML*Lower_bound_Cs_DML )
          w(i,mom(1)) = MAX(w(i,mom(1)),-rho_bound_DML*Lower_bound_Cs_DML )
       ENDDO

    CASE default
       CALL mpistop("BC not specified")
    END SELECT
  END SUBROUTINE special_bound
  !=============================================================================






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Subrutine for the total source computation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE tot_Force(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,&
     qt,w,x)
    USE mod_global_parameters

    INTEGER         , INTENT(in)    :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,&
       iwmax
    DOUBLE PRECISION, INTENT(in)    :: qdt, qtC, qt
    DOUBLE PRECISION, INTENT(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    DOUBLE PRECISION, INTENT(inout) :: w(ixImin1:ixImax1,1:nw)

    DOUBLE PRECISION :: g_CAK_CT(ixImin1:ixImax1), g_CAK_new(ixImin1:ixImax1),&
        g_RAD(ixImin1:ixImax1), g_grav_sc(ixImin1:ixImax1)
    DOUBLE PRECISION :: vel(ixImin1:ixImax1)

    INTEGER :: i,j
    INTEGER :: jxmin1,jxmax1, hxmin1,hxmax1

    DOUBLE PRECISION :: dum(ixOmin1:ixOmax1)

    DOUBLE PRECISION :: T_CT_CGS(ixImin1:ixImax1), T_CGS(ixImin1:ixImax1),&
        Cs2(ixImin1:ixImax1), pth(ixImin1:ixImax1)
    DOUBLE PRECISION :: d_CT_CGS(ixImin1:ixImax1),&
        tau_ph_cont_CT(ixImin1:ixImax1), tau_sp_new(ixImin1:ixImax1)
    DOUBLE PRECISION :: kappa_CGS(ixImin1:ixImax1),&
        kappa_cont(ixImin1:ixImax1)
    !-----------------------------------------------------------------------------
    jxmin1=ixOmin1+kr(1,1);jxmax1=ixOmax1+kr(1,1);
    hxmin1=ixOmin1-kr(1,1);hxmax1=ixOmax1-kr(1,1);

    !Collect old variables to compute the opacity and optical depth
    d_CT_CGS(ixOmin1:ixOmax1) = wCT(ixOmin1:ixOmax1,rho_) * unit_density
    g_CAK_CT(ixOmin1:ixOmax1) = wCT(ixOmin1:ixOmax1,i_g_cak)
    T_CT_CGS(ixOmin1:ixOmax1) = wCT(ixOmin1:ixOmax1,i_temper)*unit_temperature

    ! Compute krammer opacity using new density and olt temperature
    CALL kap_opal(T_CT_CGS, d_CT_CGS, kappa_CGS,ixImin1,ixImax1,ixOmin1,&
       ixOmax1)

    ! compute new photospheric optical depth using continuum kappa(t_old,denc_new)
    CALL tau_sp(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,qt,w,x,&
       kappa_CGS,tau_ph_cont_CT)

    ! CAK acceleration
    CALL CAK_Force(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,qt,&
       w,x, g_CAK_new )

    IF(qt < ecut_time) THEN
       IF(qt < 1.7) THEN
          g_CAK_new(ixOmin1:ixOmax1) = g_CAK_new(ixOmin1:ixOmax1)*EXP(-&
             tau_ph_cont_CT(ixOmin1:ixOmax1))
       ELSE
          CALL tau_ph(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,&
             qt,w,x,kappa_CGS,tau_ph_cont_CT)
          g_CAK_new(ixOmin1:ixOmax1) = g_CAK_new(ixOmin1:ixOmax1)*EXP(-&
             tau_ph_cont_CT(ixOmin1:ixOmax1))
       ENDIF
    ENDIF
    w(ixOmin1:ixOmax1,i_g_cak) = g_CAK_new(ixOmin1:ixOmax1)


    ! radiation acceleration
    g_RAD(ixOmin1:ixOmax1) = g_CAK_new(ixOmin1:ixOmax1) + &
       kappa_CGS(ixOmin1:ixOmax1)/unit_opacity_my * L_star_DML /(4 * dpi * &
       x(ixOmin1:ixOmax1,1)**2 * const_c_DML)

    ! kappa_cak = 4*r^2*g_cak/L_star
    kappa_CGS(ixOmin1:ixOmax1) = kappa_CGS(ixOmin1:ixOmax1) + 4 * dpi * &
       x(ixOmin1:ixOmax1,1)**2 * const_c_DML * g_CAK_new(ixOmin1:ixOmax1) / &
       L_star_DML * unit_opacity_my

    !save kappa
    w(ixOmin1:ixOmax1,i_kap) = kappa_CGS(ixOmin1:ixOmax1)

    ! compute new optical depth using full kappa
    CALL tau_sp(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,qt,w,x,&
       kappa_CGS, tau_sp_new)
    w(ixOmin1:ixOmax1,i_tau_sp) = tau_sp_new(ixOmin1:ixOmax1)

    ! compute temperature
    CALL Compute_Lucy_temperature(T_CGS,tau_sp_new,x,ixImin1,ixImax1,ixOmin1,&
       ixOmax1)
    w(ixOmin1:ixOmax1,i_temper) = T_CGS(ixOmin1:ixOmax1)/unit_temperature
    Temperature_test_arrey(ixOmin1:ixOmax1) = T_CGS(ixOmin1:ixOmax1)
    Lower_bound_Cs_DML = dsqrt(const_kB/(mu_molec * &
       mp_cgs)*T_CT_CGS(ixOmin1))/unit_velocity

    ! for testing -- recompute the preassure and save in the output
    Cs2(ixOmin1:ixOmax1) = const_kB/(mu_molec * mp_cgs) * &
       T_CT_CGS(ixOmin1:ixOmax1)
    pth(ixOmin1:ixOmax1) = Cs2(ixOmin1:ixOmax1) * d_CT_CGS(ixOmin1:ixOmax1) / &
       unit_pressure
    w(ixOmin1:ixOmax1,i_P_thrm) = pth(ixOmin1:ixOmax1)

    ! gravitation acceleration
    g_grav_sc(ixOmin1:ixOmax1) = G_dp_DML * M_star_DML /( x(ixOmin1:ixOmax1,&
       1)**2 )

    Gcak_test_arrey(ixOmin1:ixOmax1) = g_CAK_new(ixOmin1:ixOmax1)/g_grav_sc(&
       ixOmin1:ixOmax1)

    !w(ixO^S,i_g_cak) =  g_RAD(ixO^S) - g_grav_sc(ixO^S)
    !w = w + qdt*gsource
    w(ixOmin1:ixOmax1,mom(1)) = w(ixOmin1:ixOmax1,&
       mom(1)) + qdt * (g_RAD(ixOmin1:ixOmax1) - g_grav_sc(ixOmin1:ixOmax1)) * &
       wCT(ixOmin1:ixOmax1,rho_)

    dum = (x(jxmin1:jxmax1,1)-x(ixOmin1:ixOmax1,&
       1))/(g_RAD(ixOmin1:ixOmax1) - g_grav_sc(ixOmin1:ixOmax1))
    new_timestep = 0.3* MINVAL( ABS(dum))**half
  END SUBROUTINE tot_Force
  !=============================================================================







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Subrutine for computation of temperature  OUTPUT in CGS UNITS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Compute_Lucy_temperature(T_out_CGS, tau, x, ixImin1,ixImax1,&
     ixOmin1,ixOmax1)
    USE mod_global_parameters

    INTEGER         , INTENT(in)   :: ixImin1,ixImax1, ixOmin1,ixOmax1
    DOUBLE PRECISION, INTENT(in)   :: tau(ixImin1:ixImax1)
    DOUBLE PRECISION, INTENT(in)   :: x(ixImin1:ixImax1)
    DOUBLE PRECISION, INTENT(out)  :: T_out_CGS(ixImin1:ixImax1)

    DOUBLE PRECISION :: xx(ixImin1:ixImax1)
    DOUBLE PRECISION :: W(ixImin1:ixImax1)

    xx(ixImin1:ixImax1) = R_star_DML/x(ixImin1:ixImax1)
    ! dilution factor for lucy 1971
    W(ixOmin1:ixOmax1) = half * (1 - (1 - xx(ixOmin1:ixOmax1)**2 )**0.5  )
    W(ixImin1:ixOmin1) = 0.5d0

    ! lucy temperature structure
    T_out_CGS(1:ixOmax1) = T_star_CGS * ( (W(1:ixOmax1) + 3.0 * &
       tau(1:ixOmax1)/4.0) )**0.25

    ! set the minimal cut-off at 10000 kK
    T_out_CGS(ixOmin1:ixOmax1) =   MAX(T_out_CGS(ixOmin1:ixOmax1), 1.0d4)
  END SUBROUTINE Compute_Lucy_temperature
  !=============================================================================




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Subrutine for computation of OPAL OPACITIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE kap_opal(T_CGS, d_CGS, kappa_CGS,ixImin1,ixImax1,ixOmin1,ixOmax1)
    USE mod_global_parameters

    INTEGER         , INTENT(in)  :: ixImin1,ixImax1, ixOmin1,ixOmax1
    DOUBLE PRECISION, INTENT(out) :: kappa_CGS(ixImin1:ixImax1)
    DOUBLE PRECISION, INTENT(in)  :: T_CGS(ixImin1:ixImax1)
    DOUBLE PRECISION, INTENT(in)  :: d_CGS(ixImin1:ixImax1)

    DOUBLE PRECISION :: K_neigb(2,2)
    DOUBLE PRECISION :: R_neigb(2)
    DOUBLE PRECISION :: T_neigb(2)
    DOUBLE PRECISION :: K_inter

    INTEGER :: i , nR(2), nT(2)

    DOUBLE PRECISION :: LgR,LgT



    DO i = ixOmin1,ixOmax1
       LgT = LOG10(T_CGS(i))
       LgR = LOG10(d_CGS(i)) - 3d0*LgT + 18

       CALL Find_neighbours(Opal_R_list,41,LgR,nR)!Opal_R_list(2)
       CALL Find_neighbours(Opal_T_list,75,LgT,nT)!Opal_T_list(2)
       !
       ! PRINT*, Opal_K_list(nR(1),nT(2)),Opal_K_list(nR(2),nT(2))
       ! PRINT*, Opal_K_list(nR(1),nT(1)),Opal_K_list(nR(2),nT(1))

       K_neigb(1,1) = Opal_K_list(nR(1),nT(1))
       K_neigb(1,2) = Opal_K_list(nR(1),nT(2))
       K_neigb(2,1) = Opal_K_list(nR(2),nT(1))
       K_neigb(2,2) = Opal_K_list(nR(2),nT(2))
       T_neigb(1)  = Opal_T_list(nT(1))
       T_neigb(2)  = Opal_T_list(nT(2))
       R_neigb(1)  = Opal_R_list(nR(1))
       R_neigb(2)  = Opal_R_list(nR(2))
       ! PRINT*, "X(1,2)" , TT_Qxqq(1),TT_Qxqq(2)
       ! PRINT*, "Y(1,2)" , TT_Qyqq(1),TT_Qyqq(2)

       CALL D2_interpol(K_neigb,R_neigb,T_neigb,K_inter,LgR,LgT) !,Opal_R_list(2),Opal_T_list(2),-6.7d0,3.015d0
       !
       ! PRINT*,TT_Fqq
       ! STOP
       kappa_CGS(i) = 10d0**K_inter
    END DO
  END SUBROUTINE kap_opal
  !=============================================================================






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Subrutine for computation of Shperically modified optical depth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE tau_sp(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,qt,&
     w,x, kappa_CGS, tau_out)
    USE mod_global_parameters

    INTEGER         , INTENT(in)    :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,&
       iwmax
    DOUBLE PRECISION, INTENT(in)    :: qdt, qtC, qt
    DOUBLE PRECISION, INTENT(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    DOUBLE PRECISION, INTENT(inout) :: w(ixImin1:ixImax1,1:nw)
    DOUBLE PRECISION, INTENT(in)    :: kappa_CGS(ixImin1:ixImax1)
    DOUBLE PRECISION, INTENT(out)   :: tau_out(ixImin1:ixImax1)

    DOUBLE PRECISION :: kappa(ixImin1:ixImax1)
    DOUBLE PRECISION :: dtau(ixImin1:ixImax1),tau(ixImin1:ixImax1),tau_0
    DOUBLE PRECISION :: d_CGS(ixImin1:ixImax1),r_CGS(ixImin1:ixImax1)

    INTEGER :: i

    r_CGS(ixImin1:ixImax1) = x(ixImin1:ixImax1,1)*unit_length
    d_CGS(ixImin1:ixImax1) = wCT(ixImin1:ixImax1,rho_)*unit_density

    tau_0 = kappa_e_CGS * d_CGS(ixOmax1) * R_star_CGS**2 /(3 * r_CGS(ixOmax1))
    ! dtau(ixO^S) = -d_CGS(ixO^S) * kappa_CGS(ixO^S)*(R_star_CGS/r_CGS(ixO^S))**2
    dtau(1:ixOmax1) = -d_CGS(1:ixOmax1) * &
       kappa_CGS(1:ixOmax1)*(R_star_CGS/r_CGS(1:ixOmax1))**2

    tau(ixOmax1) = tau_0
    DO i = ixOmax1-1,ixImin1, -1
       tau(i) = half*(dtau(i) + dtau(i+1))*(r_CGS(i) - r_CGS(i+1)) + tau(i+1)
    END DO

    tau_out(ixOmin1:ixOmax1) = tau(ixOmin1:ixOmax1) !w(ixOmin1:ixOmax1,i_tau_sp) = tau(ixOmin1:ixOmax1)
  END SUBROUTINE tau_sp
  !=============================================================================



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Subrutine for computation of Shperically modified optical depth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE tau_ph(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,qt,&
     w,x, kappa_CGS, tau_out)
    USE mod_global_parameters

    INTEGER         , INTENT(in)    :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,&
       iwmax
    DOUBLE PRECISION, INTENT(in)    :: qdt, qtC, qt
    DOUBLE PRECISION, INTENT(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    DOUBLE PRECISION, INTENT(inout) :: w(ixImin1:ixImax1,1:nw)
    DOUBLE PRECISION, INTENT(in)    :: kappa_CGS(ixImin1:ixImax1)
    DOUBLE PRECISION, INTENT(out)   :: tau_out(ixImin1:ixImax1)

    DOUBLE PRECISION :: kappa(ixImin1:ixImax1)
    DOUBLE PRECISION :: dtau(ixImin1:ixImax1),tau(ixImin1:ixImax1),tau_0
    DOUBLE PRECISION :: d_CGS(ixImin1:ixImax1),r_CGS(ixImin1:ixImax1)

    INTEGER :: i

    r_CGS(ixImin1:ixImax1) = x(ixImin1:ixImax1,1)*unit_length
    d_CGS(ixImin1:ixImax1) = wCT(ixImin1:ixImax1,rho_)*unit_density

    tau_0 = kappa_e_CGS * d_CGS(ixOmax1) * r_CGS(ixOmax1)

    dtau(1:ixOmax1) = -d_CGS(1:ixOmax1) * kappa_CGS(1:ixOmax1)

    tau(ixOmax1) = tau_0
    DO i = ixOmax1-1,ixImin1, -1
       tau(i) = half*(dtau(i) + dtau(i+1))*(r_CGS(i) - r_CGS(i+1)) + tau(i+1)
    END DO

    tau_out(ixOmin1:ixOmax1) = tau(ixOmin1:ixOmax1)
  END SUBROUTINE tau_ph
  !=============================================================================




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Subrutine for computation of CAK force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CAK_Force(qdt, ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,iwmax, qtC,&
      wCT, qt, w, x, g_CAK)
    USE mod_global_parameters

    INTEGER         , INTENT(in)    :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,&
       iwmax
    DOUBLE PRECISION, INTENT(in)    :: qdt, qtC, qt
    DOUBLE PRECISION, INTENT(in)    :: wCT(ixImin1:ixImax1,1:nw),&
        x(ixImin1:ixImax1,1:ndim)
    DOUBLE PRECISION, INTENT(inout) :: w(ixImin1:ixImax1,1:nw)
    DOUBLE PRECISION, INTENT(out)   :: g_CAK(ixImin1:ixImax1)

    DOUBLE PRECISION :: vel(ixImin1:ixImax1)
    DOUBLE PRECISION :: grad_fwd(ixImin1:ixImax1),grad_bwd(ixImin1:ixImax1),&
       grad(ixImin1:ixImax1)

    DOUBLE PRECISION :: xx(ixImin1:ixImax1)
    DOUBLE PRECISION :: a_mod_cak(ixImin1:ixImax1)


    INTEGER :: i,j
    INTEGER :: jxmin1,jxmax1, hxmin1,hxmax1

    xx(ixOmin1:ixOmax1) =  1 - R_star_DML/x(ixOmin1:ixOmax1,1)
    !-------------------------------------------------------------------------
    jxmin1=ixOmin1+kr(1,1);jxmax1=ixOmax1+kr(1,1);
    hxmin1=ixOmin1-kr(1,1);hxmax1=ixOmax1-kr(1,1);
    !-------------------------------------------------------------------------

    vel(ixImin1:ixImax1) = wCT(ixImin1:ixImax1,mom(1))/wCT(ixImin1:ixImax1,&
       rho_)
    !-------------------------------------------------------------------------
    ! calculate grad_v
    ! forward difference
    grad_fwd(ixOmin1:ixOmax1) = (vel(jxmin1:jxmax1)-&
       vel(ixOmin1:ixOmax1))/(x(jxmin1:jxmax1,1)-x(ixOmin1:ixOmax1,1))
    ! backward difference
    grad_bwd(ixOmin1:ixOmax1) = (vel(ixOmin1:ixOmax1)-&
       vel(hxmin1:hxmax1))/(x(ixOmin1:ixOmax1,1)-x(hxmin1:hxmax1,1))
    ! central difference
    grad(ixOmin1:ixOmax1) = half*ABS(grad_fwd(ixOmin1:ixOmax1)+&
       grad_bwd(ixOmin1:ixOmax1))

    !-------------------------------------------------------------------------
    !calculate the modified a_cak as a function of xx = 1 - R/r
    DO i = ixOmin1,ixOmax1
       IF(xx(i) <= a_cak_x0) THEN
          a_mod_cak(i) = a_cak_max
       ELSE IF(xx(i) < a_cak_x1) THEN
          a_mod_cak(i) = a_cak_max + (a_cak_min - a_cak_max) * (a_cak_x0 - &
             xx(i)) / (a_cak_x0 - a_cak_x1)
       ELSE
          a_mod_cak(i) = a_cak_min
       ENDIF
    END DO

    ! calculate g_CAK
    g_CAK(ixOmin1:ixOmax1) = 1./(1.-a_mod_cak(ixOmin1:ixOmax1)) * kappa_e_DML &
       * L_star_DML * Q_cak / (4. * dpi * x(ixOmin1:ixOmax1,&
       1)**2 * const_c_DML) *(grad(ixOmin1:ixOmax1) / (wCT(ixOmin1:ixOmax1,&
       rho_) * const_c_DML * Q_cak * kappa_e_DML))**a_mod_cak(ixOmin1:ixOmax1)

    w(ixOmin1:ixOmax1,i_a_cak) = a_mod_cak(ixOmin1:ixOmax1)
  END SUBROUTINE CAK_Force
  !=============================================================================



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Subrutine for computation of special time step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE special_dt(w,ixImin1,ixImax1,ixOmin1,ixOmax1,dtnew,dx1,x)
    USE mod_global_parameters

    INTEGER         , INTENT(in)    :: ixImin1,ixImax1, ixOmin1,ixOmax1
    DOUBLE PRECISION, INTENT(in)    :: dx1, x(ixImin1:ixImax1,1:ndim)
    DOUBLE PRECISION, INTENT(in)    :: w(ixImin1:ixImax1,1:nw)
    DOUBLE PRECISION, INTENT(inout) :: dtnew

    IF (it .GE. 1) dtnew = new_timestep
    dtnew = MIN(dtnew,3.0d-5)

  END SUBROUTINE special_dt
  !=============================================================================

SUBROUTINE manual_csound2(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,cs2)
    USE mod_global_parameters

    INTEGER         , INTENT(in)    :: ixImin1,ixImax1, ixOmin1,ixOmax1
    DOUBLE PRECISION, INTENT(in)    :: x(ixImin1:ixImax1,1:ndim)
    DOUBLE PRECISION, INTENT(in)    :: w(ixImin1:ixImax1,1:nw)
    DOUBLE PRECISION, INTENT(out)   :: cs2(ixImin1:ixImax1)

    cs2(ixOmin1:ixOmax1) = const_kB/(mu_molec * mp_cgs) * w(ixOmin1:ixOmax1,&
       i_temper) * unit_temperature / unit_velocity**2

  END SUBROUTINE manual_csound2
  !=============================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Subrutine for computation of special thermal pressure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE p_thermal(w,x,ixImin1,ixImax1,ixOmin1,ixOmax1,pth)
    USE mod_global_parameters

    INTEGER         , INTENT(in)    :: ixImin1,ixImax1, ixOmin1,ixOmax1
    DOUBLE PRECISION, INTENT(in)    :: x(ixImin1:ixImax1,1:ndim)
    DOUBLE PRECISION, INTENT(in)    :: w(ixImin1:ixImax1,1:nw)
    DOUBLE PRECISION, INTENT(out)   :: pth(ixImin1:ixImax1)

    DOUBLE PRECISION :: Cs2(ixImin1:ixImax1), d(ixImin1:ixImax1)

    d(ixImin1:ixImax1) = w(ixImin1:ixImax1,rho_)*unit_density

    Cs2(ixOmin1:ixOmax1) = const_kB/(mu_molec * mp_cgs) * w(ixOmin1:ixOmax1,&
       i_temper) * unit_temperature

    Cs2(ixOmax1+1:ixImax1) = Cs2(ixOmax1)
    Cs2(ixImin1:ixOmin1-1) = Cs2(ixOmin1)

    pth(ixImin1:ixImax1) = Cs2(ixImin1:ixImax1) * d(ixImin1:ixImax1) / &
       unit_pressure

    ! pth(ixI^S) = zero
  END SUBROUTINE p_thermal
  !=============================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Subrutine for bilinear interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE D2_interpol(Q,Qx,Qy,F,Fx,Fy)
    DOUBLE PRECISION, INTENT(in)    :: Q(2,2)
    DOUBLE PRECISION, INTENT(in)    :: Qx(2)
    DOUBLE PRECISION, INTENT(in)    :: Qy(2)
    DOUBLE PRECISION, INTENT(out)   :: F
    DOUBLE PRECISION, INTENT(in)    :: Fx
    DOUBLE PRECISION, INTENT(in)    :: Fy

    DOUBLE PRECISION :: a0,a1,Qx1,Qx2

    INTEGER :: INTERP_METHOD

    ! PRINT*,Qx(1) - Qx(2) , Qy(1) - Qy(2)
    IF((Qx(1) == Qx(2)) .AND. (Qy(1) == Qy(2))) THEN
       INTERP_METHOD = 3
    ELSE IF(Qx(1) == Qx(2)) THEN
       INTERP_METHOD = 2
    ELSE IF(Qy(1) == Qy(2)) THEN
       INTERP_METHOD = 1
    ELSE
       INTERP_METHOD = 0
    ENDIF

    SELECT CASE (INTERP_METHOD)
    CASE(3)
       ! PRINT*,"CONST"
       F = Q(1,1)
    CASE(2)
       ! PRINT*,"CONST X"
       F = Q(1,1) + (Q(1,2) - Q(1,1))/(Qy(2) - Qy(1))*(Fy - Qy(1))
    CASE(1)
       ! PRINT*,"CONST Y"
       F = Q(1,1) + (Q(2,1) - Q(1,1))/(Qx(2) - Qx(1))*(Fx - Qx(1))
    CASE DEFAULT
       ! PRINT*,"DEAFAUL"
       Qx1 = Q(1,1) + (Q(2,1) - Q(1,1))/(Qx(2) - Qx(1))*(Fx - Qx(1))
       Qx2 = Q(1,2) + (Q(2,2) - Q(1,2))/(Qx(2) - Qx(1))*(Fx - Qx(1))
       ! PRINT*,"Qx1, Qx2", Qx1,Qx2
       F = Qx1 + (Qx2 - Qx1)/(Qy(2) - Qy(1))*(Fy - Qy(1))
    END SELECT


  END SUBROUTINE D2_interpol
  !=============================================================================



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Subrutine for bilinear interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Find_neighbours(x_list,x_list_len,x,nx)
    INTEGER         , INTENT(in)    :: x_list_len
    DOUBLE PRECISION, INTENT(in)    :: x_list(x_list_len)
    DOUBLE PRECISION, INTENT(in)    :: x
    INTEGER         , INTENT(out)   :: nx(2)

    DOUBLE PRECISION :: Diff(x_list_len)
    DOUBLE PRECISION :: Diff_up_neig_val,Diff_lo_neig_val

    Diff = x_list - x

    Diff_up_neig_val = MINVAL(Diff, MASK = Diff .GE. 0)
    Diff_lo_neig_val = MAXVAL(Diff, MASK = Diff .LE. 0)

    nx(1) = MINLOC(ABS(Diff-Diff_lo_neig_val),1)
    nx(2) = MINLOC(ABS(Diff-Diff_up_neig_val),1)

  END SUBROUTINE Find_neighbours
  !=============================================================================


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  Subrutine for computation of avarege values on d, v, T,
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Stat_analysis()
    !
    ! Compute and store hydro quantities for statistical quantities
    !
    USE mod_global_parameters

    INTEGER :: iigrid, igrid, i
    LOGICAL :: first_time = .TRUE.
    ! DOUBLE PRECISION :: d(ixG^T) v(ixG^T)

    ! if carent time is less then start time of avaraging
    IF (global_time <= avst_time_DML) THEN
       ! sat the avarege values to 0
       dvav(ixGlo1:ixGhi1)  = 0.0d0
       vav(ixGlo1:ixGhi1)   = 0.0d0
       Tav(ixGlo1:ixGhi1)   = 0.0d0
       Gav(ixGlo1:ixGhi1)   = 0.0d0
       tnorm        = 1.0d0
       CONTINUE
    ELSE
       ! PRINT*, "starting stat"
       ! Switch all grids (with ghost cells) to primitive variables
       DO iigrid = 1,igridstail; igrid = igrids(iigrid)
          CALL hd_to_primitive(ixGlo1,ixGhi1,ixGlo1,ixGhi1,ps(igrid)%w,&
             ps(igrid)%x)
          ! d(ixG^T) = ps(igrid)%w(ixG^T,rho_)
          ! v(ixG^T) = ps(igrid)%w(ixG^T,mom(1))/d(ixG^T)
          ! PRINT*, igrid, d(igrid)
       END DO

       xc(ixGlo1:ixGhi1) = ps(igrid)%x(ixGlo1:ixGhi1,1)

       IF (first_time) THEN
          !print*,'accessing first time statistical routine, initializing'
          tmp1(ixGlo1:ixGhi1) = ps(igrid)%w(ixGlo1:ixGhi1,&
             rho_) * unit_density * ps(igrid)%w(ixGlo1:ixGhi1,&
             mom(1)) * unit_velocity
          tmp2(ixGlo1:ixGhi1) = ps(igrid)%w(ixGlo1:ixGhi1,&
             mom(1)) * unit_velocity
          tmp3(ixGlo1:ixGhi1) = Temperature_test_arrey(ixGlo1:ixGhi1)
          tmp4(ixGlo1:ixGhi1) = 1.0d0

          first_time = .FALSE.
       ELSE
          tnorm = global_time - avst_time_DML

          dt = dt*unit_time
          ! momentum_avarege*t  = momentum_avarege*t + 0.5(momentum_new + momentum_old)*dt
          dvav(ixGlo1:ixGhi1)  = dvav(ixGlo1:ixGhi1) + &
             0.5d0*(ps(igrid)%w(ixGlo1:ixGhi1,&
             rho_)*unit_density *ps(igrid)%w(ixGlo1:ixGhi1,&
             mom(1))*unit_velocity + tmp1(ixGlo1:ixGhi1))*dt
          tmp1(ixGlo1:ixGhi1) = ps(igrid)%w(ixGlo1:ixGhi1,&
             rho_)*unit_density*ps(igrid)%w(ixGlo1:ixGhi1,&
             mom(1))*unit_velocity !momentum_old

          ! velocity_avarege*t  = velocity_avarege*t + 0.5(velocity_new + velocity_old)*dt
          vav(ixGlo1:ixGhi1)  = vav(ixGlo1:ixGhi1) + &
             0.5d0*(ps(igrid)%w(ixGlo1:ixGhi1,&
             mom(1))* unit_velocity + tmp2(ixGlo1:ixGhi1))*dt
          tmp2(ixGlo1:ixGhi1) = ps(igrid)%w(ixGlo1:ixGhi1,&
             mom(1)) * unit_velocity !velocity_old

          ! temperature_avarege*t  = temperature_avarege*t + 0.5(temperature_new + temperature_old)*dt
          Tav(ixGlo1:ixGhi1)  = Tav(ixGlo1:ixGhi1) + 0.5d0*( &
             Temperature_test_arrey(ixGlo1:ixGhi1) + tmp3(ixGlo1:ixGhi1))*dt
          tmp3(ixGlo1:ixGhi1) = Temperature_test_arrey(ixGlo1:ixGhi1) !temperature_old

          Gav(ixGlo1:ixGhi1)  = Gav(ixGlo1:ixGhi1) + 0.5d0*( &
             Gcak_test_arrey(ixGlo1:ixGhi1) + tmp4(ixGlo1:ixGhi1))*dt
          tmp4(ixGlo1:ixGhi1) = Gcak_test_arrey(ixGlo1:ixGhi1)

          dt = dt/unit_time

       ENDIF

       ! Switching all back to conservative variables
       DO iigrid=1,igridstail; igrid=igrids(iigrid)
          CALL hd_to_conserved(ixGlo1,ixGhi1,ixGlo1,ixGhi1,ps(igrid)%w,&
             ps(igrid)%x)
       END DO

    ENDIF

  END SUBROUTINE Stat_analysis

  !===============================================================================

  SUBROUTINE myextravar_output(ixImin1,ixImax1,ixOmin1,ixOmax1,w,x,normconv)
    !
    ! Stores statistical variables of interest in auxiliary slots of w-array
    !
    USE mod_global_parameters

    INTEGER, INTENT(in) :: ixImin1,ixImax1,ixOmin1,ixOmax1
    REAL(8), INTENT(in) :: x(ixImin1:ixImax1,1:ndim)
    REAL(8)             :: w(ixImin1:ixImax1,nw+nwauxio)
    REAL(8)             :: normconv(0:nw+nwauxio)

    ! PRINT*, SHAPE(w)
    ! PRINT*, nw
    ! stop
    ! Average density
    w(ixOmin1:ixOmax1,nw+1) = dvav(ixOmin1:ixOmax1) /(tnorm*unit_time)
    ! PRINT*, 1
    ! Average velocity
    w(ixOmin1:ixOmax1,nw+2) = vav(ixOmin1:ixOmax1) /(tnorm*unit_time)
    ! PRINT*, 2
    ! Measure of speed of clumps
    w(ixOmin1:ixOmax1,nw+3) = Tav(ixOmin1:ixOmax1) /(tnorm*unit_time)
    ! PRINT*, 3
    ! stop
    ! PRINT*, SHAPE(w)
    w(ixOmin1:ixOmax1,nw+4) = Gav(ixOmin1:ixOmax1) /(tnorm*unit_time)

  END SUBROUTINE myextravar_output

  !===============================================================================

  SUBROUTINE myextravarnames_output(varnames)
    !
    ! Additional variables to be stored in .dat file
    !
    CHARACTER(len=*) :: varnames
    varnames = '<rho*v> <v> <T> <G_cak>'

  END SUBROUTINE myextravarnames_output

END MODULE mod_usr
! a0 =  Q(1,1)*Qx(2)*Qy(2)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2))) &
!      -Q(1,2)*Qx(2)*Qy(1)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2))) &
!      -Q(2,1)*Qx(1)*Qy(2)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2))) &
!      +Q(2,2)*Qx(1)*Qy(1)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2)))
!
! a1 = -Q(1,1)*Qy(2)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2))) &
!      +Q(1,2)*Qy(1)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2))) &
!      +Q(2,1)*Qy(2)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2))) &
!      -Q(2,2)*Qy(1)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2)))
!
! a2 = -Q(1,1)*Qx(2)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2))) &
!      +Q(1,2)*Qx(2)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2))) &
!      +Q(2,1)*Qx(1)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2))) &
!      -Q(2,2)*Qx(1)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2)))
!
! a2 =  Q(1,1)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2))) &
!      -Q(1,2)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2))) &
!      -Q(2,1)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2))) &
!      +Q(2,2)/((Qx(1)-Qx(2))*(Qy(1)-Qy(2)))
!
! F = a0 + a1*Fx + a2*Fy + a3*Fx*Fy