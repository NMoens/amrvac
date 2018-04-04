module mod_srhd_phys

  implicit none
  private

  !-------------------------Comments-----------------------------------------------
  !*DM* Normally there are 2 auxiliary variables in the srhd module
  !*DM* the Lorentz factor and the pressure
  !*DM* These are used in the calculation of the proper density (d), 
  !*DM* the momentum density (rmom, previously s) and the total energy density (tau)
  !--------------------------------------------------------------------------------
  !*DM* In the previous version, srhd and srhdeos used an ideal or a Mathews EOS
  !*DM* respectively. Now we add a switch 
  !-------------------------------------------------------------------------------- 

  !> The adiabatic index
  double precision, public :: srhd_gamma = 5.d0/3.0d0
  
  !> The equation of state (options: ideal, mathews)
  character(len=std_len) :: srhd_eos = "ideal"

  integer            :: eos_type    = 1
  integer, parameter :: eos_ideal   = 1
! *DM* do the ideal case first
!  integer, parameter :: eos_mathews = 2 

!  !> Index of the energy density
!  integer, public, protected              :: e_

! *DM* Primitive variables (following the old version) 

  !> Index of the density (lab frame)
  integer, public, protected              :: d_

!  !> Index of the  density (primitive?) (*DM* we delete this one I guess...)
!  integer, public, protected              :: s0_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: rmom(:)

  !> Indices of the four velocity
  integer, allocatable, public, protected :: uvel(:)

  !> Indices of the three-velocity
  integer, allocatable, public, protected :: vvel(:)

  !> Index of the total energy density 
  integer, public, protected              :: tau_

  !> Number of tracer species
  integer, public, protected              :: srhd_n_tracer = 0

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !Auxiliary variables
  !*DM* we have to make sure that these are available everywhere needed!
  !> Index of Lorentz factor
  integer, public, protected :: lfac_
  
  !> Index of pressure
  integer, public, protected :: p_

  integer, parameter :: nvector=1                             ! No. vector vars
  integer, dimension(nvector), parameter :: iw_vector=(/ s0_ /)

  integer, parameter:: gamma_=1,neqpar=1                     ! equation parameters

  double precision :: smalltau,smallxi,minrho,minp

 !Public methods
  public :: srhd_phys_init
  public :: srhd_to_conserved
  public :: srhd_to_primitive

contains

  !> Read this module's parameters from a file
  subroutine srhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /srhd_list/ srhd_eos, srhd_n_tracer, srhd_gamma

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, hd_list, end=111)
111    close(unitpar)
    end do

    select case (srhd_eos)
    case ("ideal")
       eos_type = eos_ideal
    case ("mathews")
       eos_type = eos_mathews
    case default
       call mpistop("Invalid srhd_eos, available eos are Ideal and Mathews")
    end select

  end subroutine srhd_read_params
!=============================================================================
!*DM* 
!add the srhd write info subroutine ?
! Can do this later
!*DM* 
!---------------------------------------------------------------------
!> Initialize the module
  subroutine srhd_phys_init()
    use mod_global_parameters
    use mod_physics

    integer :: itr, idir

    call srhd_read_params(par_files)

    physics_type = "srhd"
!    phys_energy  = srhd_energy !*DM* Not necessary now, we choose between two eos !

    ! Determine flux variables
    rho_ = var_set_rho()
    d_ = rho_

    allocate(rmom(ndir))
    rmom(:) = var_set_momentum(ndir)
    !*DM* I guess we also need those now...
    allocate(uvel(ndir))
    uvel(:) = var_set_fourvelocity(ndir)
    allocate(vvel(ndir))
    vvel(:) = var_set_threevelocity(ndir)

    if (eos_type == eos_mathews) then
       ! Set index of energy variable
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if

    lfac_ = var_set_extravar("lfac", "lfac")
    p_ = var_set_extravar("pressure", "pressure")

    phys_get_dt              => srhd_get_dt
    phys_get_cmax            => srhd_get_cmax
    phys_get_cbounds         => srhd_get_cbounds
    phys_get_flux            => srhd_get_flux
    phys_add_source_geom     => srhd_add_source_geom
    phys_add_source          => srhd_add_source
    phys_to_conserved        => srhd_to_conserved
    phys_to_primitive        => srhd_to_primitive
    phys_check_params        => srhd_check_params
    phys_check_w             => srhd_check_w
    phys_write_info          => srhd_write_info
    phys_handle_small_values => srhd_handle_small_values ! Can fix later

    ! Whether diagonal ghost cells are required for the physics (TODO: true?)
    phys_req_diagonal = .false.

!*DM* derive units from basic units (TO BE ADDED FOR SRHD) ! Can fix later
!    call srhd_physical_units()

    allocate(tracer(srhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, srhd_n_tracer
       tracer(itr) = var_set_fluxvar("tr", "tr", itr, need_bc=.false.)
    end do

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector)) !*DM* Are we keeping iw in the end ??
    iw_vector(1) = rmom(1) - 1

    if (eos_type == eos_ideal) then
      minp    = max(zero,smallp)
      minrho  = max(zero,smallrho)
      smalltau= minp/(srhd_gamma-1)
      smallxi = minrho + minp*srhd_gamma/(srhd_gamma-1)
    else
      minp   = max(zero,smallp)
      minrho = max(zero,smallrho)
      ! call smallvaluesEOS() ! Check later
   end if

  end subroutine srhd_phys_init
!==============================================================================
  subroutine srhd_physical_units
    use mod_global_parameters
    double precision :: mp
!    ! Derive scaling units
!*DM* I used the proton mass here again      
    if(SI_unit) then
      mp=mp_SI
    else
      mp=mp_cgs
    end if

      unit_density=mp*unit_numberdensity
      unit_pressure=unit_density*unit_velocity**2
!      unit_temperature=unit_pressure/((2.d0+3.d0*He_abundance)*unit_numberdensity*kB)
      unit_time=unit_length/unit_velocity

  end subroutine srhd_physical_units
!==============================================================================
! Fix later, probably have to switch flag (see mod_hd_phys)
  subroutine srhd_check_w(checkprimitive,ixI^L,ixO^L,w,flag)
    use mod_global_parameters

    logical, intent(in)          :: checkprimitive
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,nw)
    logical, intent(out)         :: flag(ixG^T)
!    double precision             :: tmp(ixI^S)
    !-----------------------------------------------------------------------------

    flag(iO^S)= 0

! *DM* Deleted the useprimitiveRel option, always true

    if (checkprimitive) then
       ! check   rho>=0, p>=smallp
       flag(ixO^S) = (w(ixO^S,rho_) >= minrho).and. &
                     (w(ixO^S,pp_)  >= minp)
    else
       ! Check D>=0 and lower limit for tau
       flag(ixO^S) = (w(ixO^S,d_)    >= minrho).and. &
            (w(ixO^S,tau_)  >= smalltau)
    endif

  end subroutine srhd_check_w
!=============================================================================
  subroutine srhd_to_conserved(ixI^L,ixO^L,w,x) 

    !*DM*
    ! Transform primitive variables into conservative ones
    ! (rho,v,p) ---> (D,S,tau) **THIS IS THE OLD VERSION**
    ! (rho,v,p) ---> (D,rmom, tau) **THIS SHOULD BE THE NEW TRANSFORMATION**

    !*DM*
    use mod_global_parameters

    integer, intent(in)               :: ixI^L, ixO^L
    double precision, intent(inout)   :: w(ixI^S,nw)
    integer                         :: idir, itr

    double precision,dimension(ixG^T) :: xi

  !-----------------------------------------------------------------------------
  
  !*DM* useprimitiveRel is always true now, I deleted everything else
     ! This assumes four velocity computed in primitive (rho u p) with u=lfac*v

  ! assumes four velocity computed in primitive (rho u p) with u=lfac*v
     xi(ixO^S)=1.0d0+sum(w(ixO^S,uvel(:))**2, dim=ndim+1)
      ! fill the auxiliary variable lfac_ (Lorentz Factor) and p_ (pressure)
     w(ixO^S,lfac_)=dsqrt(xi(ixO^S))
     w(ixO^S,p_)=w(ixO^S,pp_)

   ! compute xi=Lfac w  (enthalphy w)
      xi(ixO^S)=w(ixO^S,lfac_)*(w(ixO^S,rho_)+w(ixO^S,p_)*eqpar(gamma_)/(eqpar(gamma_)-one))

    w(ixO^S,d_)=w(ixO^S,rho_)*w(ixO^S,lfac_)
    do idir=1, ndir
       w(ixO^S,rmom(idir))=xi(ixO^S)*w(ixO^S,uvel(idir))
    end do
    w(ixO^S,tau_)=xi(ixO^S)*w(ixO^S,lfac_) - w(ixO^S,p_) - w(ixO^S,d_)

! *DM* add this later
!   if(check_small_values) call srhd_handle_smallvalues(.false.,w,x,ixI^L,ixO^L,"srhd_to_conserved")

  end subroutine srhd_to_conserved
  !=============================================================================
  subroutine srhd_to_primitive(ixI^L,ixO^L,w,x)
    !*DM*
    ! Transform conservative variables into primitive ones
    ! (D,S,tau) ---> (rho,v,p) *Notation in the old MPI-AMRVAC version*
    ! (D,rmom,tau) ---> (rho,v,p) *Notation in the new version*

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    integer                         :: itr, idir
    !-----------------------------------------------------------------------------

    ! Calculate pressure and Lorentz factor from conservative variables only
    ! these are put in lfac_ and p_ auxiliaries

  do idir=1, ndir
     w(ixO^S,uvel(idir))=w(ixO^S,lfac_)*w(ixO^S,rmom(idir))/(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_));
  end do

  ! compute density
  w(ixO^S,rho_)=w(ixO^S,d_)/w(ixO^S,lfac_)
  ! fill pressure
  w(ixO^S,pp_)=w(ixO^S,p_)

  end subroutine srhd_to_primitive(ixI^L,ixO^L,w,x) 
!=============================================================================
  subroutine srhd_get_v(w,x,ixI^L,ixO^L,idim,v)

!*DM* THIS SHOULD BE THE SAME AS IN THE OLD VERSION....
    ! Calculate v_idim=m_idim/rho within ixO^L

    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L, idim
    double precision, intent(in)  :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out) :: v(ixI^S)
    !-----------------------------------------------------------------------------
    ! assuming getv FOLLOWS a getaux call for updated p_ entry and
    ! well-behaved xi=d+tau+p

    v(ixO^S) = w(ixO^S,rmom(idim))/(w(ixO^S,d_)+w(ixO^S,tau_)+w(ixO^S,p_))

  end subroutine srhd_get_v
!=============================================================================
  ! Look at mod_hd_phys
  subroutine srhd_get_cmax(w,x,ixI^L,ixO^L,idim,cmax)

!*DM* THIS IS DIFFERENT FOR SRHD/SRHDEOS
!*DM*  INCLUDE THE IDEAL CASE
    ! Calculate cmax_idim=csound+abs(v_idim) within ixO^L

    use mod_global_parameters

rhoh(ixO^S) = w(ixO^S,d_)/w(ixO^S,lfac_)+eqpar(gamma_)*w(ixO^S,p_)/(eqpar(gamma_)-one)
csound2(ixO^S)=eqpar(gamma_)*w(ixO^S,p_)/rhoh(ixO^S)
if(.not.needcmin)then
! xi(ixO^S)=1.0d0+sum(w(ixO^S,uvel(:))**2, dim=ndim+1)
   v2(ixO^S)= sum(w(ixO^S,rmom(:))**2, dim=ndim+1)/((rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_))**2.0d0)
else
   vidim(ixO^S)=(w(ixO^S,s0_+idim)/(rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_)))
endif

if(.not.needcmin)then
  if (ndir==1) then
     cmax(ixO^S)= (dsqrt(v2(ixO^S))+dsqrt(csound2(ixO^S)))/(1.0d0+dsqrt(csound2(ixO^S)*v2(ixO^S)))
  else
     vidim2(ixO^S)=(w(ixO^S,s0_+idim)/rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_)))**2.0d0

     cmax(ixO^S)=( dsqrt(vidim2(ixO^S))*(1.0d0-csound2(ixO^S)) + &
                dsqrt(csound2(ixO^S)*(1.0d0-v2(ixO^S))*( &
        1.0d0-v2(ixO^S)*csound2(ixO^S)-vidim2(ixO^S)*(1.0d0-csound2(ixO^S))) &
               ) ) / (1.0d0-v2(ixO^S)*csound2(ixO^S))
  endif
else
  if (ndir==1) then
     cmax(ixO^S)= min(1.0d0,max(0.0,(vidim(ixO^S)+dsqrt(csound2(ixO^S)))/ &
                  (1.0d0+dsqrt(csound2(ixO^S))*vidim(ixO^S))))
     cmin(ixO^S)= max(-1.0d0,min(0.0,(vidim(ixO^S)-dsqrt(csound2(ixO^S)))/ &
                  (1.0d0-dsqrt(csound2(ixO^S))*vidim(ixO^S))))
  else
     v2(ixO^S)=sum(w(ixO^S,rmom(:))**2, dim=ndim+1)/ &
             (rhoh(ixO^S)*w(ixO^S,lfac_)*w(ixO^S,lfac_))**2.0d0

     cmax(ixO^S)=min(1.0d0,max(0.0,(vidim(ixO^S)*(1.0d0-csound2(ixO^S)) + &
      dsqrt(csound2(ixO^S)*(1.0d0-v2(ixO^S))*( &
      1.0d0-v2(ixO^S)*csound2(ixO^S)-vidim(ixO^S)**2.0d0*(1.0d0-csound2(ixO^S))) &
       ) ) / (1.0d0-v2(ixO^S)*csound2(ixO^S))))

     cmin(ixO^S)=max(-1.0d0,min(0.0,(vidim(ixO^S)*(1.0d0-csound2(ixO^S)) - &
       dsqrt(csound2(ixO^S)*(1.0d0-v2(ixO^S))*( &
       1.0d0-v2(ixO^S)*csound2(ixO^S)-vidim(ixO^S)**2.0d0*(1.0d0-csound2(ixO^S))) &
                 ) ) / (1.0d0-v2(ixO^S)*csound2(ixO^S))))
  endif
endif

  end subroutine srhd_get_cmax
!=============================================================================
  ! TODO: add srhd_get_cbounds (look at mod_hd_phys)
  subroutine srhd_get_cbounds(wLC, wRC, wLp, wRp, x, ixI^L, ixO^L, idim, cmax, cmin)

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative left and right status
    double precision, intent(in)    :: wLC(ixI^S, nw), wRC(ixI^S, nw)
    ! primitive left and right status
    double precision, intent(in)    :: wLp(ixI^S, nw), wRp(ixI^S, nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(inout) :: cmax(ixI^S)
    double precision, intent(inout), optional :: cmin(ixI^S)

    double precision :: wmean(ixI^S,nw)
    double precision, dimension(ixI^S) :: umean, dmean, csoundL, csoundR, tmp1,tmp2,tmp3

 if (typeboundspeed/='cmaxmean') then
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.

      tmp1(ixO^S)=sqrt(wLp(ixO^S,rho_))
      tmp2(ixO^S)=sqrt(wRp(ixO^S,rho_))
      tmp3(ixO^S)=1.d0/(sqrt(wLp(ixO^S,rho_))+sqrt(wRp(ixO^S,rho_)))
      umean(ixO^S)=(wLp(ixO^S,mom(idim))*tmp1(ixO^S)+wRp(ixO^S,mom(idim))*tmp2(ixO^S))*tmp3(ixO^S)

!      if(srhd_energy) then
!*DM* ideal
        csoundL(ixO^S)=srhd_gamma*wLp(ixO^S,p_)/wLp(ixO^S,rho_)
        csoundR(ixO^S)=srhd_gamma*wRp(ixO^S,p_)/wRp(ixO^S,rho_)
!      else
!        csoundL(ixO^S)=hd_gamma*hd_adiab*wLp(ixO^S,rho_)**(hd_gamma-one)
!        csoundR(ixO^S)=hd_gamma*hd_adiab*wRp(ixO^S,rho_)**(hd_gamma-one)
!      end if

      dmean(ixO^S) = (tmp1(ixO^S)*csoundL(ixO^S)+tmp2(ixO^S)*csoundR(ixO^S)) * &
           tmp3(ixO^S) + 0.5d0*tmp1(ixO^S)*tmp2(ixO^S)*tmp3(ixO^S)**2 * &
           (wRp(ixO^S,mom(idim))-wLp(ixO^S,mom(idim)))**2

      dmean(ixO^S)=sqrt(dmean(ixO^S))
      if(present(cmin)) then
        cmin(ixO^S)=umean(ixO^S)-dmean(ixO^S)
        cmax(ixO^S)=umean(ixO^S)+dmean(ixO^S)
      else
        cmax(ixO^S)=abs(umean(ixO^S))+dmean(ixO^S)
      end if

 else
      wmean(ixO^S,1:nwflux)=0.5d0*(wLC(ixO^S,1:nwflux)+wRC(ixO^S,1:nwflux))
      tmp1(ixO^S)=wmean(ixO^S,mom(idim))/wmean(ixO^S,rho_)
      call hd_get_csound2(wmean,x,ixI^L,ixO^L,csoundR)
      csoundR(ixO^S) = sqrt(csoundR(ixO^S))

      if(present(cmin)) then
        cmax(ixO^S)=max(tmp1(ixO^S)+csoundR(ixO^S),zero)
        cmin(ixO^S)=min(tmp1(ixO^S)-csoundR(ixO^S),zero)
      else
        cmax(ixO^S)=abs(tmp1(ixO^S))+csoundR(ixO^S)
      end if

 end if

  end subroutine srhd_get_cbounds
!=============================================================================
 ! Calculate flux f_idim[iw]
  subroutine srhd_get_flux(wC, w, x, ixI^L, ixO^L, idim, f)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, idim
    ! conservative w
    double precision, intent(in)    :: wC(ixI^S, 1:nw)
    ! primitive w
    double precision, intent(in)    :: w(ixI^S, 1:nw)
    double precision, intent(in)    :: x(ixI^S, 1:ndim)
    double precision, intent(out)   :: f(ixI^S, nwflux)
    double precision                :: pth(ixO^S)
    integer                         :: idir, itr

!*DM* this is also not needed...
!    if (hd_energy) then
!       pth(ixO^S) = w(ixO^S,p_)
!    else
!       pth(ixO^S) = hd_adiab * w(ixO^S, rho_)**hd_gamma
!    end if

!*DM* This is also ok ?
    f(ixO^S, rho_) = w(ixO^S,mom(idim)) * w(ixO^S, rho_)

    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
      f(ixO^S, mom(idir)) = w(ixO^S,mom(idim)) * wC(ixO^S, mom(idir))
    end do

    f(ixO^S, mom(idim)) = f(ixO^S, mom(idim)) + pth(ixO^S)

!*DM* not needed ?
!    if(hd_energy) then
      ! Energy flux is v_i*e + v*p ! Check? m_i/rho*p
!      f(ixO^S, e_) = w(ixO^S,mom(idim)) * (wC(ixO^S, e_) + w(ixO^S,p_))
!    end if

    do itr = 1, srhd_n_tracer
       f(ixO^S, tracer(itr)) = w(ixO^S,mom(idim)) * w(ixO^S, tracer(itr))
    end do

  end subroutine srhd_get_flux

!=============================================================================
  subroutine srhd_con2prim(pressure,lfac,d,rmom,tau,ierror) 
  !*DM* Change s^C to rmom

    !use ieee_arithmetic
    use mod_global_parameters

    double precision:: pressure,lfac
    double precision:: d,rmom,tau  !*DM* rmom instead of s^C
    integer:: ierror

    integer:: ni,niiter
    double precision:: pcurrent,pnew
    double precision:: er,er1,ff,df,dp,vvel !vvel instead of v^C
    double precision:: pmin,lfac2inv,pLabs,pRabs,pprev
    double precision:: s2overcubeG2rh,sqrs
    double precision:: xicurrent
    double precision:: oldff1,oldff2
    double precision:: Nff
    double precision:: pleft,pright,pnewi
    integer::nit,n2it,ni2,ni3
    double precision:: h,dhdp
    !-----------------------------------------------------------------------------

    ierror=0
    ! ierror=0 : ok
    ! 
    ! ierror<>0
    ! 
    ! ierror=1 : error on entry: must have D>=minrho, tau>=smalltau
    ! ierror=2 : maxitnr reached without convergence
    ! ierror=3 : final pressure value < smallp or xi<smallxi during iteration
    ! ierror=4 : final v^2=1 hence problem as lfac=1/0
    ! ierror=5 : nonmonotonic function f?
    ! ierror=7 : stop due to strictnr violation

    if(d<minrho .or. tau<smalltau) then
       ierror=1
       return
    endif

    sqrs={s^C**2+} !*DM* what?

    if(sqrs==zero)then
       call pressureNoFlow(pressure,tau,d)
       lfac=one
       return
    endif

    ! left and right brackets for p-range
    pmin=sqrt(sqrs)/(one-dmaxvel)-tau-d
    pLabs=max(minp,pmin)
    pRabs=1.0d99
    ! start value from input
    !opedit: try to use previous value:
    if (pLabs .le. pressure .and. pressure .le. pRabs) then
       pcurrent = pressure
    else
       pcurrent = pLabs
    end if

    er1=one
    pprev=pcurrent

    ! Fudge Parameters
    oldff1=1.0d7  ! High number
    oldff2=1.0d9  ! High number bigger then oldff1
    n2it = 0
    nit  = 0

    LoopNR:  do ni=1,maxitnr
       nit = nit + 1
       !============= Controle ~1~=============!
       if(nit>maxitnr/4)then
          !print *,'ni,er,p',ni,er,pcurrent
          ! mix pressure value for convergence
          pcurrent=half*(pcurrent+pprev)
          ! relax accuracy requirement
          er1=10.*er1
          nit = nit - maxitnr/10
       endif
       !=======================================!

       niiter=ni  
       xicurrent=tau+d+pcurrent

       if(xicurrent<smallxi) then
          if(strictgetaux)then
             print*,'!--- amrvacphys/t.srhd-- con2prim ---!'
             print *,'stop: too small xi iterate:',xicurrent
             print *,'for pressure iterate p',pcurrent
             print *,'pressure bracket pLabs pRabs',pLabs,pRabs
             print *,'iteration number:',ni
             print *,'values for d,s,tau,rmom^2:',d,rmom,tau,sqrs
          endif
          ierror=3 
          return
       endif

!*DM* Check the following substitution...
!       {v^C=s^C/xicurrent\}
!       lfac2inv=one - ({v^C**2+})
       vvel(:)=rmom(:)/xicurrent
       lfac2inv=1-sum(vvel(:)**2)
       if(lfac2inv>zero) then
          lfac=one/sqrt(lfac2inv)
       else
          if(strictgetaux)then
             print*,'!--- amrvacphys/t.srhd-- con2prim ---!'
             print *,'stop: negative or zero factor 1-v^2:',lfac2inv
             print *,'for pressure iterate p',pcurrent
             print *,'absolute pressure bracket pLabs pRabs',pLabs,pRabs
             print *,'iteration number:',ni
             print *,'values for d,s,tau,s^2:',d,rmom(:),tau,sqrs
             print *,'values for v,xi:',vvel(:),xicurrent
          endif
          ierror=4
          return
       endif


       s2overcubeG2rh=sqrs/(xicurrent**3.0d0)
       !== ZM calculation done using the EOS ==!
       call FuncEnthalpy(pcurrent,lfac2inv,d,rmom(:),tau,sqrs,xicurrent,&
            s2overcubeG2rh,h,dhdp,ierror)
       !=======================================!   
       ff=-xicurrent*lfac2inv + h 
       df=- two*sqrs/(xicurrent)**2  + dhdp - lfac2inv

       if (ff*df==zero) then
          if (ff==zero) then
             exit ! zero found
          else
             if(strictgetaux)print *,'stop: df becomes zero, non-monotonic f(p)!!'
             ierror=5
             return
          endif
       else 
          pnew=pcurrent-ff/df
          if (ff*df>zero) then
             ! pressure iterate has decreased
             ! restrict to left 
             pnew=max(pnew,pLabs)
          else  ! ff*df<0
             ! pressure iterate has increased
             ! restrict to right 
             pnew=min(pnew,pRabs)
          endif
       endif


       ! handle special case where NR incorrectly believes in convergence
       if(pnew == pLabs .and. pcurrent==pnew .and. &
            abs(ff)> absaccnr .and. sqrs > zero)then
          pnewi=pnew
          ! try 2 higher pressure values to locate a sign change for f(p)
          LoopCor:  do ni2=1,2
             !=====================!
             pcurrent=pnewi*500.0d0
             xicurrent=tau+d+pcurrent
             vvel(:)=rmom(:)/xicurrent !*DM*
             lfac2inv=1.0d0-sum(vvel(:)**2)
             !{v^C=s^C/xicurrent\}
             !lfac2inv=one - ({v^C**2+})
             !=====================!

             !=====================!
             if(lfac2inv>zero)then
                lfac=one/sqrt(lfac2inv)
             else
                ierror=4
                return
             endif
             !=====================!

             s2overcubeG2rh=-sqrs/(xicurrent**3.0d0)
             !==== Calculate enthalpy and derivative ====!
             call Bisection_Enthalpy(pcurrent,lfac2inv,d,s^C,&
                  tau,sqrs,xicurrent,h,ierror)
             Nff=-xicurrent*lfac2inv + h

             !== Save old value of pressure ==!
             pnewi=pcurrent
             !================================!

             !== find the interval where is the root ==!
             if(Nff * ff <=zero)then
                pnew=pcurrent
                exit LoopCor
             endif
             !=========================================!
          enddo LoopCor

          !== No possible solution, correct all including the conservatives ==!
          if( Nff*ff>zero)then

             ! following is in accord with trick done in smallvalues
             d   = 2.0d0*(one + 10.0d0 * minrho) * minrho
             tau = 2.0d0*(one + 10.0d0 * smalltau) * smalltau
             {^C&s^C =zero;}
             pressure     = (eqpar(gamma_)-one)*tau
             lfac = one

             if(strictnr)ierror=7
             ! leave the do loop here
             return
          endif
       endif
       !===============================================!
       dp=pcurrent-pnew
       er=two*dabs(dp)/(pnew+pcurrent)
       if(((er<tolernr*er1).or.(dabs(dp)<absaccnr))) exit LoopNR
       !===============================================!

       ! For very small values of pressure, NR algorithm is not efficient to
       ! find root, use Euler algorithm to find precise value of pressure
       if((dabs(oldff2-ff) < 1.0d-8 .or. niiter >= maxitnr-maxitnr/20).and.&
            ff * oldff1 < zero    .and.  dabs(ff)>absaccnr)then

          n2it=n2it+1
          if(n2it <= 3) pcurrent=half*(pnew+pcurrent)
          if(n2it >3)then
             pright =pcurrent
             pleft=pprev
             pcurrent=half*(pleft+pright)
             Dicho:  do ni3=1,maxitnr
                !===================!
                xicurrent=tau+d+pcurrent
                vvel(:)=rmom(:)/xicurrent
                lfac2inv=1-sum(vvel(:)**2)
                !{v^C=s^C/xicurrent\}
                !lfac2inv=one - ({v^C**2+})
                if(lfac2inv>zero)then
                   lfac=one/sqrt(lfac2inv)
                else
                   ierror=4
                   return
                endif
                !===================!


                !== ZM calculation done using the EOS ==!
                call Bisection_Enthalpy(pnew,lfac2inv,d,rmom,&
                     tau,sqrs,xicurrent,h,ierror)
                Nff=-xicurrent*lfac2inv + h 
                !=======================================!
                !==== Iterate ====!
                if(ff * Nff < zero)then
                   pleft=pcurrent
                else
                   pright=pcurrent
                endif

                pcurrent=half*(pleft+pright)
                !==================!

                !=== The iteration converge ===!
                if(2.0d0*dabs(pleft-pright)/(pleft+pright) < absaccnr &
                     .or. dabs(ff)<absaccnr)then
                   pnew=pcurrent
                   exit LoopNR
                endif
                !==============================!

                !=== conserve the last value of Nff ===!
                ff=Nff
                !======================================!
             enddo    Dicho
          endif

       else
          !====== There is no problems, continue the NR iteration ======!
          pprev=pcurrent
          pcurrent=pnew
          !=============================================================!
       endif


       !=== keep the values of the 2 last ff ===!
       oldff2=oldff1
       oldff1=ff
       !========================================!
    enddo LoopNR

    if(niiter==maxitnr)then
       !print*,' ff = ',ff,' df = ',df
       !print*,'reachs maxitnr = ', niiter
       ierror=2
       return
    endif

    if(pcurrent<minp) then
       ierror=3
       return
    endif

    !--end result for pressure and lorentz factor------!
    pressure=pcurrent
    xicurrent=tau+d+pressure
!    {v^C = s^C/xicurrent\}
!    lfac2inv=one - ({v^C**2+})
    vvel(:)=rmom(:)/xicurrent
    lfac2inv=1-sum(vvel(:)**2)
    if(lfac2inv>zero) then
       lfac=one/sqrt(lfac2inv)
    else
       ierror=4
       return
    endif
    !------------------------------!

  end subroutine srhd_con2prim
!=============================================================================
  ! Fix later (first do Cartesian)
  subroutine srhd_add_geometry(qdt,ixI^L,ixO^L,wCT,w,x)

!*DM* I did nothing here... How do we handle these geometrical switches now?

    ! Add geometrical source terms to w

    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(inout) ::  w(ixI^S,1:nw)

    integer          :: iw,idir, h1x^L{^NOONED, h2x^L}
    double precision :: tmp(ixG^T)
    logical          :: angmomfix=.false.
    !-----------------------------------------------------------------------------
!!! now handled in tvdmusclf/getaux/hll variants
!!!if(typeaxial /= 'slab') call getaux(.true.,wCT,ixI^L,ixO^L,'addgeometry')

!*DM* For now, we only deal with cartesian
    select case (typeaxial)
    case ('slab')
       ! No source terms in slab symmetry
    case default
       call mpistop("Non-slab geometry not yet supported")
    end select

    ! case ('spherical')
    !    h1x^L=ixO^L-kr(1,^D); {^NOONED h2x^L=ixO^L-kr(2,^D);}
    !    do iw=1,nwflux
    !       select case(iw)
    !          ! s[s1]=((Stheta**2+Sphi**2)/xi+2*p)/r
    !       case(s1_)
    !          tmp(ixO^S)=wCT(ixO^S,p_)*x(ixO^S,1) &
    !               *(mygeo%surfaceC1(ixO^S)-mygeo%surfaceC1(h1x^S)) &
    !               /mygeo%dvolume(ixO^S){&^CE&
    !               +wCT(ixO^S,s^CE_)**2&
    !               /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_)) }
    !          {^NOONEC
    !          ! s[s2]=-(Sr*Stheta/xi)/r
    !          !       + cot(theta)*(Sphi**2/xi+p)/r
    !       case(s2_)
    !          }
    !          {^NOONED
    !          tmp(ixO^S) = +wCT(ixO^S,p_)
    !          w(ixO^S,iw)=w(ixO^S,iw) &
    !               +qdt*tmp(ixO^S)*(mygeo%surfaceC2(ixO^S)-mygeo%surfaceC2(h2x^S)) &
    !               /mygeo%dvolume(ixO^S)
    !          }
    !          {^NOONEC
    !          tmp(ixO^S)=-(wCT(ixO^S,s1_)*wCT(ixO^S,s2_)&
    !               /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_)))
    !          }
    !          {^IFTHREEC
    !          {^NOONED
    !          tmp(ixO^S)=tmp(ixO^S)+(wCT(ixO^S,s3_)**2.0&
    !               /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))) &
    !               *cos(x(ixO^S,2))/sin(x(ixO^S,2))
    !          }
    !          ! s[s3]=-(sphi*sr/xi)/r
    !          !       -cot(theta)*(stheta*sphi/xi)/r
    !       case(s3_)
    !          if(.not.angmomfix) &
    !               tmp(ixO^S)=-(wCT(ixO^S,s3_)*wCT(ixO^S,s1_)&
    !               /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))){^NOONED &
    !               -(wCT(ixO^S,s2_)*wCT(ixO^S,s3_)&
    !               /(wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))) &
    !               *cos(x(ixO^S,2))/sin(x(ixO^S,2)) }
    !          }
    !       end select
    !       ! Divide by radius and add to w
    !       if(iw==s1_{^NOONEC.or.iw==s2_}{^IFTHREEC  &
    !            .or.(iw==s3_.and..not.angmomfix)}) &
    !            w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !    end do
    ! case ('cylindrical')
    !    do iw=1,nwflux
    !       select case (iw)
    !          ! source[sr]=sphi*vphi/radius + p/radius
    !       case (sr_)
    !          w(ixO^S,sr_)=w(ixO^S,sr_)+qdt*wCT(ixO^S,p_)/x(ixO^S,1)
    !          tmp(ixO^S)=zero
    !          {^IFPHI
    !          tmp(ixO^S)=tmp(ixO^S)+wCT(ixO^S,sphi_)**2/ &
    !               (wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))
    !          ! source[sphi]=(-sphi*vr)/radius
    !       case (sphi_)
    !          tmp(ixO^S)=-wCT(ixO^S,sphi_)*wCT(ixO^S,sr_)/ &
    !               (wCT(ixO^S,tau_)+wCT(ixO^S,d_)+wCT(ixO^S,p_))
    !          }
    !       end select
    !       ! Divide by radius and add to w
    !       if (iw==sr_.or.iw==sphi_) then
    !          w(ixO^S,iw)=w(ixO^S,iw)+qdt*tmp(ixO^S)/x(ixO^S,1)
    !       end if
    !    end do

  end subroutine srhd_add_geometry
!=============================================================================
  ! Fix later
  subroutine srhd_correctaux(ixI^L,ixO^L,w,x,patchierror,subname)

    use mod_global_parameters

    integer, intent(in)         :: ixI^L, ixO^L
    integer, intent(in)         :: patchierror(ixG^T)
    character(len=*), intent(in)   :: subname
    double precision, intent(inout):: w(ixI^S,1:nw)
    double precision, intent(in)      :: x(ixI^S,1:ndim)

    integer        :: iw, kxO^L, ix^D, i
    logical        :: patchw(ixG^T)
    !-----------------------------------------------------------------------------

    {do ix^D= ixO^LIM^D\}
    ! point with local failure identified by patchierror
    ! patchierror=-1 then from smallvalues call
    ! patchierror=1  then from getaux      call
    if (patchierror(ix^D)/=0) then
       ! verify in cube with border width nflatgetaux the presence
       ! of cells where all went ok
       do i=1,nflatgetaux
          {kxOmin^D= max(ix^D-i,ixOmin^D);
          kxOmax^D= min(ix^D+i,ixOmax^D);\}
          ! in case cells are fine within smaller cube than 
          ! the userset nflatgetaux: use that smaller cube
          if (any(patchierror(kxO^S)==0)) exit
       end do
       if (any(patchierror(kxO^S)==0))then
          ! within surrounding cube, cells without problem were found
          if (useprimitive) then
             patchw(kxO^S)=(patchierror(kxO^S)/=0)
             call primitiven(ixI^L,kxO^L,w,patchw)
          end if
          ! faulty cells are corrected by averaging here
          ! only average those which were ok and replace faulty cells
          do iw = 1,nw
             w(ix^D,iw)=sum(w(kxO^S,iw),patchierror(kxO^S)==0)&
                  /count(patchierror(kxO^S)==0)
          end do
          if (useprimitive) then
             ! in addition to those switched to primitive variables
             ! above, also switch the corrected variables
             patchw(ix^D)=.false.
             call conserven(ixI^L,kxO^L,w,patchw)
          end if
       else
          ! no cells without error were found in cube of size nflatgetaux
          ! --> point of no recovery
          write(*,*)'Getaux error:',patchierror(ix^D),'ix^D=',ix^D
          !write(*,*)'New ','d=',w(ix^D,d_),'s=', &
          !        {^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
          !write(*,*)'position  ', px(saveigrid)%x(ix^D, 1:ndim)
          write(*,*)'Called from: ',subname
          if (patchierror(ix^D)<0) then
             call mpistop("---------correctaux from smallvalues-----")
          else
             call mpistop("---------correctaux from getaux----------")
          end if
       end if
    end if
    {enddo^D&\}

  end subroutine srhd_correctaux
!=============================================================================
  ! Fix later
  subroutine srhd_handle_small_values(w,x,ixI^L,ixO^L,subname)

    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    ::subname

    !!integer                         :: posvec(ndim)
    integer, dimension(ixG^T)       :: patchierror
    !-----------------------------------------------------------------------------

    if(any(w(ixO^S,d_) < minrho) .or. any(w(ixO^S,tau_) < smalltau))then
       if(strictsmall)then
          write(*,*)'smallvalues :: for tau = ',minval(w(ixO^S,tau_)), &
               'With limit=',smalltau,' For D =', minval(w(ixO^S,d_)),&
               ' With limit=',smallrho
          write(*,*)'From::  ', subname,' iteration ', it
          call mpistop("Smallvalues with strictsmall=T failed")
       else
          if(strictgetaux)then
             where(w(ixO^S,d_) < minrho .or. w(ixO^S,tau_) < smalltau)
                w(ixO^S,d_)  = 2.0*(1.0d0 + 10.0d0 * minrho)*minrho
                w(ixO^S,tau_)= 2.0*(1.0d0 + 10.0d0 * minp)*smalltau
                {^C&w(ixO^S,s^C_) =zero;}
                w(ixO^S,lfac_)=one
             end where
          else
             where(w(ixO^S,d_) < minrho .or. w(ixO^S,tau_) < smalltau)
                patchierror(ixO^S)=-1
             else where
                patchierror(ixO^S)=0
             end where
             call correctaux(ixI^L,ixO^L,w,x,patchierror,subname)
          end if
       end if ! strict
    end if

!*DM* Well, I don't think we should include the tracers in that now...
    {#IFDEF TRACER
    if(any(w(ixO^S,Dtr1_) < minrho) .or. &
         any(w(ixO^S,Dtr1_) > 10.d10 ))then
       where(w(ixO^S,Dtr1_) < minrho .or. &
            w(ixO^S,Dtr1_) > 10.d10) 
          w(ixO^S,Dtr1_) = 0.d0
       end where
    end if
    if(any(w(ixO^S,Dtr1_) .NE. w(ixO^S,Dtr1_)))then
       where(w(ixO^S,Dtr1_) .NE. w(ixO^S,Dtr1_)) 
          w(ixO^S,Dtr1_) = 0.d0
       end where
    end if\}

    end subroutine srhd_handle_small_values
 !============================================================================
  subroutine srhd_getaux(clipping,w,x,ixI^L,ixO^L,subname)

    ! Calculate auxilary variables ixO^L from non-auxiliary entries in w
    ! clipping can be set to .true. to e.g. correct unphysical pressures,
    ! densities, v>c,  etc.

    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: clipping
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S,nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    character(len=*), intent(in)    :: subname

    integer          :: err,ix^D
    double precision :: dold,tauold,{^C&sold^C_},pold,lfacold
    integer          :: patchierror(ixG^T)
    !-----------------------------------------------------------------------------

    ! artificial replacement of small D and tau values
    ! with corresponding smallrho/smallp settings,
    ! together with nullifying momenta
    call srhd_handle_small_values(w,x,ixI^L,ixO^L,subname)

    ! we compute auxiliaries p, lfac from D,S,tau
    ! put the p and lfac in the auxiliary fields lfac_ and p_
    ! on entry: p_ field may contain previous value for iteration
    ! however, for filling ghost cells, this does not exist, so we no longer
    ! use it

    {do ix^D= ixO^LIM^D\}
    dold=w(ix^D,d_)
    tauold=w(ix^D,tau_)
    pold=w(ix^D,p_)
    lfacold=w(ix^D,lfac_)
    { ^C&sold^C_=w(ix^D,s^C_);}

    call srhd_con2prim(w(ix^D,p_),w(ix^D,lfac_), &
         w(ix^D,d_),{^C&w(ix^D,s^C_)},w(ix^D,tau_),err)
    patchierror(ix^D)=err
    if (err/=0.and.strictgetaux) then
       write(*,*)'Getaux error:',err,'ix^D=',ix^D,' global timestep it=',it
       write(*,*)'Postion: ',x(ix^D,1:ndim)
       write(*,*)'start value for p=',pold
       write(*,*)'start value for lfac=',lfacold
       write(*,*)'end value for p=',w(ix^D,p_)
       write(*,*)'end value for lfac=',w(ix^D,lfac_)
       write(*,*)'exit  d=',w(ix^D,d_),'s=',{^C&w(ix^D,s^C_)},'tau=',w(ix^D,tau_)
       write(*,*)'input d=',dold,'s=',{^C&sold^C_},'tau=',tauold
       write(*,*)'Called from: ',subname
       call mpistop("problem in getaux")
    endif
    {enddo^D&\}

    if(.not.strictgetaux.and.any(patchierror(ixO^S)/=0)) &
         call srhd_correctaux(ixI^L,ixO^L,w,x,patchierror,subname)

  end subroutine srhd_getaux
!=============================================================================
end module mod_srhd_phys
