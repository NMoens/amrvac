module mod_usr
  use mod_hd
  implicit none
  double precision :: beta, eta1, eta2, f, q
  double precision :: v_inf, v_esc
  logical :: grav2, Finite_Disk
  double precision :: R_jet, h_jet, Mdot, Edot
  ! q : M_donor / M_accretor = M_primary / M_secondary
  ! f : filling factor defined w.r.t. dust condensation radius
  ! eta1 : vinf / vorb !> Not used anymore, v_inf now calculated from cak
  ! eta2 : vinf / cs @ dust condensation radius

  double precision :: alpha, Gamma, Qbar, rho0
  ! CAK parameters

  double precision :: Rd, gm1, gm2, Mdot4pi, a, Egg, csd
  double precision :: dr_smooth

  double precision :: i_r
  double precision :: i_cak, i_rch1, i_rch2, i_rch3
  double precision :: i_geff, i_cor1, i_cor2, i_cor3

contains

  subroutine usr_init()

    call set_coordinate_system('spherical_3D')

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_source        => pt_grav_source
    usr_get_dt        => get_dt_pt_grav
    usr_refine_grid   => specialrefine_grid
    usr_special_bc    => specialbound_usr
    usr_internal_bc   => jet_launch
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    call params_read(par_files)

    call hd_activate()

    i_r = var_set_extravar("r", "r")
    i_cak = var_set_extravar("CAK", "CAK")
    i_geff = var_set_extravar("geff", "geff")
    i_rch1 = var_set_extravar("rch1", "rch1")
    i_rch2 = var_set_extravar("rch2", "rch2")
    i_rch3 = var_set_extravar("rch3", "rch3")
    i_cor1 = var_set_extravar("cor1", "cor1")
    i_cor2 = var_set_extravar("cor2", "cor2")
    i_cor3 = var_set_extravar("cor3", "cor3")

  end subroutine usr_init

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /my_list/ beta, eta1, eta2, f, q, grav2

    namelist /jet_list/ R_jet, h_jet, Mdot, Edot

    namelist /cak_list/ alpha, Gamma, Qbar, rho0, Finite_Disk

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, my_list, end=111)
111    rewind(unitpar)
       read(unitpar, jet_list, end=112)
112    rewind(unitpar)
       read(unitpar, cak_list, end=113)
113    close(unitpar)
    end do

  end subroutine params_read

  subroutine initglobaldata_usr
    use mod_global_parameters
    integer :: i, pp
    double precision :: dr1, drp

    if (eta2<1.d0) call mpistop("no solution for vinf < cs")

    ! ------------------------------------------------------
    ! The geometry
    ! ------------------------------------------------------

    Egg = ( 0.49*q**(2./3.) ) / ( 0.6*q**(2./3.) + dlog(1.+q**(1./3.)) )
    Rd = f*Egg ! in units of a
    a  = 1.d0/(f*Egg) ! in units of Rd
    if (a>xprobmax1) call mpistop('Simulation space not extended to capture the orbital separation')

    ! Length normalization
    Rd   = 1.d0

    ! Convert from orb. sep. to Rd units
    R_jet=R_jet*a
    h_jet=h_jet*a

    ! ------------------------------------------------------
    ! The speeds
    ! ------------------------------------------------------

    ! Sound speed @ dust condensation radius in units of aOmega
    v_esc = dsqrt(2/(f*Egg)**2*q/(1+q))
    v_inf = dsqrt(alpha/(1-alpha)) * v_esc

    csd=v_inf/eta2

    ! ------------------------------------------------------
    ! The accretion radius, used for refinment
    ! ------------------------------------------------------

    ! The smoothing length for the graviational potential,
    ! defined as dr @ the orbital separation on the finest
    ! AMR level
    ! pp : index of cell which contains the secondary on the 1st AMR level
    if (.not. any(stretched_dim)) call mpistop('stretching needed to compute smoothing length')
    pp=int(domain_nx1*(dlog(a/xprobmin1)/dlog(xprobmax1/xprobmin1)))
    ! dr of 1st cell @ the 1st AMR level
    dr1=(xprobmax1-xprobmin1)*(1.-qstretch(1,1))/(1.-qstretch(1,1)**(dble(domain_nx1)))
    ! dr @ r=a @ the 1st AMR level
    drp=dr1*qstretch(1,1)**dble(pp-1)
    ! dr @ r=a @ the finest AMR level
    dr_smooth=drp/(2.d0**(refine_max_level-1))

    ! ------------------------------------------------------
    ! The physics
    ! ------------------------------------------------------

    ! rhos (density @ sonic point) set to 1 for normalization
    ! Mdot4pi=csd*Rd**2.d0 ! css*Rs**2.d0
    Mdot4pi=rho0*csd*Rd**2.d0 ! css*Rs**2.d0

    ! CAK Masslosrate, use this one to set density at base
    if (Finite_Disk) Mdot4pi = Mdot4pi/(1+alpha)**(1.d0/alpha)

    ! S computed @ sonic point. Should yield Mach=1 @ sonic point.
    if (hd_adiab>0.d0) hd_adiab=((v_inf/eta2)**2.d0/hd_gamma) ! css**2.d0/(hd_gamma*1.d0**(hd_gamma-1.d0))
    ! hd_adiab=0.d0

    gm1=(q/(1.d0+q))/(f*Egg)
    gm2=(1.d0/(1.d0+q))/(f*Egg)

    if (mype==0) then
      print*, 'Rd, Rr, a, hd_adiab'
      print*, Rd, Egg*a, a, hd_adiab
      print*, eta1, beta
    endif

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision :: parker(ixG^S), vbeta(ixG^S), cor(ixG^S,1:ndir)

    call compute_beta(ixG^L,ix^L,x,vbeta)

    w(ix^S,rho_)  = Mdot4pi/(x(ix^S,1)**2.d0*vbeta(ix^S))
    w(ix^S,mom(1))= Mdot4pi/x(ix^S,1)**2.d0
    w(ix^S,mom(2))= 0.d0
    w(ix^S,mom(3))= 0.d0

    w(ix^S,i_r) = x(ix^S,1)

  end subroutine initonegrid_usr

  subroutine pt_grav_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: bbb, acc(ixI^S), cent_str(ixI^S,1:ndir), cor(ixI^S,1:ndir), roche(ixI^S,1:ndir)
    double precision :: cak(ixO^S), g_eff(ixO^S)
    integer :: i

    integer :: ix^D

      ! Full acceleration
      ! call compute_beta(ixI^L,ixO^L,x,acc)
      ! acc(ixO^S) = acc(ixO^S) * beta * (1.d0/x(ixO^S,1)) * (acc(ixO^S)-eta1/eta2) / (x(ixO^S,1)-1.d0)
      ! w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt * wCT(ixO^S,rho_  ) * acc(ixO^S)

      !CAK and electron scattering:
      call get_cak_acc(ixI^L,ixO^L,x,wCT,cak)
      w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt * wCT(ixO^S,rho_  ) * cak(ixO^S)
      w(ixO^S,i_cak) = cak(ixO^S)

      call get_g_eff(ixI^L,ixO^L,x,g_eff)
      w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt * wCT(ixO^S,rho_  ) * g_eff(ixO^S)
      w(ixO^S,i_geff) = g_eff(ixO^S)

      ! ! Additional forces due to companion :
      ! ! Wobbling around center of mass + Gravity of companion
      ! call get_roche_else(ixI^L,ixO^L,x,roche)
      ! w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt * wCT(ixO^S,rho_  ) * roche(ixO^S,1)
      ! w(ixO^S,mom(2)) = w(ixO^S,mom(2)) + qdt * wCT(ixO^S,rho_  ) * roche(ixO^S,2)
      ! w(ixO^S,mom(3)) = w(ixO^S,mom(3)) + qdt * wCT(ixO^S,rho_  ) * roche(ixO^S,3)
      ! w(ixO^S,i_rch1) = roche(ixO^S,1)
      ! w(ixO^S,i_rch2) = roche(ixO^S,2)
      ! w(ixO^S,i_rch3) = roche(ixO^S,3)

      ! Cheaty Roche-lobe fix
      ! {do ix^D = ixOmin^D,ixOmax^D\}
      ! if (w(ix^D,mom(1)) .ne. w(ix^D,mom(1))) &
      !   w(ixO^S,mom(1)) = wCT(ixO^S,mom(1))
      ! ! gradv(ix^D) = dabs(gradv(ix^D))
      ! {enddo\}

      ! Compute centrifugal centered @ the star and Coriolis and add them together
      ! call get_coriolis(ixI^L,ixO^L,x,wCT,cor)
      ! w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt * wCT(ixO^S,rho_  ) * cor(ixO^S,1)
      ! w(ixO^S,mom(2)) = w(ixO^S,mom(2)) + qdt * wCT(ixO^S,rho_  ) * cor(ixO^S,2)
      ! w(ixO^S,mom(3)) = w(ixO^S,mom(3)) + qdt * wCT(ixO^S,rho_  ) * cor(ixO^S,3)
      ! w(ixO^S,i_cor1) = cor(ixO^S,1)
      ! w(ixO^S,i_cor2) = cor(ixO^S,2)
      ! w(ixO^S,i_cor3) = cor(ixO^S,3)

  end subroutine pt_grav_source

  subroutine get_dt_pt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)
    use mod_global_parameters
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
    double precision, intent(in) :: w(ixG^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: cak(ix^S), g_eff(ix^S), force(ixG^S,1:ndir), &
      cent_str(ixG^S,1:ndir), cor(ixG^S,1:ndir), roche(ixG^S,1:ndir)
    integer :: i

    dtnew=bigdouble

    ! call compute_beta(ixG^L,ix^L,x,acc)
    call get_cak_acc(ixG^L,ix^L,x,w,cak)
    call get_g_eff(ixG^L,ix^L,x,g_eff)
    ! Compute centrifugal centered @ the star and Coriolis and add them together
    !call get_coriolis(ixG^L,ix^L,x,w,cor)
    !call get_roche_else(ixG^L,ix^L,x,roche)
    !force(ix^S,1:ndim)=cor(ix^S,1:ndim)+roche(ix^S,1:ndim)
    ! force(ix^S,1)=force(ix^S,1)+cak(ix^S)+g_eff(ix^S)
    force(ix^S,1)=cak(ix^S)+g_eff(ix^S)
    dtnew=min(dtnew,minval(dsqrt(block%dx(ix^S,1)/&
    dabs(force(ix^S,1)))),&
                  minval(dsqrt(block%dx(ix^S,1)*block%dx(ix^S,2)/&
    dabs(force(ix^S,2)))),&
                  minval(dsqrt(block%dx(ix^S,1)*dsin(block%dx(ix^S,2))*block%dx(ix^S,3)/&
    dabs(force(ix^S,3)))))

  end subroutine get_dt_pt_grav


  subroutine specialbound_usr(qt,ixG^L,ixB^L,iB,w,x)
    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    integer :: ii,jj,kk
    double precision :: vbeta(ixG^S)

    select case(iB)
    case(1)
      call compute_beta(ixG^L,ixB^L,x,vbeta)
      ! w(ixB^S,rho_)  =Mdot4pi/(x(ixB^S,1)**2.d0*vbeta(ixB^S))
      ! w(ixB^S,mom(1))= Mdot4pi/x(ixB^S,1)**2.d0

      w(ixB^S,rho_)  = rho0

      do ii = ixBmax1,ixBmin1,-1
        do jj = ixGmin2,ixBmax2
          do kk = ixGmin3,ixGmax3
            w(ii,jj,kk,mom(1)) = w(ii+1,jj,kk,mom(1))
            ! w(ii,jj,kk,mom(1)) = min(w(ii,jj,kk,mom(1)), 0.99*csd*rho0)
            ! w(ii,jj,kk,mom(1)) = min(w(ii,jj,kk,mom(1)), -0.5*csd*rho0)
          end do
        end do
      end do
      w(ixB^S,mom(2))= 0.d0
      w(ixB^S,mom(3))= 0.d0

    case default
      call mpistop("special BC not defined")
    end select

  end subroutine specialbound_usr

  ! -------------------------------------------------------------------------------
  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  integer, intent(in)                :: ixI^L,ixO^L
  double precision, intent(in)       :: x(ixI^S,1:ndim)
  double precision                   :: w(ixI^S,nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)
  double precision :: parker(ixI^S), force(ixI^S,1:ndir), vbeta(ixI^S), acc(ixI^S)
  integer :: idir, jxO^L, hxO^L
  double precision :: dmdot(ixI^S), mdot(ixI^S)

   call compute_beta(ixI^L,ixO^L,x,vbeta)
   w(ixO^S,nw+1) = vbeta(ixO^S) ! eta1 * (1.d0-Rd/x(ixO^S,1))**beta
   w(ixO^S,nw+2) = dsqrt(w(ixO^S,mom(1))/w(ixO^S,rho_)) / &
                   dsqrt(hd_gamma*hd_adiab*w(ixO^S,rho_)**(hd_gamma-1.d0))
   ! w(ixO^S,nw+2) = dabs(w(ixO^S,mom(1))/w(ixO^S,rho_)) / &
   !                  dsqrt((hd_gamma-1.d0)*(w(ixO^S,e_)-0.5d0*sum(w(ixO^S,mom(:))**2.d0,dim=ndim+1)/w(ixO^S,rho_))/w(ixO^S,rho_))
   ! w(ixO^S,nw+3) = (w(ixO^S,mom(1))/w(ixO^S,rho_)) / &
   !                  dsqrt(hd_gamma*hd_adiab*w(ixO^S,rho_)**(hd_gamma-1.d0))
   ! call get_roche_else(ixI^L,ixO^L,x,w,force)
   acc(ixO^S) = beta*((1.d0-1.d0/eta2)*(1.-1.d0/x(ixO^S,1))**(beta-1.d0)/x(ixO^S,1)**2.d0)*&
    (1.d0/eta2+(1.d0-1.d0/eta2)*(1.-1.d0/x(ixO^S,1))**beta)
   w(ixO^S,nw+3) = acc(ixO^S)
   w(ixO^S,nw+4) = w(ixO^S,mom(1))/w(ixO^S,rho_)
   w(ixO^S,nw+5) = w(ixO^S,mom(2))/w(ixO^S,rho_)
   w(ixO^S,nw+6) = w(ixO^S,mom(3))/w(ixO^S,rho_)
   ! r^0.6 x temperature
   ! w(ixO^S,nw+7) = x(ixO^S,1)**0.6d0*(hd_gamma-1.d0)*(w(ixO^S,e_)-0.5d0*sum(w(ixO^S,mom(:))**2.d0,dim=ndim+1)/w(ixO^S,rho_))/w(ixO^S,rho_)

   ! gradient of mdot divided by mdot, logarithmic gradient of mdot
    w(ixO^S,nw+7) = 0.d0
    do idir=1,ndir
      hxO^L=ixO^L-kr(idir,^D);
      jxO^L=ixO^L+kr(idir,^D);
      ! delta of the mdot function rho*v*r2
      dmdot(ixO^S)=dsqrt(sum(w(jxO^S,mom(:))**2.d0,dim=ndim+1))*x(jxO^S,1)**2.d0-&
                   dsqrt(sum(w(hxO^S,mom(:))**2.d0,dim=ndim+1))*x(hxO^S,1)**2.d0
      ! geometric average of mdots between hxO, ixO & jxO
      mdot(ixO^S)=(dsqrt(sum(w(hxO^S,mom(:))**2.d0,dim=ndim+1))*x(hxO^S,1)**2.d0*&
                   dsqrt(sum(w(ixO^S,mom(:))**2.d0,dim=ndim+1))*x(ixO^S,1)**2.d0*&
                   dsqrt(sum(w(jxO^S,mom(:))**2.d0,dim=ndim+1))*x(jxO^S,1)**2.d0)**(1.d0/3.d0)
      select case(idir)
      case(1) ! log grad along r
        w(ixO^S,nw+7) = w(ixO^S,nw+7) + (dmdot(ixO^S)/mdot(ixO^S))/(x(jxO^S,1)-x(hxO^S,1))
      case(2) ! + log grad along theta
        w(ixO^S,nw+7) = w(ixO^S,nw+7) + (dmdot(ixO^S)/mdot(ixO^S))/((x(jxO^S,2)-x(hxO^S,2))*x(ixO^S,1))
      case(3) ! + log grad along phi
        w(ixO^S,nw+7) = w(ixO^S,nw+7) + (dmdot(ixO^S)/mdot(ixO^S))/((x(jxO^S,3)-x(hxO^S,3))*x(ixO^S,1)*dsin(x(ixO^S,2)))
      end select
    enddo

  end subroutine specialvar_output
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='vbeta Mach Frad vr vt vp cs2 log_grad_mdot'
  end subroutine specialvarnames_output
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
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
    double precision :: RtoB(ixG^S), xrefine

    ! Set the center of refinement in xrefine,
    ! a bit ahead of the accretor (orb. sep. minus the accretion radius typically)
    xrefine=0.95*a
    RtoB(ix^S) = dsqrt( x(ix^S,1)**2.d0 + xrefine**2.d0 - 2.d0*x(ix^S,1)*xrefine*dsin(x(ix^S,2))*dcos(x(ix^S,3)) )

    ! ! Enable max refinement in the vicinity of the C.O.
    ! ! if (any(RtoB(ix^S)<2.d0*block%dx(ix^S,1))) refine = 1
    ! if (all(RtoB(ix^S)>0.05*a)) refine = -1
    !
    ! ! if within R_jet from jet axis, refinement enabled
    ! if (all(dsqrt((x(ix^S,1)*dsin(x(ix^S,2))*dcos(x(ix^S,3))-a)**2.d0+ &
    !           (x(ix^S,1)*dsin(x(ix^S,2))*dsin(x(ix^S,3))  )**2.d0)>R_jet)) &
    !             refine = -1

    ! Enable max refinement in the vicinity of the C.O. and in the jet
    ! if (all(RtoB(ix^S)>0.05*a)) then
      if (all(dsqrt((x(ix^S,1)*dsin(x(ix^S,2))*dcos(x(ix^S,3))-a)**2.d0+ &
                    (x(ix^S,1)*dsin(x(ix^S,2))*dsin(x(ix^S,3))  )**2.d0)>R_jet)) &
                      refine = -1
    ! endif

    ! Never refine @ poles
    ! if (node(pig2_,igrid)<2**level) refine = -1

  end subroutine specialrefine_grid
  ! -------------------------------------------------------------------------------

  ! -----------------------------------------------------------------------------------
  ! Compute beta law velocity profile
  ! -----------------------------------------------------------------------------------
  subroutine compute_beta(ixI^L,ixO^L,x,speed)
  use mod_global_parameters
  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in)  :: x(ixI^S,ndim)
  double precision, intent(out) :: speed(ixI^S)
  integer :: i

  ! Simplified beta-law (ie for csd<<vinf)
  ! speed(ixO^S) = 1.d0 * (1.d0-Rd/x(ixO^S,1))**beta
  ! Full (aka modified) beta law
  ! speed(ixO^S) = 1.d0/eta2 + (1.d0-1.d0/eta2) * (1.d0-Rd/x(ixO^S,1))**beta
  ! speed(ixO^S) = eta1/eta2 + (eta1-eta1/eta2) * (1.d0-Rd/x(ixO^S,1))**beta

  speed(ixO^S) = csd + (eta2-csd)* (1.d0-Rd/x(ixO^S,1))**beta
  ! speed(ixO^S) = csd + (v_inf-csd)* (1.d0-Rd/x(ixO^S,1))**beta

  end subroutine compute_beta
  ! -----------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  ! Center @ center of masses
  ! x axis along line joining 2 bodies, from donor to accretor
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_cntfgl_str(ixI^L,ixO^L,x,F_Roche)
  use mod_global_parameters
  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in)  :: x(ixI^S,1:ndim)
  double precision, intent(out) :: F_Roche(ixI^S,1:ndir)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! centrifugal centered on star rotating @ Omega orbital so as the star does not spin
  ! in the inertial frame (but it does wobble around the center of mass, hence the need
  ! for the centrifugal term in the next subroutine)
  F_Roche(ixO^S,1)= (f*Egg)**2.d0 * (x(ixO^S,1)-x(ixO^S,1)*(dcos(x(ixO^S,2)))**2.d0)
  F_Roche(ixO^S,2)= (f*Egg)**2.d0 * x(ixO^S,1)*dcos(x(ixO^S,2))*dsin(x(ixO^S,2))
  F_Roche(ixO^S,3)= 0.d0
  end subroutine get_cntfgl_str
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  ! Rotation axis along +z, w/ y axis pointing upward in slices of the orbital plane
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_coriolis(ixI^L,ixO^L,x,wCT,F_Coriolis)
  use mod_global_parameters
  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in)  :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
  double precision, intent(out) :: F_Coriolis(ixI^S,1:ndir)
  integer :: i
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  F_Coriolis(ixO^S,1) = - 2. * (f*Egg)**1.d0 * (-wCT(ixO^S,mom(3))*dsin(x(ixO^S,2)))
  F_Coriolis(ixO^S,2) = - 2. * (f*Egg)**1.d0 * (-wCT(ixO^S,mom(3))*dcos(x(ixO^S,2)))
  F_Coriolis(ixO^S,3) = - 2. * (f*Egg)**1.d0 * ( wCT(ixO^S,mom(2))*dcos(x(ixO^S,2))+wCT(ixO^S,mom(1))*dsin(x(ixO^S,2)) )
  do i=1,ndir
    F_Coriolis(ixO^S,i) = F_Coriolis(ixO^S,i) / wCT(ixO^S,rho_  )
  enddo
  end subroutine get_coriolis
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  ! Center @ center of masses
  ! x axis along line joining 2 bodies, from donor to accretor
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_roche_else(ixI^L,ixO^L,x,F_Roche)
  use mod_global_parameters
  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in)  :: x(ixI^S,1:ndim)
  double precision, intent(out) :: F_Roche(ixI^S,1:ndim)
  double precision :: xG
  integer :: i
  double precision :: RtoB(ixI^S)

!  RtoB(ixO^S) = dsqrt( x(ixO^S,1)**2.d0 + a**2.d0 - 2.d0*x(ixO^S,1)*a*dsin(x(ixO^S,2))*dcos(x(ixO^S,3)) + 3.d0*(a*(x(3,3,4,3)-x(3,3,3,3)))**2.d0 )
  RtoB(ixO^S) = dsqrt( x(ixO^S,1)**2.d0 + a**2.d0 - 2.d0*x(ixO^S,1)*a*dsin(x(ixO^S,2))*dcos(x(ixO^S,3)) + 5.d0*dr_smooth**2.d0 )

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! For centrifugal, w.r.t CM in xG=a*(1/(1+q))
  xG=a/(1.d0+q)

  F_Roche=0.d0
  ! grav of secondary + centrifugal
  if (grav2) then
    F_Roche(ixO^S,1)= &
      ( - (gm2*(x(ixO^S,1)-a*dsin(x(ixO^S,2))*dcos(x(ixO^S,3)))) / RtoB(ixO^S)**3.d0 &
        + (f*Egg)**2.d0 * (x(ixO^S,1)-xG*dsin(x(ixO^S,2))*dcos(x(ixO^S,3))-x(ixO^S,1)*(dcos(x(ixO^S,2)))**2.d0) )
    F_Roche(ixO^S,2)= &
      ( - (gm2*(          -a*dcos(x(ixO^S,2))*dcos(x(ixO^S,3)))) / RtoB(ixO^S)**3.d0 &
        + (f*Egg)**2.d0 * (          -xG*dcos(x(ixO^S,2))*dcos(x(ixO^S,3))+x(ixO^S,1)*dcos(x(ixO^S,2))*dsin(x(ixO^S,2))  ) )
    F_Roche(ixO^S,3)= &
      ( - (gm2*(           a                 *dsin(x(ixO^S,3)))) / RtoB(ixO^S)**3.d0 &
        + (f*Egg)**2.d0 * (           xG                 *dsin(x(ixO^S,3))) )
    ! F_Roche(ixO^S,1)= - (gm2*(x(ixO^S,1)-a*dsin(x(ixO^S,2))*dcos(x(ixO^S,3)))) / RtoB(ixO^S)**3.d0
    ! F_Roche(ixO^S,2)= - (gm2*(          -a*dcos(x(ixO^S,2))*dcos(x(ixO^S,3)))) / RtoB(ixO^S)**3.d0
    ! F_Roche(ixO^S,3)= - (gm2*(           a                 *dsin(x(ixO^S,3)))) / RtoB(ixO^S)**3.d0
  else
    F_Roche(ixO^S,1)= &
      ( (f*Egg)**2.d0 * (x(ixO^S,1)-xG*dsin(x(ixO^S,2))*dcos(x(ixO^S,3))-x(ixO^S,1)*(dcos(x(ixO^S,2)))**2.d0) )
    F_Roche(ixO^S,2)= &
      ( (f*Egg)**2.d0 * (          -xG*dcos(x(ixO^S,2))*dcos(x(ixO^S,3))+x(ixO^S,1)*dcos(x(ixO^S,2))*dsin(x(ixO^S,2))  ) )
    F_Roche(ixO^S,3)= &
      ( (f*Egg)**2.d0 * (           xG                 *dsin(x(ixO^S,3))) )
  endif

  end subroutine get_roche_else
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine jet_launch(level,qt,ixI^L,ixO^L,w,x)
  integer, intent(in) :: ixI^L,ixO^L,level
  double precision, intent(in) :: qt
  double precision, intent(inout) :: w(ixI^S,1:nw)
  double precision, intent(in) :: x(ixI^S,1:ndim)
  double precision :: rho_jet, rhovz_jet
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! Define the jet launching region as a cylinder w/:
  !   - radius R_jet
  !   - height h_jet
  !   - Mdot jet
  !   - Edot jet

  rho_jet   = 0.01d0
  rhovz_jet = 0.5d0

  ! if within R_jet from jet axis
  where (dsqrt((x(ixO^S,1)*dsin(x(ixO^S,2))*dcos(x(ixO^S,3))-a)**2.d0+ &
            (x(ixO^S,1)*dsin(x(ixO^S,2))*dsin(x(ixO^S,3))  )**2.d0)<R_jet)
    ! and |z| < h_jet
    where (dabs(x(ixO^S,1)*dcos(x(ixO^S,2)))<h_jet)
      ! fill material
      w(ixO^S,rho_)=rho_jet
      w(ixO^S,mom(1))=rhovz_jet*dcos(x(ixO^S,2))
      w(ixO^S,mom(2))=rhovz_jet*(-dsin(x(ixO^S,2))) ! For upper half of jet
      w(ixO^S,mom(3))=0.d0
    end where
  end where

  end subroutine jet_launch
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine get_cak_acc(ixI^L,ixO^L,x,w,g_cak)
  use mod_global_parameters
  use mod_geometry
  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in)  :: x(ixI^S,1:ndim), w(ixI^S,1:nw)
  double precision, intent(out) :: g_cak(ixO^S)

  integer ix^D
  double precision :: vr(ixI^S), gradv(ixI^S)
  double precision :: F_fd(ixO^S), beta_fd(ixO^S)

  vr(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
  call gradient(vr,ixI^L,ixO^L,1,gradv)

  {do ix^D = ixOmin^D,ixOmax^D\}
  ! gradv(ix^D) = max(0.d0,gradv(ix^D))
  gradv(ix^D) = dabs(gradv(ix^D))
  {enddo\}

  !> Finite disk correction
  if (Finite_Disk) then
  beta_fd(ixO^S) = (1 - (vr(ixO^S)/x(ixO^S,1))/gradv(ixO^S))*(1/x(ixO^S,1))**2
  {do ix^D = ixOmin^D,ixOmax^D\}
      if (beta_fd(ix^D) .ge. 1.) then
        F_fd(ix^D) = 1./(1+alpha)
      else if (beta_fd(ix^D).lt.-1.d10) then
        F_fd(ix^D) = abs(beta)**alpha/(1+alpha)
      else if (abs(beta_fd(ix^D)).gt.1.d-3) then
        F_fd(ix^D) = (1.-(1.-beta_fd(ix^D))**(1+alpha))/(beta_fd(ix^D)*(1+alpha))
      else
        F_fd(ix^D) = 1.d0-0.5d0*alpha*beta_fd(ix^D)*(1.d0+0.333333d0*(1.-alpha)*beta_fd(ix^D))
      end if
      if (F_fd(ix^D) .lt. smalldouble) F_fd(ix^D) = one
      if (F_fd(ix^D) .gt. 5.d0) F_fd(ix^D) = one
  {enddo\}
  else
    F_fd(ix^D) = one
  endif

  g_cak(ixO^S) = (1.d0/(f*Egg)*q/(q+1))**(1-alpha) &
  *Gamma*Qbar/(1-alpha)*(1.d0/Qbar)**alpha &
  *(((1-alpha)/(alpha*Gamma))*(1-Gamma)/(Gamma*Qbar))**(1-alpha) &
  *(x(ixO^S,1)**2*vr(ixO^S)*gradv(ixO^S))**alpha*1.d0/(x(ixO^S,1)**2)

  g_cak(ixO^S) = g_cak(ixO^S)*F_fd(ixO^S)
  end subroutine get_cak_acc
  ! -------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  subroutine get_g_eff(ixI^L,ixO^L,x,grav)
  use mod_global_parameters
  use mod_geometry
  integer, intent(in) :: ixI^L, ixO^L
  double precision, intent(in)  :: x(ixI^S,1:ndim)
  double precision, intent(out) :: grav(ixO^S)

  grav(ixO^S) = -(1-Gamma)*gm1/x(ixO^S,1)**2
  end subroutine get_g_eff

end module mod_usr
