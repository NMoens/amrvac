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
    if (a>xprobmax1) call mpistop&
       ('Simulation space not extended to capture the orbital separation')

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
    if (.not. any(stretched_dim)) call mpistop(&
       'stretching needed to compute smoothing length')
    pp=int(domain_nx1*(dlog(a/xprobmin1)/dlog(xprobmax1/xprobmin1)))
    ! dr of 1st cell @ the 1st AMR level
    dr1=(xprobmax1-xprobmin1)*(1.-qstretch(1,1))/(1.-qstretch(1,&
       1)**(dble(domain_nx1)))
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
    if (hd_adiab>0.d0) hd_adiab=((v_inf/eta2)**2.d0/hd_gamma) !css**2.d0/(hd_gamma*1.d0**(hd_gamma-1.d0))
    ! hd_adiab=0.d0

    gm1=(q/(1.d0+q))/(f*Egg)
    gm2=(1.d0/(1.d0+q))/(f*Egg)

    if (mype==0) then
      print*, 'Rd, Rr, a, hd_adiab'
      print*, Rd, Egg*a, a, hd_adiab
      print*, eta1, beta
    endif

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)
    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision :: parker(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3), vbeta(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3), cor(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       1:ndir)

    call compute_beta(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,ixmin1,&
       ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,x,vbeta)

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       rho_)  = Mdot4pi/(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       1)**2.d0*vbeta(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3))
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       mom(1))= Mdot4pi/x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)**2.d0
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(2))= 0.d0
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(3))= 0.d0

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,i_r) = x(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3,1)

  end subroutine initonegrid_usr

  subroutine pt_grav_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iwmin,iwmax,qtC,&
     wCT,qt,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iwmin,iwmax
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)
    double precision :: bbb, acc(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3), cent_str(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir), cor(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir), roche(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)
    double precision :: cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
        g_eff(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    integer :: i

    integer :: ix1,ix2,ix3

      ! Full acceleration
      ! call compute_beta(ixI^L,ixO^L,x,acc)
      ! acc(ixO^S) = acc(ixO^S) * beta * (1.d0/x(ixO^S,1)) * (acc(ixO^S)-eta1/eta2) / (x(ixO^S,1)-1.d0)
      ! w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt * wCT(ixO^S,rho_  ) * acc(ixO^S)

      !CAK and electron scattering:
      call get_cak_acc(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,wCT,cak)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1)) + qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_  ) * cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         i_cak) = cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

      call get_g_eff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
         ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,g_eff)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1)) + qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_  ) * g_eff(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         i_geff) = g_eff(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

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

  subroutine get_dt_pt_grav(w,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,dtnew,dx1,dx2,dx3,x)
    use mod_global_parameters
    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in) :: dx1,dx2,dx3, x(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)
    double precision, intent(in) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: cak(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3),&
        g_eff(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3),&
        force(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndir), &
      cent_str(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndir),&
          cor(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndir),&
          roche(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndir)
    integer :: i

    dtnew=bigdouble

    ! call compute_beta(ixG^L,ix^L,x,acc)
    call get_cak_acc(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,ixmin1,&
       ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,x,w,cak)
    call get_g_eff(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,ixmin1,&
       ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,x,g_eff)
    ! Compute centrifugal centered @ the star and Coriolis and add them together
    !call get_coriolis(ixG^L,ix^L,x,w,cor)
    !call get_roche_else(ixG^L,ix^L,x,roche)
    !force(ix^S,1:ndim)=cor(ix^S,1:ndim)+roche(ix^S,1:ndim)
    ! force(ix^S,1)=force(ix^S,1)+cak(ix^S)+g_eff(ix^S)
    force(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)=cak(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3)+g_eff(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3)
    dtnew=min(dtnew,minval(dsqrt(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3,1)/&
    dabs(force(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)))),&
                  minval(dsqrt(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
                     ixmin3:ixmax3,1)*block%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
                     ixmin3:ixmax3,2)/&
    dabs(force(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)))),&
                  minval(dsqrt(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
                     ixmin3:ixmax3,1)*dsin(block%dx(ixmin1:ixmax1,&
                     ixmin2:ixmax2,ixmin3:ixmax3,2))*block%dx(ixmin1:ixmax1,&
                     ixmin2:ixmax2,ixmin3:ixmax3,3)/&
    dabs(force(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,3)))))

  end subroutine get_dt_pt_grav


  subroutine specialbound_usr(qt,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
     ixGmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3,iB,w,x)
    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)

    integer :: ii,jj,kk
    double precision :: vbeta(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3)

    select case(iB)
    case(1)
      call compute_beta(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
         ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3,x,vbeta)
      ! w(ixB^S,rho_)  =Mdot4pi/(x(ixB^S,1)**2.d0*vbeta(ixB^S))
      ! w(ixB^S,mom(1))= Mdot4pi/x(ixB^S,1)**2.d0

      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,rho_)  = rho0

      do ii = ixBmax1,ixBmin1,-1
        do jj = ixGmin2,ixBmax2
          do kk = ixGmin3,ixGmax3
            w(ii,jj,kk,mom(1)) = w(ii+1,jj,kk,mom(1))
            ! w(ii,jj,kk,mom(1)) = min(w(ii,jj,kk,mom(1)), 0.99*csd*rho0)
            ! w(ii,jj,kk,mom(1)) = min(w(ii,jj,kk,mom(1)), -0.5*csd*rho0)
          end do
        end do
      end do
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(2))= 0.d0
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(3))= 0.d0

    case default
      call mpistop("special BC not defined")
    end select

  end subroutine specialbound_usr

  ! -------------------------------------------------------------------------------
  subroutine specialvar_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
  integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
  double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)
  double precision :: parker(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
      force(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndir),&
      vbeta(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
      acc(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
  integer :: idir, jxOmin1,jxOmin2,jxOmin3,jxOmax1,jxOmax2,jxOmax3, hxOmin1,&
     hxOmin2,hxOmin3,hxOmax1,hxOmax2,hxOmax3
  double precision :: dmdot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
      mdot(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

   call compute_beta(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
      ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,vbeta)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
      nw+1) = vbeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) !eta1 * (1.d0-Rd/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))**beta
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
      nw+2) = dsqrt(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
      mom(1))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)) / &
                   dsqrt(hd_gamma*hd_adiab*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      ixOmin3:ixOmax3,rho_)**(hd_gamma-1.d0))
   ! w(ixO^S,nw+2) = dabs(w(ixO^S,mom(1))/w(ixO^S,rho_)) / &
   !                  dsqrt((hd_gamma-1.d0)*(w(ixO^S,e_)-0.5d0*sum(w(ixO^S,mom(:))**2.d0,dim=ndim+1)/w(ixO^S,rho_))/w(ixO^S,rho_))
   ! w(ixO^S,nw+3) = (w(ixO^S,mom(1))/w(ixO^S,rho_)) / &
   !                  dsqrt(hd_gamma*hd_adiab*w(ixO^S,rho_)**(hd_gamma-1.d0))
   ! call get_roche_else(ixI^L,ixO^L,x,w,force)
   acc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
      beta*((1.d0-1.d0/eta2)*(1.-1.d0/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3,1))**(beta-1.d0)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3,1)**2.d0)*&
    (1.d0/eta2+(1.d0-1.d0/eta2)*(1.-1.d0/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1))**beta)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
      nw+3) = acc(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+4) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))/w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+5) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))/w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+6) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))/w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   ! r^0.6 x temperature
   ! w(ixO^S,nw+7) = x(ixO^S,1)**0.6d0*(hd_gamma-1.d0)*(w(ixO^S,e_)-0.5d0*sum(w(ixO^S,mom(:))**2.d0,dim=ndim+1)/w(ixO^S,rho_))/w(ixO^S,rho_)

   ! gradient of mdot divided by mdot, logarithmic gradient of mdot
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+7) = 0.d0
    do idir=1,ndir
      hxOmin1=ixOmin1-kr(idir,1);hxOmin2=ixOmin2-kr(idir,2)
      hxOmin3=ixOmin3-kr(idir,3);hxOmax1=ixOmax1-kr(idir,1)
      hxOmax2=ixOmax2-kr(idir,2);hxOmax3=ixOmax3-kr(idir,3);
      jxOmin1=ixOmin1+kr(idir,1);jxOmin2=ixOmin2+kr(idir,2)
      jxOmin3=ixOmin3+kr(idir,3);jxOmax1=ixOmax1+kr(idir,1)
      jxOmax2=ixOmax2+kr(idir,2);jxOmax3=ixOmax3+kr(idir,3);
      ! delta of the mdot function rho*v*r2
      dmdot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=dsqrt(sum(w(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
         jxOmin3:jxOmax3,mom(:))**2.d0,dim=ndim+1))*x(jxOmin1:jxOmax1,&
         jxOmin2:jxOmax2,jxOmin3:jxOmax3,1)**2.d0-&
                   dsqrt(sum(w(hxOmin1:hxOmax1,hxOmin2:hxOmax2,hxOmin3:hxOmax3,&
                      mom(:))**2.d0,dim=ndim+1))*x(hxOmin1:hxOmax1,&
                      hxOmin2:hxOmax2,hxOmin3:hxOmax3,1)**2.d0
      ! geometric average of mdots between hxO, ixO & jxO
      mdot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)=(dsqrt(sum(w(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
         hxOmin3:hxOmax3,mom(:))**2.d0,dim=ndim+1))*x(hxOmin1:hxOmax1,&
         hxOmin2:hxOmax2,hxOmin3:hxOmax3,1)**2.d0*&
                   dsqrt(sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
                      mom(:))**2.d0,dim=ndim+1))*x(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2.d0*&
                   dsqrt(sum(w(jxOmin1:jxOmax1,jxOmin2:jxOmax2,jxOmin3:jxOmax3,&
                      mom(:))**2.d0,dim=ndim+1))*x(jxOmin1:jxOmax1,&
                      jxOmin2:jxOmax2,jxOmin3:jxOmax3,1)**2.d0)**(1.d0/3.d0)
      select case(idir)
      case(1) ! log grad along r
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           nw+7) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           nw+7) + (dmdot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/mdot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))/(x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,1)-x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3,1))
      case(2) ! + log grad along theta
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           nw+7) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           nw+7) + (dmdot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/mdot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))/((x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,2)-x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3,2))*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,1))
      case(3) ! + log grad along phi
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           nw+7) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           nw+7) + (dmdot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3)/mdot(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3))/((x(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
           jxOmin3:jxOmax3,3)-x(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
           hxOmin3:hxOmax3,3))*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           ixOmin3:ixOmax3,2)))
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
  subroutine specialrefine_grid(igrid,level,ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
     ixGmax2,ixGmax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,qt,w,x,refine,&
     coarsen)
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
    integer, intent(in) :: igrid, level, ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
       ixGmax2,ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in) :: qt, w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw), x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision :: RtoB(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3),&
        xrefine

    ! Set the center of refinement in xrefine,
    ! a bit ahead of the accretor (orb. sep. minus the accretion radius typically)
    xrefine=0.95*a
    RtoB(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = dsqrt( x(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3,1)**2.d0 + xrefine**2.d0 - &
       2.d0*x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       1)*xrefine*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       2))*dcos(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,3)) )

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
      if (all(dsqrt((x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         1)*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
         2))*dcos(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,3))-a)**2.d0+ &
                    (x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                       1)*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                       2))*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
                       3))  )**2.d0)>R_jet)) &
                      refine = -1
    ! endif

    ! Never refine @ poles
    ! if (node(pig2_,igrid)<2**level) refine = -1

  end subroutine specialrefine_grid
  ! -------------------------------------------------------------------------------

  ! -----------------------------------------------------------------------------------
  ! Compute beta law velocity profile
  ! -----------------------------------------------------------------------------------
  subroutine compute_beta(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,speed)
  use mod_global_parameters
  integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,ndim)
  double precision, intent(out) :: speed(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3)
  integer :: i

  ! Simplified beta-law (ie for csd<<vinf)
  ! speed(ixO^S) = 1.d0 * (1.d0-Rd/x(ixO^S,1))**beta
  ! Full (aka modified) beta law
  ! speed(ixO^S) = 1.d0/eta2 + (1.d0-1.d0/eta2) * (1.d0-Rd/x(ixO^S,1))**beta
  ! speed(ixO^S) = eta1/eta2 + (eta1-eta1/eta2) * (1.d0-Rd/x(ixO^S,1))**beta

  speed(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = csd + (eta2-csd)* &
     (1.d0-Rd/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))**beta
  ! speed(ixO^S) = csd + (v_inf-csd)* (1.d0-Rd/x(ixO^S,1))**beta

  end subroutine compute_beta
  ! -----------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  ! Center @ center of masses
  ! x axis along line joining 2 bodies, from donor to accretor
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_cntfgl_str(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,F_Roche)
  use mod_global_parameters
  integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
  double precision, intent(out) :: F_Roche(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndir)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! centrifugal centered on star rotating @ Omega orbital so as the star does not spin
  ! in the inertial frame (but it does wobble around the center of mass, hence the need
  ! for the centrifugal term in the next subroutine)
  F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)= (f*Egg)**2.d0 * (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)*(dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))**2.d0)
  F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     2)= (f*Egg)**2.d0 * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     2))*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
  F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)= 0.d0
  end subroutine get_cntfgl_str
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  ! Rotation axis along +z, w/ y axis pointing upward in slices of the orbital plane
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_coriolis(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,wCT,F_Coriolis)
  use mod_global_parameters
  integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim), wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw)
  double precision, intent(out) :: F_Coriolis(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndir)
  integer :: i
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1) = - 2. * (f*Egg)**1.d0 * (-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,mom(3))*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,2)))
  F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     2) = - 2. * (f*Egg)**1.d0 * (-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,mom(3))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,2)))
  F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     3) = - 2. * (f*Egg)**1.d0 * ( wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,mom(2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,2))+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     mom(1))*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)) )
  do i=1,ndir
    F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i) = F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i) / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_  )
  enddo
  end subroutine get_coriolis
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  ! Center @ center of masses
  ! x axis along line joining 2 bodies, from donor to accretor
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_roche_else(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,F_Roche)
  use mod_global_parameters
  integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
  double precision, intent(out) :: F_Roche(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
  double precision :: xG
  integer :: i
  double precision :: RtoB(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)

!  RtoB(ixO^S) = dsqrt( x(ixO^S,1)**2.d0 + a**2.d0 - 2.d0*x(ixO^S,1)*a*dsin(x(ixO^S,2))*dcos(x(ixO^S,3)) + 3.d0*(a*(x(3,3,4,3)-x(3,3,3,3)))**2.d0 )
  RtoB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = dsqrt( &
     x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)**2.d0 + a**2.d0 - 2.d0*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,1)*a*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,3)) + 5.d0*dr_smooth**2.d0 )

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! For centrifugal, w.r.t CM in xG=a*(1/(1+q))
  xG=a/(1.d0+q)

  F_Roche=0.d0
  ! grav of secondary + centrifugal
  if (grav2) then
    F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)= &
      ( - (gm2*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)-a*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)))) / RtoB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**3.d0 &
        + (f*Egg)**2.d0 * (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)-xG*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           3))-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           1)*(dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
           2)))**2.d0) )
    F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)= &
      ( - (gm2*(          -a*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3)))) / RtoB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3)**3.d0 &
        + (f*Egg)**2.d0 * (          -xG*dcos(x(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*dcos(x(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))+x(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)*dcos(x(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))*dsin(x(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))  ) )
    F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)= &
      ( - (gm2*(           a                 *dsin(x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)))) / RtoB(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)**3.d0 &
        + (f*Egg)**2.d0 * (           xG                 &
           *dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))) )
    ! F_Roche(ixO^S,1)= - (gm2*(x(ixO^S,1)-a*dsin(x(ixO^S,2))*dcos(x(ixO^S,3)))) / RtoB(ixO^S)**3.d0
    ! F_Roche(ixO^S,2)= - (gm2*(          -a*dcos(x(ixO^S,2))*dcos(x(ixO^S,3)))) / RtoB(ixO^S)**3.d0
    ! F_Roche(ixO^S,3)= - (gm2*(           a                 *dsin(x(ixO^S,3)))) / RtoB(ixO^S)**3.d0
  else
    F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)= &
      ( (f*Egg)**2.d0 * (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)-xG*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3))-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)*(dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)))**2.d0) )
    F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)= &
      ( (f*Egg)**2.d0 * (          -xG*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,3))+x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2))*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))  ) )
    F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)= &
      ( (f*Egg)**2.d0 * (           xG                 *dsin(x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))) )
  endif

  end subroutine get_roche_else
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine jet_launch(level,qt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
  integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,level
  double precision, intent(in) :: qt
  double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw)
  double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
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
  where (dsqrt((x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     3))-a)**2.d0+ &
            (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               2))*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
               3))  )**2.d0)<R_jet)
    ! and |z| < h_jet
    where (dabs(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       1)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))<h_jet)
      ! fill material
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)=rho_jet
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1))=rhovz_jet*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2))
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(2))=rhovz_jet*(-dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ixOmin3:ixOmax3,2))) !For upper half of jet
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))=0.d0
    end where
  end where

  end subroutine jet_launch
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine get_cak_acc(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,w,g_cak)
  use mod_global_parameters
  use mod_geometry
  integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim), w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw)
  double precision, intent(out) :: g_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3)

  integer ix1,ix2,ix3
  double precision :: vr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
      gradv(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)
  double precision :: F_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3),&
      beta_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)

  vr(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = w(ixImin1:ixImax1,&
     ixImin2:ixImax2,ixImin3:ixImax3,mom(1))/w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,rho_)
  call gradient(vr,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,1,gradv)

  do ix1 = ixOmin1,ixOmax1
  do ix2 = ixOmin2,ixOmax2
  do ix3 = ixOmin3,ixOmax3
  ! gradv(ix^D) = max(0.d0,gradv(ix^D))
  gradv(ix1,ix2,ix3) = dabs(gradv(ix1,ix2,ix3))
  enddo
  enddo
  enddo

  !> Finite disk correction
  if (Finite_Disk) then
  beta_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3) = (1 - (vr(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1))/gradv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3))*(1/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1))**2
  do ix1 = ixOmin1,ixOmax1
  do ix2 = ixOmin2,ixOmax2
  do ix3 = ixOmin3,ixOmax3
      if (beta_fd(ix1,ix2,ix3) .ge. 1.) then
        F_fd(ix1,ix2,ix3) = 1./(1+alpha)
      else if (beta_fd(ix1,ix2,ix3).lt.-1.d10) then
        F_fd(ix1,ix2,ix3) = abs(beta)**alpha/(1+alpha)
      else if (abs(beta_fd(ix1,ix2,ix3)).gt.1.d-3) then
        F_fd(ix1,ix2,ix3) = (1.-(1.-beta_fd(ix1,ix2,&
           ix3))**(1+alpha))/(beta_fd(ix1,ix2,ix3)*(1+alpha))
      else
        F_fd(ix1,ix2,ix3) = 1.d0-0.5d0*alpha*beta_fd(ix1,ix2,&
           ix3)*(1.d0+0.333333d0*(1.-alpha)*beta_fd(ix1,ix2,ix3))
      end if
      if (F_fd(ix1,ix2,ix3) .lt. smalldouble) F_fd(ix1,ix2,ix3) = one
      if (F_fd(ix1,ix2,ix3) .gt. 5.d0) F_fd(ix1,ix2,ix3) = one
  enddo
  enddo
  enddo
  else
    F_fd(ix1,ix2,ix3) = one
  endif

  g_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
     (1.d0/(f*Egg)*q/(q+1))**(1-alpha) &
  *Gamma*Qbar/(1-alpha)*(1.d0/Qbar)**alpha &
  *(((1-alpha)/(alpha*Gamma))*(1-Gamma)/(Gamma*Qbar))**(1-alpha) &
  *(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2*vr(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,ixOmin3:ixOmax3)*gradv(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3))**alpha*1.d0/(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,1)**2)

  g_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
     g_cak(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3)*F_fd(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
  end subroutine get_cak_acc
  ! -------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  subroutine get_g_eff(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
     ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,grav)
  use mod_global_parameters
  use mod_geometry
  integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
  double precision, intent(out) :: grav(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3)

  grav(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = &
     -(1-Gamma)*gm1/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2
  end subroutine get_g_eff

end module mod_usr
