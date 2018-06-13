!> Module with Roe-type Riemann solver for relativistic hydrodynamics
!*DM* This should be independent of the srhd/srhdeos choice...

module mod_srhd_roe
  use mod_srhd_phys
  use mod_physics_roe

  implicit none
  private

!  integer,parameter :: soundRW_ = 1,soundLW_=2,entropW_=3,shearW0_=3 ! waves
!  integer,parameter :: nworkroe = 3

  integer :: soundRW_ = -1
  integer :: soundLW_ = -1
  integer :: entropW_ = -1
  integer :: shearW0_ = -1

  public :: srhd_roe_init

contains

    subroutine srhd_roe_init()
    use mod_global_parameters, only: entropycoef, nw

    integer :: iw

!*DM* Copy paste from the hd module
    if (srhd_energy) then
       ! Characteristic waves
       soundRW_ = 1
       soundLW_ = 2
       entropW_ = 3
       shearW0_ = 3
       nworkroe = 3

       phys_average => srhd_average
       phys_get_eigenjump => srhd_get_eigenjump
       phys_rtimes => srhd_rtimes
    else
       ! Characteristic waves
       soundRW_ = 1
       soundLW_ = 2
       shearW0_ = 2
       nworkroe = 1

!*DM* Are those needed ? We won't use isothermal... Modify for ideal/Mathews?
       phys_average => srhd_average_iso
       phys_get_eigenjump => srhd_get_eigenjump_iso
       phys_rtimes => srhd_rtimes_iso
    end if

    allocate(entropycoef(nw))

    do iw = 1, nw
       if (iw == soundRW_ .or. iw == soundLW_) then
          ! TODO: Jannis: what's this?
          entropycoef(iw) = 0.2d0
       else
          entropycoef(iw) = -1.0d0
       end if
    end do

  end subroutine srhd_roe_init
!=============================================================================
  !*DM* Calculate the Roe average of w, assignment of variables:
  !*DM* rho -> rho, m -> v, e -> h but these are for the hd
  
  !*DM* What is the equivalent for srhd?
  !*DM* The old code assumed
  !*DM* rho -> v0, m -> v, e-> v4

  subroutine srhd_average(wL,wR,x,ix^L,idim,wroe,workroe)
    use mod_global_parameters
    integer, intent(in)             :: ix^L, idim
    double precision, intent(in)    :: wL(ixG^T, nw), wR(ixG^T, nw)
    double precision, intent(inout) :: wroe(ixG^T, nw)
    double precision, intent(inout) :: workroe(ixG^T, nworkroe)
    double precision, intent(in)    :: x(ixG^T, 1:^ND)
    integer                         :: idir

    ! call average2(wL,wR,x,ix^L,idim,wroe,workroe(ixG^T,1),workroe(ixG^T,2))
    workroe(ix^S, 1) = sqrt(wL(ix^S,rho_))
    workroe(ix^S, 2) = sqrt(wR(ix^S,rho_))

    ! The averaged density is sqrt(rhoL*rhoR)
    wroe(ix^S,rho_)  = workroe(ix^S, 1)*workroe(ix^S, 2)

    ! Now the ratio sqrt(rhoL/rhoR) is put into workroe(ix^S, 1)
    workroe(ix^S, 1) = workroe(ix^S, 1)/workroe(ix^S, 2)

    ! Roe-average velocities
    do idir = 1, ndir
   !*DM* Here I changed mom to rmom 
       wroe(ix^S,rmom(idir)) = (wL(ix^S,rmom(idir))/wL(ix^S,rho_) * workroe(ix^S,
1)+&
            wR(ix^S,rmom(idir))/wR(ix^S,rho_))/(1.0d0+workroe(ix^S, 1))
    end do

    ! Calculate enthalpyL, then enthalpyR, then Roe-average. Use tmp2 for
    ! pressure.
    call srhd_get_pthermal(wL,x,ixG^LL,ix^L, workroe(ixG^T, 2))

    wroe(ix^S,e_)    = (workroe(ix^S, 2)+wL(ix^S,e_))/wL(ix^S,rho_)

    call srhd_get_pthermal(wR,x,ixG^LL,ix^L, workroe(ixG^T, 2))

    workroe(ix^S, 2) = (workroe(ix^S, 2)+wR(ix^S,e_))/wR(ix^S,rho_)
    wroe(ix^S,e_)    = (wroe(ix^S,e_)*workroe(ix^S, 1) + workroe(ix^S,
2))/(1.0d0+workroe(ix^S, 1))
  end subroutine srhd_average
!=============================================================================
  
subroutine average(wL,wR,x,ix^L,idim,wroe,workroe)

! Calculate the Roe average of w, assignment of variables:
! rho -> v0, m -> v, e -> v4
! check thesis Eulderink pg.18

  use mod_global_parameters

  integer:: ix^L,idim,idir
  double precision, dimension(ixG^T,nw):: wL,wR,wroe
  double precision, dimension(ixG^T,nworkroe):: workroe
  double precision, dimension(ixG^T,1:ndim):: x
  !---------------------------------------------------------------------------
  call average2(wL,wR,x,ix^L,idim,wroe,workroe(ixG^T,1),workroe(ixG^T,2))

end subroutine average
!=============================================================================
subroutine average2(wL,wR,x,ix^L,idim,wroe,tmp,tmp2)

! Calculate the Roe average of w, assignment of variables:
! rho -> v0, m -> v, e -> v4
! check thesis Eulderink pg.18

  use mod_global_parameters
  
  integer:: ixG^L,ix^L,idim,idir
  double precision, dimension(ixG^T,nw):: wL,wR,wroe
  double precision, dimension(ixG^T,1:ndim):: x
  double precision, dimension(ixG^T):: tmp,tmp2,lfL,lfR
  !---------------------------------------------------------------------------
  
  call getaux(.true.,wL,x,ixG^LL,ix^L,'average2_wL')
  call getaux(.true.,wR,x,ixG^LL,ix^L,'average2_wR')
  
  ! Calculate K_L
  tmp(ix^S) =sqrt((wL(ix^S,d_)/wL(ix^S,lfac_))+ &
       srhd_gamma*wL(ix^S,p_)/(srhd_gamma-1) )
  ! Calculate K_R, K=sqrt(rho*h)
  tmp2(ix^S) =sqrt((wR(ix^S,d_)/wR(ix^S,lfac_))+ &
       srhd_gamma*wR(ix^S,p_)/(srhd_gamma-1) )

  !!! Lorentz factor
  !!lfL(ix^S)=1/sqrt(1-(^C&wL(ix^S,v^C_)**2+))
  !!lfR(ix^S)=1/sqrt(1-(^C&wR(ix^S,v^C_)**2+))
  
  ! V^0, see thesis Eulderink
  wroe(ix^S,d_)=(wR(ix^S,lfac_)*tmp2(ix^S)+tmp(ix^S)*wL(ix^S,lfac_))/&
       (tmp2(ix^S)+tmp(ix^S))
  
  ! V^i, see thesis Eulderink, equivalent to Roe-average velocities
  do idir=1,ndir
     wroe(ix^S,s0_+idir)=( (wR(ix^S,s0_+idir)/(wR(ix^S,lfac_)*tmp2(ix^S))) &
                          +(wL(ix^S,s0_+idir)/(wL(ix^S,lfac_)*tmp(ix^S))) &
                         )/(tmp2(ix^S)+tmp(ix^S))
  end do
  
  ! V^4, see thesis Eulderink
  wroe(ix^S,tau_)=(wR(ix^S,p_)/tmp2(ix^S)+wL(ix^S,p_)/tmp(ix^S))/&
       (tmp2(ix^S)+tmp(ix^S))
  
end subroutine average2
!=============================================================================
subroutine srhd_geteigenjump(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump,workroe)

! Calculate the il-th characteristic speed and the jump in the il-th 
! characteristic variable in the idim direction within ixL. 
! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

  use mod_global_parameters
  
  integer:: ix^L,il,idim
  double precision, dimension(ixG^T,nw):: wL,wR,wroe
  double precision, dimension(ixG^T)   :: smalla,a,jump
  double precision, dimension(ixG^T,1:ndim):: x
  double precision, dimension(ixG^T,nworkroe) :: workroe
  !---------------------------------------------------------------------------
  call geteigenjump2(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump, &
       workroe(ixG^T,1),workroe(ixG^T,2),workroe(ixG^T,3))

end subroutine srhd_geteigenjump
!=============================================================================
subroutine geteigenjump2(wL,wR,wroe,x,ix^L,il,idim,smalla,a,jump, &
     csound,del,dv)

! Calculate the il-th characteristic speed and the jump in the il-th 
! characteristic variable in the idim direction within ixL. 
! The eigenvalues and the L=R**(-1) matrix is calculated from wroe. 
! jump(il)=Sum_il L(il,iw)*(wR(iw)-wL(iw))

  use mod_global_parameters
  
  integer:: ix^L,il,idim,idir
  double precision, dimension(ixG^T,nw):: wL,wR,wroe
  double precision, dimension(ixG^T)   :: smalla,a,jump,tmp,tmp2
  double precision, dimension(ixG^T,1:ndim):: x
  double precision, dimension(ixG^T)   :: csound,del,dv,del0,cp,e,k,y2
  !!save dpperc2,dvperc
  !!common /roe/ csound
  !---------------------------------------------------------------------------

  if(il==1)then
     !Square of sound speed: s^2=0.5*gam*v4*(1+v0^2-v^2)-0.5(gam-1)(1-v0^2+v^2)
     csound(ix^S)=0.5d0*srhd_gamma*wroe(ix^S,tau_)*(1.0d0+ &
     wroe(ix^S,d_)*wroe(ix^S,d_)-(^C&wroe(ix^S,s^C_)**2+))-0.5d0* &
     (srhd_gamma-1.0d0)*(1.0d0-wroe(ix^S,d_)*wroe(ix^S,d_)+(^C&wroe(ix^S,s^C_)**2+))
     
     ! Make sure that csound**2 is positive
     !csound(ix^S)=max(srhd_gamma*smalldouble/wroe(ix^S,d_),csound(ix^S))
     
     ! Calculate uR-uL
     del(ix^S)=wR(ix^S,d_)-wL(ix^S,d_)
     dv(ix^S)=wR(ix^S,s0_+idim)-wL(ix^S,s0_+idim)
     del0(ix^S)=wR(ix^S,tau_)-wL(ix^S,tau_)
     
     !Now get the correct sound speed
     csound(ix^S)=sqrt(csound(ix^S))
  endif

  !Some help variables
  cp(ix^S)=1.0d0+srhd_gamma*wroe(ix^S,tau_)/(srhd_gamma-1.0d0)
  e(ix^S)=wroe(ix^S,d_)*wroe(ix^S,d_)-wroe(ix^S,s0_+idim)*wroe(ix^S,s0_+idim)
  k(ix^S)=wroe(ix^S,d_)*(del0(ix^S)+del(ix^S))-wroe(ix^S,s0_+idim)*dv(ix^S)
  y2(ix^S)=(1.0d0-srhd_gamma*wroe(ix^S,tau_))*e(ix^S)+csound(ix^S)*csound(ix^S)

  select case(il)
  case(soundRW_)
     !lambda+=lambda2=((1-g*v4)*v0*v1+s*y)/((1-g*v4)*v0*v0+s^2)
     a(ix^S)=((1.0d0-srhd_gamma*wroe(ix^S,tau_))*wroe(ix^S,d_)*wroe(ix^S,s0_+idim)+&
          csound(ix^S)*sqrt(y2(ix^S)))/((1.0d0-srhd_gamma*wroe(ix^S,tau_))*&
          wroe(ix^S,d_)*wroe(ix^S,d_)+csound(ix^S)*csound(ix^S))
     !alp2=(s^2*k-s*y*(v0*dv-v1*(del0+del)+(g-1)*e*(del+cp*(-v0*(del0+del)+v1*dv)))/(-2*e*s^2)
     jump(ix^S)=(csound(ix^S)*csound(ix^S)*k(ix^S)-csound(ix^S)*sqrt(y2(ix^S))*&
          (wroe(ix^S,d_)*dv(ix^S)-wroe(ix^S,s0_+idim)*(del0(ix^S)+&
          del(ix^S)))+(srhd_gamma-1.0d0)*e(ix^S)*(del(ix^S)+&
          cp(ix^S)*(-wroe(ix^S,d_)*(del0(ix^S)+del(ix^S))+&
          (^C&wroe(ix^S,s^C_)*(wR(ix^S,s^C_)-wL(ix^S,s^C_))+)  )))/&
          (-2*e(ix^S)*csound(ix^S)*csound(ix^S))
  case(soundLW_)
     !lambda-=lambda1=((1-g*v4)*v0*v1-s*y)/((1-g*v4)*v0*v0+s^2)
     a(ix^S)=((1.0d0-srhd_gamma*wroe(ix^S,tau_))*wroe(ix^S,d_)*wroe(ix^S,s0_+idim)-&
          csound(ix^S)*sqrt(y2(ix^S)))/((1.0d0-srhd_gamma*wroe(ix^S,tau_))*&
          wroe(ix^S,d_)*wroe(ix^S,d_)+csound(ix^S)*csound(ix^S))
     !alp1=(s^2*k+s*y*(v0*dv-v1*(del0+del)+(g-1)*e*(del+cp*(-v0*(del0+del)+v1*dv)))/(-2*e*s^2)
     jump(ix^S)=(csound(ix^S)*csound(ix^S)*k(ix^S)+csound(ix^S)*sqrt(y2(ix^S))*&
          (wroe(ix^S,d_)*dv(ix^S)-wroe(ix^S,s0_+idim)*(del0(ix^S)+&
          del(ix^S)))+(srhd_gamma-1.0d0)*e(ix^S)*(del(ix^S)+&
          cp(ix^S)*(-wroe(ix^S,d_)*(del0(ix^S)+del(ix^S))+&
          (^C&wroe(ix^S,s^C_)*(wR(ix^S,s^C_)-wL(ix^S,s^C_))+) )))/&
          (-2.0d0*e(ix^S)*csound(ix^S)*csound(ix^S))
  case(entropW_)
     !lambda0=lambda3=v1/v0
     a(ix^S)=wroe(ix^S,s0_+idim)/wroe(ix^S,d_)
     !alp3=(2*s^2*k+(g-1)*e*(del+cp*(-v0*(del0+del)+v1*dv)))/(e*s^2)
     jump(ix^S)=(2.0d0*csound(ix^S)*csound(ix^S)*k(ix^S)+(srhd_gamma-1.0d0)*e(ix^S)*&
          (del(ix^S)+cp(ix^S)*(-wroe(ix^S,d_)*(del0(ix^S)+del(ix^S))+&
          (^C&wroe(ix^S,s^C_)*(wR(ix^S,s^C_)-wL(ix^S,s^C_))+) )))/&
          (e(ix^S)*csound(ix^S)*csound(ix^S))
  case default
     !Determine the direction of the shear wave
     idir=il-shearW0_; if(idir>=idim)idir=idir+1
     a(ix^S)=wroe(ix^S,s0_+idim)/wroe(ix^S,d_)
     !alp4_5=del2_3-kv2_3/e
     jump(ix^S)=wR(ix^S,s0_+idir)-wL(ix^S,s0_+idir)-k(ix^S)*wroe(ix^S,s0_+idir)/e(ix^S)
  end select

  ! Calculate "smalla" or modify "a" based on the "typeentropy" switch
  ! Put left and right eigenvalues, if needed, into tmp and tmp2
  ! OK, since subroutines getpthermal and entropyfix do not use tmp and tmp2

  select case(typeentropy(il))
  case('yee')
     ! Based on Yee JCP 68,151 eq 3.23
     smalla(ix^S)=entropycoef(il)
  case('harten','powell')
     call getaux(.true.,wL,x,ixG^LL,ix^L,'geteigenjump2_wL')
     call getaux(.true.,wR,x,ixG^LL,ix^L,'geteigenjump2_wR')
     ! Based on Harten & Hyman JCP 50, 235 and Zeeuw & Powell JCP 104,56
     select case(il)
     case(soundRW_)
        tmp(ix^S)=wL(ix^S,s0_+idim)/wL(ix^S,d_)&
             + sqrt(srhd_gamma*wL(ix^S,p_)/wL(ix^S,d_))
        tmp2(ix^S)=wR(ix^S,s0_+idim)/wR(ix^S,d_)&
             + sqrt(srhd_gamma*wR(ix^S,p_)/wR(ix^S,d_))
     case(soundLW_)
        tmp(ix^S)=wL(ix^S,s0_+idim)/wL(ix^S,d_)&
             - sqrt(srhd_gamma*wL(ix^S,p_)/wL(ix^S,d_))
        tmp2(ix^S)=wR(ix^S,s0_+idim)/wR(ix^S,d_)&
             - sqrt(srhd_gamma*wR(ix^S,p_)/wR(ix^S,d_))
     case default
        tmp(ix^S) =wL(ix^S,s0_+idim)/wL(ix^S,d_)
        tmp2(ix^S)=wR(ix^S,s0_+idim)/wR(ix^S,d_)
     end select
  end select
  
  call entropyfix(ix^L,il,tmp,tmp2,a,smalla)
  
end subroutine geteigenjump2
!=============================================================================
subroutine srhd_rtimes(q,wroe,ix^L,iw,il,idim,rq,workroe)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

  use mod_global_parameters
  
  integer::          ix^L,iw,il,idim
  double precision:: wroe(ixG^T,nw)
  double precision, dimension(ixG^T):: q,rq
  double precision, dimension(ixG^T,nworkroe):: workroe
!-----------------------------------------------------------------------------

  call rtimes2(q,wroe,ix^L,iw,il,idim,rq,workroe(ixG^T,1))

end subroutine srhd_rtimes
!=============================================================================
subroutine rtimes2(q,wroe,ix^L,iw,il,idim,rq,csound)

! Multiply q by R(il,iw), where R is the right eigenvalue matrix at wroe

  use mod_global_parameters

  integer::          ix^L,iw,il,idim,idir
  double precision:: wroe(ixG^T,nw)
  double precision, dimension(ixG^T):: q,rq,csound,cm,cp,e,y
  logical:: shearwave
  !!common /roe/ csound
!-----------------------------------------------------------------------------

  shearwave=il>shearW0_
  if(shearwave)then
     ! Direction of shearwave increases with il plus idir==idim is jumped over
     idir=il-shearW0_; if(idir>=idim)idir=idir+1
  endif

  cm(ix^S)=1.0d0-srhd_gamma*wroe(ix^S,tau_)/(srhd_gamma-1)
  cp(ix^S)=1.0d0+srhd_gamma*wroe(ix^S,tau_)/(srhd_gamma-1)
  e(ix^S)=wroe(ix^S,d_)*wroe(ix^S,d_)-&
       wroe(ix^S,s0_+idim)*wroe(ix^S,s0_+idim)
  y(ix^S)=sqrt((1.0d0-srhd_gamma*wroe(ix^S,tau_))*e(ix^S)+&
       csound(ix^S)*csound(ix^S))
  
  select case(iw)
  case(d_)
     select case(il)
     case(soundRW_)
        rq(ix^S)=q(ix^S)*cm(ix^S)
     case(soundLW_)
        rq(ix^S)=q(ix^S)*cm(ix^S)
     case(entropW_)
        rq(ix^S)=q(ix^S)*(cm(ix^S)+csound(ix^S)*csound(ix^S)&
             /(srhd_gamma-1.0d0))
     case default
        rq(ix^S)=-q(ix^S)*wroe(ix^S,s0_+idir)*cp(ix^S)
     end select
  case(tau_)
     select case(il)
     case(soundRW_)
        rq(ix^S)=q(ix^S)*(wroe(ix^S,d_)+csound(ix^S)*&
             wroe(ix^S,s0_+idim)/y(ix^S)-cm(ix^S))
     case(soundLW_)
        rq(ix^S)=q(ix^S)*(wroe(ix^S,d_)-csound(ix^S)*&
             wroe(ix^S,s0_+idim)/y(ix^S)-cm(ix^S))
     case(entropW_)
        rq(ix^S)=q(ix^S)*(wroe(ix^S,d_)-cm(ix^S)-&
             csound(ix^S)*csound(ix^S)/(srhd_gamma-1))
     case default
        rq(ix^S)=q(ix^S)*wroe(ix^S,s0_+idir)*cp(ix^S)
     end select
  case default
     if(iw==s0_+idim)then
        select case(il)
        case(soundRW_)
           rq(ix^S)=q(ix^S)*(wroe(ix^S,s0_+idim)+&
                csound(ix^S)*wroe(ix^S,d_)/y(ix^S))
        case(soundLW_)
           rq(ix^S)=q(ix^S)*(wroe(ix^S,s0_+idim)-&
                csound(ix^S)*wroe(ix^S,d_)/y(ix^S))
        case(entropW_)
           rq(ix^S)=q(ix^S)*wroe(ix^S,s0_+idim)
        case default
           rq(ix^S)=zero
        end select
     else
        if(shearwave)then
           if(iw==s0_+idir)then
              rq(ix^S)=q(ix^S)
           else
              rq(ix^S)=zero
           endif
        else
           rq(ix^S)=q(ix^S)*wroe(ix^S,iw)
        endif
     endif
  end select
  
end subroutine rtimes2
!=============================================================================
! *DM* Are we going to add here the relevant functions for the ideal gas ?
! *DM* (the same way that isothermal is included in the hd_module?


end module mod_srhd_roe
