!program to compute initial condition for time-dependet
!(AMRVAC) 2D FLD simulations of a radiation-dominated,
!massive star atmosphere
!Using beta-law type wind as upper boundary for
!hydrostatic structure with given opacity law

!Supply file init_indat with:
!0.5,50.,20.,0.1,100,0.03,200.0,'constant'
!(gamma0,mass,rstar,Yhe,ny_amrvac,dy_amrvac,kap_law)
!ny_amrvac = number of y-points in amrvac sim
!dy_amrvac = y step-size in amrvac, in units of Hp
!IF you want to change overlying wind condtion,
!change in program directly
!---------------------------------------------------

module global

!Basic constants
  double precision, parameter :: gnewt=6.6740000d-08,rsun=6.9599000d10,&
       lsun=3.8268000d33,msun=1.9891000d33,sigb=5.6704000d-05,&
       cc=2.9979246d10,year=31557000.d0,pi=3.1415926536d0,&
       ar=4.*sigb/cc,kb=1.3806503d-16,mh=1.6735344d-24,&
       sigth=6.6524585d-25

!Basic variables for structure
double precision, dimension(:), allocatable :: pg,tg,er,rho,fr,pr,&
     gammar,tr,tau,mc,yarr,kap,t_ana
double precision :: pg0,tg0,er0,rho0,fr0,pr0,gamma0,tr0,tau0,mc0,y0,kap0
double precision :: mu,mstar,lstar,g,rstar,dy,teff,yhe,dm,tr_dum
double precision :: rhob,sigma_b,n,afw,bfw,tlim,hp,hpl,mdot,vinf,dy_vac,tau_max
integer :: ny,ny_vac
character(8) :: kap_law

!The below just for odeint routines
INTEGER :: KMAX,KOUNT
double precision :: DXSAV, XP(200),YP(10,200)

end module global

!------------------------
!------------------------

PROGRAM for_nico_init

  use global
  implicit none

  call set_boundary
  call calc_structure
  print*,'el acabose !'
  return

END PROGRAM for_nico_init

!------------------------
!------------------------

subroutine set_boundary

  use global
  implicit none

  double precision :: geff,a20
  double precision :: frac,tt,bb,v0

  !Just try simple read-in file with various input params.
  !IMPORTANT: In this version, integrating outside in. SO
  !stellar input radius is defined at outermost point,
  !NOT innermost core radius, where AMRVAC simulation
  !presumably will start from.
  !Thus set: g(amrvac) = g(lower boundary), e.g. using value
  !given in init_params_amrvac output-file
  open(1,file='init_indat')
  read(1,*) gamma0,mstar,rstar,yhe,ny_vac,dy_vac,tau_max,kap_law
  close(1)
  mstar = mstar*msun
  rstar = rstar*rsun
!  tau_max = 100.0
!  gamma0 = 0.3
!  ny_vac = 100
!  dy_vac = 0.03
  kap0 = sigth/mh * (1.+2.*yhe)/(1.+4.*yhe)
  mu = (1.+4*yhe)/(2.+3.*yhe)
  !mean molecular weight, H+He, compl. ionization
  !approx. optical depth at lower boundary
  !Start at gamma0=X as base, set other
  !parameters accordingly
  lstar = gamma0*4.*pi*gnewt*mstar*cc/kap0
  g = gnewt*mstar/rstar**2
  teff = (lstar/(4.*pi*rstar**2*sigb))**0.25
  geff = g*(1.-gamma0)
  !  a20 = kb*tg0/(mu*mh)
  a20 = kb*teff/(mu*mh)
  !Trying now to integrate outside in
  !with overlying wind condition, using grey T-structure
  mdot = 1.d-6 * msun/year
  vinf = 2000.*1d5
  tt = kap0*mdot/(4.*pi*vinf*rstar)
  !NOTE: If you want different connecting point for
  !setting upper boundary, change velocity here
  v0 = 0.5*sqrt(a20)
  !Will assume beta=1 law for now
  bb = 1.-v0/vinf
  mc0 = - tt*log(1.-bb)/bb
  tau0 = - tt*kap0*(bb*(2.+bb)+2.*log(1.-bb))/(2.*bb**3)
  tr0 = teff*(0.5+0.75*tau0)**0.25
  rho0 = mdot/(4.*pi*rstar**2*v0)
  pg0 = rho0*tr0*kb/(mu*mh)
  y0 = rstar !0.0
  tlim = 0.1 !minimum temp.
  er0 = ar*tr0**4.
  tg0 = tr0
  pr0 = er0/3.
  !
  hp = a20/geff
  frac = 0.3 !0.03
  dy = hp * frac
  !For now, ny just maximum number of integration points
  ny = 10000

  !allocate basic variables
  allocate(pg(ny),tg(ny),er(ny),rho(ny),fr(ny),pr(ny),&
     gammar(ny),tr(ny),tau(ny),mc(ny),yarr(ny),kap(ny),t_ana(ny))
  pg(1) = pg0
  tg(1) = tg0
  er(1) = er0
  rho(1) = rho0
  pr(1) = pr0
  tr(1) = tr0
  yarr(1) = y0
  tau(1) = tau0
  mc(1) = mc0

end subroutine set_boundary

subroutine calc_structure

use global
implicit none
external derivs_hydro_diff
external derivs_hydro_diff_mc
external rkqc_hydro

double precision :: x1,x2,h1,eps,xx,y1,y2,yy
integer :: nok,nbad,i,j,narr
double precision, dimension(4) :: ystart_mc
double precision, dimension(3) :: ystart

double precision, dimension(ny) :: derdtau
double precision :: kap_dum,rho_dum,dm0

double precision, dimension(:), allocatable :: rho_vac,pg_vac,er_vac,&
       tg_vac,y_vac,tau_vac

nok = 0
nbad = 0
eps = 1.d-5
!Runge-Kutta integrate through simulation space

rhob = 0.5d0*1.d-9 !rho0
n = 1.5d0
sigma_b = 2.d-2
afw = 7.500317d30
bfw = -4.5854081
!parameters for opacity laws
!the kramer fastwind-like ones are set from a
!ZetPup like converged full CMF hydro model

x1 = mc0
ystart_mc =[pg0,tr0,tau0,y0]

dm = mc0/10.
do i=2,ny

   !Need to set dm here !
   dm = dm*1.005
   x2 = x1 + dm
   h1 = (x2-x1)/10.
   print*,i-1,rho(i-1),tg(i-1),yarr(i-1)/rsun,tau(i-1)
   call ODEINT_HYDRO(YSTART_mc,4,X1,X2,EPS,H1,0.0,NOK,NBAD,&
        derivs_hydro_diff_mc,RKQC_HYDRO)
   pg(i) = ystart_mc(1)
   tr(i) = ystart_mc(2)
   tau(i) = ystart_mc(3)
   yarr(i) = ystart_mc(4)
   er(i) = ar*tr(i)**4.
   tg(i)=tr(i)
   rho(i) = pg(i)*mu*mh/(kb*tg(i))
   mc(i) = x2
!   if (tg(i).gt.200000. .and. tau(i).gt.100.) then
!   if (tg(i).gt.200000. .or. tau(i).gt.100.) then
   if (tau(i).gt.tau_max) then
      print*,i,rho(i),tg(i),yarr(i)/rsun,tau(i)
      exit
   endif
   x1 = x2
enddo
narr = i

print*,'opacity law:  ',kap_law
if (kap_law.eq.'constant') then
   gammar=gamma0
else if (kap_law.eq.'kramerfw') then
   gammar = gamma0*(1.+afw*rho*tr**(bfw))
else if (kap_law.eq.'rho_bump') then
   gammar = gamma0*(1. + n*exp(-1./sigma_b*(log(rho/rhob))**2.))
else
   print*,'not implemented opacity law !'
   stop
endif

print*,'lower boundary:'
print*,'tau,mc=',tau(narr),mc(narr)
print*,'temp,rho=',tr(narr),rho(narr)
!Test vs. pp-version of temperature
!Need modified optical dept hscale for this test !
t_ana = 0.75*teff**4.*(tau-tau0) + tr0**4.
t_ana = t_ana**0.25
x2 = 0.0
do i=1,narr
   x1 = abs(tr(i)/t_ana(i)-1.)
   x2 = max(x1,x2)
enddo
print*,'max err t_sph vs. t_pp =',x2
!Flux comservation test
call deriv_3p(derdtau(1:narr),tau(1:narr),er(1:narr),narr)
derdtau = derdtau*yarr**2.
print*,'max err (sph) flux conservation:',abs(minval(derdtau(1:narr))/maxval(derdtau(1:narr))-1.)
print*,'---------------'
!Scale height at lower boundary
x1 = kb*tg(narr)/(mu*mh)
x2 = gnewt*mstar/yarr(narr)**2 * (1.-gammar(narr))
hpl = x1/x2

  !print some basic stuff to output param file
  open(1,file='init_params')
  write(1,*) 'mass              =',mstar/msun
  write(1,*) 'lum               =',lstar/lsun
  write(1,*) 'gamma_base        =',gamma0
  write(1,*) 'mdot, vinf        =',mdot/msun*year,vinf/1.d5
  write(1,*) 'r_uppperb         =',rstar/rsun
  write(1,*) 'r_lowerb          =',yarr(narr)/rsun
  write(1,*) 'log g_eff(r_up)   =',log10(g*(1.-gamma0))
  write(1,*) 'log g(r_up)       =',log10(g)
  write(1,*) 'log g(r_low)      =',log10(gnewt*mstar/yarr(narr)**2)
  write(1,*) 'Tg0=Tr0           =',tg0
  write(1,*) 'T_lowerb          =',tg(narr)
  write(1,*) 'Teff(upperb)      =',teff
  write(1,*) 'hp(upb)           =',hp/rstar
  write(1,*) 'hp(lowb)          =',hpl/yarr(narr)
  write(1,*) 'mu                =',mu
  write(1,*) 'kap0              =',kap0
  write(1,*) 'tau0              =',tau0
  write(1,*) 'tau_max           =',tau(narr)
  close(1)

open(1,file='init_struc_total')
do i=1,narr
   write(1,11) rho(i),tr(i),pg(i),tau(i),mc(i),gammar(i),er(i),yarr(i)
enddo
11 FORMAT(8d18.8)
close(1)

!Finally interpolate structure upon amrvac grid
allocate(rho_vac(ny_vac),pg_vac(ny_vac),er_vac(ny_vac),&
     tg_vac(ny_vac),y_vac(ny_vac),tau_vac(ny_vac))

y_vac(1) = yarr(narr)
rho_vac(1) = rho(narr)
pg_vac(1) = pg(narr)
er_vac(1) = er(narr)
tg_vac(1) = tg(narr)
tau_vac(1) = tau(narr)

do i=2,ny_vac
   y_vac(i) = y_vac(i-1) + dy_vac*hpl
   xx = y_vac(i)
   if (y_vac(i).gt.yarr(1)) then
      print*,'range in hydro sim too big, reduce dy_vac or ny_vac'
      print*,ny_vac,dy_vac,y_vac(i)/yarr(narr)
      stop
   endif
   !
   do j=2,narr
      if (yarr(j).lt.xx .and. yarr(j-1).ge.xx) exit
   enddo
   !
   x1 = yarr(j)
   x2 = yarr(j-1)
   !rho
   y1 = rho(j)
   y2 = rho(j-1)
   call int_log_or_lin(xx,x1,x2,y1,y2,yy)
   rho_vac(i) = yy
   !pg
   y1 = pg(j)
   y2 = pg(j-1)
   call int_log_or_lin(xx,x1,x2,y1,y2,yy)
   pg_vac(i) = yy
   !er
   y1 = er(j)
   y2 = er(j-1)
   call int_log_or_lin(xx,x1,x2,y1,y2,yy)
   er_vac(i) = yy
   !tg
   y1 = tg(j)
   y2 = tg(j-1)
   call int_log_or_lin(xx,x1,x2,y1,y2,yy)
   tg_vac(i) = yy
   !tau
   y1 = tau(j)
   y2 = tau(j-1)
   call int_log_or_lin(xx,x1,x2,y1,y2,yy)
   tau_vac(i) = yy
   !
enddo

open(1,file='init_struc_amrvac')
do i=1,ny_vac
   write(1,12) i, y_vac(i),rho_vac(i),tg_vac(i),pg_vac(i),er_vac(i),tau_vac(i)
enddo
12 FORMAT(1i4,1e20.10,5e18.8)
close(1)
!
open(1,file='init_params_amrvac')
write(1,*) 'mass',mstar/msun
write(1,*) 'lum',lstar/lsun
write(1,*) 'gamma_base',gamma0
write(1,*) 'r_core=rc',y_vac(1)
write(1,*) 'r_up_b',y_vac(ny_vac)
write(1,*) 'rc,r_up:rsun',y_vac(1)/rsun,y_vac(ny_vac)/rsun
write(1,*) 'log_g(rc)',log10(gnewt*mstar/y_vac(1)**2)
write(1,*) 'Tgc=Trc',tg_vac(1)
write(1,*) 'rhoc',rho_vac(1)
write(1,*) 'hp(rc):rc',hpl/y_vac(1)
write(1,*) 'mu',mu
write(1,*) 'kap0',kap0
write(1,*) 'tauc',tau_vac(1)
write(1,*) 'tau_up',tau_vac(ny_vac)
close(1)

end subroutine calc_structure



!--------------------------------------
!--------------------------------------
!--------------------------------------
!--------------------------------------
!--------------------------------------
!--------------------------------------
!various subroutines below; th first one
!is the one computing derivatives for
!Runge-Kutta ODE solver


subroutine derivs_hydro_diff_mc(x,y,dydx)

!NOTE; HERE

! Y(1) --- gas pressure
! Y(2) --- radiation temperature
! Y(3) --- optical depth proxy
! Y(4) --- yarr

!AND X IS double precisionLY MASS_COLUMN MC !

  use global
  implicit none
  double precision :: x
  double precision, dimension(4) :: y,dydx
  double precision :: gam, temp, rho_dum,kap_dum,dwdr,dtaudr,dtdr

  if (Y(2).le.tlim*teff) Y(2)=tlim*teff
  !Potentially jump start to tmin
  !(never actually used in this set-up, since we're
  !not computing outer wind, where this might happen)
  rho_dum = y(1)*mu*mh/(kb*y(2))
  !density dummy variable

  !Set potential variation of opacity (and thus gamma)
  if (kap_law.eq.'constant') then
     kap_dum = kap0
  else if (kap_law.eq.'kramerfw') then
     kap_dum = kap0*(1.+afw*rho_dum*temp**(bfw))
  else if (kap_law.eq.'rho_bump') then
     kap_dum = kap0*(1. + n*exp(-1./sigma_b*(log(rho_dum/rhob))**2.))
  else
     print*,'not implemented opacity law !'
     stop
  endif
  gam = gamma0*(kap_dum/kap0)

  !gas pressure gradient
  !(accouting for spherical gravity, take away if wanted)
  dydx(1) = g*(1.-gam) * (rstar/y(4))**2

  !temperature gradient (now spherical fastwind-like)
  !Here we could include variation of Eddington factor
  !3, for e.g. Lucy's grey model, FLD, etc.
  if (Y(2).gt.tlim*teff) then
     dydx(2) = 3.*kap_dum*teff**4./(16.*y(2)**3.) &
          * (rstar/y(4))**2
  else
     dydx(2) = 0.0
  endif

  !tau
  dydx(3) = kap_dum

  !radius = yarr
  if (rho_dum.le.0.0) then
     dydx(4) = 0.0
  else
     dydx(4) = -1./rho_dum
  endif

end subroutine derivs_hydro_diff_mc



!-------------------------------------
!-------------------------------------
!-------------------------------------
!-------------------------------------


SUBROUTINE ODEINT_HYDRO(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS, &
&                  RKQC)
!JS-CHANGE: NMAX = 4 now

USE global, ONLY: KMAX,KOUNT,DXSAV,XP,YP

IMPLICIT NONE
!
!     .. parameters ..
!integer(i4b), parameter :: MAXSTP=10000,NMAX=1
integer, parameter :: MAXSTP=30000,NMAX=4
double precision, parameter :: TWO=2.D0,ZERO=0.D0,TINY=1.D-3
!INTEGER(i4b) ::  MAXSTP,NMAX
!double precision(dp) ::  TWO,ZERO,TINY
!PARAMETER (MAXSTP=10000,NMAX=1,TWO=2.D0,ZERO=0.D0,TINY=1.D-30)
!     ..
!     .. scalar arguments ..
double precision ::  EPS,H1,HMIN,X1,X2
INTEGER ::  NBAD,NOK,NVAR
!     ..
!     .. array arguments ..
double precision ::  YSTART(NVAR)
!     ..
!     .. subroutine arguments ..
EXTERNAL DERIVS,RKQC
!     ..
!     .. local scalars ..
double precision ::  H,HDID,HNEXT,X,XSAV
INTEGER ::  I,NSTP
!     ..
!     .. local arrays ..
double precision ::  DYDX(NMAX),Y(NMAX),YSCAL(NMAX)
!     ..
!     .. intrinsic functions ..
INTRINSIC ABS,SIGN
!     ..

X = X1
H = SIGN(H1,X2-X1)
NOK = 0
NBAD = 0
KOUNT = 0

DO I = 1,NVAR
     Y(I) = YSTART(I)
END DO

XSAV = X - DXSAV*TWO

ITLOOP: DO NSTP = 1,MAXSTP

     CALL DERIVS(X,Y,DYDX)

     DO I = 1,NVAR
          YSCAL(I) = ABS(Y(I)) + ABS(H*DYDX(I)) + TINY
     END DO

     IF (KMAX.GT.0) THEN

          IF (ABS(X-XSAV).GT.ABS(DXSAV)) THEN
               IF (KOUNT.LT.KMAX-1) THEN
                    KOUNT = KOUNT + 1
                    XP(KOUNT) = X
                    DO I = 1,NVAR
                         YP(I,KOUNT) = Y(I)
                    END DO

                    XSAV = X
               END IF
          END IF
     END IF

     IF ((X+H-X2)* (X+H-X1).GT.ZERO) H = X2 - X

     CALL RKQC(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
     !JS-NEW : BUG ?? (calling old rkqc routine...)
     !JS-NOTE: No bug, since determined from calling of odeint_hydro
     !(like derivs)

     IF (HDID.EQ.H) THEN
          NOK = NOK + 1
     ELSE
          NBAD = NBAD + 1
     END IF

     IF ((X-X2)* (X2-X1).GE.ZERO) THEN
          DO I = 1,NVAR
               YSTART(I) = Y(I)
          END DO

          IF (KMAX.NE.0) THEN
               KOUNT = KOUNT + 1
               XP(KOUNT) = X
               DO I = 1,NVAR
                    YP(I,KOUNT) = Y(I)
               END DO
          END IF

          RETURN
     END IF

     !IF (ABS(HNEXT).LT.HMIN) STOP 'STEPSIZE SMALLER THAN MINIMUM.'

     H = HNEXT

END DO ITLOOP

STOP 'TOO MANY STEPS, HYDROW_ODE'

END SUBROUTINE ODEINT_HYDRO


SUBROUTINE RKQC_HYDRO(Y,DYDX,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)

IMPLICIT NONE
!
!     .. parameters ..
double precision, parameter :: FCOR=.0666666667D0,ONE=1.D0,SAFETY=0.9D0, &
&          ERRCON=6.D-4
!JS-NEW, nmax =4
integer, parameter :: NMAX=4
!INTEGER(i4b) ::  NMAX
!double precision(dp) ::  FCOR,ONE,SAFETY,ERRCON
!PARAMETER (NMAX=1,FCOR=.0666666667D0,ONE=1.D0,SAFETY=0.9D0, &
!&          ERRCON=6.D-4)
!     ..
!     .. scalar arguments ..
double precision ::  EPS,HDID,HNEXT,HTRY,X
INTEGER ::  N
!     ..
!     .. array arguments ..
double precision ::  DYDX(N),Y(N),YSCAL(N)
!     ..
!     .. subroutine arguments ..
EXTERNAL DERIVS
!     ..
!     .. local scalars ..
double precision ::  ERRMAX,H,HH,PGROW,PSHRNK,XSAV
INTEGER ::  I
!     ..
!     .. local arrays ..
double precision ::  DYSAV(NMAX),YSAV(NMAX),YTEMP(NMAX)
!     ..
!     .. external subroutines ..
EXTERNAL RK4_HYDRO
!     ..
!     .. intrinsic functions ..
INTRINSIC ABS,MAX
!     ..

PGROW = -0.20D0
PSHRNK = -0.25D0
XSAV = X

DO I = 1,N
     YSAV(I) = Y(I)
     DYSAV(I) = DYDX(I)
END DO

H = HTRY

   20 CONTINUE

HH = 0.5D0*H

CALL RK4_HYDRO(YSAV,DYSAV,N,XSAV,HH,YTEMP,DERIVS)

X = XSAV + HH

CALL DERIVS(X,YTEMP,DYDX)
CALL RK4_HYDRO(YTEMP,DYDX,N,X,HH,Y,DERIVS)

X = XSAV + H

!IF (X.EQ.XSAV) PRINT *,'STEPSIZE NOT SIGNIFICANT IN RKQC.'
IF (X.EQ.XSAV) THEN
!IF (ABS(X-XSAV) THEN
   PRINT *,'STEPSIZE NOT SIGNIFICANT IN RKQC.'
   PRINT*,X,XSAV,H
   PRINT*,'-------------------'
ENDIF

CALL RK4_HYDRO(YSAV,DYSAV,N,XSAV,H,YTEMP,DERIVS)
ERRMAX = 0.D0

DO I = 1,N
     YTEMP(I) = Y(I) - YTEMP(I)
     ERRMAX = MAX(ERRMAX,ABS(YTEMP(I)/YSCAL(I)))
END DO

ERRMAX = ERRMAX/EPS

IF (ERRMAX.GT.ONE) THEN
     H = SAFETY*H* (ERRMAX**PSHRNK)
     GO TO 20
ELSE
     HDID = H
     IF (ERRMAX.GT.ERRCON) THEN
          HNEXT = SAFETY*H* (ERRMAX**PGROW)
     ELSE
          HNEXT = 4.D0*H
     END IF
END IF

DO I = 1,N
     Y(I) = Y(I) + YTEMP(I)*FCOR
END DO

RETURN
END SUBROUTINE RKQC_HYDRO


SUBROUTINE RK4_HYDRO(Y,DYDX,N,X,H,YOUT,DERIVS)

IMPLICIT NONE
!
!     .. parameters ..
INTEGER, parameter ::  NMAX=4
!PARAMETER (NMAX=4)
!     ..
!     .. scalar arguments ..
double precision ::  H,X
INTEGER ::  N
!     ..
!     .. array arguments ..
double precision ::  DYDX(N),Y(N),YOUT(N)
!     ..
!     .. subroutine arguments ..
EXTERNAL DERIVS
!     ..
!     .. local scalars ..
double precision ::  H6,HH,XH
INTEGER ::  I
!     ..
!     .. local arrays ..
double precision ::  DYM(NMAX),DYT(NMAX),YT(NMAX)
!     ..

HH = H*0.5D0
H6 = H/6.D0
XH = X + HH

DO I = 1,N
     YT(I) = Y(I) + HH*DYDX(I)
END DO

CALL DERIVS(XH,YT,DYT)

DO I = 1,N
     YT(I) = Y(I) + HH*DYT(I)
END DO

CALL DERIVS(XH,YT,DYM)


DO I = 1,N
     YT(I) = Y(I) + H*DYM(I)
     DYM(I) = DYT(I) + DYM(I)
END DO


CALL DERIVS(X+H,YT,DYT)


DO I = 1,N
     YOUT(I) = Y(I) + H6* (DYDX(I)+DYT(I)+2.D0*DYM(I))
END DO

RETURN

END SUBROUTINE RK4_HYDRO

SUBROUTINE DERIV_3P(D,X,Y,N)

! Perform numerical differentiation using 3-point, Lagrangian interpolation.
! df/dx = y0*(2x-x1-x2)/(x01*x02)+y1*(2x-x0-x2)/(x10*x12)+y2*(2x-x0-x1)/(x20*x21)
! Where: x01 = x0-x1, x02 = x0-x2, x12 = x1-x2, etc.

! taken from IDL; note that cshift (f90) and shift(idl) use different sign conventions


IMPLICIT NONE

INTEGER :: N, N2
double precision, DIMENSION (N) :: X, Y, D

double precision, DIMENSION(N) :: X12, X01, X02

IF (N.LT.3) STOP 'N LT 3 IN DERIV_3P'

X12 = X - CSHIFT(X,1)     !x1 - x2
X01 = CSHIFT(X,-1) - X    !x0 - x1
X02 = CSHIFT(X,-1) - CSHIFT(X,1) !x0 - x2

!Middle points
D = CSHIFT(Y,-1) * (X12 / (X01*X02)) + Y * (1./X12 - 1./X01) - &
&   CSHIFT(Y,1) * (X01 / (X02 * X12))

! Formulae for the first and last points:
D(1) =  Y(1) * (X01(2)+X02(2))/(X01(2)*X02(2)) - &
&       Y(2) * X02(2)/(X01(2)*X12(2)) + &
&       Y(3) * X01(2)/(X02(2)*X12(2))

N2 = N-1
D(N) = -Y(N-2) * X12(N2)/(X01(N2)*X02(N2)) + &
&       Y(N-1) * X02(N2)/(X01(N2)*X12(N2)) - &
&         Y(N) * (X02(N2)+X12(N2)) / (X02(N2)*X12(N2))

RETURN
END SUBROUTINE DERIV_3P

subroutine int_log_or_lin(xx,x1,x2,y1,y2,yy)

!Performs logarithmic inteporlation if allowed,
!linear otherwise

implicit none

double precision, intent(in) :: xx,x1,x2,y1,y2
!double precision, intent(inout) :: yy
double precision :: yy
!
double precision :: x_int,xx1,xx2,yy1,yy2,factor,ydum,yerr1,yerr2

if (y1.gt.0.0d0 .and. y2.gt.0.0d0 .and. x1.gt.0.0d0 .and. x2.gt.0.0d0) then
   !
   x_int = log(xx)
   xx1 = log(x1)
   xx2 = log(x2)
   !
   yy1 = log(y1)
   yy2 = log(y2)
   factor = (x_int-xx1)/(xx2-xx1)
   ydum = yy1 + (yy2-yy1) * factor
   yy = exp(ydum)
   !
else
   x_int = xx
   xx1 = x1
   xx2 = x2
   !
   yy1 = y1
   yy2 = y2
   factor = (x_int-xx1)/(xx2-xx1)
   ydum = yy1 + (yy2-yy1) * factor
   yy = ydum
   !
endif
!
if (yy.gt.max(y1,y2).or.yy.lt.min(y1,y2)) then
   yerr1 = 0.d0
   if (yy.gt.max(y1,y2)) yerr1 = abs(1.-max(y1,y2)/yy)
   if (yy.lt.min(y1,y2)) yerr1 = abs(1.-min(y1,y2)/yy)
   if (yerr1.gt.1.d-5) then
      print*,'WARNING, WARNING -- extrapolation in int_log_or_lin'
      print*,xx,x1,x2,y1,y2,yy
      print*,yerr1
      print*,'-------------------'
   endif
endif

Return
end subroutine int_log_or_lin
