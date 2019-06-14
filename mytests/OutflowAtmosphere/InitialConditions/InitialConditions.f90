!program to compute initial condition for time-dependent
!(AMRVAC) 2D FLD simulations of a radiation-dominated,
!massive star atmosphere with a photon-tired
!wind outflow
!
!NOTE: THESE INITIAL CONDITIONS ASSUME PURE
!DIFFUSION LIMIT FROM OWOCKI+ 2017
!
!--> THUS: For sims with gamma0*m >> 1
!the relaxed state will BE DIFFERENT
!(since initial conditions neglect radiative
!enthalpy term)


module global

!Basic constants
  double precision, parameter :: gnewt=6.6740000d-08,rsun=6.9599000d10,&
       lsun=3.8268000d33,msun=1.9891000d33,sigb=5.6704000d-05,&
       cc=2.9979246d10,year=31557000.d0,pi=3.1415926536d0,&
       ar=4.*sigb/cc,kb=1.3806503d-16,mh=1.6735344d-24,&
       sigth=6.6524585d-25

!Basic variables for structure
double precision, dimension(:), allocatable :: pg,tg,er,rho,fr,pr,&
     gammar,tr,tau,mc,yarr,kap,t_ana,w,v,pscale,rarr
double precision :: pg0,tg0,er0,rho0,fr0,pr0,gamma0,gammab,tr0,tau0,mc0,y0,kap0,lum
double precision :: mu,mstar,lstar,g,rstar,dy,teff,yhe,dm,tr_dum,vesc
double precision :: rhob,sigma_b,n,afw,bfw,tlim,hp,hpl,mdot,vinf,dy_vac,tau_max,mm
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

  !----------------------
  open(1,file='init_indat_tired')
  read(1,*) gamma0,mm,rstar,vesc,kap0
  close(1)
  rstar = rstar*rsun
  vesc = vesc * 1.d5
  mstar = vesc**2.*0.5*rstar/gnewt
  !This will set the scaling radius !
  lum = gamma0*4.*pi*gnewt*mstar*cc/kap0
  print*,'Input params:'
  print*,'tiring param m=',mm
  print*,'gamma0 > 1 =',gamma0
  print*,'Scaling params:'
  print*,'mass,lum,rad=',mstar/msun,lum/lsun,rstar/rsun
  mdot = mm*rstar*lum/(gnewt*mstar)
  print*,'mdot=',mdot * year/msun
  print*,'-------- NOW CALC STARTS !'
  ny = 1000
  !allocate basic variables
  allocate(er(ny),rho(ny),yarr(ny),&
       w(ny),v(ny),pscale(ny),rarr(ny))

end subroutine set_boundary

subroutine calc_structure

use global
implicit none
external derivs_hydro_diff
external derivs_hydro_diff_mc
external derivs_hydro_diff_tired
external rkqc_hydro

double precision :: x1,x2,h1,eps,xx,y1,y2,yy
integer :: nok,nbad,i,j,narr
double precision, dimension(4) :: ystart_mc
double precision, dimension(1) :: ystart

double precision, dimension(ny) :: derdtau
double precision :: kap_dum,rho_dum,dm0,wdum,last_it

double precision, dimension(:), allocatable :: rho_vac,er_vac,y_vac,v_vac

nok = 0
nbad = 0
eps = 1.d-5
!Runge-Kutta integrate through simulation space


x1 = 1.d0 - 1.d0/1000.d0
wdum = (1.-exp(-mm*gamma0*x1) )/mm - x1
ystart = -x1 * ( -(1.-x1)**2/sqrt(wdum) *mm*gamma0*(1.-mm*(wdum+x1)) )
yarr(1) = x1
pscale(1) = ystart(1)
w(1) = wdum

!NOTE: dm dummy for dx, where x = 1-Rstar/r !
dm = 0.002d0
last_it = 0.d0
do i=2,ny
   !Need to set dm here !
   !dm = dm*1.005
   x2 = x1 - dm
   wdum = (1.-exp(-mm*gamma0*x2) )/mm - x2
   if (x2.le.0.0d0) then
      x2 = 0.0d0+dm*0.1d0
      last_it = 1.d0
      if (x1.lt.x2) stop 'x1>x2'
   endif
   h1 = (x2-x1)/10.d0
   call ODEINT_HYDRO(YSTART,1,X1,X2,EPS,H1,0.d0,NOK,NBAD,&
        derivs_hydro_diff_tired,RKQC_HYDRO)
   wdum = (1.d0-dexp(-mm*gamma0*x2) )/mm - x2
!   print*,ystart,x2,wdum
   yarr(i) = x2
   pscale(i) = ystart(1)
   w(i) = wdum
   if (last_it.gt.0.d5) exit
   x1 = x2
enddo
narr = i

do i=1,narr
   v(i) = sqrt(w(i)) * vesc
   rarr(i) = rstar/(1.d0-yarr(i))
   rho(i) = mdot/(4.d0*pi*rarr(i)**2.*v(i))
   er(i) = pscale(i)*lum*3.d0/(4.d0*pi*rstar**2.*vesc)
   print*,rarr(i)/rstar,v(i)/1.d5,rho(i),(er(i)/ar)**0.25d0
enddo

open(1,file='init_struc_tired_scaled')
do i=1,narr
   write(1,22) i,yarr(i),pscale(i),w(i)
enddo
22 FORMAT(i8,3d18.8)
close(1)

open(1,file='init_struc_tired_cgs')
write(1,*)  'index, radius, density, velocity, radiation energy density'
do i=1,narr
   write(1,23) i,rarr(i),rho(i),v(i),er(i)
enddo
23 FORMAT(i8,4d18.8)
close(1)

ny_vac = 256+2*2
dy_vac = (15.d0-1.d0)/ny_vac

!Finally interpolate structure upon radially uniform amrvac grid
allocate(rho_vac(ny_vac),v_vac(ny_vac),er_vac(ny_vac),y_vac(ny_vac))

y_vac(1) = rarr(narr)
rho_vac(1) = rho(narr)
er_vac(1) = er(narr)
v_vac(1)= v(narr)

do i=2,ny_vac
   y_vac(i) = y_vac(i-1) + dy_vac*rstar
   xx = y_vac(i)
   !
   do j=2,narr
      if (rarr(j).lt.xx .and. rarr(j-1).ge.xx) exit
   enddo
   !
   x1 = rarr(j)
   x2 = rarr(j-1)
   !
   y1 = rho(j)
   y2 = rho(j-1)
   call int_log_or_lin(xx,x1,x2,y1,y2,yy)
   rho_vac(i) = yy
   !rho
   y1 = er(j)
   y2 = er(j-1)
   call int_log_or_lin(xx,x1,x2,y1,y2,yy)
   er_vac(i) = yy
   !er
   y1 = v(j)
   y2 = v(j-1)
   call int_log_or_lin(xx,x1,x2,y1,y2,yy)
   v_vac(i) = yy
   !v
enddo

print*,'AMRVAC structure:'
open(1,file='init_struc_amrvac')
write(1,*) 'n, dr =',ny_vac, dy_vac
write(1,*) 'index, radius, velocity, density, radiaiton energy density'
do i=1,ny_vac
   write(1,112) i, y_vac(i),v_vac(i),rho_vac(i),er_vac(i)
   print*,y_vac(i)/rstar,v_vac(i)/1.d5,rho_vac(i),(er_vac(i)/ar)**0.25
enddo
112 FORMAT(i8,4d18.8)
close(1)
!
open(1,file='init_params_amrvac')
write(1,*) 'gamma0',gamma0
write(1,*) 'TiringParamM=',mm
write(1,*) 'mass',mstar/msun
write(1,*) 'lum',lum/lsun
write(1,*) 'rstar',rstar/rsun
write(1,*) 'r_core=rc/rsun',y_vac(1)/rsun
write(1,*) 'MassLossRate=',mdot * year/msun
write(1,*) 'rho0=',rho_vac(1)
close(1)


end subroutine calc_structure



!--------------------------------------
!--------------------------------------
!--------------------------------------
!--------------------------------------
!--------------------------------------
!--------------------------------------
!various subroutines below; the first one
!is the one computing derivatives for
!Runge-Kutta ODE solver


subroutine derivs_hydro_diff_tired(x,y,dydx)

!NOTE; HERE

! Y(1) --- gas pressure like variable

!AND X is 1-rstar/r

  use global
  implicit none
  double precision :: x
  double precision, dimension(4) :: y,dydx
  double precision :: gam, temp, rho_dum, wdum

  wdum = (1.-exp(-mm*gamma0*x) )/mm - x
  if (wdum.lt.0.0d0) stop 'neg w !'

  dydx(1) = - (1.d0-x)**2/dsqrt(wdum) *mm*gamma0 * (1.d0-mm*(wdum+x) )

end subroutine derivs_hydro_diff_tired



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
