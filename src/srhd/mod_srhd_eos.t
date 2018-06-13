!###########################################################################
! module amrvacphys - srhdeos  
! This module is developed using the paper Meliani et al 2004

!*DM* I assume it means the "Relativistic Parker winds with 
!*DM* variable effective polytropic index"

!*DM  The srhd_eos uses the Mathews equation of state
!*DM* The srhd uses an ideal law
!===========================================================================
module mod_srhd_eos

contains

subroutine srhd_enthalpy(w,ixI^L,ixO^L,patchw,rhoh)

!================== IMPORTANT ==================!
!This subroutine is used only with primitive variables on input w
!===============================================!

use mod_global_parameters

integer:: ixI^L,ixO^L
double precision, dimension(ixI^S,1:nw):: w
double precision, dimension(ixG^T)  :: E_Th,E,rhoh
logical,          dimension(ixG^T)  :: patchw
!--------------------------------------! 

where(.not.patchw(ixO^S))
   !== thermal energy in polytropic definition ==!
   E_Th(ixO^S) = w(ixO^S,pp_)/(srhd_gamma-1.0d0) 
   ! internal energy
   E(ixO^S) = (E_Th(ixO^S) + dsqrt(E_th(ixO^S)**2.0d0+w(ixO^S,rho_)**2.0d0))
   ! puposely writing rho/E/rho instead of rho^2/E for numerics
   rhoh(ixO^S) = 0.5d0*((srhd_gamma+1.0d0)*E(ixO^S)-&
	    (srhd_gamma-1.0d0)*w(ixO^S,rho_)*(w(ixO^S,rho_)/E(ixO^S)))
end where

return
end subroutine srhd_enthalpy
!===========================================================================
subroutine Einternal(w,ixI^L,ixO^L,patchw,varconserve,E)

!================== IMPORTANT ==================!
! if varconserve=.true.
!The subroutine is used only with conserve variables on input w
! if varconserve=.false.
!The subroutine is used only with primitive variables on input w
!===============================================!

use mod_global_parameters

integer:: ixI^L,ixO^L
double precision, dimension(ixI^S,1:nw):: w
logical, intent(in)                 :: varconserve
double precision, dimension(ixG^T)  :: E_Th,E
logical,          dimension(ixG^T)  :: patchw
!--------------------------------------!
if (varconserve) then
where(.not.patchw(ixO^S))
   !== thermal energy in polytropic definition ==!
   E_Th(ixO^S) = w(ixO^S,p_)/(srhd_gamma-1.0d0)
   ! internal energy
   E(ixO^S) = (E_Th(ixO^S) + dsqrt(E_th(ixO^S)**2.0d0+(w(ixO^S,d_)/w(ixO^S,lfac_))**2.0d0))&
                     -(w(ixO^S,d_)/w(ixO^S,lfac_))
end where
else
   !== thermal energy in polytropic definition ==!
   E_Th(ixO^S) = w(ixO^S,pp_)/(srhd_gamma-1.0d0)
   ! internal energy
   E(ixO^S) = (E_Th(ixO^S) + dsqrt(E_th(ixO^S)**2.0d0+w(ixO^S,rho_)**2.0d0))&
                     -w(ixO^S,rho_)

end if

return
end subroutine Einternal
!=============================================================================
subroutine getcsound2(w,ixI^L,ixO^L,vacconserve,rhoh,csound2)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!===============================================!

use mod_global_parameters

integer:: ixI^L,ixO^L
double precision                      :: w(ixI^S,1:nw)
logical, intent(in)                   :: vacconserve
double precision, dimension(ixG^T)    :: rho
double precision, dimension(ixG^T)    :: E_th,E
double precision, dimension(ixG^T)    :: rhoh,csound2
!_______________________________________________________________!
if (vacconserve) then
  !== the thermal energy in polytropic definition ==!
  E_th(ixO^S) = (1.0d0/(srhd_gamma-1.0d0) * w(ixO^S,p_))
  rho(ixO^S) = (w(ixO^S,d_)/w(ixO^S,lfac_))
  E(ixO^S) = (E_th(ixO^S) + dsqrt(E_th(ixO^S)**2.0d0+rho(ixO^S)**2.0d0))
  rhoh(ixO^S) =  0.5d0*((srhd_gamma+1.0d0)*E(ixO^S)-&
   	    (srhd_gamma-1.0d0)*rho(ixO^S)*(rho(ixO^S)/E(ixO^S)))

  !====== The sound speed ======!
  csound2(ixO^S)=0.5d0*w(ixO^S,p_)/rhoh(ixO^S)&
                 *((srhd_gamma+1.0d0)&
                 +(srhd_gamma-1.0d0)*(rho(ixO^S)/E(ixO^S))**2.0d0)
else
  !== the thermal energy in polytropic definition ==!
  E_th(ixO^S) = (1.0d0/(srhd_gamma-1.0d0) * w(ixO^S,p_))
  E(ixO^S) = (E_th(ixO^S) + dsqrt(E_th(ixO^S)**2.0d0+w(ixO^S,rho_)**2.0d0))
  rhoh(ixO^S) =  0.5d0*((srhd_gamma+1.0d0)*E(ixO^S)-&
   	    (srhd_gamma-1.0d0)*w(ixO^S,rho_)*(w(ixO^S,rho_)/E(ixO^S)))

  !====== The sound speed ======!
  csound2(ixO^S)=0.5d0*w(ixO^S,p_)/rhoh(ixO^S)&
                 *((srhd_gamma+1.0d0)&
                 +(srhd_gamma-1.0d0)*(w(ixO^S,rho_)/E(ixO^S))**2.0d0)
end if
end subroutine getcsound2
!=============================================================================
subroutine FuncEnthalpy(pcurrent,lfac2inv,d,s^C,tau,sqrs,xicurrent,dv2d2p,h,dhdp,ierror)

use mod_global_parameters
  
integer:: ixI^L,ixO^L
double precision:: pcurrent,lfac2inv,d,s^C,tau,sqrs,xicurrent
double precision:: h,dhdp
integer::ierror

double precision:: rho
double precision:: E_th,E
double precision:: dv2d2p,dE_thdp,dEdp
!_______________________________________________________________!
rho=d*sqrt(lfac2inv)
!== The Classical definition of the thermal energy ==!
E_th = pcurrent/(srhd_gamma-1.0d0)
E = (E_th + dsqrt(E_th**2.0+rho**2.0d0))
!== Enthalpy ==!
h = 0.5d0 *((srhd_gamma+1.0d0)*E-(srhd_gamma-1.0d0)*rho/(E/rho))

!=== Derivative of thermal energy ===!
dE_thdp = (1.0d0/(srhd_gamma-1.0d0))

!=== Derivative of internal energy ===!
dEdp = dE_thdp * (1.0d0+E_th/dsqrt(E_th**2.0d0+rho**2.0d0))&
	+  d**2.0d0*dv2d2p/dsqrt(E_th**2.0d0+rho**2.0d0)

!====== Derivative of Enthalpy ======!
dhdp = 0.5d0*((srhd_gamma+1.0d0)*dEdp + &
 (srhd_gamma-1.0d0)*(rho**2.0d0/E)*(-2.0d0*dv2d2p/lfac2inv+dEdp/E)) 

return
end subroutine FuncEnthalpy
!=============================================================================!
subroutine smallvaluesEOS

use mod_global_parameters

double precision::Lsmallrho,Lsmallp,LsmallE
!_______________________________________________________________!
Lsmallrho=(1.0d0 + 10.0d0 * minrho) * minrho
Lsmallp=(1.0d0 + 10.0d0 * minp) * minp
LsmallE=Lsmallp/(srhd_gamma-1.0d0)+&
	dsqrt((Lsmallp/(srhd_gamma-1.0d0))**2.0d0+Lsmallrho**2.0d0)

smalltau=0.5d0*((srhd_gamma+1.0d0)*LsmallE-&
    (srhd_gamma-1.0d0)*Lsmallrho**2.0d0/LsmallE)-Lsmallp-Lsmallrho
! may need to replace by smallxi above

smallxi=0.5d0*((srhd_gamma+1.0d0)*LsmallE-&
   (srhd_gamma-1.0d0)*Lsmallrho**2.0d0/LsmallE)

end subroutine smallvaluesEOS
!=============================================================================!
subroutine Bisection_Enthalpy(pcurrent,lfac2inv,d,s^C,tau,sqrs,xicurrent,h,ierror)

use mod_global_parameters
  
integer:: ixI^L,ixO^L
double precision:: pcurrent,lfac2inv,d,s^C,tau,sqrs,xicurrent
double precision:: h
integer::ierror

double precision:: rho
double precision:: E_th,E
!_______________________________________________________________!
rho=d*sqrt(lfac2inv)
E_th = (1.0d0/(srhd_gamma-1.0d0) * pcurrent)
E = (E_th + dsqrt(E_th**2.0d0+rho**2.0d0))
!== Enthalpy ==!
h = 0.5d0 *((srhd_gamma+1.0d0) * E-(srhd_gamma-1.0d0) * rho*(rho/E))

return
end subroutine Bisection_Enthalpy
!=============================================================================
subroutine pressureNoFlow(pressure,tau,d)

use mod_global_parameters

double precision:: pressure,tau,d
!_______________________________________________________________!

pressure=0.5d0*(srhd_gamma-1.0d0)*(tau+d-d**2/(tau+d))

end subroutine pressureNoFlow
!=============================================================================
subroutine Calcule_Geffect(w,ixI^L,ixO^L,varconserve,Geff)

!================== IMPORTANT ==================!
!This subroutine is used with conserved variables in w when varconserve=T
!This subroutine is used with primitive variables in w when varconserve=F
!   both cases assume updated auxiliary variables p_ en lfac_
!===============================================!
use mod_global_parameters

integer:: ixI^L,ixO^L
double precision:: w(ixI^S,1:nw)
logical,intent(in)         :: varconserve
double precision, dimension(ixG^T)    :: Geff
!-----------------------------------
! assume we have the pressure in w(ixO^S,p_)
! and the Lorentz factor in lfac_ and the conserved variables
if (varconserve) then
  Geff(ixO^S) = srhd_gamma-0.5d0*(srhd_gamma-1.0d0) *          &
       (1.0d0-((w(ixO^S,d_)/w(ixO^S,lfac_))/(w(ixO^S, p_)/(srhd_gamma-1.0d0)+  &
       dsqrt((w(ixO^S,p_)/(srhd_gamma-1.0d0))**2.0d0+&
       (w(ixO^S,d_)/w(ixO^S,lfac_))**2.0d0)))**2.0d0)
else
  Geff(ixO^S) = srhd_gamma-0.5d0*(srhd_gamma-1.0d0) *          &
       (1.0d0-(w(ixO^S,rho_)/(w(ixO^S, p_)/(srhd_gamma-1.0d0)+  &
       dsqrt((w(ixO^S,p_)/(srhd_gamma-1.0d0))**2.0d0+&
       w(ixO^S,rho_)**2.0d0)))**2.0d0)
end if
end subroutine Calcule_Geffect

end module mod_srhd_eos
