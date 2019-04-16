!> Module with finite volume methods for fluxes
module mod_finite_volume
  implicit none
  private

  public :: finite_volume
  public :: hancock
  public :: reconstruct_LR

contains

  !> The non-conservative Hancock predictor for TVDLF
  !>
  !> on entry:
  !> input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBnghostcells
  !> one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2
  !> on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2
  subroutine hancock(qdt,ixI^L,ixO^L,idims^LIM,qtC,sCT,qt,snew,dx^D,x)
    use mod_physics
    use mod_global_parameters
    use mod_source, only: addsource2

    integer, intent(in) :: ixI^L, ixO^L, idims^LIM
    double precision, intent(in) :: qdt, qtC, qt, dx^D, x(ixI^S,1:ndim)
    type(state) :: sCT, snew

    double precision, dimension(ixI^S,1:nw) :: wprim, wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision :: fLC(ixI^S, nwflux), fRC(ixI^S, nwflux)
    double precision :: dxinv(1:ndim)
    integer :: idims, iw, ix^L, hxO^L

    associate(wCT=>sCT%w,wnew=>snew%w)
    ! Expand limits in each idims direction in which fluxes are added
    ix^L=ixO^L;
    do idims= idims^LIM
       ix^L=ix^L^LADDkr(idims,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
         call mpistop("Error in Hancock: Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    ^D&dxinv(^D)=-qdt/dx^D;
    do idims= idims^LIM
       block%iw0=idims
       ! Calculate w_j+g_j/2 and w_j-g_j/2
       ! First copy all variables, then upwind wLC and wRC.
       ! wLC is to the left of ixO, wRC is to the right of wCT.
       hxO^L=ixO^L-kr(idims,^D);

       wRp(hxO^S,1:nwflux)=wprim(ixO^S,1:nwflux)
       wLp(ixO^S,1:nwflux)=wprim(ixO^S,1:nwflux)

       call reconstruct_LR(ixI^L,ixO^L,hxO^L,idims,wprim,wprim,wLC,wRC,wLp,wRp,x,.false.)

       ! Calculate the fLC and fRC fluxes
       call phys_get_flux(wRC,wRp,x,ixI^L,hxO^L,idims,fRC)
       call phys_get_flux(wLC,wLp,x,ixI^L,ixO^L,idims,fLC)

       ! Advect w(iw)
       do iw=1,nwflux
          if (slab) then
             wnew(ixO^S,iw)=wnew(ixO^S,iw)+dxinv(idims)* &
                  (fLC(ixO^S, iw)-fRC(hxO^S, iw))
          else
             wnew(ixO^S,iw)=wnew(ixO^S,iw)-qdt/block%dvolume(ixO^S) &
                  *(block%surfaceC(ixO^S,idims)*fLC(ixO^S, iw) &
                  -block%surfaceC(hxO^S,idims)*fRC(hxO^S, iw))
          end if
       end do
    end do ! next idims
    block%iw0=0

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)
    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,wnew,x,ixI^L,ixO^L,'finite_volume')
    end associate
  end subroutine hancock

  !> finite volume method
  subroutine finite_volume(method,qdt,ixI^L,ixO^L,idims^LIM, &
       qtC,sCT,qt,snew,sold,fC,fE,dx^D,x)

    use mod_physics
    use mod_global_parameters
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2
    use mod_usr_methods
    use mod_constrained_transport

    character(len=*), intent(in)                          :: method
    double precision, intent(in)                          :: qdt, qtC, qt, dx^D
    integer, intent(in)                                   :: ixI^L, ixO^L, idims^LIM
    double precision, dimension(ixI^S,1:ndim), intent(in) ::  x
    type(state)                                           :: sCT, snew, sold
    double precision, dimension(ixI^S,1:nwflux,1:ndim)    :: fC
    double precision, dimension(ixI^S,1:ndir)             :: fE

    ! primitive w at cell center
    double precision, dimension(ixI^S,1:nw) :: wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision, dimension(ixI^S,1:nwflux) :: fLC, fRC
    double precision, dimension(ixI^S)      :: cmaxC
    double precision, dimension(ixI^S)      :: cminC
    double precision, dimension(ixO^S)      :: inv_volume
    double precision, dimension(1:ndim)     :: dxinv
    double precision, dimension(ixI^S,1:ndir,2) :: vbarC                                                                                      
    double precision, dimension(ixI^S,1:ndir,2) :: vbarLC,vbarRC                                                                              
    double precision, dimension(ixI^S,ndim)   :: cbarmin,cbarmax                                                                              
    integer                                 :: idimE,idimN
    integer, dimension(ixI^S)               :: patchf
    integer :: idims, iw, ix^L, hxO^L, ixC^L, ixCR^L, kxC^L, kxR^L

    associate(wCT=>sCT%w, wnew=>snew%w, wold=>sold%w)
    staggered : associate(wCTs=>sCT%ws, wnews=>snew%ws, wolds=>sold%ws)
    if (idimsmax>idimsmin .and. typelimited=='original')&
         call mpistop("Error in fv: Unsplit dim. and original is limited")

    fC=0.d0

    ! The flux calculation contracts by one in the idims direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ix^L=ixO^L;
    do idims= idims^LIM
       ix^L=ix^L^LADD2*kr(idims,^D);
    end do
    if (ixI^L^LTix^L|.or.|.or.) &
         call mpistop("Error in fv : Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)

    ^D&dxinv(^D)=-qdt/dx^D;
    do idims= idims^LIM
       ! use interface value of w0 at idims
       block%iw0=idims

       hxO^L=ixO^L-kr(idims,^D);

       kxCmin^D=ixImin^D; kxCmax^D=ixImax^D-kr(idims,^D);
       kxR^L=kxC^L+kr(idims,^D);

       if(stagger_grid) then
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax^D=ixOmax^D; ixCmin^D=hxOmin^D;
       else
         ixCmax^D=ixOmax^D+nghostcells-nghostcells*kr(idims,^D); ixCmin^D=hxOmin^D-nghostcells+nghostcells*kr(idims,^D);
       end if

       ! wRp and wLp are defined at the same locations, and will correspond to
       ! the left and right reconstructed values at a cell face. Their indexing
       ! is similar to cell-centered values, but in direction idims they are
       ! shifted half a cell towards the 'lower' direction.
       wRp(kxC^S,1:nw)=wprim(kxR^S,1:nw)
       wLp(kxC^S,1:nw)=wprim(kxC^S,1:nw)

       ! Determine stencil size
       {ixCRmin^D = max(ixCmin^D - phys_wider_stencil,ixGlo^D)\}
       {ixCRmax^D = min(ixCmax^D + phys_wider_stencil,ixGhi^D)\}

       ! apply limited reconstruction for left and right status at cell interfaces
       select case (typelimited)
       case ('previous')
         call reconstruct_LR(ixI^L,ixCR^L,ixCR^L,idims,wold,wprim,wLC,wRC,wLp,wRp,x,.true.)
       case ('predictor')
         call reconstruct_LR(ixI^L,ixCR^L,ixCR^L,idims,wprim,wprim,wLC,wRC,wLp,wRp,x,.false.)

         if(stagger_grid) then
           wLC(ixCR^S,iw_mag(idims))=wCTs(ixCR^S,idims)
           wRC(ixCR^S,iw_mag(idims))=wCTs(ixCR^S,idims)
           wLp(ixCR^S,iw_mag(idims))=wCTs(ixCR^S,idims)
           wRp(ixCR^S,iw_mag(idims))=wCTs(ixCR^S,idims)
         end if
       case default
         call mpistop("Error in reconstruction: no such base for limiter")
       end select

       ! special modification of left and right status before flux evaluation
       call phys_modify_wLR(wLp, wRp, ixI^L, ixC^L, idims)

       ! evaluate physical fluxes according to reconstructed status
       call phys_get_flux(wLC,wLp,x,ixI^L,ixC^L,idims,fLC)
       call phys_get_flux(wRC,wRp,x,ixI^L,ixC^L,idims,fRC)

       ! estimating bounds for the minimum and maximum signal velocities
       if(method=='tvdlf'.or.method=='tvdmu') then
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixC^L,idims,cmaxC)
       else
         call phys_get_cbounds(wLC,wRC,wLp,wRp,x,ixI^L,ixC^L,idims,cmaxC,cminC)
       end if

       if(stagger_grid) then
         ! Store magnitude of characteristics
         cbarmin(ixC^S,idims)=max(-cminC(ixC^S),zero)
         cbarmax(ixC^S,idims)=max( cmaxC(ixC^S),zero)

         idimN=mod(idims,ndir)+1 ! 'Next' direction
         idimE=mod(idims+1,ndir)+1 ! Electric field direction
         ! Store velocities
         call phys_get_v_idim(wLp,x,ixI^L,ixC^L,idimN,vbarLC(ixI^S,idims,1))
         call phys_get_v_idim(wRp,x,ixI^L,ixC^L,idimN,vbarRC(ixI^S,idims,1))
         vbarC(ixC^S,idims,1)=(cbarmax(ixC^S,idims)*vbarLC(ixC^S,idims,1) &
              +cbarmin(ixC^S,idims)*vbarRC(ixC^S,idims,1))&
              /(cbarmax(ixC^S,idims) + cbarmin(ixC^S,idims))

         call phys_get_v_idim(wLp,x,ixI^L,ixC^L,idimE,vbarLC(ixI^S,idims,2))
         call phys_get_v_idim(wRp,x,ixI^L,ixC^L,idimE,vbarRC(ixI^S,idims,2))
         vbarC(ixC^S,idims,2)=(cbarmax(ixC^S,idims)*vbarLC(ixC^S,idims,2) &
              +cbarmin(ixC^S,idims)*vbarRC(ixC^S,idims,2))&
              /(cbarmax(ixC^S,idims) + cbarmin(ixC^S,idims))
       end if

       ! use approximate Riemann solver to get flux at interfaces
       select case(method)
       case('tvdmu')
         call get_Riemann_flux_tvdmu()
       case('tvdlf')
         call get_Riemann_flux_tvdlf()
       case('hll')
         call get_Riemann_flux_hll()
       case('hllc','hllcd')
         call get_Riemann_flux_hllc()
       case('hlld')
         call get_Riemann_flux_hlld()
       case default
         call mpistop('unkown Riemann flux')
       end select

       if(associated(usr_set_flux)) call usr_set_flux(ixI^L,ixC^L,idims,fC)

    end do ! Next idims
    block%iw0=0

    if(stagger_grid) call update_faces_uct2(ixI^L,ixO^L,qdt,vbarC,cbarmin,cbarmax,fE,snew)

    do idims= idims^LIM
       hxO^L=ixO^L-kr(idims,^D);

       ! Multiply the fluxes by -dt/dx since Flux fixing expects this
       if (slab) then
          fC(ixI^S,1:nwflux,idims)=dxinv(idims)*fC(ixI^S,1:nwflux,idims)
          wnew(ixO^S,1:nwflux)=wnew(ixO^S,1:nwflux) &
               + (fC(ixO^S,1:nwflux,idims)-fC(hxO^S,1:nwflux,idims))
       else
          fC(ixI^S,1:nwflux,idims)=-qdt*fC(ixI^S,1:nwflux,idims)
          if (.not. angmomfix) then ! default case
            inv_volume = 1.0d0/block%dvolume(ixO^S)
            do iw=1,nwflux
              wnew(ixO^S,iw)=wnew(ixO^S,iw) + (fC(ixO^S,iw,idims)-fC(hxO^S,iw,idims)) * &
                  inv_volume
            enddo
          else
            ! If angular momentum conserving way to solve the equations,
            ! some fluxes additions need to be treated specifically
            call phys_angmomfix(fC,x,wnew,ixI^L,ixO^L,idims)
          endif
       end if

       ! For the MUSCL scheme apply the characteristic based limiter
       if (method=='tvdmu') &
            call tvdlimit2(method,qdt,ixI^L,ixC^L,ixO^L,idims,wLC,wRC,wnew,x,fC,dx^D)

    end do ! Next idims

    if(stagger_grid) call faces2centers(ixO^L,snew)

    if (.not.slab.and.idimsmin==1) &
         call phys_add_source_geom(qdt,ixI^L,ixO^L,wCT,wnew,x)

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), &
         ixI^L,ixO^L,1,nw,qtC,wCT,qt,wnew,x,.false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,wnew,x,ixI^L,ixO^L,'finite_volume')

  end associate staggered
  end associate
  contains

    subroutine get_Riemann_flux_tvdmu()

      do iw=1,nwflux
         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixC^S, iw)=half*(fLC(ixC^S, iw)+fRC(ixC^S, iw))
         if (slab) then
           fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
         else
           fC(ixC^S,iw,idims)=block%surfaceC(ixC^S,idims)*fLC(ixC^S, iw)
         end if
      end do
    end subroutine get_Riemann_flux_tvdmu

    subroutine get_Riemann_flux_tvdlf()
      double precision :: fac(ixC^S)

      fac = -0.5d0*tvdlfeps*cmaxC(ixC^S)

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=1,nwflux

         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixC^S, iw)=0.5d0*(fLC(ixC^S, iw)+fRC(ixC^S, iw))

         ! Add TVDLF dissipation to the flux
         if (flux_type(idims, iw) /= flux_no_dissipation) then
            fLC(ixC^S, iw)=fLC(ixC^S, iw) + fac*(wRC(ixC^S,iw)-wLC(ixC^S,iw))
         end if

         if (slab) then
           fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
         else
           fC(ixC^S,iw,idims)=block%surfaceC(ixC^S,idims)*fLC(ixC^S, iw)
         end if

      end do ! Next iw
    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll()

      double precision :: fac(ixC^S), div(ixC^S)

      where(cminC(ixC^S) >= zero)
        patchf(ixC^S) = -2
      elsewhere(cmaxC(ixC^S) <= zero)
        patchf(ixC^S) =  2
      elsewhere
        patchf(ixC^S) =  1
      endwhere

      fac = tvdlfeps*cminC(ixC^S)*cmaxC(ixC^S)
      div = 1/(cmaxC(ixC^S)-cminC(ixC^S))

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=1,nwflux
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixC^S, iw) = half*(fLC(ixC^S, iw) + fRC(ixC^S, iw) &
                 -tvdlfeps*max(cmaxC(ixC^S), dabs(cminC(ixC^S))) * &
                 (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
         else
            where(patchf(ixC^S)==1)
               ! Add hll dissipation to the flux
               fLC(ixC^S, iw) = (cmaxC(ixC^S)*fLC(ixC^S, iw)-cminC(ixC^S) * fRC(ixC^S, iw) &
                    +fac*(wRC(ixC^S,iw)-wLC(ixC^S,iw))) * div
            elsewhere(patchf(ixC^S)== 2)
               fLC(ixC^S, iw)=fRC(ixC^S, iw)
            elsewhere(patchf(ixC^S)==-2)
               fLC(ixC^S, iw)=fLC(ixC^S, iw)
            endwhere
         endif

         if (slab) then
           fC(ixC^S,iw,idims)=fLC(ixC^S, iw)
         else
           fC(ixC^S,iw,idims)=block%surfaceC(ixC^S,idims)*fLC(ixC^S, iw)
         end if

      end do ! Next iw
    end subroutine get_Riemann_flux_hll

    subroutine get_Riemann_flux_hllc()
      implicit none
      double precision, dimension(ixI^S,1:nwflux)     :: whll, Fhll, fCD
      double precision, dimension(ixI^S)              :: lambdaCD

      patchf(ixC^S) =  1
      where(cminC(ixC^S) >= zero)
         patchf(ixC^S) = -2
      elsewhere(cmaxC(ixC^S) <= zero)
         patchf(ixC^S) =  2
      endwhere
      ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4
      if(method=='hllcd') &
           call phys_diffuse_hllcd(ixI^L,ixC^L,idims,wLC,wRC,fLC,fRC,patchf)

      !---- calculate speed lambda at CD ----!
      if(any(patchf(ixC^S)==1)) &
           call phys_get_lCD(wLC,wRC,fLC,fRC,cminC,cmaxC,idims,ixI^L,ixC^L, &
           whll,Fhll,lambdaCD,patchf)

      ! now patchf may be -1 or 1 due to phys_get_lCD
      if(any(abs(patchf(ixC^S))== 1))then
         !======== flux at intermediate state ========!
         call phys_get_wCD(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,&
              cminC,cmaxC,ixI^L,ixC^L,idims,fCD)
      endif ! Calculate the CD flux

      do iw=1,nwflux
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixC^S,iw) = 0.5d0 * (fLC(ixC^S,iw) + fRC(ixC^S,iw) - tvdlfeps * &
                 max(cmaxC(ixC^S), abs(cminC(ixC^S))) * &
                 (wRC(ixC^S,iw) - wLC(ixC^S,iw)))
         else
            where(patchf(ixC^S)==-2)
               fLC(ixC^S,iw)=fLC(ixC^S,iw)
            elsewhere(abs(patchf(ixC^S))==1)
               fLC(ixC^S,iw)=fCD(ixC^S,iw)
            elsewhere(patchf(ixC^S)==2)
               fLC(ixC^S,iw)=fRC(ixC^S,iw)
            elsewhere(patchf(ixC^S)==3)
               ! fallback option, reducing to HLL flux
               fLC(ixC^S,iw)=Fhll(ixC^S,iw)
            elsewhere(patchf(ixC^S)==4)
               ! fallback option, reducing to TVDLF flux
               fLC(ixC^S,iw) = half*((fLC(ixC^S,iw)+fRC(ixC^S,iw)) &
                    -tvdlfeps * max(cmaxC(ixC^S), dabs(cminC(ixC^S))) * &
                    (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
            endwhere
         end if

         if (slab) then
           fC(ixC^S,iw,idims)=fLC(ixC^S,iw)
         else
           fC(ixC^S,iw,idims)=block%surfaceC(ixC^S,idims)*fLC(ixC^S,iw)
         end if

      end do ! Next iw
    end subroutine get_Riemann_flux_hllc

    !> HLLD Riemann flux from Miyoshi 2005 JCP, 208, 315
    subroutine get_Riemann_flux_hlld()
      use mod_mhd_phys
      implicit none
      double precision, dimension(ixI^S,1:nwflux) :: w1R,w1L,f1R,f1L
      double precision, dimension(ixI^S,1:nwflux) :: w2R,w2L
      double precision, dimension(ixI^S) :: sm,s1R,s1L,suR,suL,Bx
      double precision, dimension(ixI^S) :: pts,ptR,ptL,signBx,r1L,r1R,tmp
      double precision, dimension(ixI^S,ndir) :: vRC, vLC
      integer :: ip1,ip2,ip3,idir

      f1R=0.d0
      f1L=0.d0
      ip1=idims
      ip3=3
      vRC(ixC^S,:)=wRp(ixC^S,mom(:))
      vLC(ixC^S,:)=wLp(ixC^S,mom(:))
      ! estimate normal magnetic field at cell interfaces
      Bx(ixC^S)=0.5d0*(wRC(ixC^S,mag(ip1))+wLC(ixC^S,mag(ip1)))
      suR(ixC^S)=(cmaxC(ixC^S)-vRC(ixC^S,ip1))*wRC(ixC^S,rho_)
      suL(ixC^S)=(cminC(ixC^S)-vLC(ixC^S,ip1))*wLC(ixC^S,rho_)
      ptR(ixC^S)=wRp(ixC^S,e_)+0.5d0*sum(wRC(ixC^S,mag(:))**2,dim=ndim+1)
      ptL(ixC^S)=wLp(ixC^S,e_)+0.5d0*sum(wLC(ixC^S,mag(:))**2,dim=ndim+1)
      ! equation (38)
      sm(ixC^S)=(suR(ixC^S)*vRC(ixC^S,ip1)-suL(ixC^S)*vLC(ixC^S,ip1)-&
                 ptR(ixC^S)+ptL(ixC^S))/(suR(ixC^S)-suL(ixC^S))
      ! equation (39)
      w1R(ixC^S,mom(ip1))=sm(ixC^S)
      w1L(ixC^S,mom(ip1))=sm(ixC^S)
      w2R(ixC^S,mom(ip1))=sm(ixC^S)
      w2L(ixC^S,mom(ip1))=sm(ixC^S)
      w1R(ixC^S,mag(ip1))=Bx(ixC^S)
      w1L(ixC^S,mag(ip1))=Bx(ixC^S)
      w2R(ixC^S,mag(ip1))=Bx(ixC^S)
      w2L(ixC^S,mag(ip1))=Bx(ixC^S)
      ! equation (41)
      pts(ixC^S)=(suR(ixC^S)*ptL(ixC^S)-suL(ixC^S)*ptR(ixC^S)+suR(ixC^S)*suL(ixC^S)*&
                 (vRC(ixC^S,ip1)-vLC(ixC^S,ip1)))/(suR(ixC^S)-suL(ixC^S))
      ! equation (43)
      w1R(ixC^S,rho_)=suR(ixC^S)/(cmaxC(ixC^S)-sm(ixC^S))
      w1L(ixC^S,rho_)=suL(ixC^S)/(cminC(ixC^S)-sm(ixC^S))
      ! equation (44) ~ (47)
      ip2=mod(ip1+1,ndir)
      if(ip2==0) ip2=ndir
      r1R(ixC^S)=suR(ixC^S)*(cmaxC(ixC^S)-sm(ixC^S))-Bx(ixC^S)**2
      where(r1R(ixC^S)/=0.d0)
        r1R(ixC^S)=1.d0/r1R(ixC^S)
      endwhere
      r1L(ixC^S)=suL(ixC^S)*(cminC(ixC^S)-sm(ixC^S))-Bx(ixC^S)**2
      where(r1L(ixC^S)/=0.d0)
        r1L(ixC^S)=1.d0/r1L(ixC^S)
      endwhere
      w1R(ixC^S,mom(ip2))=vRC(ixC^S,ip2)-Bx(ixC^S)*wRC(ixC^S,mag(ip2))*&
        (sm(ixC^S)-vRC(ixC^S,ip1))*r1R(ixC^S)
      w1L(ixC^S,mom(ip2))=vLC(ixC^S,ip2)-Bx(ixC^S)*wLC(ixC^S,mag(ip2))*&
        (sm(ixC^S)-vLC(ixC^S,ip1))*r1L(ixC^S)
      w1R(ixC^S,mag(ip2))=(suR(ixC^S)*(cmaxC(ixC^S)-vRC(ixC^S,ip1))-Bx(ixC^S)**2)*r1R(ixC^S)
      w1L(ixC^S,mag(ip2))=(suL(ixC^S)*(cminC(ixC^S)-vLC(ixC^S,ip1))-Bx(ixC^S)**2)*r1L(ixC^S)
      if(ndir==3) then
        ip3=mod(ip1+2,ndir)
        if(ip3==0) ip3=ndir
        w1R(ixC^S,mom(ip3))=vRC(ixC^S,ip3)-Bx(ixC^S)*wRC(ixC^S,mag(ip3))*&
          (sm(ixC^S)-vRC(ixC^S,ip1))*r1R(ixC^S)
        w1L(ixC^S,mom(ip3))=vLC(ixC^S,ip3)-Bx(ixC^S)*wLC(ixC^S,mag(ip3))*&
          (sm(ixC^S)-vLC(ixC^S,ip1))*r1L(ixC^S)
        w1R(ixC^S,mag(ip3))=wRC(ixC^S,mag(ip3))*w1R(ixC^S,mag(ip2))
        w1L(ixC^S,mag(ip3))=wLC(ixC^S,mag(ip3))*w1L(ixC^S,mag(ip2))
      end if
      w1R(ixC^S,mag(ip2))=wRC(ixC^S,mag(ip2))*w1R(ixC^S,mag(ip2))
      w1L(ixC^S,mag(ip2))=wLC(ixC^S,mag(ip2))*w1L(ixC^S,mag(ip2))
      ! equation (48)
      if(mhd_energy) then
        w1R(ixC^S,e_)=((cmaxC(ixC^S)-vRC(ixC^S,ip1))*wRC(ixC^S,e_)-ptR(ixC^S)*vRC(ixC^S,ip1)+&
          pts(ixC^S)*sm(ixC^S)+Bx(ixC^S)*(sum(vRC(ixC^S,:)*wRC(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w1R(ixC^S,mom(:))*w1R(ixC^S,mag(:)),dim=ndim+1)))/(cmaxC(ixC^S)-sm(ixC^S))
        w1L(ixC^S,e_)=((cminC(ixC^S)-vLC(ixC^S,ip1))*wLC(ixC^S,e_)-ptL(ixC^S)*vLC(ixC^S,ip1)+&
          pts(ixC^S)*sm(ixC^S)+Bx(ixC^S)*(sum(vLC(ixC^S,:)*wLC(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w1L(ixC^S,mom(:))*w1L(ixC^S,mag(:)),dim=ndim+1)))/(cminC(ixC^S)-sm(ixC^S))
      end if
      ! equation (49)
      w2R(ixC^S,rho_)=w1R(ixC^S,rho_)
      w2L(ixC^S,rho_)=w1L(ixC^S,rho_)
      r1R(ixC^S)=sqrt(w1R(ixC^S,rho_))
      r1L(ixC^S)=sqrt(w1L(ixC^S,rho_))
      tmp(ixC^S)=1.d0/(r1R(ixC^S)+r1L(ixC^S))
      signBx(ixC^S)=sign(1.d0,Bx(ixC^S))
      ! equation (51)
      s1R(ixC^S)=sm(ixC^S)+abs(Bx(ixC^S))/r1R(ixC^S)
      s1L(ixC^S)=sm(ixC^S)-abs(Bx(ixC^S))/r1L(ixC^S)
      ! equation (59)
      w2R(ixC^S,mom(ip2))=(r1L(ixC^S)*w1L(ixC^S,mom(ip2))+r1R(ixC^S)*w1R(ixC^S,mom(ip2))+&
          (w1R(ixC^S,mag(ip2))-w1L(ixC^S,mag(ip2)))*signBx(ixC^S))*tmp(ixC^S)
      w2L(ixC^S,mom(ip2))=w2R(ixC^S,mom(ip2))
      ! equation (61)
      w2R(ixC^S,mag(ip2))=(r1L(ixC^S)*w1R(ixC^S,mag(ip2))+r1R(ixC^S)*w1L(ixC^S,mag(ip2))+&
          r1L(ixC^S)*r1R(ixC^S)*(w1R(ixC^S,mom(ip2))-w1L(ixC^S,mom(ip2)))*signBx(ixC^S))*tmp(ixC^S)
      w2L(ixC^S,mag(ip2))=w2R(ixC^S,mag(ip2))
      if(ndir==3) then
        ! equation (60)
        w2R(ixC^S,mom(ip3))=(r1L(ixC^S)*w1L(ixC^S,mom(ip3))+r1R(ixC^S)*w1R(ixC^S,mom(ip3))+&
            (w1R(ixC^S,mag(ip3))-w1L(ixC^S,mag(ip3)))*signBx(ixC^S))*tmp(ixC^S)
        w2L(ixC^S,mom(ip3))=w2R(ixC^S,mom(ip3))
        ! equation (62)
        w2R(ixC^S,mag(ip3))=(r1L(ixC^S)*w1R(ixC^S,mag(ip3))+r1R(ixC^S)*w1L(ixC^S,mag(ip3))+&
            r1L(ixC^S)*r1R(ixC^S)*(w1R(ixC^S,mom(ip3))-w1L(ixC^S,mom(ip3)))*signBx(ixC^S))*tmp(ixC^S)
        w2L(ixC^S,mag(ip3))=w2R(ixC^S,mag(ip3))
      end if
      ! equation (63)
      if(mhd_energy) then
        w2R(ixC^S,e_)=w1R(ixC^S,e_)+r1R(ixC^S)*(sum(w1R(ixC^S,mom(:))*w1R(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w2R(ixC^S,mom(:))*w2R(ixC^S,mag(:)),dim=ndim+1))*signBx(ixC^S)
        w2L(ixC^S,e_)=w1L(ixC^S,e_)-r1L(ixC^S)*(sum(w1L(ixC^S,mom(:))*w1L(ixC^S,mag(:)),dim=ndim+1)-&
          sum(w2L(ixC^S,mom(:))*w2L(ixC^S,mag(:)),dim=ndim+1))*signBx(ixC^S)
      end if
      ! convert velocity to momentum
      do idir=1,ndir
        w1R(ixC^S,mom(idir))=w1R(ixC^S,mom(idir))*w1R(ixC^S,rho_)
        w1L(ixC^S,mom(idir))=w1L(ixC^S,mom(idir))*w1L(ixC^S,rho_)
        w2R(ixC^S,mom(idir))=w2R(ixC^S,mom(idir))*w2R(ixC^S,rho_)
        w2L(ixC^S,mom(idir))=w2L(ixC^S,mom(idir))*w2L(ixC^S,rho_)
      end do
      ! get fluxes of intermedate states
      do iw=1,nwflux
        if (flux_type(idims, iw) == flux_tvdlf) then
          fC(ixC^S,iw,ip1)=0.5d0*(fLC(ixC^S,iw) + fRC(ixC^S,iw))
          !fC(ixC^S,iw,ip1) = 0.5d0 * (fLC(ixC^S,iw) + fRC(ixC^S,iw) - tvdlfeps * &
          !     max(cmaxC(ixC^S), abs(cminC(ixC^S))) * &
          !     (wRC(ixC^S,iw)-wLC(ixC^S,iw)))
          cycle
        end if
        f1L(ixC^S,iw)=fLC(ixC^S,iw)+cminC(ixC^S)*(w1L(ixC^S,iw)-wLC(ixC^S,iw))
        f1R(ixC^S,iw)=fRC(ixC^S,iw)+cmaxC(ixC^S)*(w1R(ixC^S,iw)-wRC(ixC^S,iw))
        where(cminC(ixC^S)>0.d0)
          fC(ixC^S,iw,ip1)=fLC(ixC^S,iw)
        else where(s1L(ixC^S)>=0.d0)
          fC(ixC^S,iw,ip1)=f1L(ixC^S,iw)
        else where(sm(ixC^S)>=0.d0)
          fC(ixC^S,iw,ip1)=f1L(ixC^S,iw)+s1L(ixC^S)*(w2L(ixC^S,iw)-w1L(ixC^S,iw))
        else where(s1R(ixC^S)>=0.d0)
          fC(ixC^S,iw,ip1)=f1R(ixC^S,iw)+s1R(ixC^S)*(w2R(ixC^S,iw)-w1R(ixC^S,iw))
        else where(cmaxC(ixC^S)>=0.d0)
          fC(ixC^S,iw,ip1)=f1R(ixC^S,iw)
        else where(cmaxC(ixC^S)<0.d0)
          fC(ixC^S,iw,ip1)=fRC(ixC^S,iw)
        end where
        if(.not.slab) then
          fC(ixC^S,iw,ip1)=block%surfaceC(ixC^S,ip1)*fC(ixC^S,iw,ip1)
        end if
      end do

    end subroutine get_Riemann_flux_hlld

  end subroutine finite_volume

  !> Determine the upwinded wLC(ixL) and wRC(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR(ixI^L,ixL^L,ixR^L,idims,w,wCT,wLC,wRC,wLp,wRp,x,needprim)
    use mod_physics
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: ixI^L, ixL^L, ixR^L, idims
    logical, intent(in) :: needprim
    double precision, dimension(ixI^S,1:nw) :: w, wCT
    ! left and right constructed status in conservative form
    double precision, dimension(ixI^S,1:nw) :: wLC, wRC
    ! left and right constructed status in primitive form
    double precision, dimension(ixI^S,1:nw) :: wLp, wRp
    double precision, dimension(ixI^S,1:ndim) :: x

    integer            :: jxR^L, ixC^L, jxC^L, iw
    double precision   :: ldw(ixI^S), rdw(ixI^S), dwC(ixI^S)
    ! integer            :: flagL(ixI^S), flagR(ixI^S)

    ! Transform w,wL,wR to primitive variables
    if (needprim) then
       call phys_to_primitive(ixI^L,ixI^L,w,x)
    end if

    if (typelimiter == limiter_mp5) then
       call MP5limiter(ixI^L,ixL^L,idims,w,wLp,wRp)
    else if (typelimiter == limiter_ppm) then
       call PPMlimiter(ixI^L,ixM^LL,idims,w,wCT,wLp,wRp)
    else
       jxR^L=ixR^L+kr(idims,^D);
       ixCmax^D=jxRmax^D; ixCmin^D=ixLmin^D-kr(idims,^D);
       jxC^L=ixC^L+kr(idims,^D);

       do iw=1,nwflux
          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=dlog10(w(ixCmin^D:jxCmax^D,iw))
             wLp(ixL^S,iw)=dlog10(wLp(ixL^S,iw))
             wRp(ixR^S,iw)=dlog10(wRp(ixR^S,iw))
          end if

          dwC(ixC^S)=w(jxC^S,iw)-w(ixC^S,iw)

          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixI^L,ixC^L,idims,typelimiter,ldw,rdw)
          wLp(ixL^S,iw)=wLp(ixL^S,iw)+half*ldw(ixL^S)
          wRp(ixR^S,iw)=wRp(ixR^S,iw)-half*rdw(jxR^S)

          if (loglimit(iw)) then
             w(ixCmin^D:jxCmax^D,iw)=10.0d0**w(ixCmin^D:jxCmax^D,iw)
             wLp(ixL^S,iw)=10.0d0**wLp(ixL^S,iw)
             wRp(ixR^S,iw)=10.0d0**wRp(ixR^S,iw)
          end if
       end do

       !! TODO: does this actually help? if not, remove
       !call phys_check_w(.true., ixI^L, ixL^L, wLtmp, flagL)
       !call phys_check_w(.true., ixI^L, ixR^L, wRtmp, flagR)

       !do iw=1,nwflux
       !   where (flagL(ixL^S) == 0 .and. flagR(ixR^S) == 0)
       !      wLC(ixL^S,iw)=wLtmp(ixL^S,iw)
       !      wRC(ixR^S,iw)=wRtmp(ixR^S,iw)
       !   end where

       !   ! Elsewhere, we still need to convert back when using loglimit
       !   if (loglimit(iw)) then
       !      where (flagL(ixL^S) /= 0 .or. flagR(ixR^S) /= 0)
       !         wLC(ixL^S,iw)=10.0d0**wLC(ixL^S,iw)
       !         wRC(ixR^S,iw)=10.0d0**wRC(ixR^S,iw)
       !      end where
       !   end if
       !enddo
    endif

    ! Transform w,wL,wR back to conservative variables
    if(needprim)then
       call phys_to_conserved(ixI^L,ixI^L,w,x)
    endif
    wLC(ixL^S,1:nw)=wLp(ixL^S,1:nw)
    wRC(ixR^S,1:nw)=wRp(ixR^S,1:nw)
    call phys_to_conserved(ixI^L,ixL^L,wLC,x)
    call phys_to_conserved(ixI^L,ixR^L,wRC,x)

  end subroutine reconstruct_LR

end module mod_finite_volume
