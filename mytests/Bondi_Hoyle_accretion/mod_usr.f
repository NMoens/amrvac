module mod_usr
  use mod_hd
  implicit none
  double precision :: mach, rho0, vel, Racc, gm, pini

contains

  subroutine usr_init()
    call set_coordinate_system('Cartesian_3D')

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_source        => pt_grav_source
    usr_get_dt        => get_dt_pt_grav
    usr_special_bc    => specialbound_usr
    ! usr_internal_bc   => intern_bc

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr

    hd_gamma=5.d0/3.d0
    mach=4.d0
    vel=one
    gm=half
    rho0=one
    Racc=one
    pini=((1/hd_gamma)*rho0*vel**two)/mach**two

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

    ! initialize one grid

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)

    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)   =  rho0/1.d5
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(1)) =  0.d0 !-rho0*vel*dcos(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2))
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(2)) =  0.d0 !rho0*vel*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2))
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(3)) =  0.d0 !rho0*vel*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2))
    w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       e_)     =  (pini/1.d5)/(hd_gamma-one) !+half*rho0*vel**two

  end subroutine initonegrid_usr

  subroutine pt_grav_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iwmin,iwmax,qtC,&
     wCT,qt,w,x)

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iwmin,iwmax
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    double precision :: radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3)
    integer :: ii

    radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3) = dsqrt(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,1)**2 + x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,2)**2 + x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ixOmin3:ixOmax3,3)**2)

    do ii = 1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(ii))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(ii)) - qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)    * (gm*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         ii)/radius(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)**3.d0)

      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_ )=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_ ) - qdt * wCT(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(ii))  * (gm*x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3,ii)/radius(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,ixOmin3:ixOmax3)**3.d0)
    enddo

    ! if (x(1,3,1,1)<xprobmin1) print*, it, qdt, qtC

  end subroutine pt_grav_source

  subroutine get_dt_pt_grav(w,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,dtnew,dx1,dx2,dx3,x)

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in) :: dx1,dx2,dx3, x(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)
    double precision, intent(in) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: radius(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)

    radius(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = dsqrt(x(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3,1)**2 + x(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3,2)**2 + x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       3)**2)

    dtnew=bigdouble
    dtnew=min(dtnew,minval(dsqrt(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3,1)/(gm/radius(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3)**two))))

  end subroutine get_dt_pt_grav

  subroutine specialbound_usr(qt,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
     ixGmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision :: rho(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3)

    select case (iB)

    case(1)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,rho_) = rho0
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(1)) = rho0*vel
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(2)) = zero
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(3)) = zero
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         e_) = pini/(hd_gamma-one)+half*rho0*vel**two
    case default
      call mpistop("not defined this spec_bound")
    end select

  end subroutine specialbound_usr

  ! subroutine intern_bc(level,qt,ixI^L,ixO^L,w,x)
  !
  !   integer, intent(in) :: ixI^L,ixO^L,level
  !   double precision, intent(in) :: qt
  !   double precision, intent(inout) :: w(ixI^S,1:nw)
  !   double precision, intent(in) :: x(ixI^S,1:ndim)
  !   double precision :: R1, R2, rho(ixI^S), zeta(ixI^S)
  !   double precision :: radius(ixO^S)
  !
  !   radius(ixO^S) = dsqrt(x(ixO^S,1)**2 + x(ixO^S,2)**2 + x(ixO^S,3)**2)
  !
  !   R1=half*Racc ! xprobmax1/(335./20.)
  !   R2=one *Racc ! xprobmax1/(335./35.)
  !   zeta(ixO^S)=half*radius(ixO^S)*dsin(x(ixO^S,2))*&
  !        ( one + dsqrt(one+4.d0*(gm/vel**two)*&
  !        ((one-dcos(x(ixO^S,2)))/(radius(ixO^S)*dsin(x(ixO^S,2))**two))) )
  !   rho(ixO^S)=rho0*&
  !        (zeta(ixO^S)**two/&
  !        (x(ixO^S,1)*dsin(x(ixO^S,2))*(two*zeta(ixO^S)-x(ixO^S,1)*dsin(x(ixO^S,2)))))
  !   ! Is it ahead the shock? THEN FIX! No need to compare to B-K solution
  !   where ( x(ixO^S,2)<dpi/two .AND. x(ixO^S,1)>R2/(1+((R2-R1)/R1)*dcos(x(ixO^S,2))) )
  !      w(ixO^S,rho_)= rho(ixO^S)
  !      w(ixO^S,mom(1)) = -rho0*&
  !           (zeta(ixO^S)**two/&
  !           (x(ixO^S,1)*dsin(x(ixO^S,2))*(two*zeta(ixO^S)-x(ixO^S,1)*dsin(x(ixO^S,2)))))*&
  !           dsqrt(vel**two+(two*gm)/x(ixO^S,1)-&
  !           ((zeta(ixO^S)*vel)/x(ixO^S,1))**two)
  !      w(ixO^S,mom(2)) = rho0*&
  !           (zeta(ixO^S)**two/&
  !           (x(ixO^S,1)*dsin(x(ixO^S,2))*(two*zeta(ixO^S)-x(ixO^S,1)*dsin(x(ixO^S,2)))))*&
  !           ((zeta(ixO^S)*vel)/x(ixO^S,1))
  !      w(ixO^S,e_)  = ((rho0*vel**two)/hd_gamma)*&
  !           (one/(mach**two*(hd_gamma-one))+half)*&
  !           (zeta(ixO^S)**two/&
  !           (x(ixO^S,1)*dsin(x(ixO^S,2))*(two*zeta(ixO^S)-x(ixO^S,1)*dsin(x(ixO^S,2)))))+&
  !           (rho0/hd_gamma)*&
  !           (zeta(ixO^S)**two/&
  !           (x(ixO^S,1)*dsin(x(ixO^S,2))*(two*zeta(ixO^S)-x(ixO^S,1)*dsin(x(ixO^S,2)))))*&
  !           (half*(hd_gamma-one)* (dsqrt(vel**two+(two*gm)/x(ixO^S,1)-&
  !           ((zeta(ixO^S)*vel)/x(ixO^S,1))**two)**two+&
  !           ((zeta(ixO^S)*vel)/x(ixO^S,1))**two) + gm/x(ixO^S,1))
  !   endwhere
  !
  ! end subroutine intern_bc

end module mod_usr
