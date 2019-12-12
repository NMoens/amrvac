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

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    w(ix^S,rho_)   =  rho0/1.d5
    w(ix^S,mom(1)) =  0.d0 ! -rho0*vel*dcos(x(ix^S,2))
    w(ix^S,mom(2)) =  0.d0 ! rho0*vel*dsin(x(ix^S,2))
    w(ix^S,mom(3)) =  0.d0 ! rho0*vel*dsin(x(ix^S,2))
    w(ix^S,e_)     =  (pini/1.d5)/(hd_gamma-one) ! +half*rho0*vel**two

  end subroutine initonegrid_usr

  subroutine pt_grav_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: radius(ixO^S)
    integer :: ii

    radius(ixO^S) = dsqrt(x(ixO^S,1)**2 + x(ixO^S,2)**2 + x(ixO^S,3)**2)

    do ii = 1,ndir
      w(ixO^S,mom(ii))=w(ixO^S,mom(ii)) - &
           qdt * wCT(ixO^S,rho_)    * (gm*x(ixO^S,ii)/radius(ixO^S)**3.d0)

      w(ixO^S,e_ )=w(ixO^S,e_ ) - &
           qdt * wCT(ixO^S,mom(ii))  * (gm*x(ixO^S,ii)/radius(ixO^S)**3.d0)
    enddo

    ! if (x(1,3,1,1)<xprobmin1) print*, it, qdt, qtC

  end subroutine pt_grav_source

  subroutine get_dt_pt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
    double precision, intent(in) :: w(ixG^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: radius(ix^S)

    radius(ix^S) = dsqrt(x(ix^S,1)**2 + x(ix^S,2)**2 + x(ix^S,3)**2)

    dtnew=bigdouble
    dtnew=min(dtnew,minval(dsqrt(block%dx(ix^S,1)/(gm/radius(ix^S)**two))))

  end subroutine get_dt_pt_grav

  subroutine specialbound_usr(qt,ixG^L,ixB^L,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision :: rho(ixG^S)

    select case (iB)

    case(1)
      w(ixB^S,rho_) = rho0
      w(ixB^S,mom(1)) = rho0*vel
      w(ixB^S,mom(2)) = zero
      w(ixB^S,mom(3)) = zero
      w(ixB^S,e_) = pini/(hd_gamma-one)+half*rho0*vel**two
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
