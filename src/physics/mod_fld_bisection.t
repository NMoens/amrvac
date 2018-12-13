module mod_fld_bisection
  use mod_fld
  implicit none

contains

  subroutine Energy_interaction(w, x, ixI^L, ixO^L)
    use mod_global_parameters
    use mod_physics, only: phys_get_pthermal

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: fld_kappa(ixO^S)
    double precision :: rad_pressure(ixO^S)
    double precision :: temperature(ixI^S), div_v(ixI^S), vel(ixI^S,1:ndim)
    double precision :: a1(ixO^S), a2(ixO^S), a3(ixO^S)
    double precision :: c0(ixO^S), c1(ixO^S)
    double precision :: e_gas(ixO^S), E_rad(ixO^S)

    integer :: i,j,idir

    !> Calculate the radiative flux using the FLD Approximation
    call fld_get_radpress(w,x,ixI^L,ixO^L,rad_pressure)

    !> Get pressure
    call phys_get_pthermal(w,x,ixI^L,ixO^L,temperature)

    !> calc Temperature as p/rho
    temperature(ixO^S)=temperature(ixO^S)/w(ixO^S,iw_rho)

    !> set Opacity
    call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)

    !> calc photon tiring term
    do idir=1,ndim
      vel(ixI^S,idir)= w(ixI^S,iw_mom(idir))/w(ixI^S,iw_rho)
    enddo

    call divvector(vel,ixI^L,ixO^L,div_v)

    e_gas(ixO^S) = w(ixO^S,iw_e)
    E_rad(ixO^S) = w(ixO^S,iw_r_e)

    !> Calculate coefficients for polynomial
    a1(ixO^S) = 4*fld_kappa(ixO^S)*w(ixO^S,iw_rho)*fld_sigma_0*(temperature(ixO^S)/e_gas(ixO^S))**4.d0*dt
    a2(ixO^S) = fld_speedofligt_0*fld_kappa(ixO^S)*w(ixO^S,iw_rho)*dt
    a3(ixO^S) = div_v(ixO^S)*rad_pressure(ixO^S)/E_rad(ixO^S)*dt

    c0(ixO^S) = ((one + a1(ixO^S) + a3(ixO^S))*e_gas(ixO^S) + a2(ixO^S)*E_rad(ixO^S))/(a1(ixO^S)*(one + a3(ixO^S)))
    c1(ixO^S) = (one + a1(ixO^S) + a3(ixO^S))/(a1(ixO^S)*(one + a3(ixO^S)))

    !> Loop over every cell for bisection method
    do i = ixOmin1,ixOmax1
    do j =  ixOmin2,ixOmax2
      call Bisection_method(e_gas(i,j), E_rad(i,j), c0(i,j), c1(i,j))
    enddo
    enddo

    !> Update gas-energy in w
    w(ixO^S,iw_e) = e_gas(ixO^S)

    !> Calculate new radiation energy
    !> Get new pressure
    call phys_get_pthermal(w,x,ixI^L,ixO^L,temperature)

    !> calc new Temperature as p/rho
    temperature(ixO^S)=(temperature(ixO^S)/w(ixO^S,iw_rho))

    !> Update a1
    a1(ixO^S) = 4*fld_kappa(ixO^S)*w(ixO^S,iw_rho)*fld_sigma_0*(temperature(ixO^S)/e_gas(ixO^S))**4.d0*dt

    !> advance E_rad
    E_rad(ixO^S) = (a1*e_gas(ixO^S)**4.d0 + E_rad(ixO^S))/(one + a2 + a3)

    !> Update rad-energy in w
    w(ixO^S,iw_r_e) = E_rad(ixO^S)
  end subroutine Energy_interaction


  subroutine Bisection_method(e_gas, E_rad, c0, c1)
    use mod_global_parameters

    double precision, intent(in)    :: c0, c1
    double precision, intent(in)    :: E_rad
    double precision, intent(inout) :: e_gas

    double precision :: bisect_a, bisect_b, bisect_c

    bisect_a = zero
    bisect_b = max(abs(c0/c1),abs(c0)**(1.d0/4.d0))

    do while (abs(Polynomial_Bisection(bisect_b, c0, c1)-Polynomial_Bisection(bisect_a, c0, c1))&
       .ge. fld_bisect_tol*min(e_gas,E_rad))
      bisect_c = (bisect_a + bisect_b)/two

      if (Polynomial_Bisection(bisect_a, c0, c1)*&
      Polynomial_Bisection(bisect_b, c0, c1) .le. zero) then

        if (Polynomial_Bisection(bisect_a, c0, c1)*&
        Polynomial_Bisection(bisect_c, c0, c1) .le. zero) then
          bisect_b = bisect_c
        elseif (Polynomial_Bisection(bisect_b, c0, c1)*&
        Polynomial_Bisection(bisect_c, c0, c1) .le. zero) then
          bisect_a = bisect_c
        else
          call mpistop("Problem with fld bisection method")
        endif
      else
        bisect_a = e_gas
        bisect_b = e_gas
        print*, "IGNORING ENERGY GAS-RAD EXCHANGE ", c0, c1
      endif
    enddo

    e_gas = (bisect_a + bisect_b)/two
  end subroutine Bisection_method


  function Polynomial_Bisection(e_gas, c0, c1) result(pol_result)
    use mod_global_parameters

    double precision, intent(in) :: e_gas
    double precision, intent(in) :: c0, c1
    double precision :: pol_result

    pol_result = e_gas**4.d0 + c1*e_gas - c0
  end function Polynomial_Bisection


end module mod_fld_bisection
