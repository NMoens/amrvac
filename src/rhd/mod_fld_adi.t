module mod_fld_adi
  use mod_fld
  implicit none

contains

  subroutine Evolve_E_rad(w, x, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision :: E_new(ixI^S), E_old(ixI^S), ADI_Error
    double precision :: frac_grid
    integer :: w_max, frac_dt
    logical :: converged

    integer :: i

    E_new(ixI^S) = w(ixI^S,iw_r_e)

    converged = .false.
    ADI_Error = bigdouble
    w_max = 1
    frac_grid = two
    frac_dt = 1

    do while (converged .eqv. .false.)

      !> Check if solution converged
      if (ADI_Error .lt. fld_adi_tol) then
        !> If converged in former loop, break loop
        converged = .true.
      else
        !> Reset E_new
        E_old(ixI^S) = w(ixI^S,iw_r_e)
        E_new(ixI^S) = w(ixI^S,iw_r_e)

        !> If no convergence, adapt pseudostepping
        w_max = 2*w_max
        frac_grid = 2*frac_grid
      endif

      !> Evolve using ADI
      if (converged .eqv. .false.) then

        call Evolve_ADI(w, x, E_new, E_old, w_max, frac_grid, ixI^L, ixO^L)
        call Error_check_ADI(w, x, E_new, E_old, ixI^L, ixO^L, ADI_Error) !> SHOULD THIS BE DONE EVERY ITERATION???
        if (ADI_Error .lt. fld_adi_tol) then
          converged = .true.
        endif
      endif

      !> If adjusting pseudostep doesn't work, divide the actual timestep in smaller parts
      if (w_max .gt. fld_maxdw) then
        if (converged .eqv. .false.) then
          !> use a smaller timestep than the hydrodynamical one
          call half_timestep_ADI(w, x, E_new, E_old, ixI^L, ixO^L, converged)
          call Error_check_ADI(w, x, E_new, E_old, ixI^L, ixO^L, ADI_Error)
          if (ADI_Error .lt. fld_adi_tol) then
            converged = .true.
          endif
        endif
      endif
    enddo

    w(ixO^S,iw_r_e) = E_new(ixO^S)
  end subroutine Evolve_E_rad


  subroutine half_timestep_ADI(w, x, E_new, E_old, ixI^L, ixO^L, converged)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim)
    double precision, intent(in) :: E_old(ixI^S)
    double precision, intent(out) :: E_new(ixI^S)
    logical, intent(inout) :: converged
    double precision :: frac_grid
    double precision :: E_loc(ixI^S)
    double precision :: saved_dt, ADI_Error
    integer :: i,  w_max, frac_dt

    saved_dt = dt
    ADI_Error = bigdouble
    frac_dt = 1
    5231 frac_dt = 2*frac_dt
    w_max = 1
    frac_grid = two

    if (frac_dt .gt. fld_max_fracdt) call mpistop("No convergence after halving timestep N times")
    dt = dt/frac_dt

    E_loc = E_old

    do i = 1,frac_dt
      !---------------------------------------------------------------
      do while (converged .eqv. .false.)
        !> Check if solution converged
        if (ADI_Error .lt. fld_adi_tol) then
          !> If converged in former loop, break loop
          converged = .true.
          goto 7895
        else
          !> If no convergence, adapt pseudostepping
          w_max = 2*w_max
          frac_grid = 2*frac_grid
        endif

        !> Evolve using ADI
        call Evolve_ADI(w, x, E_new, E_loc, w_max, frac_grid, ixI^L, ixO^L)
        call Error_check_ADI(w, x, E_new, E_loc, ixI^L, ixO^L, ADI_Error) !> SHOULD THIS BE DONE EVERY ITERATION???

        !> If adjusting pseudostep doesn't work, divide the actual timestep in smaller parts
        if (w_max .gt. fld_maxdw) goto 5231

        if (ADI_Error .lt. fld_adi_tol) then
          converged = .true.
        endif

      enddo
      !---------------------------------------------------------------
      7895 E_loc = E_new
    enddo

    dt = saved_dt
  end subroutine half_timestep_ADI


  subroutine Error_check_ADI(w, x, E_new, E_old, ixI^L, ixO^L, ADI_Error)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S, 1:ndim), w(ixI^S, 1:nw)
    double precision, intent(in) :: E_new(ixI^S), E_old(ixI^S)
    double precision, intent(out) :: ADI_Error
    double precision :: LHS(ixO^S), RHS(ixO^S), D(ixI^S,1:ndim)
    integer :: jx1^L, hx1^L,jx2^L, hx2^L

    integer :: i

    jx1^L=ixO^L+kr(1,^D);
    hx1^L=ixO^L-kr(1,^D);
    jx2^L=ixO^L+kr(2,^D);
    hx2^L=ixO^L-kr(2,^D);

    call fld_get_diffcoef(w, x, ixI^L, ixO^L, D)

    !> LHS = dx^2/dt * (E_new - E_old)
    LHS(ixO^S) = (x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))*&
    (x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))/dt*&
    (E_new(ixO^S) - E_old(ixO^S))

    !> RHS = D1(E_+ - E) - D1(E - E_-) + D2(E_+ - E) - D2(E - E_-)
    RHS(ixO^S) = &
      D(jx1^S,1)*(E_new(jx1^S) - E_new(ixO^S)) &
    - D(ixO^S,1)*(E_new(ixO^S) - E_new(hx1^S)) &
    + D(jx2^S,2)*(E_new(jx2^S) - E_new(ixO^S)) &
    - D(ixO^S,2)*(E_new(ixO^S) - E_new(hx2^S))

    ADI_Error = maxval(abs((RHS-LHS)/(E_old/dt))) !> Try mean value or smtn
    !ADI_Error = sum(abs((RHS-LHS)/(E_old/dt)))/((ixOmax1-ixOmin1)*(ixOmax2-ixOmin2))
  end subroutine Error_check_ADI


  subroutine Evolve_ADI(w, x, E_new, E_old, w_max, frac_grid, ixI^L, ixO^L)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, w_max
    double precision, intent(in) :: w(ixI^S, 1:nw), x(ixI^S, 1:ndim), frac_grid
    double precision, intent(in) :: E_old(ixI^S)
    double precision, intent(out):: E_new(ixI^S)
    double precision :: E_m(ixI^S), E_n(ixI^S)
    double precision :: diag1(ixImax1,ixImax2),sub1(ixImax1,ixImax2),sup1(ixImax1,ixImax2),bvec1(ixImax1,ixImax2)
    double precision :: diag2(ixImax2,ixImax1),sub2(ixImax2,ixImax1),sup2(ixImax2,ixImax1),bvec2(ixImax2,ixImax1)
    double precision :: Evec1(ixImin1:ixImax1), Evec2(ixImin2:ixImax2)
    double precision :: dw, w0, w1
    integer :: m, j, i

    w0 = (x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))*(x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2))/frac_grid
    w1 = (x(ixOmax1,ixOmin2,1)-x(ixOmin1,ixOmin2,1))*(x(ixOmin1,ixOmax2,2)-x(ixOmin1,ixOmin2,2))/frac_grid !4.d0

    E_m = E_old

    do m = 1,w_max
      E_n = E_old

      !> Set pseudotimestep
      dw = w0*(w1/w0)**((m-one)/(w_max-one))

      call ADI_boundary_conditions(ixI^L,ixO^L,E_m,w,x)

      !> Setup matrix and vector for sweeping in direction 1
      call make_matrix(x,w,dw,E_m,E_n,1,ixImax1,ixI^L, ixO^L,diag1,sub1,sup1,bvec1,diag2,sub2,sup2,bvec2)
      do j = ixImin2,ixImax2
        Evec1(ixImin1:ixImax1) = E_m(ixImin1:ixImax1,j)
        call solve_tridiag(ixOmin1,ixOmax1,ixImin1,ixImax1,diag1(:,j),sub1(:,j),sup1(:,j),bvec1(:,j),Evec1)
        !E_m(ixOmin1:ixOmax1,j) = Evec1(ixOmin1:ixOmax1)
        E_m(iximin1:ixImax1,j) = Evec1(ixImin1:ixImax1)
      enddo

      call ADI_boundary_conditions(ixI^L,ixO^L,E_m,w,x)

      !> Setup matrix and vector for sweeping in direction 2
      call make_matrix(x,w,dw,E_m,E_n,2,ixImax2,ixI^L,ixO^L,diag1,sub1,sup1,bvec1,diag2,sub2,sup2,bvec2)
      do j = ixImin1,ixImax1
        Evec2(ixImin2:ixImax2) = E_m(j,ixImin2:ixImax2)
        call solve_tridiag(ixOmin2,ixOmax2,ixImin2,ixImax2,diag2(:,j),sub2(:,j),sup2(:,j),bvec2(:,j),Evec2)
        !E_m(j,ixOmin2:ixOmax2) = Evec2(ixOmin2:ixOmax2)
        E_m(j,ixImin2:ixImax2) = Evec2(ixImin2:ixImax2)
      enddo

      call ADI_boundary_conditions(ixI^L,ixO^L,E_m,w,x)

    enddo
    E_new = E_m
  end subroutine Evolve_ADI


  subroutine fld_get_diffcoef(w, x, ixI^L, ixO^L, D)
    use mod_global_parameters

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(out):: D(ixI^S,1:ndim)
    double precision :: fld_kappa(ixO^S)
    double precision :: fld_lambda(ixO^S), fld_R(ixO^S)
    double precision :: D_center(ixI^S)
    integer :: idir,i,j

    if (fld_diff_testcase) then
      ! D = unit_length/unit_velocity
      D(ixI^S,1) = x(ixI^S,2)/maxval(x(ixI^S,2))*unit_length/unit_velocity
      D(ixI^S,2) = x(ixI^S,2)/maxval(x(ixI^S,2))*unit_length/unit_velocity
    else
      !> calculate lambda
      call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, fld_lambda, fld_R)

      !> set Opacity
      call fld_get_opacity(w, x, ixI^L, ixO^L, fld_kappa)

      !> calculate diffusion coefficient
      D_center(ixO^S) = fld_speedofligt_0*fld_lambda(ixO^S)/(fld_kappa(ixO^S)*w(ixO^S,iw_rho))

      !> Extrapolate lambda to ghostcells
      !> Edges
      !> To calculate the diffusion coefficient at the ghostcells, copy lambda from grid, but use correct kappa and rho
      do i = 0,nghostcells-1
        D_center(ixImin1+i,:) = D_center(ixImin1+nghostcells,:)
        D_center(ixImax1-i,:) = D_center(ixImax1-nghostcells,:)
        D_center(:,ixImin2+i) = D_center(:,ixImin2+nghostcells)
        D_center(:,ixImax2-i) = D_center(:,ixImax2-nghostcells)
      end do

      !call Diff_boundary_conditions(ixI^L,ixO^L,D)

      !> Corners
      do i = 0,nghostcells-1
        do j = 0, nghostcells-1
          D_center(ixImin1+i,ixImax2-j) = D_center(ixImin1+nghostcells,ixImax2-nghostcells)
          D_center(ixImax1-i,ixImax2-j) = D_center(ixImax1-nghostcells,ixImax2-nghostcells)
          D_center(ixImin1+i,ixImin2+j) = D_center(ixImin1+nghostcells,ixImin2+nghostcells)
          D_center(ixImax1-i,ixImin2+j) = D_center(ixImax1-nghostcells,ixImin2+nghostcells)
        end do
      end do

      !> Go from cell center to cell face
      do i = ixImin1+1, ixImax1
      do j = ixImin2+1, ixImax2
         ! D(i,j,1) = (D_center(i,j) + D_center(i-1,j))/two
         ! D(i,j,2) = (D_center(i,j) + D_center(i,j-1))/two
        D(i,j,1) = (2*D_center(i,j) + 2*D_center(i-1,j)&
                  + D_center(i,j+1) + D_center(i-1,j+1)&
                  + D_center(i,j-1) + D_center(i-1,j-1))/8.d0
        D(i,j,2) = (2*D_center(i,j) + 2*D_center(i,j-1)&
                  + D_center(i+1,j) + D_center(i+1,j-1)&
                  + D_center(i-1,j) + D_center(i-1,j-1))/8.d0
      enddo
      enddo
      D(ixImin1,:,1) = D_center(ixImin1,:)
      D(:,ixImin2,1) = D_center(:,ixImin2)
      D(ixImin1,:,2) = D_center(ixImin1,:)
      D(:,ixImin2,2) = D_center(:,ixImin2)

      !D(:,ixImax2-2,:) = D(:,ixImax2-3,:)

    endif
  end subroutine fld_get_diffcoef


  subroutine make_matrix(x,w,dw,E_m,E_n,sweepdir,ixImax,ixI^L,ixO^L,diag1,sub1,sup1,bvec1,diag2,sub2,sup2,bvec2)
    use mod_global_parameters

    integer, intent(in) :: sweepdir, ixImax
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S, 1:nw), dw
    double precision, intent(in) :: x(ixI^S, 1:ndim)
    double precision, intent(in) :: E_n(ixI^S), E_m(ixI^S)
    double precision, intent(out):: diag1(ixImin1:ixImax1,ixImin2:ixImax2),sub1(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out):: sup1(ixImin1:ixImax1,ixImin2:ixImax2),bvec1(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out):: diag2(ixImin2:ixImax2,ixImin1:ixImax1),sub2(ixImin2:ixImax2,ixImin1:ixImax1)
    double precision, intent(out):: sup2(ixImin2:ixImax2,ixImin1:ixImax1),bvec2(ixImin2:ixImax2,ixImin1:ixImax1)
    double precision :: D(ixI^S,1:ndim), h, beta(0:ixImax), delta_x
    integer :: idir,i,j

    call fld_get_diffcoef(w, x, ixI^L, ixO^L, D)

    !calculate h
    if (sweepdir == 1) then
      delta_x = x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1)
      !delta_x = x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2)
    elseif (sweepdir == 2) then
      !delta_x = x(ixOmin1+1,ixOmin2,1)-x(ixOmin1,ixOmin2,1)
      delta_x = x(ixOmin1,ixOmin2+1,2)-x(ixOmin1,ixOmin2,2)
    endif
    h = dw/(two*delta_x**two)

    !> Matrix depends on sweepingdirection
    if (sweepdir == 1) then
      !calculate matrix for sweeping in 1-direction
      do j = ixImin2,ixImax2
       !calculate beta
       do i = ixImin1,ixImax1-1
         beta(i) = one + dw/(two*dt) + h*(D(i+1,j,1)+D(i,j,1))
       enddo

       do i = ixImin1,ixImax1-1
         diag1(i,j) = beta(i)
         sub1(i+1,j) = -h*D(i+1,j,1)
         sup1(i,j) = -h*D(i+1,j,1)
         bvec1(i,j) = (one - h*(D(i,j+1,2)+D(i,j,2)))*E_m(i,j) &
         + h*D(i,j+1,2)*E_m(i,j+1) + h*D(i,j,2)*E_m(i,j-1) + dw/(two*dt)*E_n(i,j)
       enddo

       !> Boundary conditions on matrix
       sub1(ixImin1,j) = zero
       sup1(ixImax1,j) = zero
       diag1(ixImin1,j) = beta(ixImin1) - h*D(ixImin1,j,1)
       diag1(ixImax1,j) = beta(ixImax1-1) - h*D(ixImax1,j,1)
       bvec1(ixImax1,j) = (one - h*(D(ixImax1,j+1,2)+D(ixImax1,j,2)))*E_m(ixImax1,j) &
       + h*D(ixImax1,j+1,2)*E_m(ixImax1,j+1) + h*D(ixImax1,j,2)*E_m(ixImax1,j-1) + dw/(two*dt)*E_n(ixImax1,j)

      enddo

    elseif ( sweepdir == 2 ) then
      !calculate matrix for sweeping in 2-direction
      do j = ixImin1,ixImax1
       !calculate beta
       do i = ixImin2,ixImax2-1
         beta(i) = one + dw/(two*dt) + h*(D(j,i+1,2)+D(j,i,2))
       enddo

       do i = ixImin2,ixImax2-1
         diag2(i,j) = beta(i)
         sub2(i+1,j) = -h*D(j,i+1,2)
         sup2(i,j) = -h*D(j,i+1,2)
         bvec2(i,j) = (one - h*(D(j+1,i,1)+D(j,i,1)))*E_m(j,i) &
         + h*D(j+1,i,1)*E_m(j+1,i) + h*D(j,i,1)*E_m(j-1,i) + dw/(two*dt)*E_n(j,i)
       enddo

       !> Boundary conditions on matrix
       sub2(ixImin2,j) = zero
       sup2(ixImax2,j) = zero
       diag2(ixImin2,j) = beta(ixImin2) - h*D(j,ixImin2,2)
       diag2(ixImax2,j) = beta(ixImax2-1) - h*D(j,ixImax2,2)
       bvec2(ixImax2,j) = (one - h*(D(j+1,ixImax2,1)+D(j,ixImax2,1)))*E_m(j,ixImax2) &
       + h*D(j+1,ixImax2,1)*E_m(j+1,ixImax2) + h*D(j,ixImax2,1)*E_m(j-1,ixImax2) + dw/(two*dt)*E_n(j,ixImax2)
      enddo

    else
      call mpistop("sweepdirection unknown")
    endif
  end subroutine make_matrix


  subroutine solve_tridiag(ixOmin,ixOmax,ixImin,ixImax,diag,sub,sup,bvec,Evec)
    use mod_global_parameters
    implicit none

    integer, intent(in) :: ixOmin,ixOmax,ixImin,ixImax
    double precision, intent(in) :: diag(ixImin:ixImax), bvec(ixImin:ixImax)
    double precision, intent(in) :: sub(ixImin:ixImax), sup(ixImin:ixImax)
    double precision, intent(out) :: Evec(ixImin:ixImax)
    double precision :: cp(ixImin:ixImax), dp(ixImin:ixImax)
    integer :: i

    ! initialize c-prime and d-prime
    cp(ixImin) = sup(ixImin)/diag(ixImin)
    dp(ixImin) = bvec(ixImin)/diag(ixImin)

    ! solve for vectors c-prime and d-prime
    do i = ixImin+1 ,ixImax-1
      cp(i) = sup(i)/(diag(i)-cp(i-1)*sub(i))
      dp(i) = (bvec(i)-dp(i-1)*sub(i))/(diag(i)-cp(i-1)*sub(i))
    enddo
    dp(ixImax) = (bvec(ixImax)-dp(ixImax-1)*sub(ixImax))/(diag(ixImax)-cp(ixImax-1)*sub(ixImax))

    ! initialize x
    Evec(ixImax-1) = dp(ixImax-1)

    ! solve for x from the vectors c-prime and d-prime
    do i = ixImax-2, ixImin+1, -1
      Evec(i) = dp(i)-cp(i)*Evec(i+1)
    end do
  end subroutine solve_tridiag


  subroutine ADI_boundary_conditions(ixI^L,ixO^L,E_m,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(inout) :: E_m(ixI^S)
    integer g, h

    select case (fld_bound_min2)
    case('periodic')
      E_m(:,ixImin2:ixImin2+1) = E_m(:,ixImax2-3:ixImax2-2)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      E_m(:,ixOmin2-1) = 2.d0*E_m(:,ixOmin2) - E_m(:,ixOmin2+1)
      E_m(:,ixImin2) = 2.d0*E_m(:,ixOmin2-1) - E_m(:,ixOmin2)
    case('fixed')
      E_m(:,ixImin2:ixOmin2-1) = w(:,ixImin2:ixOmin2-1,iw_r_e)
    case default
      call mpistop("ADI boundary not defined")
    end select

    select case (fld_bound_max2)
    case('periodic')
      E_m(:,ixImax2-1:ixImax2) = E_m(:,ixImin2+2:ixImin2+3)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      ! E_m(:,ixOmax2+1) = 2.d0*E_m(:,ixOmax2) - E_m(:,ixOmax2-1)
      ! E_m(:,ixImax2) = 2.d0*E_m(:,ixOmax2+1) - E_m(:,ixOmax2)

      E_m(:,ixOmax2+1) = max(abs(w(:, ixOmax2 - 1, iw_rho)/w(:, ixOmax2, iw_rho)&
      *(w(:, ixOmax2, r_e) - w(:, ixOmax2-2, r_e)) + w(:, ixOmax2-1, r_e)),zero)
      E_m(:,ixOmax2+2) = max(abs(w(:, ixOmax2, iw_rho)/w(:, ixOmax2+1, iw_rho)&
      *(w(:, ixOmax2+1, r_e) - w(:, ixOmax2-1, r_e)) + w(:, ixOmax2, r_e)),zero)

      ! E_m(:,ixOmax2+1) = max(abs(x(:,ixOmax2+1,2)- x(:,ixOmax2,2))**two/two*(2*w(:,ixOmax2,r_e) - w(:,ixOmax2-1,r_e)),zero)
      ! E_m(:,ixImax2) = max(abs(x(:,ixOmax2+1,2)- x(:,ixOmax2,2))**two/two*(2*w(:,ixImax2-1,r_e) - w(:,ixImax2-2,r_e)),zero)

    case('fixed')
      E_m(:,ixImax2:ixOmax2+1) = w(:,ixImax2:ixOmax2+1,iw_r_e)
    case default
      call mpistop("ADI boundary not defined")
    end select

    select case (fld_bound_min1)
    case('periodic')
      E_m(ixImin1:ixImin1+1,:) = E_m(ixImax1-3:ixImax1-2,:)
      !E_m(ixImin1:ixImin1+1,:) = w(ixImax1-3:ixImax1-2,:,r_e)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      E_m(ixOmin1-1,:) = 2.d0*E_m(ixOmin1,:) - E_m(ixOmin1+1,:)
      E_m(ixImin1,:) = 2.d0*E_m(ixOmin1-1,:) - E_m(ixOmin1,:)
    case('fixed')
      E_m(ixImin1:ixOmin1-1,:) = w(ixImin1:ixOmin1-1,:,iw_r_e)
    case default
      call mpistop("ADI boundary not defined")
    end select

    select case (fld_bound_max1)
    case('periodic')
      E_m(ixImax1-1:ixImax1,:) = E_m(ixImin1+2:ixImin1+3,:)
      !E_m(ixImax1-1:ixImax1,:) = w(ixImin1+2:ixImin1+3,:,r_e)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      E_m(ixOmax1+1,:) = 2.d0*E_m(ixOmax1,:) - E_m(ixOmax1-1,:)
      E_m(ixImax1,:) = 2.d0*E_m(ixOmax1+1,:) - E_m(ixOmax1,:)
    case('fixed')
      E_m(ixImax1:ixOmax1+1,:) = w(ixImax1:ixOmax1+1,:,iw_r_e)
    case default
      call mpistop("ADI boundary not defined")
    end select

    !Corners
    ! do g = 0,nghostcells-1
    !   do h = 0, nghostcells-1
    !     E_m(ixImin1+g,ixImax2-h) = w(ixImin1+nghostcells,ixImax2-nghostcells,iw_r_e)
    !     E_m(ixImax1-g,ixImax2-h) = w(ixImax1-nghostcells,ixImax2-nghostcells,iw_r_e)
    !     E_m(ixImin1+g,ixImin2+h) = w(ixImin1+nghostcells,ixImin2+nghostcells,iw_r_e)
    !     E_m(ixImax1-g,ixImin2+h) = w(ixImax1-nghostcells,ixImin2+nghostcells,iw_r_e)
    !   end do
    ! end do
  end subroutine ADI_boundary_conditions


  subroutine Diff_boundary_conditions(ixI^L,ixO^L,D)
    use mod_global_parameters

    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(inout) :: D(ixI^S)
    double precision :: Dmn1(ixI^S),Dmx1(ixI^S),Dmn2(ixI^S),Dmx2(ixI^S)

    select case (fld_bound_min1)
    case('periodic')
      D(ixImin1:ixOmin1-1,:) = D(ixOmax1-1:ixOmax1,:)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      D(ixOmin1-1,:) = 2.d0*D(ixOmin1,:) - D(ixOmin1+1,:)
      D(ixImin1,:) = 2.d0*D(ixOmin1-1,:) - D(ixOmin1,:)
    case('fixed')
      if (it==0) then
        Dmn1 = D
      endif
      D(ixImin1:ixOmin1-1,:) = Dmn1(ixImin1:ixOmin1-1,:)
    case default
      call mpistop("ADI boundary not defined")
    end select

    select case (fld_bound_max1)
    case('periodic')
      D(ixImax1:ixOmax1+1,:) = D(ixOmin1+1:ixOmin1,:)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      D(ixOmax1+1,:) = 2.d0*D(ixOmax1,:) - D(ixOmax1-1,:)
      D(ixImax1,:) = 2.d0*D(ixOmax1+1,:) - D(ixOmax1,:)
    case('fixed')
      if (it==0) then
        Dmx1 = D
      endif
      D(ixImax1:ixOmax1+1,:) = Dmx1(ixImax1:ixOmax1+1,:)
    case default
      call mpistop("ADI boundary not defined")
    end select

    select case (fld_bound_min2)
    case('periodic')
      D(:,ixImin2:ixOmin2-1) = D(:,ixOmax2-1:ixOmax2)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      D(:,ixOmin2-1) = 2.d0*D(:,ixOmin2) - D(:,ixOmin2+1)
      D(:,ixImin2) = 2.d0*D(:,ixOmin2-1) - D(:,ixOmin2)
    case('fixed')
      if (it==0) then
        Dmn2 = D
      endif
      D(:,ixImin2:ixOmin2-1) = Dmn2(:,ixImin2:ixOmin2-1)
      ! D(:,ixImin2+1) = Dmn2(:,ixOmin2)
      ! D(:,ixImin2) = Dmn2(:,ixOmin2)
    case default
      call mpistop("ADI boundary not defined")
    end select

    select case (fld_bound_max2)
    case('periodic')
      D(:,ixImax2:ixOmax2+1) = D(:,ixOmin2+1:ixOmin2)
    case('cont')
      if (nghostcells .ne. 2) call mpistop("continious ADI boundary conditions not defined for more than 2 ghostcells")
      D(:,ixOmax2+1) = 2.d0*D(:,ixOmax2) - D(:,ixOmax2-1)
      D(:,ixImax2) = 2.d0*D(:,ixOmax2+1) - D(:,ixOmax2)
    case('fixed')
      if (it==0) then
        Dmx2 = D
      endif
      D(:,ixImax2:ixOmax2+1) = Dmx2(:,ixImax2:ixOmax2+1)
    case default
      call mpistop("ADI boundary not defined")
    end select
  end subroutine Diff_boundary_conditions


end module mod_fld_adi
