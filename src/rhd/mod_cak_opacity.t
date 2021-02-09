module mod_cak_opacity
    implicit NONE

    !> min and max indices for R,T-range in opacity table
    integer, parameter :: rmin = 2
    integer, parameter :: rmax = 21
    integer, parameter :: tmin = 2
    integer, parameter :: tmax = 21

    !> The opacity tables are read once and stored globally in Kappa_vals
    double precision, public :: alpha_vals(2:21,2:21)
    double precision, public :: Qbar_vals(2:21,2:21)
    double precision, public :: Q0_vals(2:21,2:21)
    double precision, public :: kappa_e_vals(2:21,2:21)

    double precision, public :: Log_D_list(2:21)
    double precision, public :: Log_T_list(2:21)

    character(*), parameter, public :: AMRVAC_DIR = "/lhome/nicolasm/amrvac/" ! use call getenv("AMRVAC_DIR", AMRVAC_DIR)
    character(*), parameter, public :: fileplace = AMRVAC_DIR//"src/rhd/CAK_tables/"

    public :: init_cak
    public :: set_cak_opacity

  contains

!> This routine is called when the fld radiation module is initialised.
!> Here, the tables for different He Abndcs are read and interpolated
subroutine init_cak()
  ! use mod_global_parameters

  call read_table(Log_D_list, Log_T_list, alpha_vals, "alg_TD")
  call read_table(Log_D_list, Log_T_list, Qbar_vals, "Q_TD")
  call read_table(Log_D_list, Log_T_list, Q0_vals, "Q0g_TD")
  call read_table(Log_D_list, Log_T_list, kappa_e_vals, "K_TD")

  print*, "Read Luka's tables"

end subroutine init_cak

!> This subroutine calculates the opacity for
!> a given temperature-density structure.
!> The opacities are read from a table that has the initialised metalicity
subroutine set_cak_opacity(rho,temp, gradv,kappa_cak)
  double precision, intent(in) :: rho, temp, gradv
  double precision, intent(out) :: kappa_cak

  double precision, PARAMETER :: const_c     = 2.99792458d10   ! cm s^-1           ; Speed of light

  double precision :: D_input, T_input
  double precision :: alpha_output, Qbar_output, Q0_output, kappa_e_output
  double precision :: alpha, Qbar, Q0, kappa_e

  double precision :: tau, M_t

  D_input = dlog10(rho)
  T_input = dlog10(temp)

  ! print*, 'input D and T'
  ! print*, D_input, T_input

  D_input = min(-10.d0-1.d-5, D_input)
  D_input = max(-20.d0+1.d-5, D_input)
  T_input = min(4.7d0-1.d-5, T_input)
  T_input = max(3.7d0+1.d-5, T_input)

  call get_val(alpha_vals, Log_D_list, Log_T_list, D_input, T_input, alpha_output)
  call get_val(Qbar_vals, Log_D_list, Log_T_list, D_input, T_input, Qbar_output)
  call get_val(Q0_vals, Log_D_list, Log_T_list, D_input, T_input, Q0_output)
  call get_val(kappa_e_vals, Log_D_list, Log_T_list, D_input, T_input, kappa_e_output)

  ! !> If the outcome is 9.999, look right in the table
  ! do while (K_output .gt. 9.0d0)
  !     ! print*, 'R,T datapoint out of opal table'
  !     D_input = D_input + 0.5
  !     call get_val(Kappa_vals, Log_D_list, Log_T_list, D_input, T_input, K_output)
  ! enddo
  !
  ! !> If the outcome is NaN, look left in the table
  ! do while (K_output .eq. 0.0d0)
  !     ! print*, 'R,T datapoint out of opal table'
  !     D_input = D_input - 0.5d0
  !     call get_val(Kappa_vals, Log_D_list, Log_T_list, D_input, T_input, K_output)
  ! enddo

  alpha = alpha_output
  Qbar = Qbar_output
  Q0 = Q0_output
  kappa_e = kappa_e_output

  ! print*, alpha, Qbar, Q0, kappa_e

  tau = Q0*kappa_e*rho*const_c/gradv
  M_t = Qbar/(1-alpha)*((1+tau)**(1-alpha) - 1)/tau
  kappa_cak = kappa_e*M_t

end subroutine set_cak_opacity

!> This routine reads out values and arguments from a table
subroutine read_table(D, T, K, filename)
    !> This routine reads in the the values for log kappa, and the values for log T and log R on the x and y axis

    double precision, intent(out) :: K(2:21,2:21), D(2:21), T(2:21)
    character(*), intent(in) :: filename

    character :: dum
    integer :: row, col

    OPEN(1,status = 'old', FILE=fileplace//filename)

    !> Read logT
    READ(1,*) dum,T(2:20)

    !> Read T and K
    do row = 2,21 !> NOT READING ENTIRE TABLE
      ! READ(1,'(f4.2,19f7.3)') D(row), K(row,2:20)
      READ(1,*) D(row), K(row,2:20)
    enddo

    CLOSE(1)

end subroutine read_table

!>This subroutine looks in the table for the four couples (T,R)
!surrounding a given input for T and R
subroutine get_val(Kappa_vals, Log_D_list, Log_T_list, D, T, K)

    double precision, intent(in) :: Kappa_vals(2:21,2:21)
    double precision, intent(in) :: Log_D_list(2:21)
    double precision, intent(in) :: Log_T_list(2:21)

    double precision, intent(in) :: D, T
    double precision, intent(out) :: K

    integer :: low_r_index, up_r_index
    integer :: low_t_index, up_t_index

    if (D .gt. maxval(Log_D_list)) then
        ! print*, 'Extrapolating in logR'
        low_r_index = 20
        up_r_index = 21
    elseif (D .lt. minval(Log_D_list)) then
        ! print*, 'Extrapolating in logR'
        low_r_index = 2
        up_r_index = 3
    else
        call get_low_up_index(D, Log_D_list, 2, 21, low_r_index, up_r_index)
    endif

    if (T .gt. maxval(Log_T_list)) then
        ! print*, 'Extrapolating in logT'
        low_t_index = 20
        up_t_index = 21
    elseif ( T .lt. minval(Log_T_list)) then
        ! print*, 'Extrapolating in logT'
        low_t_index = 2
        up_t_index = 3
    else
        call get_low_up_index(T, Log_T_list, 2, 21, low_t_index, up_t_index)
    endif

    call interpolate_KRT(low_r_index, up_r_index, low_t_index, up_t_index, Log_D_list, Log_T_list, Kappa_vals, D, T, K)

end subroutine get_val


!> this subroutine finds the indexes in R and T arrays of the two values surrounding the input R and T
subroutine get_low_up_index(x, x_list, imin, imax, low_i, up_i)

    integer, intent(in) :: imin, imax
    double precision, intent(in) :: x
    double precision, intent(in) :: x_list(imin:imax)

    integer, intent(out) :: low_i, up_i

    double precision :: low_val, up_val
    double precision :: res(imin:imax)

    res = x_list - x

    up_val = minval(res, MASK = res .ge. 0) + x
    low_val = maxval(res, MASK = res .le. 0) + x

    up_i = minloc(abs(x_list - up_val),1) + imin -1
    low_i = minloc(abs(x_list - low_val),1) + imin -1


    if (up_i .eq. low_i) low_i = low_i - 1

end subroutine get_low_up_index

!> This subroutine does a bilinear interpolation in the R,T-plane
subroutine interpolate_KRT(low_r, up_r, low_t, up_t, Log_D_list, Log_T_list, Kappa_vals, D, T, k_interp)

    integer, intent(in) :: low_r, up_r, low_t, up_t
    double precision, intent(in) :: Kappa_vals(2:21,2:21)
    double precision, intent(in) :: Log_D_list(2:21)
    double precision, intent(in) :: Log_T_list(2:21)
    double precision, intent(in) :: D,T
    double precision, intent(out) :: k_interp

    double precision :: r1,r2,t1,t2
    double precision :: k1, k2, k3, k4
    double precision :: ka, kb

    !Cool ascii drawing of interpolation scheme: first interpolate twice in the T coord to get
    !ka and kb, then interpolate in the R coord to get ki

!   r_1    R        r_2
!     |                |
!     |                |
! ----k1--------ka-----k2----- t_1
!     |          |     |
!     |          |     |
!   T |          |     |
!     |          |     |
!     |          ki    |
!     |          |     |
! ----k3--------kb-----k4----- t_2
!     |                |
!     |                |


    r1 = Log_D_list(low_r)
    r2 = Log_D_list(up_r)
    t1 = Log_T_list(low_t)
    t2 = Log_T_list(up_t)

    ! k1 = Kappa_vals(low_t, low_r)
    ! k2 = Kappa_vals(low_t, up_r)
    ! k3 = Kappa_vals(up_t, low_r)
    ! k4 = Kappa_vals(up_t, up_r)

    k1 = Kappa_vals(low_r, low_t)
    k2 = Kappa_vals(low_r, up_t)
    k3 = Kappa_vals(up_r, low_t)
    k4 = Kappa_vals(up_r, up_t)

    ! print*, 'surounding indexes'
    ! print*, low_r, up_r, low_t, up_t
    ! print*, 'surounding input values'
    ! print*, r1,r2,t1,t2
    ! print*, 'surounding table values'
    ! print*, k1, k2, k3, k4

    call interpolate1D(r1,r2,D,k1,k2,ka)
    call interpolate1D(r1,r2,D,k3,k4,kb)
    call interpolate1D(t1,t2,T,ka,kb,k_interp)

    ! print*, 'interpolated value'
    ! print*, ka, kb, k_interp

end subroutine interpolate_KRT


!> Interpolation in one dimension
subroutine interpolate1D(x1, x2, x, y1, y2, y)

    double precision, intent(in) :: x, x1, x2
    double precision, intent(in) :: y1, y2
    double precision, intent(out) :: y

    y = y1 + (x-x1)*(y2-y1)/(x2-x1)

end subroutine interpolate1D

end module mod_cak_opacity

! !> Interpolation on logarithmic scale
! subroutine log_interpolate1D(x1, x2, x, y1, y2, y)
!
!    double precision, intent(in) :: x, x1, x2
!    double precision, intent(in) :: y1, y2
!    double precision, intent(out) :: y
!
!    double precision :: expx, expx1, expx2
!    double precision :: expy1, expy2
!
!    expx = 10**x
!    expx1 = 10**x1
!    expx2 = 10**x2
!
!    expy1 = 10**y1
!    expy2 = 10**y2
!
!    y = expy1 + (expx-expx1)*(expy2-expy1)/(expx2-expx1)
!    y = log10(y)
!
! end subroutine log_interpolate1D
