!> This module reads in opacities from opal tables.

module mod_opacity
    implicit NONE

    integer, parameter :: rmin = 2
    integer, parameter :: rmax = 20
    integer, parameter :: tmin = 7
    integer, parameter :: tmax = 76

    real, public :: Kappa_vals(7:76,2:20)
    real, public :: Kappa_vals1(7:76,2:20)
    real, public :: Kappa_vals2(7:76,2:20)
    real, public :: Log_R_list(2:20)
    real, public :: Log_T_list(7:76)

    public :: init_opal
    public :: set_opal_opacity


        !> This is a fortran program to read in one of the tables found in the GN93hz Opal opacities files.
        ! This only works for fixed metalicity. Part of the table has to be copied to a file called 'My_table'
        !
        ! Update: One can now interpolate between two tables as well
        !
        !
        ! integer, parameter :: rmin = 2
        ! integer, parameter :: rmax = 20
        ! integer, parameter :: tmin = 7
        ! integer, parameter :: tmax = 76
        !
        ! integer :: i
        !
        ! real :: Kappa_vals(7:76,2:20)
        ! real :: Kappa_vals1(7:76,2:20)
        ! real :: Kappa_vals2(7:76,2:20)
        ! real :: Log_R_list(2:20)
        ! real :: Log_T_list(7:76)
        !
        ! real :: R_input, T_input
        ! real :: K_output
        !
        ! call read_table(Log_R_list, Log_T_list, Kappa_vals1,'My_table1')
        ! call read_table(Log_R_list, Log_T_list, Kappa_vals2,'My_table2')
        !
        ! call interpolate_two_tables(1.000,9.990, 9.995, Kappa_vals1, Kappa_vals2, Kappa_vals)
        !
        ! !print*, Kappa_vals2(7:76,2:20)
        !
        !
        ! do i = 7,76
        !     print*, Kappa_vals1(i,:)
        ! enddo
        !
        ! R_input = 0.0
        ! T_input = 8.0
        !
        ! call get_kappa(Kappa_vals, Log_R_list, Log_T_list, R_input, T_input, K_output)
        !
        ! !> If the outcome is 9.999, look right in the table
        ! do while (K_output .gt. 9.0)
        !     print*, 'K > 9', K_output, R_input, T_input
        !     R_input = R_input + 0.5
        !     call get_kappa(Kappa_vals, Log_R_list, Log_T_list, R_input, T_input, K_output)
        ! enddo
        !
        ! !> If the outcome is NaN, look left in the table
        ! do while (K_output .eq. 0.0)
        !     print*, 'K = NaN', K_output, R_input, T_input
        !     R_input = R_input - 0.5
        !     call get_kappa(Kappa_vals, Log_R_list, Log_T_list, R_input, T_input, K_output)
        ! enddo
        !
        ! print*, K_output

subroutine init_opal(He_abundance)

  double precision, intent(in) :: He_abundance

  call read_table(Log_R_list, Log_T_list, Kappa_vals1,'My_table1')
  call read_table(Log_R_list, Log_T_list, Kappa_vals2,'My_table2')

  call interpolate_two_tables(1.000,9.990, He_abundance, Kappa_vals1, Kappa_vals2, Kappa_vals)

end subroutine init_opal()


subroutine set_opal_opacity(rho,temp,kappa)
  double precision, intent(in) :: rho, temp
  double precision, intent(out) :: kappa

  double precision :: R_input, T_input, K_output

  R = rho
  T = temp

  call get_kappa(Kappa_vals, Log_R_list, Log_T_list, R_input, T_input, K_output)
  
  !> If the outcome is 9.999, look right in the table
  do while (K_output .gt. 9.0)
      print*, 'K > 9', K_output, R_input, T_input
      R_input = R_input + 0.5
      call get_kappa(Kappa_vals, Log_R_list, Log_T_list, R_input, T_input, K_output)
  enddo

  !> If the outcome is NaN, look left in the table
  do while (K_output .eq. 0.0)
      print*, 'K = NaN', K_output, R_input, T_input
      R_input = R_input - 0.5
      call get_kappa(Kappa_vals, Log_R_list, Log_T_list, R_input, T_input, K_output)
  enddo

end subroutine set_opal_opacity

subroutine read_table(R, T, K, filename)
    !> This routine reads in the the values for log kappa, and the values for log T and log R on the x and y axis

    real, intent(out) :: K(7:76,2:20), R(2:20), T(7:76)
    character(len=9), intent(in) :: filename

    character :: dum
    integer :: row, col

    OPEN(1,FILE=filename)

    !> Skip first 4 lines
    do row = 1,4
        READ(1,*)
    enddo

    !> Read R
    READ(1,*) dum,R(2:20)

    READ(1,*)

    !> Read T and K
    do row = 7,76 !> NOT READING ENTIRE TABLE
        READ(1,'(f4.2,19f7.3)') T(row), K(row,2:20)
    enddo

    CLOSE(1)

end subroutine read_table



subroutine interpolate_two_tables(Y1, Y2, Y_in, K1, K2, K_interp)
    !> This subroutine creates a new table for a given metalicity,
    ! by interpolating to known tables at every point in the R,T plane
    real, intent(in) :: K1(7:76,2:20), K2(7:76,2:20)
    real, intent(in) :: Y1, Y2, Y_in
    real, intent(out) :: K_interp(7:76,2:20)

    integer row, colum

    do colum=2,20
    do row=7,76
        call interpolate1D(Y1,Y2,Y_in,K1(row,colum),K2(row,colum),K_interp(row,colum))
    enddo
    enddo

end subroutine interpolate_two_tables



subroutine get_kappa(Kappa_vals, Log_R_list, Log_T_list, R, T, K)

    !>This subroutine looks in the table for the four couples (T,R)
    !surrounding a given input for T and R

    real, intent(in) :: Kappa_vals(7:76,2:20)
    real, intent(in) :: Log_R_list(2:20)
    real, intent(in) :: Log_T_list(7:76)

    real, intent(in) :: R, T
    real, intent(out) :: K

    integer :: low_r_index, up_r_index
    integer :: low_t_index, up_t_index

    if (R .ge. maxval(Log_R_list)) then
        print*, 'Extrapolating in logR'
        low_r_index = 20
        up_r_index = 20
    elseif (R .le. minval(Log_R_list)) then
        print*, 'Extrapolating in logR'
        low_r_index = 2
        up_r_index = 2
    else
        call get_low_up_index(R, Log_R_list, 2, 20, low_r_index, up_r_index)
    endif

    if (T .ge. maxval(Log_T_list)) then
        print*, 'Extrapolating in logT'
        low_t_index = 76
        up_t_index = 76
    elseif ( T .le. minval(Log_T_list)) then
        print*, 'Extrapolating in logT'
        low_t_index = 7
        up_t_index = 7
    else
        call get_low_up_index(T, Log_T_list, 7, 76, low_t_index, up_t_index)
    endif

    call interpolate_KRT(low_r_index, up_r_index, low_t_index, up_t_index, Log_R_list, Log_T_list, Kappa_vals, R, T, K)

end subroutine get_kappa



subroutine get_low_up_index(x, x_list, imin, imax, low_i, up_i)

    !> this subroutine finds the indexes in R and T arrays of the two values surrounding the input R and T

    integer, intent(in) :: imin, imax
    real, intent(in) :: x
    real, intent(in) :: x_list(imin:imax)

    integer, intent(out) :: low_i, up_i

    real :: low_val, up_val
    real :: res(imin:imax)

    res = x_list - x

    up_val = minval(res, MASK = res .ge. 0) + x
    low_val = maxval(res, MASK = res .le. 0) + x

    up_i = minloc(abs(x_list - up_val),1) + imin -1
    low_i = minloc(abs(x_list - low_val),1) + imin -1


    if (up_i .eq. low_i) low_i = low_i - 1

end subroutine get_low_up_index


subroutine interpolate_KRT(low_r, up_r, low_t, up_t, Log_R_list, Log_T_list, Kappa_vals, R, T, k_interp)

    !> This subroutine does a bilinear interpolation in the R,T-plane

    integer, intent(in) :: low_r, up_r, low_t, up_t
    real, intent(in) :: Kappa_vals(7:76,2:20)
    real, intent(in) :: Log_R_list(2:20)
    real, intent(in) :: Log_T_list(7:76)
    real, intent(in) :: R,T
    real, intent(out) :: k_interp

    real :: r1,r2,t1,t2
    real :: k1, k2, k3, k4
    real :: ka, kb

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


    r1 = Log_R_list(low_r)
    r2 = Log_R_list(up_r)
    t1 = Log_T_list(low_t)
    t2 = Log_T_list(up_t)

    k1 = Kappa_vals(low_t, low_r)
    k2 = Kappa_vals(low_t, up_r)
    k3 = Kappa_vals(up_t, low_r)
    k4 = Kappa_vals(up_t, up_r)

    call interpolate1D(r1,r2,R,k1,k2,ka)
    call interpolate1D(r1,r2,R,k3,k4,kb)
    call interpolate1D(t1,t2,T,ka,kb,k_interp)

!     print*, '----------------------------------------------------------'
!     print*, '----------------------------------------------------------'
!     print*, '----------------------------------------------------------'
!     print*, r1, R, r2
!     print*, t1, T, t2
!     print*, '----------------------------------------------------------'
!     print*, k1, ka, k2
!     print*, k3, kb, k4
!     print*, '----------------------------------------------------------'
!     print*, ka, k_interp, kb

end subroutine interpolate_KRT



subroutine interpolate1D(x1, x2, x, y1, y2, y)

    real, intent(in) :: x, x1, x2
    real, intent(in) :: y1, y2
    real, intent(out) :: y

    y = y1 + (x-x1)*(y2-y1)/(x2-x1)

end subroutine interpolate1D

end module mod_opacity

!subroutine log_interpolate1D(x1, x2, x, y1, y2, y)
!
!    real, intent(in) :: x, x1, x2
!    real, intent(in) :: y1, y2
!    real, intent(out) :: y
!
!    real :: expx, expx1, expx2
!    real :: expy1, expy2
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
!end subroutine log_interpolate1D
