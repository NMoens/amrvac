!> This module reads in Rosseland-mean opacities from OPAL tables. Table opacity
!> values are given in base 10 logarithm and are a function of mass density (R)
!> and temperature (T), which are also both given in base 10 logarithm.
!>
!> The table is read once on the root rank, broadcast, and stored in a 2D
!> lookup table (mod_lookup_table, LT2_t) on a regular (logT, logR) grid.
!> OPAL tables have a non-uniform logT spacing (0.05/0.1/0.2), but all nodes
!> lie on the finest spacing, so the resampling onto the regular grid is exact
!> at the table nodes and reproduces bilinear interpolation in between.
!>
!> Missing table entries (blank fields in the bottom-right corner of the table,
!> or 9.999) are filled once at initialisation by extending each logT row at
!> constant value along logR. Look-ups outside the table range are clamped to
!> the table edge (default behaviour of LT2_get_col).

module mod_opal_opacity

  use mod_lookup_table
  use mod_comm_lib, only: mpistop

  implicit none
  private

  !> Size of an OPAL table: number of logT rows and logR columns
  integer, parameter :: nT_opal = 70, nR_opal = 19

  !> Field widths in the OPAL table rows: logT (f4.2) followed by nR_opal values (f7.3)
  integer, parameter :: wT_opal = 4, wK_opal = 7

  !> Lookup table, coordinates (logT, logR), one column: log10(kappa)
  type(LT2_t) :: opal_lt

  public :: init_opal_table
  public :: set_opal_opacity

contains

  !> Read the OPAL table on the root rank, broadcast it and build the lookup table.
  !> Called when the radiation module (FLD, AFLD, M1) is initialised.
  subroutine init_opal_table(tablename)
    use mod_global_parameters, only: mype, icomm, ierrmpi
    use mpi

    character(len=*), intent(in) :: tablename

    character(len=:), allocatable :: AMRVAC_DIR, path_table_dir
    double precision :: logT_list(nT_opal), logR_list(nR_opal), logK(nT_opal, nR_opal)
    double precision :: dlogT, x_min(2), x_max(2)
    integer          :: lenv, n_reg(2)

    if (mype == 0) then
      call get_environment_variable(name="AMRVAC_DIR", length=lenv)
      allocate(character(len=lenv) :: AMRVAC_DIR)
      call get_environment_variable(name="AMRVAC_DIR", value=AMRVAC_DIR)

      ! Absolute path to OPAL table location
      path_table_dir = trim(AMRVAC_DIR)//"/src/rhd/OPAL_tables/"//trim(tablename)
      print *, 'OPAL opacity table: ', path_table_dir

      call read_table(logR_list, logT_list, logK, path_table_dir)
    end if

    call MPI_BCAST(logT_list, nT_opal, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)
    call MPI_BCAST(logR_list, nR_opal, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)
    call MPI_BCAST(logK, nT_opal*nR_opal, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)

    !> Regular grid: logT at the finest table spacing, logR at the (uniform) table spacing
    dlogT = minval(logT_list(2:nT_opal) - logT_list(1:nT_opal-1))
    if (any(abs(modulo((logT_list - logT_list(1))/dlogT + 0.5d0, 1.d0) - 0.5d0) > 1.d-6)) then
      call mpistop("init_opal_table: logT nodes do not lie on a regular grid")
    end if

    x_min    = [logT_list(1), logR_list(1)]
    x_max    = [logT_list(nT_opal), logR_list(nR_opal)]
    n_reg(1) = nint((x_max(1) - x_min(1))/dlogT) + 1
    n_reg(2) = nR_opal

    opal_lt = LT2_create(x_min, x_max, n_reg, 1)
    call LT2_set_col(opal_lt, 1, logT_list, logR_list, logK)

  end subroutine init_opal_table

  !> Rosseland-mean opacity (cgs) for a given density and temperature (cgs).
  !> Look-ups outside the table are clamped to the table edge.
  subroutine set_opal_opacity(rho, temp, kappa)

    double precision, intent(in)  :: rho, temp
    double precision, intent(out) :: kappa

    double precision :: logR_in, logT_in

    ! Convert input density and temperature to OPAL table convention
    logR_in = log10(rho/(temp*1d-6)**3)
    logT_in = log10(temp)

    kappa = 10d0**LT2_get_col(opal_lt, 1, logT_in, logR_in)

  end subroutine set_opal_opacity

  !> Read an OPAL table: logR list, logT list and log10(kappa)(logT, logR).
  !> Missing entries (blank fields, or 9.999) are filled along each logT row:
  !> 9.999 entries take the first valid value to their right, entries after the
  !> last valid value of a row take that last valid value.
  subroutine read_table(R, T, K, filename)

    character(*), intent(in)      :: filename
    double precision, intent(out) :: R(nR_opal), T(nT_opal), K(nT_opal, nR_opal)

    character(len=256) :: line
    character          :: dum
    logical            :: valid(nR_opal)
    integer            :: row, col, n_fields, unit_opal

    open(newunit=unit_opal, status='old', file=trim(filename))

    ! Skip first 4 lines
    do row = 1, 4
      read(unit_opal,*)
    end do

    ! Read logR
    read(unit_opal,*) dum, R(1:nR_opal)
    read(unit_opal,*)

    ! Read logT and log(kappa); rows may contain fewer than nR_opal values
    do row = 1, nT_opal
      read(unit_opal,'(a)') line
      read(line,'(f4.2)') T(row)
      n_fields = min(nR_opal, max(0, (len_trim(line) - wT_opal + wK_opal - 1)/wK_opal))
      do col = 1, nR_opal
        if (col <= n_fields) then
          read(line(wT_opal + (col-1)*wK_opal + 1 : wT_opal + col*wK_opal), *) K(row,col)
          valid(col) = (K(row,col) < 9.d0)
        else
          K(row,col) = 0.d0
          valid(col) = .false.
        end if
      end do

      if (.not. any(valid)) then
        call mpistop("read_table (OPAL): row without valid entries")
      end if

      ! 9.999 entries: take the first valid value to the right
      do col = nR_opal-1, 1, -1
        if (.not. valid(col) .and. col < n_fields .and. valid(col+1)) then
          K(row,col) = K(row,col+1)
          valid(col) = .true.
        end if
      end do

      ! Remaining entries (end of row): take the last valid value to the left
      do col = 2, nR_opal
        if (.not. valid(col) .and. valid(col-1)) then
          K(row,col) = K(row,col-1)
          valid(col) = .true.
        end if
      end do
    end do

    close(unit_opal)

  end subroutine read_table

end module mod_opal_opacity
