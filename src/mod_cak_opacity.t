!> This module reads in the line-force (CAK) parameters alpha, Qbar, Q0 and the
!> electron-scattering opacity kappa_e, tabulated as a function of log10(density)
!> and log10(temperature) on a regular grid.
!>
!> The four tables are read once on the root rank, broadcast, and stored as four
!> columns of a single 2D lookup table (mod_lookup_table, LT2_t), so that a
!> look-up locates the (logD, logT) cell once for all four quantities.
!> Look-ups outside the table range are clamped to the table edge.

module mod_cak_opacity

  use mod_lookup_table
  use mod_comm_lib, only: mpistop

  implicit none
  private

  !> Size of the CAK tables: number of logD rows and logT columns
  integer, parameter :: nD_cak = 20, nT_cak = 20

  !> Column indices in the lookup table
  integer, parameter :: i_alpha = 1, i_Qbar = 2, i_Q0 = 3, i_kappa_e = 4, n_cak = 4

  !> Lookup table, coordinates (logD, logT), columns alpha, Qbar, Q0, kappa_e
  type(LT2_t) :: cak_lt

  public :: init_cak
  public :: set_cak_opacity

contains

  !> Read the CAK tables on the root rank, broadcast them and build the lookup table.
  !> Called when the radiation module is initialised.
  subroutine init_cak(tablename)
    use mod_global_parameters, only: mype, icomm, ierrmpi
    use mpi

    character(len=*), intent(in) :: tablename

    character(len=*), parameter :: files(n_cak) = ['al_TD', 'Qb_TD', 'Q0_TD', 'Ke_TD']
    character(len=:), allocatable :: AMRVAC_DIR, path_table_dir
    double precision :: data(nD_cak, nT_cak, n_cak)
    double precision :: logD(nD_cak), logT(nT_cak), logD_i(nD_cak), logT_i(nT_cak)
    integer          :: lenv, icol

    if (mype == 0) then
      call get_environment_variable(name="AMRVAC_DIR", length=lenv)
      allocate(character(len=lenv) :: AMRVAC_DIR)
      call get_environment_variable(name="AMRVAC_DIR", value=AMRVAC_DIR)

      path_table_dir = trim(AMRVAC_DIR)//"/src/rhd/CAK_tables/"//trim(tablename)//"/"
      print *, 'CAK tables: ', path_table_dir

      do icol = 1, n_cak
        call read_table(logD_i, logT_i, data(:,:,icol), path_table_dir//trim(files(icol)))
        if (icol == 1) then
          logD = logD_i
          logT = logT_i
        else if (any(abs(logD_i - logD) > 1.d-6) .or. any(abs(logT_i - logT) > 1.d-6)) then
          call mpistop("init_cak: CAK tables have different (logD, logT) axes")
        end if
      end do

      call check_uniform(logD, "logD")
      call check_uniform(logT, "logT")
    end if

    call MPI_BCAST(logD, nD_cak, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)
    call MPI_BCAST(logT, nT_cak, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)
    call MPI_BCAST(data, nD_cak*nT_cak*n_cak, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)

    cak_lt = LT2_create_from_data([logD(1), logT(1)], [logD(nD_cak), logT(nT_cak)], data)

  end subroutine init_cak

  !> Line-force parameters and electron-scattering opacity for a given density
  !> and temperature (cgs). gradv is currently unused.
  subroutine set_cak_opacity(rho, temp, gradv, alpha_output, Qbar_output, Q0_output, kappa_e_output)

    double precision, intent(in)  :: rho, temp, gradv
    double precision, intent(out) :: alpha_output, Qbar_output, Q0_output, kappa_e_output

    type(LT2_loc_t) :: loc

    loc = LT2_get_loc(cak_lt, log10(rho), log10(temp))

    alpha_output   = LT2_get_col_at_loc(cak_lt, i_alpha, loc)
    Qbar_output    = LT2_get_col_at_loc(cak_lt, i_Qbar, loc)
    Q0_output      = LT2_get_col_at_loc(cak_lt, i_Q0, loc)
    kappa_e_output = LT2_get_col_at_loc(cak_lt, i_kappa_e, loc)

  end subroutine set_cak_opacity

  !> Read one CAK table: header with a dummy value and logT, then rows of logD and values
  subroutine read_table(D, T, K, filename)

    character(*), intent(in)      :: filename
    double precision, intent(out) :: D(nD_cak), T(nT_cak), K(nD_cak, nT_cak)

    double precision :: dum
    integer          :: row, unit_cak

    open(newunit=unit_cak, status='old', file=trim(filename))

    read(unit_cak,*) dum, T(1:nT_cak)
    do row = 1, nD_cak
      read(unit_cak,*) D(row), K(row,1:nT_cak)
    end do

    close(unit_cak)

  end subroutine read_table

  !> The LT2 lookup assumes a uniform grid
  subroutine check_uniform(x, name)
    double precision, intent(in) :: x(:)
    character(len=*), intent(in) :: name
    double precision :: dx

    dx = (x(size(x)) - x(1))/(size(x) - 1)
    if (any(abs(x(2:) - x(:size(x)-1) - dx) > 1.d-6*abs(dx))) then
      call mpistop("init_cak: non-uniform "//name//" axis in CAK tables")
    end if
  end subroutine check_uniform

end module mod_cak_opacity
