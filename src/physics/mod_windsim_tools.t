module mod_windsim_tools
implicit none

  !> average photosphere
  double precision, public :: R_phot
  integer, public :: i_phot

  !> Grid on which to compute photosphere
  double precision, allocatable :: radial_grid(:)
  double precision, allocatable :: dr(:)

  !> Opacity value to use to compute photosphere
  double precision :: wst_kappa

  !> index of density
  integer :: rho_

  !> Number of points on which to evaluate the photospheric surface
  integer :: np1

  public :: init_windsim_tools
  public :: wst_get_photosphere

contains

  subroutine init_windsim_tools(i_rho,kappa_in)
    use mod_physics

    double precision, intent(in) :: kappa_in
    integer, intent(in) :: i_rho

    wst_kappa = kappa_in
    rho_ = i_rho

    !> Point to allocation subroutine
    phys_before_main_loop => allocate_photosphere_grid

  end subroutine init_windsim_tools

  subroutine allocate_photosphere_grid
    use mod_global_parameters

    integer :: i_r

    ! print*, "Grid info: xmin, xmax, nx, qs"
    ! print*, xprobmin1, xprobmax1, domain_nx1, qstretch_baselevel(1)

    np1 = domain_nx1

    allocate(radial_grid(1:np1))
    allocate(dr(1:np1))

    !> Recalculate grid at cell bottom
    dr(1) = (xprobmax1 - xprobmin1)*(1-qstretch_baselevel(1))/(1-qstretch_baselevel(1)**(np1))

    !> Cell center r coordinate
    radial_grid(1) = xprobmin1 !+ 0.5d0*dr(1)
    do i_r = 2, np1
      dr(i_r) = dr(i_r - 1)*qstretch_baselevel(1)
      ! radial_grid(i_r) = radial_grid(i_r - 1) + 0.5d0*(dr(i_r - 1) +  dr(i_r))
      radial_grid(i_r) = radial_grid(i_r - 1) + dr(i_r - 1)
    end do

    !> initial calculation of photosphere on initial conditions
    call wst_get_photosphere(0,0.d0)

  end subroutine allocate_photosphere_grid

  !> This subroutine resamples the density to level one and takes the lateral average.
  !> Then, it computes the optical depth and find the cell with optical depth closest to 2/3
  !> It then stores the radial location of that cell in public variable R_phot
  subroutine wst_get_photosphere(iit,qt)
    use mod_global_parameters

    integer, intent(in)          :: iit
    double precision, intent(in) :: qt

    double precision :: dtau(np1), tau(np1)
    double precision :: dtau_S(np1), dtau_R(np1)
    double precision :: rho(np1)

    double precision :: tau_ph = 2.d0/3

    integer :: iigrid, igrid
    integer :: ii, jj, il, ih, jj_blk, nc

    if (.not. convert) then

      rho(:) = 0.d0

      !> loop over all grids for this proc
      do iigrid=1,igridstail; igrid=igrids(iigrid);
        block=>ps(igrid)

          !> cell widht in cells on lvl 1: On lvl 1: 1, lvl 2: 2, lvl 3: 4
          nc = 2**(node(plevel_,igrid)-1)

          !> Map the block index to the global index
          !> jj_blk = number of cells on level one between xprobmin1 and the beginning of this block
          jj_blk = (node(pig1_,igrid)-1)*block_nx1/nc

          ! print*, "jj_blk", jj_blk

          !> For all cells in the current block, average density, opacity, volume, ... over lateral direction.
          !> Take into account possible amr
          do ii = 1,block_nx1/nc
            !> relate global index jj to local index ii on grid
            il = nghostcells + 1 + (ii-1)*nc !> il and ih are the lower and upper cell indexes that fit within equivalent cell jj on lvl one
            ih = nghostcells + 1 + (ii)*nc -1
            jj = jj_blk + ii !> jj is the coordinate on level 1

            ! print*, "il, ih, jj", il, ih, jj

            !> Density
            {^IFONED
            rho(jj) = rho(jj) &
            + sum(block%w(il:ih,rho_)*block%dvolume(il:ih))/sum(block%dvolume(il:ih))
            ! + sum(block%w(il:ih,rho_))/nc !> average value on lvl 1
            }

            {^IFTWOD
            rho(jj) = rho(jj) &
            + sum(block%w(il:ih,ixMlo2:ixMhi2,rho_)*block%dvolume(il:ih,ixMlo2:ixMhi2))/sum(block%dvolume(il:ih,ixMlo2:ixMhi2)) !> average value on lvl 1
            ! + sum(block%w(il:ih,ixMlo2:ixMhi2,rho_))/domain_nx2/nc**2 !> average value on lvl 1
            }

            {^IFTHREED
            rho(jj) = rho(jj) &
            + sum(block%w(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,rho_)*block%dvolume(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3))/sum(block%dvolume(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3)) !> average value on lvl 1
            ! + sum(block%w(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,rho_))/(domain_nx2*domain_nx3)/nc**3 !> average value on lvl 1
            }

          end do
      end do

      !> Calculate increment of optical depth
      ! dtau(1:np1) = rho(1:np1)*unit_density*wst_kappa*unit_opacity*dr(1:np1)*unit_length
      dtau(1:np1) = rho(1:np1)*wst_kappa*dr(1:np1)

      ! print*, "dtau"
      ! print*, dtau

      !> Communicate optical depth additions over processors
      dtau_S = dtau
      ! call mpi_allreduce(dtau_S,dtau_R,np1,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
      call mpi_reduce(dtau_S,dtau_R,np1,MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
      dtau = dtau_R

      !> Calculate optical depth radius
      if (mype == 0) then
        !> Optical depth integral
        tau(np1) = dtau(np1)
        do ii = np1-1, 1, -1
          tau(ii) = tau(ii+1) + dtau(ii)
        end do

        ! print*, "ii, r, rho, tau, dr, dtau"
        ! do ii = np1, 1, -1
        !   print*, ii, radial_grid(ii), rho(ii), tau(ii), dr(ii), dtau(ii)
        ! end do

        !> cell closest to 2/3
        i_phot = minloc(abs(tau(1:np1) - tau_ph),DIM=1)
        R_phot = radial_grid(i_phot)

        !> With continuous interpolation
        if ((tau(i_phot) - tau_ph) * (tau(i_phot + 1) - tau_ph) .le. 0) then
          R_phot = radial_grid(i_phot)  &
                 + (radial_grid(i_phot + 1) - radial_grid(i_phot)) &
                 / (log10(tau(i_phot+1)) - log10(tau(i_phot))) &
                 * (log10(tau_ph) - log10(tau(i_phot)))
        else if ((tau(i_phot) - tau_ph) * (tau(i_phot - 1) - tau_ph) .le. 0) then
          R_phot = radial_grid(i_phot)  &
                 + (radial_grid(i_phot - 1) - radial_grid(i_phot)) &
                 / (log10(tau(i_phot - 1)) - log10(tau(i_phot))) &
                 * (log10(tau_ph) - log10(tau(i_phot)))
        else
          stop
        end if



      end if

      !> Broadcast photospheric radius
      call mpi_bcast(R_phot, 1, MPI_DOUBLE_PRECISION, 0, icomm, ierrmpi)

    end if

  end subroutine wst_get_photosphere

end module mod_windsim_tools
