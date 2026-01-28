module mod_substep

implicit none

contains

  subroutine advance_substep(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters

    double precision :: dt_substep, dt_factor_substep
    integer :: n_substeps, it_ss

    !> Get number of substeps
    call substep_get_dt(dt_substep)
    dt_factor_substep = dt

    !> Loop over number of substeps
    do it_ss = 1, n_substeps

      !> Advect substep
      call ss_advect1(flux_method,one,idim^LIM,global_time,ps,global_time,ps1)

      !> Implicit update substep
      call ss_global_implicit_update(one,dt,global_time+dt,ps,ps1)

      !> Call sourceterms that need to be substepped
      call ss_source()
    end do

  end subroutine advance_substep

  !> Integrate all grids by one partial step
  subroutine ss_advect1(method,dtfactor,idim^LIM,qtC,psa,qt,psb)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_fix_conserve
    use mod_physics

    integer, intent(in) :: idim^LIM
    type(state), target :: psa(max_blocks) !< Compute fluxes based on this state
    type(state), target :: psb(max_blocks) !< Update solution on this state
    double precision, intent(in) :: dtfactor !< Advance over dtfactor * dt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: qt
    integer, intent(in) :: method(nlevelshi)

    ! cell face flux
    double precision :: fC(ixG^T,1:nwflux,1:ndim)
    ! cell edge flux
    double precision :: fE(ixG^T,sdim:3)
    double precision :: qdt
    integer :: iigrid, igrid


    qdt=dtfactor*dt
    ! opedit: Just advance the active grids:
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

      call advect1_grid(method(block%level),qdt,dtfactor,ixG^LL,idim^LIM,&
        qtC,psa(igrid),qt,psb(igrid),fC,fE,rnode(rpdx1_:rnodehi,igrid),ps(igrid)%x)

      ! opedit: Obviously, flux is stored only for active grids.
      ! but we know in fix_conserve wether there is a passive neighbor
      ! but we know in conserve_fix wether there is a passive neighbor
      ! via neighbor_active(i^D,igrid) thus we skip the correction for those.
      ! This violates strict conservation when the active/passive interface
      ! coincides with a coarse/fine interface.
      if (fix_conserve_global .and. fix_conserve_at_step) then
        call store_flux(igrid,fC,idim^LIM,nwflux)
      end if

    end do
    !$OMP END PARALLEL DO

    if (fix_conserve_global .and. fix_conserve_at_step) then
      call recvflux(idim^LIM)
      call sendflux(idim^LIM)
      call fix_conserve(psb,idim^LIM,1,nwflux)
    end if

    ! For all grids: fill ghost cells
    call getbc(qt+qdt,qdt,psb,iwstart,nwgc,phys_req_diagonal)

  end subroutine advect1

end module mod_substep
