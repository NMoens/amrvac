module mod_substep

implicit none

contains

  subroutine advance_substep(dtfactor,qdt,qtC,psa,psb)

    double precision :: dt_substep, dt_factor_substep
    integer :: n_substeps, it_ss

    !> Get number of substeps
    call get_dt_substep(dt_substep)
    dt_factor_substep = dt

    !> Loop over number of substeps
    do it_ss = 1, n_substeps

      call advect_substep()

      call source_substep()

    end do



  end subroutine advance_substep

  subroutine get_dt_substep

  end subroutine get_dt_substep

  subroutine advect_substep

    qdt=dtfactor*dt

    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
      block=>ps(igrid)
      ^D&dxlevel(^D)=rnode(rpdx^D_,igrid);

      !> advect1_grid
      ixO^L=ixI^L^LSUBnghostcells;
      call tvdlimit(method,qdt,ixI^L,ixO^L,idim^LIM,sCT,qt+qdt,s,fC,dxs,x)

    end do
    !$OMP END PARALLEL DO

  end subroutine advect_substep

end module mod_substep
