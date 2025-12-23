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
      call advect1(flux_method,one,idim^LIM,global_time,ps1,global_time,ps)

      !> Reset values that are NOT to be substepped
      !> THIS IS NOT IDEAL AND SHOULD BE DONE OTHERWISE
      ps() = ps1()

      !> Call sourceterms that need to be substepped
      call substep_source()

    end do

  end subroutine advance_substep

end module mod_substep
