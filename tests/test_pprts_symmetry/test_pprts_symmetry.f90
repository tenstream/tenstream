@test(npes =[1])
subroutine test_pprts_symmetry_ex1(this)

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, zero, one, pi

    use m_pprts, only : init_pprts, set_optical_properties, t_solver_3_6, t_solver_8_10, &
      solve_pprts, set_angles, pprts_get_result, destroy_pprts
    use m_tenstream_options, only: read_commandline_options

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this


    integer(iintegers) :: numnodes, comm
    integer(iintegers) :: myid

    integer(iintegers),parameter :: nxp=9,nyp=9,nv=100
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=0, phi1=180, theta0=60
    real(ireals),parameter :: albedo=0.1, dz=dx
    real(ireals),parameter :: incSolar = 1000
    real(ireals),parameter :: atolerance = 1
    real(ireals) :: dz1d(nv)

    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g,B
    real(ireals),allocatable,dimension(:,:,:) :: fdir0,fdn0,fup0,fdiv0

    real(ireals),allocatable,dimension(:,:,:) :: fdir1,fdn1,fup1,fdiv1

    integer(iintegers) :: i,j,k

    !type(t_solver_3_6)  :: solver
    type(t_solver_8_10) :: solver

    dz1d = dz

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    !!!!!!!!!!!!!!!!!!!!  calculation for phi0 = 0  !!!!!!!!!!!!!!!!!!!!!!!
    call init_pprts(comm, nv, nxp, nyp, dx,dy, phi0, theta0, solver, dz1d)

    allocate(fdir0 (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fdn0  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fup0  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fdiv0 (solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))

    allocate(fdir1 (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fdn1  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fup1  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(fdiv1 (solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))


    allocate(kabs(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(ksca(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(g   (solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))

    kabs = 1._ireals/nv/dz
    ksca = 1._ireals/nv/dz
    g    = zero

    kabs(nv/2,nxp/2,1:nyp) = 1/dz
    ksca(nv/2,nxp/2,1:nyp) = 1/dz
    g   (nv/2,nxp/2,1:nyp) = .9

    call set_optical_properties(solver, albedo, kabs, ksca, g)!, B )
    call solve_pprts(solver, incSolar, opt_solution_uid=100)

    call set_angles(solver, phi0, theta0)
    call solve_pprts(solver, incSolar, opt_solution_uid=101)

    call pprts_get_result(solver, fdir0,fdn0,fup0,fdiv0, opt_solution_uid=100)
    call pprts_get_result(solver, fdir1,fdn1,fup1,fdiv1, opt_solution_uid=101)

    if(myid.eq.0) then
        print *,'kabs:', kabs(:,1,1)
        print *,'fdn:',  fdn1 (:,1,1)
        print *,'fup:',  fup1 (:,1,1)
        print *,'fdiv:', fdiv1(:,1,1)
    endif

    do j=lbound(fdiv0,3), ubound(fdiv0,3)
      do i=lbound(fdiv0,2), ubound(fdiv0,2)
        do k=lbound(fdiv0,1), ubound(fdiv0,1)
          !print *,k,i,j,'::', fdiv0(k,i,ubound(fdiv0,2)-j+lbound(fdiv0,2)), fdiv1(k,i,j)
          @assertEqual(fdiv0(k,i,ubound(fdiv0,2)-j+lbound(fdiv0,2)), fdiv1(k,i,j), atolerance*1e-3, 'divergence not symmetric for azimuth')
        enddo
      enddo
    enddo

    do j=lbound(fdir0,3), ubound(fdir0,3)
      do i=lbound(fdir0,2), ubound(fdir0,2)
        do k=lbound(fdir0,1), ubound(fdir0,1)
          @assertEqual(fdir0(k,i,ubound(fdir0,2)-j+lbound(fdir0,2)), fdir1(k,i,j), atolerance, 'Edirradiation not symmetric for azimuth')
          @assertEqual(fdn0(k,i,ubound(fdn0,2)-j+lbound(fdn0,2)), fdn1(k,i,j), atolerance,     'Edn radiation not symmetric for azimuth')
          @assertEqual(fup0(k,i,ubound(fup0,2)-j+lbound(fup0,2)), fup1(k,i,j), atolerance,     'Eup radiation not symmetric for azimuth')
        enddo
      enddo
    enddo


    ! Check that surface emission is the one that we stick in
    !@assertEqual(B(ubound(B,1),1,1)*pi, fup (ubound(fup,1),1,1), atolerance, 'Surface Emission not correct')

    call destroy_pprts(solver, .True.)
end subroutine
