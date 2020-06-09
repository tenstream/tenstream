module test_pprts_symmetry

  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, irealLUT, zero, one, pi, mpiint

#include "petsc/finclude/petsc.h"
  use petsc

  use m_pprts_base, only : t_solver, t_solver_3_10, t_solver_8_10, destroy_pprts
  use m_pprts, only : init_pprts, set_optical_properties, &
    solve_pprts, set_angles, pprts_get_result_toZero
  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: itoa

  use m_optprop, only: t_optprop, t_optprop_cube, &
    t_optprop_3_10, t_optprop_8_10
  use pfunit_mod

  implicit none

  class(t_optprop_cube),allocatable :: OPP
  type(t_solver_3_10) :: solver_3_10
  type(t_solver_8_10) :: solver_8_10

contains
  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    continue
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    if(allocated(OPP)) then
      call OPP%destroy()
      deallocate(OPP)
    endif

    if(solver_8_10%linitialized) then
      call destroy_pprts(solver_8_10, lfinalizepetsc=.True.)
    endif
    if(solver_3_10%linitialized) then
      call destroy_pprts(solver_3_10, lfinalizepetsc=.True.)
    endif
  end subroutine teardown

  @test(npes = [1])
  subroutine test_pprts_symmetry_ex2(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm, ierr

    real(irealLUT), allocatable :: dir2dir(:), dir2diff(:)
    integer(iintegers) :: i

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    allocate (t_optprop_3_10 :: OPP)
    call this_test(OPP)
    deallocate(OPP)

    allocate (t_optprop_8_10 :: OPP)
    call this_test(OPP)
    deallocate(OPP)

    contains
      subroutine this_test(OPP)
        class(t_optprop_cube)  :: OPP

        call init_mpi_data_parameters(comm)

        call read_commandline_options(comm)

        call OPP%init(comm)

        allocate(dir2dir(OPP%OPP_LUT%dir_streams**2))
        allocate(dir2diff(OPP%OPP_LUT%diff_streams*OPP%OPP_LUT%dir_streams))

        do i=1,ubound(dir2dir,1)
          dir2dir(i) = real(i, kind=irealLUT)
        enddo

        do i=1,ubound(dir2diff,1)
          dir2diff(i) = real(i, kind=irealLUT)
        enddo

        call OPP%dir2dir_coeff_symmetry(dir2dir, .True., .True.)
        call OPP%dir2dir_coeff_symmetry(dir2dir, .True., .True.)

        do i=1,ubound(dir2dir,1)
          @mpiassertEqual(i,dir2dir(i), 'Coeff dir2dir not equal after switching two time north-south and east-west')
        end do

        call OPP%dir2diff_coeff_symmetry(dir2diff, .True., .True.)
        call OPP%dir2diff_coeff_symmetry(dir2diff, .True., .True.)

        do i=1,ubound(dir2diff,1)
          @mpiassertEqual(i,dir2diff(i), 'Coeff dir2diff not equal after switching two time north-south and east-west')
        end do

        deallocate(dir2dir)
        deallocate(dir2diff)
        call OPP%destroy()

        call PetscFinalize(ierr)
      end subroutine
    end subroutine

  @test(npes =[1])
  subroutine test_pprts_symmetry_ex1(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm

    integer(iintegers),parameter :: nxp=5,nyp=5,nv=5
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=10, theta0=60
    real(ireals),parameter :: albedo=0., dz=dx
    real(ireals),parameter :: incSolar=1000
    real(ireals),parameter :: atolerance=.1
    real(ireals) :: dz1d(nv)

    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir0,fdn0,fup0,fdiv0
    real(ireals),allocatable,dimension(:,:,:) :: fdir1,fdn1,fup1,fdiv1
    real(ireals),allocatable,dimension(:,:,:) :: fdir2,fdn2,fup2,fdiv2
    real(ireals),allocatable,dimension(:,:,:) :: fdir3,fdn3,fup3,fdiv3

    integer(iintegers) :: j
    integer(iintegers) :: cx, cy      ! global indices of cloud

    dz1d = dz

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call this_test(solver_3_10)

    contains
      subroutine this_test(solver)
      class(t_solver), intent(inout) :: solver
    call init_pprts(comm, nv, nxp, nyp, dx,dy, phi0, theta0, solver, dz1d)

    allocate(kabs(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(ksca(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(g   (solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))

    kabs = 1._ireals/nv/dz
    ksca = 1._ireals/nv/dz
    g    = zero

    cx = int(real(nxp)/2)+1
    cy = int(real(nyp)/2)+1

    if(cx.le.(solver%C_one%xe+1) .and. cx.gt.solver%C_one%xs) then
      if(cy.le.(solver%C_one%ye+1) .and. cy.gt.solver%C_one%ys) then
        kabs(2,  cx-solver%C_one%xs,  cy-solver%C_one%ys) = 1/dz
        ksca(2,  cx-solver%C_one%xs,  cy-solver%C_one%ys) = 1/dz
        g   (2,  cx-solver%C_one%xs,  cy-solver%C_one%ys) = .9
      endif
    endif

    call set_optical_properties(solver, albedo, kabs, ksca, g)
    call set_angles(solver, 10._ireals, theta0)
    call solve_pprts(solver, incSolar, opt_solution_uid=10_iintegers)

    call set_angles(solver, 190._ireals, theta0)
    call solve_pprts(solver, incSolar, opt_solution_uid=190_iintegers)

    call pprts_get_result_toZero(solver, fdn0, fup0, fdiv0, fdir0, opt_solution_uid=10_iintegers)
    call pprts_get_result_toZero(solver, fdn1, fup1, fdiv1, fdir1, opt_solution_uid=190_iintegers)

    call set_angles(solver, 100._ireals, theta0)
    call solve_pprts(solver, incSolar, opt_solution_uid=100_iintegers)

    call set_angles(solver, 280._ireals, theta0)
    call solve_pprts(solver, incSolar, opt_solution_uid=280_iintegers)

    call pprts_get_result_toZero(solver, fdn2, fup2, fdiv2, fdir2, opt_solution_uid=100_iintegers)
    call pprts_get_result_toZero(solver, fdn3, fup3, fdiv3, fdir3, opt_solution_uid=280_iintegers)


    if(myid.eq.0) then
      do j=lbound(fdir0,3), ubound(fdir0,3)
        print *, j, 'edir0', fdir0(4, :, j)
      enddo
      do j=lbound(fdir0,3), ubound(fdir0,3)
        print *, j, 'edir1', fdir1(4, :, j)
      enddo
      do j=lbound(fdir0,3), ubound(fdir0,3)
        print *, j, 'edn0', fdn0(4, :, j)
      enddo
      do j=lbound(fdir0,3), ubound(fdir0,3)
        print *, j, 'edn1', fdn1(4, :, j)
      enddo
      fdir1(:,:,:) = fdir1(:,nxp:1:-1,nyp:1:-1)
      fdn1 (:,:,:) = fdn1 (:,nxp:1:-1,nyp:1:-1)
      fup1 (:,:,:) = fup1 (:,nxp:1:-1,nyp:1:-1)
      fdiv1(:,:,:) = fdiv1(:,nxp:1:-1,nyp:1:-1)
      @assertEqual(fdiv0, fdiv1, atolerance, '10 -> 190: divergence not symmetric for azimuth')
      @assertEqual(fdir0, fdir1, atolerance, '10 -> 190: Edirradiation not symmetric for azimuth ')
      @assertEqual(fdn0 , fdn1 , atolerance, '10 -> 190: Edn radiation not symmetric for azimuth ')
      @assertEqual(fup0 , fup1 , atolerance, '10 -> 190: Eup radiation not symmetric for azimuth ')

      fdir3(:,:,:) = fdir3(:,nxp:1:-1,nyp:1:-1)
      fdn3 (:,:,:) = fdn3 (:,nxp:1:-1,nyp:1:-1)
      fup3 (:,:,:) = fup3 (:,nxp:1:-1,nyp:1:-1)
      fdiv3(:,:,:) = fdiv3(:,nxp:1:-1,nyp:1:-1)
      @assertEqual(fdiv2, fdiv3, atolerance, '10 -> 190: divergence not symmetric for azimuth')
      @assertEqual(fdir2, fdir3, atolerance, '10 -> 190: Edirradiation not symmetric for azimuth ')
      @assertEqual(fdn2 , fdn3 , atolerance, '10 -> 190: Edn radiation not symmetric for azimuth ')
      @assertEqual(fup2 , fup3 , atolerance, '10 -> 190: Eup radiation not symmetric for azimuth ')
    endif
    end subroutine
  end subroutine
end module
