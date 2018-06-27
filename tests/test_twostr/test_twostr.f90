module test_twostr

#include "petsc/finclude/petsc.h"
  use petsc


  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, zero, one, pi, mpiint

  use m_tenstream_options, only: read_commandline_options

  use m_twostream, only: delta_eddington_twostream, petsc_delta_eddington_twostream, adding_delta_eddington_twostream

  use pfunit_mod

  implicit none

contains

  @test(npes = [1])
  subroutine test_twostr_ex1(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: myid, numnodes, comm, ierr
    integer(iintegers), parameter :: ke1=10, ke=ke1-1
    real(ireals), dimension(ke) :: dtau, w0, g
    real(ireals), dimension(ke1) :: S, Edn, Eup
    real(ireals), dimension(ke1) :: pS, pEdn, pEup
    real(ireals) :: mu0, incSolar, albedo
    integer(iintegers) :: k

    integer(iintegers), parameter :: iter = 1000
    real(ireals) :: starttime, now

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    mu0 = .5; incSolar = 1e2; albedo = .1
    dtau = one/ke
    w0 = .5
    g = .5

    call cpu_time(starttime)
    do k=1,iter
      call delta_eddington_twostream(dtau, w0, g, mu0, incSolar, albedo, S, Edn, Eup)
    enddo
    call cpu_time(now)
    print *,'Time for delta_eddington_twostream:', now-starttime

    @assertEqual(incSolar, S(1))
    @assertEqual(exp(-sum(dtau)/mu0)*incSolar, S(ke1), sqrt(epsilon(one)))

    @assertEqual(zero, Edn(1))

    @assertEqual((S(ke1)+Edn(ke1))*albedo, Eup(ke1))


    PETSC_COMM_WORLD = comm
    call PetscInitialize(PETSC_NULL_CHARACTER ,ierr)
    call init_mpi_data_parameters(comm)

    call read_commandline_options()

    call petsc_delta_eddington_twostream(dtau, w0, g, mu0, incSolar, albedo, pS, pEdn, pEup)

    do k=1,size(S)
      @assertEqual(S(k), pS(k), sqrt(epsilon(one)))
      @assertEqual(Edn(k), pEdn(k), sqrt(epsilon(one)))
      @assertEqual(Eup(k), pEup(k), sqrt(epsilon(one)))
    enddo
    call PetscFinalize(ierr)


    call cpu_time(starttime)
    do k=1,iter
      call adding_delta_eddington_twostream(dtau, w0, g, mu0, incSolar, albedo, pS, pEdn, pEup)
    enddo
    call cpu_time(now)
    print *,'Time for adding_delta_eddington_twostream:', now-starttime

    do k=1,size(S)
      @assertEqual(S(k), pS(k), sqrt(epsilon(one)))
      @assertEqual(Edn(k), pEdn(k), sqrt(epsilon(one)))
      @assertEqual(Eup(k), pEup(k), sqrt(epsilon(one)))
    enddo
  end subroutine


end module
