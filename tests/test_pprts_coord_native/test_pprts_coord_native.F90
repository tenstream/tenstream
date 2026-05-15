module test_pprts_coord_native

  use mpi, only: MPI_PROC_NULL, MPI_INTEGER, MPI_SUM, MPI_Sendrecv, MPI_STATUS_SIZE
  use m_data_parameters, only: iintegers, mpiint, i0, i1
  use m_pprts_base, only: t_coord, setup_coord_native

  use pfunit_mod

  implicit none

contains

  @test(npes=[1])
  subroutine test_coord_1rank(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid
    type(t_coord), allocatable :: C
    integer(iintegers), parameter :: Nz = 5, Nx = 8, Ny = 6, dof = 4
    integer(iintegers) :: null_rank

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    call setup_coord_native(comm, Nz, Nx, Ny, dof, C)

    ! Global dimensions
    @assertEqual(Nz, C%glob_zm)
    @assertEqual(Nx, C%glob_xm)
    @assertEqual(Ny, C%glob_ym)

    ! z: not decomposed
    @assertEqual(i0, C%zs)
    @assertEqual(Nz - i1, C%ze)
    @assertEqual(Nz, C%zm)
    @assertEqual(i0, C%gzs)
    @assertEqual(Nz - i1, C%gze)
    @assertEqual(Nz, C%gzm)

    ! x,y: single rank owns full extent
    @assertEqual(i0, C%xs)
    @assertEqual(Nx - i1, C%xe)
    @assertEqual(Nx, C%xm)
    @assertEqual(i0, C%ys)
    @assertEqual(Ny - i1, C%ye)
    @assertEqual(Ny, C%ym)

    ! ghost bounds: one cell wider on each side
    @assertEqual(-i1, C%gxs)
    @assertEqual(Nx, C%gxe)
    @assertEqual(Nx + 2*i1, C%gxm)
    @assertEqual(-i1, C%gys)
    @assertEqual(Ny, C%gye)
    @assertEqual(Ny + 2*i1, C%gym)

    ! self-neighbor
    @assertEqual(int(myid, iintegers), C%neighbors(13))

    ! with periodic topology, 1-rank wraps to itself on all 4 sides
    ! west=10, east=16, south=4, north=22
    @assertEqual(int(myid, iintegers), C%neighbors(10))
    @assertEqual(int(myid, iintegers), C%neighbors(16))
    @assertEqual(int(myid, iintegers), C%neighbors(4))
    @assertEqual(int(myid, iintegers), C%neighbors(22))

    @assertEqual(dof, C%dof)
    @assertEqual(3_iintegers, C%dim)

    deallocate (C)
  end subroutine

  @test(npes=[4])
  subroutine test_coord_4rank(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: comm, myid, mpierr
    type(t_coord), allocatable :: C
    integer(iintegers), parameter :: Nz = 3, Nx = 10, Ny = 8, dof = 2
    integer(iintegers) :: null_rank
    integer(iintegers) :: my_cells, total_cells
    integer(mpiint) :: east_nbr, west_nbr, north_nbr, south_nbr
    integer(mpiint) :: recv_from_east, recv_from_west
    integer(mpiint) :: recv_from_north, recv_from_south
    integer(mpiint) :: status(MPI_STATUS_SIZE)

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    call setup_coord_native(comm, Nz, Nx, Ny, dof, C)

    null_rank = int(MPI_PROC_NULL, iintegers)

    ! Self-neighbor
    @assertEqual(int(myid, iintegers), C%neighbors(13))

    ! Local range consistency
    @assertEqual(C%xs + C%xm - i1, C%xe)
    @assertEqual(C%ys + C%ym - i1, C%ye)
    @assertTrue(C%xs >= i0)
    @assertTrue(C%xe < Nx)
    @assertTrue(C%ys >= i0)
    @assertTrue(C%ye < Ny)

    ! Ghost bounds
    @assertEqual(C%xs - i1, C%gxs)
    @assertEqual(C%xe + i1, C%gxe)
    @assertEqual(C%gxe - C%gxs + i1, C%gxm)
    @assertEqual(C%ys - i1, C%gys)
    @assertEqual(C%ye + i1, C%gye)
    @assertEqual(C%gye - C%gys + i1, C%gym)

    ! Total coverage: sum of (xm*ym) over all ranks must equal Nx*Ny
    my_cells = C%xm * C%ym
    call mpi_reduce(my_cells, total_cells, 1, MPI_INTEGER, MPI_SUM, 0_mpiint, comm, mpierr)
    if (myid == 0) then
      @assertEqual(Nx * Ny, total_cells)
    end if

    ! Neighbor consistency: if I say my east is rank E, then rank E says its west is me.
    east_nbr = int(C%neighbors(16), mpiint)
    west_nbr = int(C%neighbors(10), mpiint)
    north_nbr = int(C%neighbors(22), mpiint)
    south_nbr = int(C%neighbors(4), mpiint)

    recv_from_east = int(MPI_PROC_NULL, mpiint)
    recv_from_west = int(MPI_PROC_NULL, mpiint)
    recv_from_north = int(MPI_PROC_NULL, mpiint)
    recv_from_south = int(MPI_PROC_NULL, mpiint)

    ! Send myid eastward, receive from west: recv should be west_nbr's own rank
    call MPI_Sendrecv( &
      myid, 1, MPI_INTEGER, east_nbr, 0, &
      recv_from_west, 1, MPI_INTEGER, west_nbr, 0, &
      comm, status, mpierr)
    if (west_nbr /= int(MPI_PROC_NULL, mpiint)) then
      @assertEqual(west_nbr, recv_from_west)
    end if

    ! Send myid westward, receive from east: recv should be east_nbr's own rank
    call MPI_Sendrecv( &
      myid, 1, MPI_INTEGER, west_nbr, 1, &
      recv_from_east, 1, MPI_INTEGER, east_nbr, 1, &
      comm, status, mpierr)
    if (east_nbr /= int(MPI_PROC_NULL, mpiint)) then
      @assertEqual(east_nbr, recv_from_east)
    end if

    ! Send myid northward, receive from south
    call MPI_Sendrecv( &
      myid, 1, MPI_INTEGER, north_nbr, 2, &
      recv_from_south, 1, MPI_INTEGER, south_nbr, 2, &
      comm, status, mpierr)
    if (south_nbr /= int(MPI_PROC_NULL, mpiint)) then
      @assertEqual(south_nbr, recv_from_south)
    end if

    ! Send myid southward, receive from north
    call MPI_Sendrecv( &
      myid, 1, MPI_INTEGER, south_nbr, 3, &
      recv_from_north, 1, MPI_INTEGER, north_nbr, 3, &
      comm, status, mpierr)
    if (north_nbr /= int(MPI_PROC_NULL, mpiint)) then
      @assertEqual(north_nbr, recv_from_north)
    end if

    deallocate (C)
  end subroutine

end module
