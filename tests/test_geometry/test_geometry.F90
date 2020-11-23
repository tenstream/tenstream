module test_geometry
  use m_data_parameters, only :     &
    mpiint, iintegers, ireals, ireal_dp, &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_boxmc_geometry, only : &
    & setup_default_unit_cube_geometry, &
    & rand_pnt_on_plane

  use pfunit_mod
  implicit none

contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
  end subroutine teardown

  @test(npes =[1])
  subroutine test_pnt_on_plane(this)
    class (MpiTestMethod), intent(inout) :: this
    real(ireal_dp), parameter :: dx=100, dy=dx*2, dz=dx*3
    real(ireal_dp), allocatable :: vertices(:)
    integer(iintegers) :: k

    call setup_default_unit_cube_geometry(dx, dy, dz, vertices)

    do k=0,8-1
      print *,'Vertex', k, '::', vertices(1+3*k:3*k+3)
    enddo

    associate( &
        & A => vertices( 1: 3), &
        & B => vertices( 4: 6), &
        & C => vertices( 7: 9), &
        & D => vertices(10:12), &
        & E => vertices(13:15), &
        & F => vertices(16:18), &
        & G => vertices(19:21), &
        & H => vertices(22:24)  )

      call check_equal_distribution_on_plane(A,B,D,C, 1, dx, 2, dy) ! bot plane
      call check_equal_distribution_on_plane(E,G,H,F, 1, dx, 2, dy) ! top plane
      call check_equal_distribution_on_plane(A,E,F,B, 1, dx, 3, dz) ! back plane
      call check_equal_distribution_on_plane(G,C,D,H, 1, dx, 3, dz) ! front plane
      call check_equal_distribution_on_plane(A,C,G,E, 2, dy, 3, dz) ! left plane
      call check_equal_distribution_on_plane(B,F,H,D, 2, dy, 3, dz) ! right plane

    end associate

    contains
      subroutine check_equal_distribution_on_plane(A,B,C,D, pi, dx, pj, dy)
        real(ireal_dp), intent(in), dimension(3) :: A,B,C,D ! should have right hand side winding order
        integer(iintegers), intent(in) :: pi, pj ! dim idx of corresponing plane
        real(ireal_dp), intent(in) :: dx, dy ! edge length along first and second dim
        integer(iintegers), parameter :: Nsample=1000000
        integer(iintegers), parameter :: Nsub=10
        integer(iintegers) :: subAreas(Nsub, Nsub)
        integer(iintegers) :: k, i, j
        real(ireal_dp), dimension(3) :: pnt, normal, U, V

        subAreas(:,:) = -1 * Nsample / Nsub**2
        do k = 1, Nsample
          call rand_pnt_on_plane(A, B, C, D, pnt, normal, U, V)
          i = 1 + int(pnt(pi) / dx * Nsub)
          j = 1 + int(pnt(pj) / dy * Nsub)
          if (i.gt.Nsub .or. j.gt.Nsub) print *,'pnt', pnt, i,j
          subAreas(i,j) = subAreas(i,j)+1
        enddo
        do j = 1, Nsub
          print *,'Samples:', j, ':', subAreas(:,j)
        enddo
        print *,''
        @assertGreaterThan(3*sqrt(real(Nsample/Nsub, ireal_dp)), abs(subAreas), 'should have all area patches sampled equally')
      end subroutine
  end subroutine

end module
