@test(npes =[1,2,4])
subroutine test_mpi_functions_dp(this)

    use m_data_parameters, only: ireal_dp, iintegers, mpiint, init_mpi_data_parameters
    use m_helper_functions_dp, only : imp_bcast

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid, i

    real(ireal_dp),allocatable :: bcast_1d_arr(:)
    real(ireal_dp),allocatable :: bcast_2d_arr(:,:)
    real(ireal_dp),allocatable :: bcast_3d_arr(:,:,:)
    real(ireal_dp),allocatable :: bcast_5d_arr(:,:,:,:,:)
    real(ireal_dp) :: bcast_real_scalar

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    ! Check if scalar and array bcasts work
    if(myid.eq.0) then
      bcast_real_scalar = 1234._ireal_dp
      allocate(bcast_1d_arr(2), source=1234._ireal_dp)
      allocate(bcast_2d_arr(2,1), source=1234._ireal_dp)
      allocate(bcast_3d_arr(2,1,1), source=1234._ireal_dp)
      allocate(bcast_5d_arr(2,1,1,1,1), source=1234._ireal_dp)
    else
      bcast_real_scalar = -1
    endif

    call imp_bcast(comm, bcast_1d_arr, 0_mpiint)
    @assertTrue(allocated(bcast_1d_arr), 'real array broadcast wrong, allocation failed')
    @assertEqual(1234, bcast_1d_arr(1), 'real array broadcast wrong')
    @assertEqual(1234, bcast_1d_arr(2), 'real array broadcast wrong')

    call imp_bcast(comm, bcast_2d_arr, 0_mpiint)
    @assertTrue(allocated(bcast_2d_arr), 'real array broadcast wrong, allocation failed')
    @assertEqual(1234, bcast_2d_arr(1,1), 'real array broadcast wrong')
    @assertEqual(1234, bcast_2d_arr(2,1), 'real array broadcast wrong')

    call imp_bcast(comm, bcast_3d_arr, 0_mpiint)
    @assertTrue(allocated(bcast_3d_arr), 'real array broadcast wrong, allocation failed')
    @assertEqual(1234, bcast_3d_arr(1,1,1), 'real array broadcast wrong')
    @assertEqual(1234, bcast_3d_arr(2,1,1), 'real array broadcast wrong')

    call imp_bcast(comm, bcast_5d_arr, 0_mpiint)
    @assertTrue(allocated(bcast_5d_arr), 'real array broadcast wrong, allocation failed')
    @assertEqual(1234, bcast_5d_arr(1,1,1,1,1), 'real array broadcast wrong')
    @assertEqual(1234, bcast_5d_arr(2,1,1,1,1), 'real array broadcast wrong')


    ! Scalar Bcasts:
    call imp_bcast(comm, bcast_real_scalar, 0_mpiint)
    @assertEqual(1234, bcast_real_scalar, 'real scalar broadcast wrong')

end subroutine

@test(npes =[1])
subroutine test_triangle_functions_dp()

    use m_data_parameters, only: ireal_dp, iintegers, mpiint, init_mpi_data_parameters, zero, one
    use m_helper_functions_dp, only : compute_normal_3d, hit_plane, pnt_in_triangle, norm

    use pfunit_mod

    implicit none

    real(ireal_dp),parameter :: A(2) = [0,0]
    real(ireal_dp),parameter :: B(2) = [2,0]
    real(ireal_dp),parameter :: C(2) = [1,1]
    real(ireal_dp) :: P(2), distance

    real(ireal_dp) :: normal(3)
    real(ireal_dp) :: normal_target(3) = [0,0,1]

    normal = compute_normal_3d([A(1),A(2),zero], [B(1),B(2),zero], [C(1),C(2),zero])
    @assertEqual([zero,zero,one], normal, '3D normal not as expected')

    normal = compute_normal_3d([A(1),A(2),zero], [C(1),C(2),zero], [B(1),B(2),zero])
    @assertEqual([zero,zero,one], -normal, '3D normal not as expected')

    @assertEqual(one, norm(normal), 'returned normal is not normed to one')

    ! Check if we can determine if a point is in a triangle
    @assertTrue(pnt_in_triangle(A,B,C, A), 'pnt_in_triangle wrong for edge case in A')
    @assertTrue(pnt_in_triangle(A,B,C, B), 'pnt_in_triangle wrong for edge case in B')
    @assertTrue(pnt_in_triangle(A,B,C, C), 'pnt_in_triangle wrong for edge case in C')
    @assertTrue(pnt_in_triangle(A,B,C, [.5_ireal_dp, .5_ireal_dp]), 'pnt_in_triangle wrong for center of triangle')

    @assertTrue(pnt_in_triangle(A,B,C, A+(C-A)/2), 'pnt_in_triangle wrong for edge case on line between A and C')
    @assertTrue(pnt_in_triangle(A,B,C, A+(B-A)/2), 'pnt_in_triangle wrong for edge case on line between A and B')
    @assertTrue(pnt_in_triangle(A,B,C, C+(B-C)/2), 'pnt_in_triangle wrong for edge case on line between B and C')

    @assertFalse(pnt_in_triangle(A,B,C, A-[one,one ]), 'pnt_in_triangle wrong for outside case 1')
    @assertFalse(pnt_in_triangle(A,B,C, B+[one,zero]), 'pnt_in_triangle wrong for outside case 2')
    @assertFalse(pnt_in_triangle(A,B,C, C+[one,one] ), 'pnt_in_triangle wrong for outside case 3')


    ! vector from C to pnt halfway between (AB):
    P = A+(B-A)/2 - C
    @assertTrue(pnt_in_triangle(A,B,C, C+P), 'pnt_in_triangle wrong for edge case on line between A and B')
    @assertFalse(pnt_in_triangle(A,B,C, C+(one+epsilon(one))*P), 'pnt_in_triangle wrong for edge case epsilon after line between A and B')


    ! Check if distance caluclations are OK
    distance = hit_plane([A(1),A(2),one], [zero,zero,-one], [A(1),A(2),zero], normal)
    @assertEqual(one, distance, 'distance calculation not correct 1')

    distance = hit_plane([A(1),A(2),one], [zero,zero,+one], [A(1),A(2),zero], normal)
    @assertEqual(-one, distance, 'distance calculation not correct 2')


    distance = hit_plane([A(1),A(2),one], [zero,zero,-one], [C(1),C(2),zero], normal)
    @assertEqual(one, distance, 'distance calculation not correct 3')

    distance = hit_plane([A(1),A(2),one], [zero,zero,+one], [C(1),C(2),zero], normal)
    @assertEqual(-one, distance, 'distance calculation not correct 4')
end subroutine
