module test_intersection

  use m_data_parameters, only: &
    & ireals, ireal_dp, &
    & iintegers, mpiint, zero, one

  use m_helper_functions, only : &
    & compute_normal_3d, &
    & distance_to_edge, &
    & determine_normal_direction

  use m_intersection, only : &
    & hit_plane, &
    & pnt_in_triangle, &
    & triangle_intersection, &
    & line_intersection_3d

  use pfunit_mod

  implicit none

contains

  !@test(npes =[1])
  subroutine test_triangle_functions_dp(this)
    class (MpiTestMethod), intent(inout) :: this

    real(ireal_dp),parameter :: zero=0, one=1, dx = 100
    real(ireal_dp),parameter :: A(2) = [zero, zero]
    real(ireal_dp),parameter :: B(2) = [dx, zero]
    real(ireal_dp),parameter :: C(2) = [dx/2.,sqrt(dx**2 - (dx/2)**2)]
    real(ireal_dp) :: P(2), distance

    real(ireal_dp) :: normal(3)

    ! Tests determining the distance of a point to a 2D line/edge
    @assertEqual(zero, distance_to_edge(A,B,A), 'from point on line, the distance to same line should be zero distance_to_edge1')
    @assertEqual(zero, distance_to_edge(A,B,[dx/2,zero]), 'from point on line, the distance to same line should be zero distance_to_edge2')
    @assertEqual(one, distance_to_edge(A,B,[dx/2,one]), 'here point line <-> distance should be one distance_to_edge3')
    @assertEqual(sqrt(epsilon(dx)), distance_to_edge(A,B,[dx/2,sqrt(epsilon(dx))]), 'here point line <-> distance should be different distance_to_edge4')
    @assertEqual(epsilon(dx), distance_to_edge(A,B,[dx/2,epsilon(dx)]), 'here point line <-> distance should be different distance_to_edge5')
    @assertEqual(epsilon(dx), distance_to_edge(A,B,[dx/2,-epsilon(dx)]), 'here point line <-> distance should be different distance_to_edge6')
    @assertEqual(one, distance_to_edge(A,B,B+[zero,one]), 'here point line <-> distance should be one test distance_to_edge7')
    @assertEqual(zero, distance_to_edge(A,B,[-one,zero]), 'from point on line, the distance to same line should be one distance_to_edge8')
    @assertEqual(zero, distance_to_edge(A,B,[dx+one,zero]), 'from point on line, the distance to same line should be one distance_to_edge9')

    ! Compute normals
    normal = compute_normal_3d([A(1),A(2),zero], [B(1),B(2),zero], [C(1),C(2),zero])
    @assertEqual([zero,zero,one], normal, '3D normal not as expected')

    normal = compute_normal_3d([A(1),A(2),zero], [C(1),C(2),zero], [B(1),B(2),zero])
    @assertEqual([zero,zero,one], -normal, '3D normal not as expected')

    @assertEqual(one, norm2(normal), 'returned normal is not normed to one')

    ! Check if we can determine if a point is in a triangle
    @assertTrue(pnt_in_triangle(A,B,C, A), 'pnt_in_triangle wrong for edge case in A')
    @assertTrue(pnt_in_triangle(A,B,C, B), 'pnt_in_triangle wrong for edge case in B')
    @assertTrue(pnt_in_triangle(A,B,C, C), 'pnt_in_triangle wrong for edge case in C')
    @assertTrue(pnt_in_triangle(A,B,C, [.5_ireal_dp, .5_ireal_dp]), 'pnt_in_triangle wrong for center of triangle')

    @assertTrue(pnt_in_triangle(A,B,C, A+(C-A)/2), 'pnt_in_triangle wrong for edge case on line between A and C')
    @assertTrue(pnt_in_triangle(A,B,C, A+(B-A)/2), 'pnt_in_triangle wrong for edge case on line between A and B')
    @assertTrue(pnt_in_triangle(A,B,C, C+(B-C)/2), 'pnt_in_triangle wrong for edge case on line between B and C')

    @assertFalse(pnt_in_triangle(A,B,C, A-[one,one ]), 'pnt_in_triangle wrong for outside case 1')
    @assertFalse(pnt_in_triangle(A,B,C, B+[one,zero]), 'pnt_in_triangle wrong for outside case 22')
    @assertFalse(pnt_in_triangle(A,B,C, C+[one,one] ), 'pnt_in_triangle wrong for outside case 3')

    @assertTrue(pnt_in_triangle(A,B,C, [0.38475394248962402_ireal_dp, zero]), 'custom edge case point should be in triangle! case 1')

    @assertFalse(pnt_in_triangle(A,B,C, [dx, dx]), 'point (dx,dy) should not be in triangle!')


    ! vector from C to pnt halfway between (AB):
    P = A+(B-A)/2 - C
    @assertTrue(pnt_in_triangle(A,B,C, C+P), 'pnt_in_triangle wrong for edge case on line between A and B')
    @assertFalse(pnt_in_triangle(A,B,C, C+(one+sqrt(epsilon(one)))*P), 'pnt_in_triangle wrong for edge case epsilon after line between A and B')

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

  !@test(npes =[1])
  subroutine test_triangle_functions_sp(this)
    class (MpiTestMethod), intent(inout) :: this

    real(ireals),parameter :: zero=0, one=1, dx = 100
    real(ireals),parameter :: A(2) = [zero, zero]
    real(ireals),parameter :: B(2) = [dx, zero]
    real(ireals),parameter :: C(2) = [dx/2.,sqrt(dx**2 - (dx/2)**2)]
    real(ireals) :: P(2), distance

    real(ireals), dimension(3) :: normal, new_loc, center_face, center_cell
    integer(iintegers) :: normal_direction

    ! Tests determining the distance of a point to a 2D line/edge
    @assertEqual(zero, distance_to_edge(A,B,A), 'from point on line, the distance to same line should be zero distance_to_edge1')
    @assertEqual(zero, distance_to_edge(A,B,[dx/2,zero]), 'from point on line, the distance to same line should be zero distance_to_edge2')
    @assertEqual(one, distance_to_edge(A,B,[dx/2,one]), 'here point line <-> distance should be one distance_to_edge3')
    @assertEqual(sqrt(epsilon(dx)), distance_to_edge(A,B,[dx/2,sqrt(epsilon(dx))]), 'here point line <-> distance should be different distance_to_edge4')
    @assertEqual(epsilon(dx), distance_to_edge(A,B,[dx/2,epsilon(dx)]), 'here point line <-> distance should be different distance_to_edge5')
    @assertEqual(epsilon(dx), distance_to_edge(A,B,[dx/2,-epsilon(dx)]), 'here point line <-> distance should be different distance_to_edge6')
    @assertEqual(one, distance_to_edge(A,B,B+[zero,one]), 'here point line <-> distance should be one test distance_to_edge7')


    ! Checks if points lie in a triangle
    new_loc = [0.38475394248962402_ireals, zero, zero]
    @assertTrue(pnt_in_triangle(A,B,C, [new_loc(1), new_loc(2)]), 'custom edge case point should be in triangle!')

    normal = compute_normal_3d([A(1),A(2),zero], [B(1),B(2),zero], [C(1),C(2),zero])
    @assertEqual([zero,zero,one], normal, 10*epsilon(normal), '3D normal not as expected')

    normal = compute_normal_3d([A(1),A(2),zero], [C(1),C(2),zero], [B(1),B(2),zero])
    @assertEqual([zero,zero,one], -normal, 10*epsilon(normal), '3D normal not as expected')

    @assertEqual(one, norm2(normal), 10*epsilon(normal), 'returned normal is not normed to one')

    ! Check if we can determine if a point is in a triangle
    @assertTrue(pnt_in_triangle(A,B,C, A), 'pnt_in_triangle wrong for edge case in A')
    @assertTrue(pnt_in_triangle(A,B,C, B), 'pnt_in_triangle wrong for edge case in B')
    @assertTrue(pnt_in_triangle(A,B,C, C), 'pnt_in_triangle wrong for edge case in C')
    @assertTrue(pnt_in_triangle(A,B,C, [one/2, one/2]), 'pnt_in_triangle wrong for center of triangle')

    @assertTrue(pnt_in_triangle(A,B,C, A+(C-A)/2), 'pnt_in_triangle wrong for edge case on line between A and C')
    @assertTrue(pnt_in_triangle(A,B,C, A+(B-A)/2), 'pnt_in_triangle wrong for edge case on line between A and B')
    @assertTrue(pnt_in_triangle(A,B,C, C+(B-C)/2), 'pnt_in_triangle wrong for edge case on line between B and C')

    @assertFalse(pnt_in_triangle(A,B,C, A-[one,one ]), 'pnt_in_triangle wrong for outside case 1')
    @assertFalse(pnt_in_triangle(A,B,C, B+[one,zero]), 'pnt_in_triangle wrong for outside case 2')
    @assertFalse(pnt_in_triangle(A,B,C, C+[one,one] ), 'pnt_in_triangle wrong for outside case 3')


    ! vector from C to pnt halfway between (AB):
    P = A+(B-A)/2 - C
    @assertTrue(pnt_in_triangle(A,B,C, C+P), 'pnt_in_triangle wrong for edge case on line between A and B')
    @assertFalse(pnt_in_triangle(A,B,C, C+(one+sqrt(epsilon(one)))*P), 'pnt_in_triangle wrong for edge case epsilon after line between A and B')

    ! Check if distance caluclations are OK
    distance = hit_plane([A(1),A(2),one], [zero,zero,-one], [A(1),A(2),zero], normal)
    @assertEqual(one, distance, 'distance calculation not correct 1')

    distance = hit_plane([A(1),A(2),one], [zero,zero,+one], [A(1),A(2),zero], normal)
    @assertEqual(-one, distance, 'distance calculation not correct 2')


    distance = hit_plane([A(1),A(2),one], [zero,zero,-one], [C(1),C(2),zero], normal)
    @assertEqual(one, distance, 'distance calculation not correct 3')

    distance = hit_plane([A(1),A(2),one], [zero,zero,+one], [C(1),C(2),zero], normal)
    @assertEqual(-one, distance, 'distance calculation not correct 4')

    ! test routines that determine direction of normal
    normal = [zero, zero, one]
    center_face = [zero, zero, zero]
    center_cell = [zero, zero, one]
    normal_direction = determine_normal_direction(normal, center_face, center_cell)
    @assertEqual(1_iintegers, normal_direction, 'direction of normal not towards cell center (case 1)')

    center_cell = [zero, zero, -one]
    normal_direction = determine_normal_direction(normal, center_face, center_cell)
    @assertEqual(-1_iintegers, normal_direction, 'direction of normal not towards cell center (case 2)')
  end subroutine

  !@test(npes =[1])
  subroutine test_triangle_intersection(this)
    class (MpiTestMethod), intent(inout) :: this

    real(ireal_dp),parameter :: zero=0, one=1, dx = 100, atol=100*epsilon(atol)
    real(ireal_dp),parameter :: rng(2) = [0._ireal_dp, huge(rng)]
    real(ireal_dp) :: A(3)
    real(ireal_dp) :: B(3)
    real(ireal_dp) :: C(3)
    real(ireal_dp) :: direction(3)
    real(ireal_dp) :: origin(3)
    real(ireal_dp) :: hit(4)
    logical :: lhit

    A = [zero, zero, zero]
    B = [  dx, zero, zero]
    C = [dx/2, zero, sqrt(dx**2 - (dx/2)**2)]

    direction = [zero, one, zero]

    origin    = [dx/2,-one, dx/2]
    call triangle_intersection(origin, direction, A, B, C, rng, lhit, hit)
    @assertTrue(lhit)
    @assertEqual(origin(1), hit(1), atol)
    @assertEqual(zero,      hit(2), atol)
    @assertEqual(origin(3), hit(3), atol)
    print *,origin,'=>',hit

    origin    = [dx/2,one, dx/2]
    call triangle_intersection(origin, direction, A, B, C, rng, lhit, hit)
    @assertFalse(lhit)

    origin    = [dx/2,-epsilon(dx), dx/2]
    call triangle_intersection(origin, direction, A, B, C, rng, lhit, hit)
    @assertTrue(lhit)
    @assertEqual(origin(1), hit(1), atol)
    @assertEqual(zero,      hit(2), atol)
    @assertEqual(origin(3), hit(3), atol)

    origin    = [dx/2,epsilon(dx), dx/2]
    call triangle_intersection(origin, direction, A, B, C, rng, lhit, hit)
    @assertFalse(lhit)

    origin    = A + [zero,-epsilon(dx), zero]
    call triangle_intersection(origin, direction, A, B, C, rng, lhit, hit)
    @assertTrue(lhit)
    @assertEqual(origin(1), hit(1), atol)
    @assertEqual(zero,      hit(2), atol)
    @assertEqual(origin(3), hit(3), atol)

    origin    = C + [zero,-epsilon(dx), zero]
    call triangle_intersection(origin, direction, A, B, C, rng, lhit, hit)
    @assertTrue(lhit)
    @assertEqual(origin(1), hit(1), atol)
    @assertEqual(zero,      hit(2), atol)
    @assertEqual(origin(3), hit(3), atol)

    origin    = C + [zero,-epsilon(dx), zero]
    call triangle_intersection(origin, direction, B, A, C, rng, lhit, hit)
    @assertTrue(lhit)
    @assertEqual(origin(1), hit(1), atol)
    @assertEqual(zero,      hit(2), atol)
    @assertEqual(origin(3), hit(3), atol)

    origin    = [dx/2,epsilon(dx), dx/2]
    call triangle_intersection(origin, direction, B, A, C, rng, lhit, hit)
    @assertFalse(lhit)
  end subroutine

  @test(npes =[1])
  subroutine test_line_intersection_3d(this)
    class (MpiTestMethod), intent(inout) :: this

    real(ireal_dp) :: origin1(3), origin2(3), direction1(3), direction2(3), coeff1, coeff2, atol
    integer(mpiint) :: ierr

    atol = 1e-3_ireal_dp
    origin1 = [zero, zero, zero]
    origin2 = [one,  one,  one ]

    direction1 = [zero, one, zero]
    direction2 = [zero, one, zero]

    print *, 'test running'
    call line_intersection_3d(origin1, direction1, origin2, direction2, coeff1, coeff2, ierr)
    @assertEqual(coeff1, -huge(coeff1), atol)
    @assertTrue(.false.)
    print *, 'ierr', ierr

  end subroutine
end module
