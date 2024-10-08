module m_LUT_param_phi
  use iso_fortran_env, only: real32, real64
  use m_data_parameters, only: irealLUT, mpiint, ireal_dp, ireal_params, ireals
  use m_data_parameters, only: pi => pi_ireal_params
  use m_helper_functions, only: angle_between_two_vec, rad2deg, approx, is_between, CHKERR, ftoa, solve_quadratic, &
                                normalize_vec
  use m_boxmc_geometry, only: wedge_halfspaces
  implicit none

  interface LUT_wedge_aspect_zx
    module procedure LUT_wedge_aspect_zx_r32, LUT_wedge_aspect_zx_r64
  end interface
  interface LUT_wedge_dz
    module procedure LUT_wedge_dz_r32, LUT_wedge_dz_r64
  end interface

  logical, parameter :: ldebug = .false.

contains
!> @brief reverse transformation from azimuth(radians) to param_phi[-2,2]
  function param_phi_from_azimuth(phi, wedge_C) result(param_phi)
    real(ireal_params), intent(in) :: phi, wedge_C(:)
    real(ireal_params) :: param_phi

    real(ireal_params) :: alpha, beta
    real(ireal_params) :: lb, rb, x1, x2! bounds of local spline

    associate (pA => [0._ireal_params, 0._ireal_params], pB => [1._ireal_params, 0._ireal_params], pC => wedge_C)
      alpha = angle_between_two_vec(pB - pA, pC - pA)
      beta = angle_between_two_vec(pA - pB, pC - pB)

      if (phi .gt. pi / 2 - alpha) then ! range [ . , -1]
        x1 = -2; x2 = -1
        lb = pi / 2 - alpha / 2
        rb = pi / 2 - alpha
      elseif (phi .lt. beta - pi / 2) then ! between [ 1, .]
        x1 = 1; x2 = 2
        lb = beta - pi / 2
        rb = beta / 2 - pi / 2
      else ! range [-1, 1]
        x1 = -1; x2 = 1
        lb = pi / 2 - alpha
        rb = beta - pi / 2
      end if
      param_phi = (x2 - x1) / (rb - lb) * (phi - lb) + x1
    end associate
  end function

!> @brief translate from the parameterized azimuth [-2,2] to azimuth [radians]
!> @details we use local splines so that it ensures we have mapped:
!>   alpha and beta are the inner angles of a wedge triangle between AB and AC or BA and BC respectively
!>   90 - alpha/2 -> -2
!>   90 - alpha   -> -1
!>   beta   - 90  ->  1
!>   beta/2 - 90  ->  2
  function azimuth_from_param_phi(param_phi, wedge_C) result(phi)
    real(ireal_params), intent(in) :: param_phi, wedge_C(:)
    real(ireal_params) :: phi

    real(ireal_params) :: alpha, beta
    real(ireal_params) :: lb, rb, x1, x2! bounds of local spline

    alpha = angle_between_two_vec([1._ireal_params, 0._ireal_params], wedge_C)
    beta = angle_between_two_vec([-1._ireal_params, 0._ireal_params], wedge_C - [1._ireal_params, 0._ireal_params])

    if (param_phi .lt. -1._ireal_params) then ! range [ . , -1]
      x1 = -2; x2 = -1
      lb = pi / 2 - alpha / 2
      rb = pi / 2 - alpha
    elseif (param_phi .gt. 1._ireal_params) then ! between [ 1, .]
      x1 = 1; x2 = 2
      lb = beta - pi / 2
      rb = beta / 2 - pi / 2
    else ! range [-1, 1]
      x1 = -1; x2 = 1
      lb = pi / 2 - alpha
      rb = beta - pi / 2
    end if
    phi = (rb - lb) / (x2 - x1) * (param_phi - x1) + lb
  end function

  subroutine phi_crit(side_normal, theta, phic, ierr)
    ! solve euqation system that arises for condition
    ! sunvec_crit = [sin(phi_crit) * sin(theta), cos(phi_crit) * sin(theta), -cos(theta)]
    ! left_side_face_normal .dot. sunvec_crit === 0
    ! i.e. solve( 0 = n1 * sin ( phi ) * sin(theta) + n2 * cos(phi) * sin(theta) - n3 * cos(theta), phi)

    real(ireal_params), intent(in) :: side_normal(3)
    real(ireal_params), intent(in) :: theta ! sun zenith in [rad]
    real(ireal_params), intent(out) :: phic
    integer(mpiint), intent(out) :: ierr

    real(ireal_params) :: a, b, c, z, sol(2)
    real(ireal_params) :: x(2), y(2)
    !real(ireal_params) :: st, ct

    if (ldebug) then
      if (.not. approx(norm2(side_normal), 1._ireal_params, 10 * epsilon(side_normal))) then
        call CHKERR(1_mpiint, 'side_normal needs to be normed: '// &
                    ftoa(side_normal)//' ( '//ftoa(norm2(side_normal))//' )')
      end if
    end if
    z = -cos(theta)
    a = side_normal(1)**2 + side_normal(2)**2
    if (abs(side_normal(1)) .gt. abs(side_normal(2))) then
      b = 2 * side_normal(2) * side_normal(3) * z
      c = (side_normal(1)**2 + side_normal(3)**2) * z**2 - side_normal(1)**2

      call solve_quadratic(a, b, c, y, ierr)
      if (ierr .ne. 0) return

      x = (-side_normal(2) * y - side_normal(3) * z) / side_normal(1)
    else
      b = 2 * side_normal(1) * side_normal(3) * z
      c = (side_normal(2)**2 + side_normal(3)**2) * z**2 - side_normal(2)**2

      call solve_quadratic(a, b, c, x, ierr)
      if (ierr .ne. 0) return

      y = (-side_normal(1) * x - side_normal(3) * z) / side_normal(2)
    end if

    sol = atan2(x, y)
    if (abs(sol(1)) .ge. (abs(sol(2)))) then
      phic = sol(2)
    else
      phic = sol(1)
    end if

    !associate( n1=>side_normal(1), n2=>side_normal(2), n3=>side_normal(3) )
    !  st = sin(theta)
    !  ct = cos(theta)
    !  if(approx(n2,0._ireal_params) .and. approx(n3,0._ireal_params) .and. .not. approx(n1*st,0._ireal_params)) then
    !    phic = 0._ireal_params
    !    return
    !  endif
    !  if(approx(n2,0._ireal_params) .and. &
    !    approx(ct, 0._ireal_params) .and. &
    !    .not. approx(n3, 0._ireal_params) .and. &
    !    .not. approx(n1*st, 0._ireal_params) ) then
    !    phic = 0._ireal_params
    !    return
    !  endif
    !  if(.not.approx(n2,0._ireal_params) .and. approx(st, -n3*ct/n2)) then
    !    phic = pi_ireal_params
    !    return
    !  endif
    !  z = n2*st+n3*ct
    !  if(.not.approx(z,0._ireal_params)) then
    !    b = n1**2*st**2 + n2**2*st**2 - n3**2*ct**2
    !    if(b.lt.0._ireal_params) ierr = 1
    !    a = n1*st / z
    !    c = n2*st + n3*ct

    !    if(approx(n1**2*st**2 - n1*st*sqrt(b) + n2**2*st**2 + n2*n3*st*ct, 0._ireal_params)) ierr = 2
    !    sol(1) = 2*atan(a - sqrt(b)/c)
    !    if(approx(n1**2*st**2 + n1*st*sqrt(b) + n2**2*st**2 + n2*n3*st*ct, 0._ireal_params)) ierr = 3
    !    sol(2) = 2*atan(sqrt(b)/c + a)
    !    !print *,'z',z,'a,b,c', a,b,c,'st', st,ct, 'sol', sol, ':', a - sqrt(b)/c, sqrt(b)/c + a
    !    if(abs(sol(1)).ge.(abs(sol(2)))) then
    !      phic = sol(2)
    !    else
    !      phic = sol(1)
    !    endif
    !    return
    !  endif
    !end associate
    !ierr = 1
  end subroutine

  function theta_crit(side_normal, phi) result(thetac)
    ! solve equation system that arises for condition
    ! sunvec_crit = [sin(phi) * sin(theta_crit), cos(phi) * sin(theta_crit), -cos(theta_crit)]
    ! base_face_normal .dot. sunvec_crit === 0
    ! i.e. solve( 0 = n1 * sin ( phi ) * sin(theta) + n2 * cos(phi) * sin(theta) - n3 * cos(theta), theta)

    real(ireal_params), intent(in) :: side_normal(3)
    real(ireal_params), intent(in) :: phi ! sun azimuth in [rad]
    real(ireal_params) :: thetac

    !real(ireal_params) :: discr

    associate (n1 => side_normal(1), n2 => side_normal(2), n3 => side_normal(3))
      !discr  = ( n1*sin(phi) + n2*cos(phi) )**2 + n3**2
      !thetac = asin(n3 / sqrt(discr))
      thetac = atan(n3 / (n1 * sin(phi) + n2 * cos(phi)))
      !thetac = asin(n3 / sqrt( ( n1*sin(phi) + n2*cos(phi) )**2 + n3**2 ) )
    end associate

    !if(ldebug .and. discr.le.0._ireal_params) then
    !  print *,'discr', discr, side_normal, norm2(side_normal), '->', side_normal(3) / discr, 'thetac', thetac
    !  call CHKERR(1_mpiint, 'no solution for theta crit')
    !endif
  end function

  subroutine iterative_phi_theta_from_param_phi_and_param_theta(wedge_coords3d, &
                                                                param_phi, param_theta, &
                                                                phi, theta, ierr)
    real(ireal_params), intent(in), dimension(:) :: wedge_coords3d ! dim 18
    real(ireal_params), intent(in) :: param_phi, param_theta
    real(ireal_params), intent(out) :: phi, theta
    integer(mpiint), intent(out) :: ierr
    real(ireal_params) :: phic3, phic4, thetac, phie3, phie4
    real(ireal_params) :: last_phi, last_theta, new_phi
    real(ireal_params), dimension(3) :: n2, n3, n4
    integer(mpiint) :: iter
    integer(mpiint), parameter :: Niter = 10000
    real(ireal_params), parameter :: eps = 1e-8

    call determine_end_phi_points_plus_geometry(wedge_coords3d, n2, n3, n4, phie3, phie4)

    phi = 0._ireal_params
    theta = 0._ireal_params
    last_phi = phi
    last_theta = theta

    do iter = 1, Niter
      thetac = theta_crit(n2, real(phi, ireal_params))
      theta = theta_from_param_theta(param_theta, thetac)

      call determine_critical_phi_points( &
        n3, n4, theta, phie3, phie4, &
        phic3, phic4)

      call phi_from_param_phi(param_phi, phic3, phie3, phic4, phie4, new_phi, ierr); call CHKERR(ierr)
      phi = phi + (new_phi - phi) * (1._ireal_params - real(iter - 1, ireal_params) / real(Niter, ireal_params))**2
      !print *,iter, 'param_phi', param_phi, 'param_theta', param_theta, ':', new_phi, &
      !  '=>', phi, theta, 'thetac', thetac, &
      !  'phie/c',rad2deg([phie3, phie4, phic3, phic4])
      if (approx(last_phi, phi, eps) .and. approx(last_theta, theta, eps)) then
        ierr = 0
        return
      end if
      last_phi = phi
      last_theta = theta
    end do
    ierr = iter
  end subroutine

  subroutine determine_end_phi_points_plus_geometry(wedge_coords3d, &
                                                    n2, n3, n4, phie3, phie4)
    real(ireal_params), intent(in), dimension(:) :: wedge_coords3d ! dim 18
    real(ireal_params), intent(out) :: n2(:), n3(:), n4(:), phie3, phie4

    real(ireal_dp) :: origins(3, 5), normals(3, 5)
    real(ireal_params) :: alpha, beta
    integer(mpiint) :: ierr

    call wedge_halfspaces(real(wedge_coords3d, ireal_dp), origins, normals)

    call normalize_vec(real(normals(:, 2), ireal_params), n2, ierr); call CHKERR(ierr)
    call normalize_vec(real(normals(:, 3), ireal_params), n3, ierr); call CHKERR(ierr)
    call normalize_vec(real(normals(:, 4), ireal_params), n4, ierr); call CHKERR(ierr)

    associate (A => wedge_coords3d(1:3), B => wedge_coords3d(4:6), C => wedge_coords3d(7:9))
      alpha = angle_between_two_vec(B - A, C - A)
      beta = angle_between_two_vec(A - B, C - B)
    end associate
    phie3 = pi / 2 - alpha / 2
    phie4 = beta / 2 - pi / 2
  end subroutine

  subroutine determine_critical_phi_points( &
    n3, n4, theta, phie3, phie4, &
    phic3, phic4)
    real(ireal_params), intent(in) :: theta
    real(ireal_params), intent(in) :: n3(:), n4(:), phie3, phie4
    real(ireal_params), intent(out) :: phic3, phic4
    integer(mpiint) :: ierr

    call phi_crit(n3, theta, phic3, ierr)
    if (ierr .ne. 0) phic3 = phie3

    call phi_crit(n4, theta, phic4, ierr)
    if (ierr .ne. 0) phic4 = phie4

    if (phic3 .lt. phie4) phic3 = phie3 ! use phie3 if there is no solution
    if (phic4 .gt. phie3) phic4 = phie4 ! use phie4 if there is no solution

    phic3 = min(phie3, phic3)
    phic4 = max(phie4, phic4)
  end subroutine

  function theta_from_param_theta(param_theta, thetac) result(theta)
    real(ireal_params), intent(in) :: param_theta
    real(ireal_params), intent(in) :: thetac
    real(ireal_params) :: theta
    real(ireal_params) :: lb, rb, x1, x2! bounds of local spline

    if (param_theta .le. 0._ireal_params) then ! range [ . , 0]
      x1 = -1; x2 = -0
      lb = 0._ireal_params
      rb = thetac
    else ! range [0, .]
      x1 = 0; x2 = 1
      lb = thetac
      rb = pi / 2
    end if
    theta = (rb - lb) / (x2 - x1) * (param_theta - x1) + lb
  end function

  subroutine phi_from_param_phi(param_phi, phic3, phie3, phic4, phie4, phi, ierr)
    real(ireal_params), intent(in) :: param_phi
    real(ireal_params), intent(in) :: phic3, phic4, phie3, phie4
    real(ireal_params), intent(out) :: phi
    integer(mpiint), intent(out) :: ierr
    real(ireal_params) :: lb, rb, x1, x2! bounds of local spline

    if (param_phi .lt. -1._ireal_params) then ! range [ . , -1]
      x1 = -2; x2 = -1
      lb = phie3
      rb = phic3
    elseif (param_phi .gt. 1._ireal_params) then ! between [ 1, .]
      x1 = 1; x2 = 2
      lb = phic4
      rb = phie4
    else ! range [-1, 1]
      x1 = -1; x2 = 1
      lb = phic3
      rb = phic4
    end if
    phi = (rb - lb) / (x2 - x1) * (param_phi - x1) + lb

    ierr = 0
  end subroutine

  subroutine param_phi_param_theta_from_phi_and_theta_withnormals(n2, n3, n4, Cx, Cy, phi, theta, &
                                                                  param_phi, param_theta, ierr)
    real(ireal_params), dimension(:), intent(in) :: n2, n3, n4 !normals of base, left, right face, dim 3
    real(ireal_params), intent(in) :: Cx, Cy     ! Cx, Cy are local coords of C-point
    real(ireal_params), intent(in) :: phi, theta ! phi, theta are local azimuth and zenith
    real(ireal_params), intent(out) :: param_phi, param_theta
    integer(mpiint), intent(out) :: ierr
    real(ireal_params) :: alpha, beta
    real(ireal_params) :: phic3, phic4, thetac, phie3, phie4

    if (approx(Cx, 0._ireal_params)) then
      alpha = 0
    else
      alpha = atan(Cy / Cx)
    end if
    if (approx(Cx, 1._ireal_params)) then
      beta = 0
    else
      beta = atan(Cy / (1._ireal_params - Cx))
    end if
    phie3 = pi / 2 - alpha / 2
    phie4 = beta / 2 - pi / 2

    thetac = theta_crit(n2, phi)
    param_theta = param_theta_from_theta(theta, thetac)

    call determine_critical_phi_points(n3, n4, theta, phie3, phie4, phic3, phic4)

    call param_phi_from_phi(phi, phic3, phie3, phic4, phie4, param_phi, ierr); call CHKERR(ierr)
  end subroutine

  subroutine param_phi_param_theta_from_phi_and_theta_withcoords(wedge_coords3d, phi, theta, param_phi, param_theta, ierr)
    real(ireal_params), intent(in), dimension(:) :: wedge_coords3d ! dim 18
    real(ireal_params), intent(in) :: phi, theta
    real(ireal_params), intent(out) :: param_phi, param_theta
    integer(mpiint), intent(out) :: ierr
    real(ireal_params) :: phic3, phic4, thetac, phie3, phie4
    real(ireal_params), dimension(3) :: n2, n3, n4

    call determine_end_phi_points_plus_geometry(wedge_coords3d, n2, n3, n4, phie3, phie4)

    thetac = theta_crit(n2, phi)
    param_theta = param_theta_from_theta(theta, thetac)

    call determine_critical_phi_points(n3, n4, theta, phie3, phie4, phic3, phic4)

    call param_phi_from_phi(phi, phic3, phie3, phic4, phie4, param_phi, ierr); call CHKERR(ierr)

    !print *,'phi', phi, 'theta', theta, thetac, '=>', param_phi, param_theta, ':', phie3, phie4, phic3, phic4
  end subroutine

  function param_theta_from_theta(theta, thetac) result(param_theta)
    real(ireal_params), intent(in) :: theta, thetac
    real(ireal_params) :: param_theta
    real(ireal_params) :: lb, rb, x1, x2! bounds of local spline

    if (approx(theta, thetac, epsilon(theta) * 10)) then
      param_theta = 0
      return
    end if

    if (theta .lt. thetac) then ! range [ . , 0]
      x1 = -1; x2 = -0
      lb = 0._ireal_params
      rb = thetac
    else ! range [0, .]
      x1 = 0; x2 = 1
      lb = thetac
      rb = pi / 2
    end if
    param_theta = real((theta - lb) / (rb - lb) * (x2 - x1) + x1, ireal_params)
  end function

!> @brief reverse transformation from azimuth, phi(radians) to param_phi[-2,2]
  subroutine param_phi_from_phi(phi, phic3, phie3, phic4, phie4, param_phi, ierr)
    real(ireal_params), intent(in) :: phi, phic3, phic4, phie3, phie4
    real(ireal_params), intent(out) :: param_phi
    integer(mpiint), intent(out) :: ierr

    real(ireal_params) :: lb, rb, x1, x2! bounds of local spline
    real(ireal_params) :: eps = 10 * epsilon(1._ireals)

    if (phi .gt. phic3) then ! range [ . , -1]
      x1 = -2; x2 = -1
      lb = phie3
      rb = phic3
    elseif (phi .lt. phic4) then ! between [ 1, .]
      x1 = 1; x2 = 2
      lb = phic4
      rb = phie4
    else ! range [-1, 1]
      x1 = -1; x2 = 1
      lb = phic3
      rb = phic4
    end if
    if (abs(rb - lb) .lt. epsilon(rb)) then
      !print *,'tiny rb-lb:', rb - lb, ':', rb, lb, &
      !  'phi:', phi, new_line(''), &
      !  ':', phic3         , phie3         , phie4         , phic4         , new_line(''), &
      !  ':', rad2deg(phic3), rad2deg(phie3), rad2deg(phie4), rad2deg(phic4), new_line(''), &
      !  ':', x1, x2
      param_phi = (x1 + x2) / 2
    else
      param_phi = (phi - lb) / (rb - lb) * (x2 - x1) + x1
    end if

    ierr = 0
    if (ldebug) then
      if (.not. is_between(param_phi, &
                           -2._ireal_params - eps, &
                           2._ireal_params + eps)) &
        call CHKERR(1_mpiint, 'not in range '// &
                    '('//ftoa(-2._ireal_params - eps)// &
                    ' '//ftoa(2._ireal_params + eps)//') '// &
                    'param_phi= '//ftoa(param_phi)//new_line('')// &
                    'phi      = '//ftoa(phi)//'('//ftoa(rad2deg(phi))//'deg)'//new_line('')// &
                    'phic3    = '//ftoa(phic3)//'('//ftoa(rad2deg(phic3))//'deg)'//new_line('')// &
                    'phie3    = '//ftoa(phie3)//'('//ftoa(rad2deg(phie3))//'deg)'//new_line('')// &
                    'phic4    = '//ftoa(phic4)//'('//ftoa(rad2deg(phic4))//'deg)'//new_line('')// &
                    'phie4    = '//ftoa(phie4)//'('//ftoa(rad2deg(phie4))//'deg)'//new_line('')// &
                    'x1       = '//ftoa(x1)//new_line('')// &
                    'x2       = '//ftoa(x2)//new_line('')// &
                    'lb       = '//ftoa(lb)//'('//ftoa(rad2deg(lb))//'deg)'//new_line('')// &
                    'rb       = '//ftoa(rb)//'('//ftoa(rad2deg(rb))//'deg)'//new_line(''))
    end if
  end subroutine

  real(real32) function LUT_wedge_aspect_zx_r32(Atop, dz) result(aspect_zx)
    real(real32), intent(in) :: Atop, dz
    aspect_zx = dz / sqrt(Atop * 2)
  end function
  real(real64) function LUT_wedge_aspect_zx_r64(Atop, dz) result(aspect_zx)
    real(real64), intent(in) :: Atop, dz
    aspect_zx = dz / sqrt(Atop * 2)
  end function
  real(real32) function LUT_wedge_dz_r32(Atop, aspect_zx) result(dz)
    real(real32), intent(in) :: Atop, aspect_zx
    dz = sqrt(Atop * 2) * aspect_zx
  end function
  real(real64) function LUT_wedge_dz_r64(Atop, aspect_zx) result(dz)
    real(real64), intent(in) :: Atop, aspect_zx
    dz = sqrt(Atop * 2) * aspect_zx
  end function
end module
