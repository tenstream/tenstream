module m_LUT_param_phi
  use m_data_parameters, only : irealLUT, mpiint, ireal_dp
  use m_data_parameters, only: pi=>pi_irealLUT
  use m_helper_functions, only : angle_between_two_vec, rad2deg, approx, CHKERR
  use m_boxmc_geometry, only: wedge_halfspaces
  implicit none

  contains
!> @brief reverse transformation from azimuth(radians) to param_phi[-2,2]
function param_phi_from_azimuth(phi, wedge_C) result (param_phi)
  real(irealLUT), intent(in) :: phi, wedge_C(:)
  real(irealLUT) :: param_phi

  real(irealLUT) :: alpha, beta
  real(irealLUT) :: lb, rb, x1, x2! bounds of local spline

  associate( pA => [0._irealLUT, 0._irealLUT], pB => [1._irealLUT, 0._irealLUT], pC => wedge_C )
    alpha = angle_between_two_vec(pB-pA, pC-pA)
    beta  = angle_between_two_vec(pA-pB, pC-pB)

    if(phi.gt.pi/2-alpha) then ! range [ . , -1]
      x1 = -2; x2 = -1
      lb = pi/2 - alpha / 2
      rb = pi/2 - alpha
    elseif (phi.lt.beta-pi/2) then ! between [ 1, .]
      x1 = 1; x2 = 2
      lb = beta - pi/2
      rb = beta/2 - pi/2
    else ! range [-1, 1]
      x1 = -1; x2 = 1
      lb = pi/2 - alpha
      rb = beta - pi/2
    endif
    param_phi = (x2-x1) / (rb - lb) * (phi - lb) + x1
  end associate
end function

!> @brief translate from the parameterized azimuth [-2,2] to azimuth [radians]
!> @details we use local splines so that it ensures we have mapped:
!>   alpha and beta are the inner angles of a wedge triangle between AB and AC or BA and BC respectively
!>   90 - alpha/2 -> -2
!>   90 - alpha   -> -1
!>   beta   - 90  ->  1
!>   beta/2 - 90  ->  2
function azimuth_from_param_phi(param_phi, wedge_C) result (phi)
  real(irealLUT), intent(in) :: param_phi, wedge_C(:)
  real(irealLUT) :: phi

  real(irealLUT) :: alpha, beta
  real(irealLUT) :: lb, rb, x1, x2! bounds of local spline

  alpha = angle_between_two_vec([ 1._irealLUT, 0._irealLUT], wedge_C)
  beta  = angle_between_two_vec([-1._irealLUT, 0._irealLUT], wedge_C - [1._irealLUT, 0._irealLUT])

  if(param_phi.lt.-1._irealLUT) then ! range [ . , -1]
    x1 = -2; x2 = -1
    lb = pi/2 - alpha / 2
    rb = pi/2 - alpha
  elseif (param_phi.gt.1._irealLUT) then ! between [ 1, .]
    x1 = 1; x2 = 2
    lb = beta - pi/2
    rb = beta/2 - pi/2
  else ! range [-1, 1]
    x1 = -1; x2 = 1
    lb = pi/2 - alpha
    rb = beta - pi/2
  endif
  phi = (rb - lb) / (x2-x1) * (param_phi - x1) + lb
end function

  subroutine phi_crit(side_normal, theta, phic, ierr)
    ! solve euqation system that arises for condition
    ! sunvec_crit = [sin(phi_crit) * sin(theta), cos(phi_crit) * sin(theta), -cos(theta)]
    ! left_side_face_normal .dot. sunvec_crit === 0
    ! i.e. solve( 0 = n1 * sin ( phi ) * sin(theta) + n2 * cos(phi) * sin(theta) - n3 * cos(theta), phi)

    use m_helper_functions, only: solve_quadratic
    real(irealLUT), intent(in) :: side_normal(3)
    real(irealLUT), intent(in) :: theta ! sun zenith in [rad]
    real(irealLUT), intent(out) :: phic
    integer(mpiint), intent(out) :: ierr

    real(irealLUT) :: a, b, c, z, sol(2)
    real(irealLUT) :: x(2), y(2)
    !real(irealLUT) :: st, ct

    z = -cos(theta)
    a = side_normal(1)**2 + side_normal(2)**2
    b = 2*side_normal(2)*side_normal(3) * z
    c = (side_normal(1)**2 + side_normal(3)**2) * z**2 - side_normal(1)**2

    call solve_quadratic(a,b,c, y, ierr)
    if(ierr.ne.0) return

    x = (-side_normal(2)*y - side_normal(3) * z) / side_normal(1)

    sol = atan2(x,y)
    if(abs(sol(1)).ge.(abs(sol(2)))) then
      phic = sol(2)
    else
      phic = sol(1)
    endif

    !associate( n1=>side_normal(1), n2=>side_normal(2), n3=>side_normal(3) )
    !  st = sin(theta)
    !  ct = cos(theta)
    !  if(approx(n2,0._irealLUT) .and. approx(n3,0._irealLUT) .and. .not. approx(n1*st,0._irealLUT)) then
    !    phic = 0._irealLUT
    !    return
    !  endif
    !  if(approx(n2,0._irealLUT) .and. &
    !    approx(ct, 0._irealLUT) .and. &
    !    .not. approx(n3, 0._irealLUT) .and. &
    !    .not. approx(n1*st, 0._irealLUT) ) then
    !    phic = 0._irealLUT
    !    return
    !  endif
    !  if(.not.approx(n2,0._irealLUT) .and. approx(st, -n3*ct/n2)) then
    !    phic = pi_irealLUT
    !    return
    !  endif
    !  z = n2*st+n3*ct
    !  if(.not.approx(z,0._irealLUT)) then
    !    b = n1**2*st**2 + n2**2*st**2 - n3**2*ct**2
    !    if(b.lt.0._irealLUT) ierr = 1
    !    a = n1*st / z
    !    c = n2*st + n3*ct

    !    if(approx(n1**2*st**2 - n1*st*sqrt(b) + n2**2*st**2 + n2*n3*st*ct, 0._irealLUT)) ierr = 2
    !    sol(1) = 2*atan(a - sqrt(b)/c)
    !    if(approx(n1**2*st**2 + n1*st*sqrt(b) + n2**2*st**2 + n2*n3*st*ct, 0._irealLUT)) ierr = 3
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
    ! solve euqation system that arises for condition
    ! sunvec_crit = [sin(phi) * sin(theta_crit), cos(phi) * sin(theta_crit), -cos(theta_crit)]
    ! base_face_normal .dot. sunvec_crit === 0
    ! i.e. solve( 0 = n1 * sin ( phi ) * sin(theta) + n2 * cos(phi) * sin(theta) - n3 * cos(theta), theta)

    real(irealLUT), intent(in) :: side_normal(3)
    real(irealLUT), intent(in) :: phi ! sun azimuth in [rad]
    real(irealLUT) :: thetac

    real(irealLUT) :: discr, frac

    !discr = side_normal(1)*sin(phi) + side_normal(2)*cos(phi)

    !frac = side_normal(3) / discr
    !thetac = atan(frac)

    associate( n1=>side_normal(1), n2=>side_normal(2), n3=>side_normal(3) )
      discr = (n1*sin(phi) + n2*cos(phi))**2 + n3**2
      frac = n3 / sqrt(discr)
      thetac = asin( frac )
    end associate

    if(discr.le.0._irealLUT) then
      print *,'discr', discr, side_normal, norm2(side_normal), '->', side_normal(3) / discr, 'thetac', thetac
      call CHKERR(1_mpiint, 'no solution for theta crit')
    endif
  end function

  subroutine iterative_phi_theta_from_param_phi_and_param_theta(wedge_coords3d, param_phi, param_theta, phi, theta, ierr)
    real(irealLUT), intent(in), dimension(:) :: wedge_coords3d ! dim 18
    real(irealLUT), intent(in) :: param_phi, param_theta
    real(irealLUT), intent(out) :: phi, theta
    integer(mpiint), intent(out) :: ierr
    real(irealLUT) :: phic3, phic4, thetac, phie3, phie4
    real(irealLUT) :: last_phi, last_theta, new_phi
    real(irealLUT), dimension(3) :: n2, n3, n4
    integer(mpiint) :: iter
    integer(mpiint), parameter :: Niter=100
    real(irealLUT), parameter :: eps=epsilon(phi)*1

    call determine_end_phi_points_plus_geometry(wedge_coords3d, n2, n3, n4, phie3, phie4)

    phi = 0._irealLUT
    theta = 0._irealLUT
    last_phi   = phi
    last_theta = theta

    do iter=1,Niter
      thetac = theta_crit(n2, phi)
      theta  = theta_from_param_theta(param_theta, thetac)

      call determine_critical_phi_points(&
        n3, n4, theta, phie3, phie4, &
        phic3, phic4)

      call phi_from_param_phi(param_phi, phic3, phie3, phic4, phie4, new_phi, ierr); call CHKERR(ierr)
      phi = phi + (new_phi - phi) * (1._irealLUT - real(iter-1, irealLUT)/real(Niter, irealLUT))**2
      !print *,iter, 'param_phi', param_phi, 'param_theta', param_theta, ':', new_phi, '=>', phi, theta, 'thetac', thetac
      if(approx(last_phi,phi,eps) .and. approx(last_theta, theta, eps)) then
        ierr = 0
        return
      endif
      last_phi = phi
      last_theta = theta
    enddo
    ierr = iter
  end subroutine

  subroutine determine_end_phi_points_plus_geometry(wedge_coords3d, &
      n2, n3, n4, phie3, phie4)
    real(irealLUT), intent(in), dimension(:) :: wedge_coords3d ! dim 18
    real(irealLUT), intent(out)  :: n2(:), n3(:), n4(:), phie3, phie4

    real(ireal_dp) :: origins(3,5), normals(3,5)
    real(irealLUT) :: alpha, beta

    call wedge_halfspaces(real(wedge_coords3d, ireal_dp), origins, normals)

    n2 = real(normals(:,2), irealLUT)
    n3 = real(normals(:,3), irealLUT)
    n4 = real(normals(:,4), irealLUT)

    associate(A => wedge_coords3d(1:3), B => wedge_coords3d(4:6), C => wedge_coords3d(7:9))
      alpha = angle_between_two_vec(B-A, C-A)
      beta  = angle_between_two_vec(A-B, C-B)
    end associate
    phie3 = pi/2 - alpha / 2
    phie4 = beta/2 - pi/2
  end subroutine

  subroutine determine_critical_phi_points(&
      n3, n4, theta, phie3, phie4, &
      phic3, phic4)
    real(irealLUT), intent(in)  :: n3(:), n4(:), theta, phie3, phie4
    real(irealLUT), intent(out) :: phic3, phic4
    integer(mpiint) :: ierr

    call phi_crit(n3, theta, phic3, ierr)
    if(ierr.ne.0) phic3 = phie3

    call phi_crit(n4, theta, phic4, ierr)
    if(ierr.ne.0) phic4 = phie4

    if(phic3.lt.phie4) phic3 = phie3 ! use phie3 if there is no solution
    if(phic4.gt.phie3) phic4 = phie4 ! use phie4 if there is no solution

    phic3 = min(phie3, phic3)
    phic4 = max(phie4, phic4)
  end subroutine

  function theta_from_param_theta(param_theta, thetac) result(theta)
    real(irealLUT), intent(in) :: param_theta, thetac
    real(irealLUT) :: theta
    real(irealLUT) :: lb, rb, x1, x2! bounds of local spline

    if(param_theta.le.0._irealLUT) then ! range [ . , 0]
      x1 = -1; x2 = -0
      lb = 0._irealLUT
      rb = thetac
    else ! range [0, .]
      x1 = 0; x2 = 1
      lb = thetac
      rb = pi/2
    endif
    theta = (rb - lb) / (x2-x1) * (param_theta - x1) + lb
  end function

  subroutine phi_from_param_phi(param_phi, phic3, phie3, phic4, phie4, phi, ierr)
    real(irealLUT), intent(in)  :: param_phi, phic3, phic4, phie3, phie4
    real(irealLUT), intent(out) :: phi
    integer(mpiint), intent(out) :: ierr
    real(irealLUT) :: lb, rb, x1, x2! bounds of local spline

    if(param_phi.lt.-1._irealLUT) then ! range [ . , -1]
      x1 = -2; x2 = -1
      lb = phie3
      rb = phic3
    elseif (param_phi.gt.1._irealLUT) then ! between [ 1, .]
      x1 = 1; x2 = 2
      lb = phic4
      rb = phie4
    else ! range [-1, 1]
      x1 = -1; x2 = 1
      lb = phic3
      rb = phic4
    endif
    phi = (rb - lb) / (x2-x1) * (param_phi - x1) + lb

    ierr = 0
  end subroutine

  subroutine param_phi_param_theta_from_phi_and_theta_withnormals(n2, n3, n4, Cx, Cy, phi, theta, &
      param_phi, param_theta, ierr)
    real(irealLUT), dimension(:), intent(in) :: n2, n3, n4 !normals of base, left, right face, dim 3
    real(irealLUT), intent(in)  :: Cx, Cy     ! Cx, Cy are local coords of C-point
    real(irealLUT), intent(in)  :: phi, theta ! phi, theta are local azimuth and zenith
    real(irealLUT), intent(out) :: param_phi, param_theta
    integer(mpiint), intent(out) :: ierr
    real(irealLUT) :: alpha, beta
    real(irealLUT) :: phic3, phic4, thetac, phie3, phie4

    alpha = atan(Cy/Cx)
    beta  = atan(Cy/(1._irealLUT - Cx))
    phie3 = pi/2 - alpha / 2
    phie4 = beta/2 - pi/2

    thetac = theta_crit(n2, phi)
    param_theta = param_theta_from_theta(theta, thetac)

    call determine_critical_phi_points(n3, n4, theta, phie3, phie4, phic3, phic4)

    call param_phi_from_phi(phi, phic3, phie3, phic4, phie4, param_phi, ierr); call CHKERR(ierr)
  end subroutine

  subroutine param_phi_param_theta_from_phi_and_theta_withcoords(wedge_coords3d, phi, theta, param_phi, param_theta, ierr)
    real(irealLUT), intent(in), dimension(:) :: wedge_coords3d ! dim 18
    real(irealLUT), intent(in)  :: phi, theta
    real(irealLUT), intent(out) :: param_phi, param_theta
    integer(mpiint), intent(out) :: ierr
    real(irealLUT) :: phic3, phic4, thetac, phie3, phie4
    real(irealLUT), dimension(3) :: n2, n3, n4

    call determine_end_phi_points_plus_geometry(wedge_coords3d, n2, n3, n4, phie3, phie4)

    thetac = theta_crit(n2, phi)
    param_theta = param_theta_from_theta(theta, thetac)

    call determine_critical_phi_points(n3, n4, theta, phie3, phie4, phic3, phic4)

    call param_phi_from_phi(phi, phic3, phie3, phic4, phie4, param_phi, ierr); call CHKERR(ierr)

    !print *,'phi', phi, 'theta', theta, '=>', param_phi, param_theta
  end subroutine

  function param_theta_from_theta(theta, thetac) result(param_theta)
    real(irealLUT), intent(in) :: theta, thetac
    real(irealLUT) :: param_theta
    real(irealLUT) :: lb, rb, x1, x2! bounds of local spline

    if(approx(theta, thetac, epsilon(theta)*10)) then
      param_theta = 0
      return
    endif

    if(theta.lt.thetac) then ! range [ . , 0]
      x1 = -1; x2 = -0
      lb = 0._irealLUT
      rb = thetac
    else ! range [0, .]
      x1 = 0; x2 = 1
      lb = thetac
      rb = pi/2
    endif
    param_theta = (theta - lb) / (rb - lb) * (x2-x1) + x1
  end function

!> @brief reverse transformation from azimuth, phi(radians) to param_phi[-2,2]
  subroutine param_phi_from_phi(phi, phic3, phie3, phic4, phie4, param_phi, ierr)
    real(irealLUT), intent(in)  :: phi, phic3, phic4, phie3, phie4
    real(irealLUT), intent(out) :: param_phi
    integer(mpiint), intent(out) :: ierr

    real(irealLUT) :: lb, rb, x1, x2! bounds of local spline

    if(phi.gt.phic3) then ! range [ . , -1]
      x1 = -2; x2 = -1
      lb = phie3
      rb = phic3
    elseif (phi.lt.phic4) then ! between [ 1, .]
      x1 = 1; x2 = 2
      lb = phic4
      rb = phie4
    else ! range [-1, 1]
      x1 = -1; x2 = 1
      lb = phic3
      rb = phic4
    endif
    param_phi = (phi - lb) / (rb - lb) * (x2-x1) + x1

    ierr = 0
  end subroutine

end module
