@test(npes =[1])
subroutine test_eddington(this)

    use m_eddington
    use m_data_parameters, only: ireals, default_str_len

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    real(ireals), parameter :: tol=1e-6
    real(ireals) :: inp(4) ! tau, omega0, g, mu0
    real(ireals) :: out(5) ! a11, a12, a13, a23, a33
    real(ireals) :: targ(5)! target values, computed with double precision at LMU_MIM

    print *, 'Checking some Eddington coefficient edge cases'

    inp = [ 0.4824516550E-01, 0.5542391539, 0.4550637007, 1.00000000000]
    targ= [ 9.526110e-01, 5.290933e-03, 4.103254e-03, 2.144391e-02, 9.529001e-01]
    call calc(inp, out)
    call check_out(inp, out, targ)

    inp = [0.2018013448, 0.3797843754, 0.4556422830, 1.00000000000]
    targ= [0.77855106460434387, 8.2984956871573792E-009, 9.7582138200468863E-003, 5.1463626191620593E-002, 0.81725726426147882]
    call calc(inp, out)
    call check_out(inp, out, targ)

    inp = [ 0.1731484532, 0.6180083156, 0.4121485054, 1.00000000000 ]
    targ= [0.85001352568220123, 2.6158207726886692E-002, 1.8420479735352102E-002, 7.3447050483338217E-002, 0.84101275451228208]
    call calc(inp, out)
    call check_out(inp, out, targ)

    inp = [ 0.1931012775E-04, 0.4384377003, -0.0000000000E+00, 0.7070999742 ]
    targ= [0.99997467389811490, 3.6386575613389649E-006, 5.9864742811154891E-006, 5.9864869794782480E-006, 0.99997269146543133]
    call calc(inp, out)
    call check_out(inp, out, targ)

    inp = [ 4.895462513 , 0.3104626103E-05 , -0.0000000000E+00 , 0.49999997019767761  ] ! resonance case
    targ= [2.0667880575476368E-004, 0.0000000000000000, 7.7215617940279102E-007, 1.6439929979903860E-009, 5.5957078888272847E-005]
    call calc(inp, out)
    call check_out(inp, out, targ)

    inp = [ 3.2662250689525390E-011 , 0.99999171495417127 , 0.0000000000000000  , 0.17364817766693041 ]
    targ= [0.99999999997550049, 2.4499513529008254E-011, 9.4056105052780989E-011, 9.4060580516466263E-011, 0.99999999981190557]
    call calc(inp, out)
    call check_out(inp, out, targ)

    inp = [ 2.9317851124478626E-012, 1.0000000000000000  , 0.0000000000000000 , 0.17364817766693041 ]
    targ= [0.99999999999779732, 2.2026824808563106E-012, 8.4443208651779663E-012, 8.4384068595938508E-012, 0.99999999998311651]
    call calc(inp, out)
    call check_out(inp, out, targ)

    inp = [1.93321303E-10, 0.999984443 , 2.22044605E-16, 0.17364817766693041 ]
    targ= [0.99999999985500665, 1.4499335065920604E-010, 5.5664139018583103E-010, 5.5664738411040951E-010, 0.99999999888670699]
    call calc(inp, out)
    call check_out(inp, out, targ)

    inp = [7.89528581E-11, 0.999988437  , 2.22044605E-16 , 0.17364817766693041 ]
    targ= [0.99999999994077626, 5.9209526170889148E-011, 2.2734076817292897E-010, 2.2733741157516502E-010, 0.99999999954532859]
    call calc(inp, out)
    call check_out(inp, out, targ)

    inp = [ 1.3865453490508738E-011 , 0.99999499320987351 , 0.0000000000000000 , 0.17364817766693041 ]
    targ= [0.99999999998959765, 1.0402345651527867E-011, 3.9931125946960987E-011, 3.9926650483275707E-011, 0.99999999992015198]
    call calc(inp, out)
    call check_out(inp, out, targ)

    inp = [ 113.59224626216431 , 2.7005225550306174E-008 , 0.0000000000000000  , 0.17364817766693041 ]
    targ= [0.0, 0.0000000000000000, 9.6352077372196252E-009, 4.1335735289420085E-024, 0.0]
    call calc(inp, out)
    call check_out(inp, out, targ)

    inp = [ 1.4503301E-02, 1.5233955E-12, 1.1920928E-07, 1.000000 ]
    targ= [   0.974944890, 0.00000000, 8.45933545E-10, 8.45891468E-10, 0.985601366]
    call calc(inp, out)
    call check_out(inp, out, targ)

  contains
    subroutine calc(inp, out)
        real(ireals) :: inp(4)
        real(ireals) :: out(5)

!        print *, 'inp ::', inp
        call eddington_coeff_zdun(inp(1), inp(2), inp(3), inp(4), out(1), out(2), out(3), out(4), out(5))
!        call eddington_coeff_bm(inp(1), inp(2), inp(3), inp(4), out(1), out(2), out(3), out(4), out(5))
!        print *, 'zdun ::', out
!        print *, ''
    end subroutine
    subroutine check_out(inp, out, targ)
        real(ireals), intent(in) :: inp(:)
        real(ireals), intent(in) :: out(:), targ(:)
        character(default_str_len) :: msg

        print *, ''
        print *, '-----------------------'
        print *, 'inp    ::', inp
        print *, 'out    ::', out
        print *, 'target ::', targ
        print *, ''
        print *, 'diff   ::', targ-out
        print *, '-----------------------'
        print *, ''


        write(msg, *) "Eddington coefficient not as "
        @assertEqual(targ , out, tol, msg)

    end subroutine
end subroutine

!@test
!subroutine test_eddington_zdun()
!
!   use m_eddington
!   use m_data_parameters, only: ireals, iintegers, zero, one, pi
!
!   use pfunit_mod
!
!
!   implicit none
!
!   real(ireals), parameter :: tol=1e-6
!
!   real(ireals) :: tau, w0, g, mu0
!   real(ireals) :: a11, a12, a13, a23, a33
!
!
!   tau = zero
!   w0  = one
!   g   = one
!   mu0 = one
!
!   call eddington_coeff_zdun(tau, w0, g, mu0, a11, a12, a13, a23, a33)
!   print *, a11, a12, a13, a23, a33
!
!!   @assertAll([one, zero, zero], [a11, a12, a13], tol)
!
!   @assertEqual(one , a11, tol)
!   @assertEqual(zero, a12, tol)
!   @assertEqual(zero, a13, tol)
!   @assertEqual(zero, a23, tol)
!   @assertEqual(one , a33, tol)
!end subroutine
