module test_boxmc_3_30
  use m_boxmc, only : t_boxmc, t_boxmc_3_30
  use m_data_parameters, only :     &
    mpiint, iintegers, ireals, ireal_dp,     &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_helper_functions, only: toStr, cstr, colored_str_by_range, deg2rad
  use m_optprop_parameters, only : stddev_atol
  use m_boxmc_geometry, only : setup_default_unit_cube_geometry

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi,theta,dx,dy,dz
  real(ireals), target :: S(30),T(3), S_target(30), T_target(3)
  real(ireals) :: S_tol(30),T_tol(3)
  real(ireal_dp), allocatable :: vertices(:)

  type(t_boxmc_3_30) :: bmc

  integer(mpiint) :: myid,mpierr,numnodes,comm
  character(len=120) :: msg

  real(ireal_dp),parameter :: sigma = 3 ! normal test range for coefficients

  real(ireal_dp),parameter :: atol=1e-3, rtol=1e-2
contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()


    call init_mpi_data_parameters(comm)

    call bmc%init(comm)

    call mpi_comm_size(comm, numnodes, mpierr)
    if(myid.eq.0) print *,numnodes,'Testing Box-MonteCarlo model with tolerances atol/rtol :: ',atol,rtol

    phi   =  0
    theta =  0

    dx = 100
    dy = dx
    dz = 100

    call setup_default_unit_cube_geometry(dx, dy, dz, vertices)

    S_target = zero
    T_target = zero
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    if(myid.eq.0) print *,'Finishing boxmc tests module'
  end subroutine teardown

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src3(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the side center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right and adjacent side vertical streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from bot face towards northern directions
    src = 3
    S_target = 0
    pS(2,1,3) = c1

    pS(2,5,2) = c2
    pS(1,5,2) = c2

    pS(1,2,2) = c3
    pS(2,2,2) = c3
    pS(2,3,3) = c3
    pS(2,5,3) = c3

    pS(2,2,3) = c4

    pS(1,2,1) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src4(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the side center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right and adjacent side vertical streams
    real(ireals), parameter :: c4=.2776 ! adjacent side downward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from top face towards northern directions
    src = 4
    S_target = 0
    pS(2,1,3) = c1
    pS(2,5,2) = c2
    pS(1,5,2) = c2

    pS(1,4,2) = c3
    pS(2,4,2) = c3
    pS(2,3,3) = c3
    pS(2,5,3) = c3

    pS(2,4,3) = c4

    pS(2,2,1) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src5(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the adjacent center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right and adjacent side vertical streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from bot face towards western directions
    src = 5
    S_target = 0
    pS(1,1,2) = c1

    pS(2,3,3) = c2
    pS(1,3,3) = c2

    pS(1,2,3) = c3
    pS(2,2,3) = c3
    pS(1,3,2) = c3
    pS(1,5,2) = c3

    pS(1,2,2) = c4

    pS(1,3,1) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src6(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the adjacent center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right (vertical) and adjacent side (horizontal) streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from top face towards western directions
    src = 6
    S_target = 0
    pS(1,1,2) = c1

    pS(1,3,3) = c2
    pS(2,3,3) = c2

    pS(1,4,3) = c3
    pS(2,4,3) = c3

    pS(1,3,2) = c3
    pS(1,5,2) = c3

    pS(1,4,2) = c4

    pS(2,3,1) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src7(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the adjacent center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right (vertical) and adjacent side (horizontal) streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from bot face towards southern directions
    src = 7
    S_target = 0
    pS(1,1,3) = c1

    pS(1,3,2) = c2
    pS(2,3,2) = c2

    pS(1,2,2) = c3
    pS(2,2,2) = c3

    pS(1,3,3) = c3
    pS(1,5,3) = c3

    pS(1,2,3) = c4

    pS(1,4,1) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src8(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the adjacent center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right (vertical) and adjacent side (horizontal) streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from top face towards southern directions
    src = 8
    S_target = 0
    pS(1,1,3) = c1

    pS(1,3,2) = c2
    pS(2,3,2) = c2

    pS(1,4,2) = c3
    pS(2,4,2) = c3

    pS(1,3,3) = c3
    pS(1,5,3) = c3

    pS(1,4,3) = c4

    pS(2,4,1) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src9(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the adjacent center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right (vertical) and adjacent side (horizontal) streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from bot face towards eastern directions
    src = 9
    S_target = 0
    pS(2,1,2) = c1

    pS(1,5,3) = c2
    pS(2,5,3) = c2

    pS(1,2,3) = c3
    pS(2,2,3) = c3

    pS(2,3,2) = c3
    pS(2,5,2) = c3

    pS(2,2,2) = c4

    pS(1,5,1) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src10(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the adjacent center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right (vertical) and adjacent side (horizontal) streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from top face towards eastern directions
    src = 10
    S_target = 0
    pS(2,1,2) = c1

    pS(1,5,3) = c2
    pS(2,5,3) = c2

    pS(1,4,3) = c3
    pS(2,4,3) = c3

    pS(2,3,2) = c3
    pS(2,5,2) = c3

    pS(2,4,2) = c4

    pS(2,5,1) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src13(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the adjacent center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right (vertical) and adjacent side (horizontal) streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from right face towards top directions
    src = 13
    S_target = 0
    pS(1,1,1) = c1

    pS(1,2,3) = c2
    pS(2,2,3) = c2

    pS(1,3,3) = c3
    pS(2,3,3) = c3

    pS(1,2,1) = c3
    pS(1,4,1) = c3

    pS(1,3,1) = c4

    pS(1,2,2) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src14(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the adjacent center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right (vertical) and adjacent side (horizontal) streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from left face towards top directions
    src = 14
    S_target = 0
    pS(1,1,1) = c1

    pS(1,2,3) = c2
    pS(2,2,3) = c2

    pS(1,5,3) = c3
    pS(2,5,3) = c3

    pS(1,2,1) = c3
    pS(1,4,1) = c3

    pS(1,5,1) = c4

    pS(2,2,2) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src19(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the adjacent center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right (vertical) and adjacent side (horizontal) streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from right face towards forward directions
    src = 19
    S_target = 0
    pS(2,1,3) = c1

    pS(1,2,1) = c2
    pS(2,2,1) = c2

    pS(1,3,1) = c3
    pS(2,3,1) = c3

    pS(2,2,3) = c3
    pS(2,4,3) = c3

    pS(2,3,3) = c4

    pS(1,5,2) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src20(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the adjacent center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right (vertical) and adjacent side (horizontal) streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from left face towards forward directions
    src = 20
    S_target = 0
    pS(2,1,3) = c1

    pS(1,2,1) = c2
    pS(2,2,1) = c2

    pS(1,5,1) = c3
    pS(2,5,1) = c3

    pS(2,2,3) = c3
    pS(2,4,3) = c3

    pS(2,5,3) = c4

    pS(2,5,2) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src27(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the adjacent center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right (vertical) and adjacent side (horizontal) streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from forward face towards bot directions
    src = 27
    S_target = 0
    pS(2,1,1) = c1

    pS(1,4,2) = c2
    pS(2,4,2) = c2

    pS(1,3,2) = c3
    pS(2,3,2) = c3

    pS(2,3,1) = c3
    pS(2,5,1) = c3

    pS(2,4,1) = c4

    pS(1,4,3) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_src30(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=.4347 ! to the adjacent center stream
    real(ireals), parameter :: c2=.0814 ! left and right side along stream
    real(ireals), parameter :: c3=.0281 ! left and right (vertical) and adjacent side (horizontal) streams
    real(ireals), parameter :: c4=.2776 ! adjacent side upward stream
    real(ireals), parameter :: c5=.0125 ! opposite side

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 0._ireal_dp, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation from backward face towards right directions
    src = 30
    S_target = 0
    pS(2,1,2) = c1

    pS(1,5,1) = c2
    pS(2,5,1) = c2

    pS(1,2,1) = c3
    pS(2,2,1) = c3

    pS(2,2,2) = c3
    pS(2,4,2) = c3

    pS(2,5,2) = c4

    pS(2,5,3) = c5

    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine


  @test(npes =[1])
  subroutine test_boxmc_select_cases_diffuse_mid(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)
    real(ireals), parameter :: c1=0.4407 ! straight
    real(ireals), parameter :: c2=0.1398 ! side


    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [0._ireal_dp, 1._ireal_dp/dz, 1._ireal_dp ]

    theta = 0; phi = 0;
    T_target = zero

    ! should send diffuse radiation into vertical angle stream and a bit to the sides
    src = 1
    S_target = 0
    pS(1,1,1) = c1
    pS(1,2,2) = c2
    pS(2,2,2) = c2
    pS(1,2,3) = c2
    pS(2,2,3) = c2
    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 2
    S_target = 0
    pS(2,1,1) = c1
    pS(1,4,2) = c2
    pS(2,4,2) = c2
    pS(1,4,3) = c2
    pS(2,4,3) = c2
    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)


    ! should send diffuse radiation into horizontal angle stream and a bit to the sides
    src = 11
    S_target = 0
    pS(1,1,2) = c1
    pS(1,3,1) = c2
    pS(2,3,1) = c2
    pS(1,3,3) = c2
    pS(2,3,3) = c2
    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 12
    S_target = 0
    pS(2,1,2) = c1
    pS(1,5,1) = c2
    pS(2,5,1) = c2
    pS(1,5,3) = c2
    pS(2,5,3) = c2
    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    ! should send diffuse radiation into horizontal angle stream and a bit to the sides
    src = 21
    S_target = 0
    pS(1,1,3) = c1
    pS(1,4,1) = c2
    pS(2,4,1) = c2
    pS(1,3,2) = c2
    pS(2,3,2) = c2
    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)

    src = 22
    S_target = 0
    pS(2,1,3) = c1
    pS(1,2,1) = c2
    pS(2,2,1) = c2
    pS(1,5,2) = c2
    pS(2,5,2) = c2
    call bmc%get_coeff(comm,bg,src,.False.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_diffuse',src
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  @test(npes =[1])
  subroutine test_boxmc_select_cases_direct_srctopface(this)
  class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), pointer :: pS(:,:,:) ! dim (2 directions[up/down], 5 streams, 3 axis)

    pS(1:2, 1:5, 1:3) => S_target(1:30)

    bg  = [1._ireal_dp/dz, 1e-1_ireal_dp, 1._ireal_dp ]

    ! should send diffuse radiation into vertical angle stream
    theta = 0
    phi = 0; src = 1
    T_target = [real(ireals) :: exp(-sum(bg(1:2))*dz), 0, 0]
    S_target = 0; pS(2,1,1) = .367879
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta

    call check(S_target,T_target, S,T, msg=msg)

    theta = 180
    phi = 0; src = 1
    T_target = [real(ireals) :: exp(-sum(bg(1:2))*dz), 0, 0]
    S_target = 0; pS(1,1,1) = .367879
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    ! a bit towards north
    theta = 45
    phi = 0; src = 1
    T_target = [real(ireals) :: 0, 0, 6.42878E-02]
    S_target = 0
    pS(2,4,3) = 4.70484E-01
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    theta = 180-45
    phi = 0; src = 1
    T_target = [real(ireals) :: 0, 0, 6.42878E-02]
    S_target = 0
    pS(2,2,3) = 4.70484E-01
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    ! a bit towards south
    theta = 45
    phi = 180; src = 1
    T_target = [real(ireals) :: 0, 0, 6.42878E-02]
    S_target = 0
    pS(1,4,3) = 4.70484E-01
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    theta = 180-45
    phi = 180; src = 1
    T_target = [real(ireals) :: 0, 0, 6.42878E-02]
    S_target = 0
    pS(1,2,3) = 4.70484E-01
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    ! a bit towards east
    theta = 45
    phi = 90; src = 1
    T_target = [real(ireals) :: 0, 6.42878E-02, 0]
    S_target = 0
    pS(2,4,2) = 4.70484E-01
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    theta = 180-45
    phi = 90; src = 1
    T_target = [real(ireals) :: 0, 6.42878E-02, 0]
    S_target = 0
    pS(2,2,2) = 4.70484E-01
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    ! a bit towards west
    theta = 45
    phi = 270; src = 1
    T_target = [real(ireals) :: 0, 6.42878E-02, 0]
    S_target = 0
    pS(1,4,2) = 4.70484E-01
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)

    theta = 180-45
    phi = 270; src = 1
    T_target = [real(ireals) :: 0, 6.42878E-02, 0]
    S_target = 0
    pS(1,2,2) = 4.70484E-01
    call bmc%get_coeff(comm,bg,src,.True.,phi,theta,vertices,S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    write(msg,*) ' test_boxmc_select_cases_direct_srctopface',src,phi,theta
    call check(S_target,T_target, S,T, msg=msg)
  end subroutine

  subroutine check(S_target,T_target, S,T, msg)
    real(ireals),intent(in),target,dimension(:) :: S_target,T_target, S,T
    character(len=*),optional :: msg

    integer(iintegers) :: i
    real(ireals) :: pS(5,6), pS_target(5,6) ! dim (5 streams, 6 sides)
    character(default_str_len) :: local_msgS, local_msgT
    real(ireals), parameter :: test_atol = real(atol, ireals) * real(sigma, ireals)

    real(ireals), parameter :: color_limits(5) = [0., 1e-5, 1e-3, .1, 1.]
    character(len=*), parameter :: colors(4) = [character(len=10):: 'black', 'blue', 'green', 'red']

    real(ireals), parameter :: diff_color_limits(8) = [-1., -.1, -1e-3, -1e-5, 1e-5, 1e-3, .1, 1.]
    character(len=*), parameter :: diff_colors(7) = [character(len=10)::'red', 'green', 'blue', 'black', 'blue', 'green', 'red']

    if(myid.eq.0) then
      print*,''

      if( present(msg) ) then
        write(local_msgS,*) trim(msg),':: Diffuse boxmc coefficient not as '
        write(local_msgT,*) trim(msg),':: Direct  boxmc coefficient not as '
        print *,msg
      else
        write(local_msgS,*) 'Diffuse boxmc coefficient not as '
        write(local_msgT,*) 'Direct  boxmc coefficient not as '
      endif

      do i=1,3
        pS(:,2*(i-1)+1) = S(2*(i-1)*5+1 : 2*i*5 : 2)
        pS(:,2*(i-1)+2) = S(2*(i-1)*5+2 : 2*i*5 : 2)
        pS_target(:,2*(i-1)+1) = S_target(2*(i-1)*5+1 : 2*i*5 : 2)
        pS_target(:,2*(i-1)+2) = S_target(2*(i-1)*5+2 : 2*i*5 : 2)
      enddo

      print*,'---------------------'
      do i=1,6
        write(*, FMT='( " diffuse :: ",I0," :: ", A)' ) &
          & i, colored_str_by_range(pS(:,i), color_limits, colors)
      enddo
      print *,cstr('---', 'blue')
      do i=1,6
        write(*, FMT='( " target  :: ",I0," :: ", A)' ) &
          & i, colored_str_by_range(pS_target(:,i), color_limits, colors)
      enddo
      print *,cstr('---', 'blue')
      do i=1,6
        write(*, FMT='( " diff side  ",I0," :: ", A )' ) &
          & i, colored_str_by_range(pS_target(:,i) - pS(:,i), diff_color_limits, diff_colors)
      enddo
      print*,''
      write(*, FMT='( " direct  ::  :: ", 8(es12.5) )' ) T
      write(*, FMT='( " target  ::  :: ", 8(es12.5) )' ) T_target
      write(*, FMT='( " diff    ::  :: ", 8(es12.5) )' ) T_target-T
      print*,'---------------------'
      print*,''

      @assertEqual(S_target, S, test_atol, local_msgS )
      @assertLessThanOrEqual   (zero, S)
      @assertGreaterThanOrEqual(one , S)

      @assertEqual(T_target, T, test_atol, local_msgT )
      @assertLessThanOrEqual   (zero, T)
      @assertGreaterThanOrEqual(one , T)
    endif
  end subroutine

end module
