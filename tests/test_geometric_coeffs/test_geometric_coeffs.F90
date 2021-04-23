module test_geometric_coeffs
  use m_boxmc, only : t_boxmc, t_boxmc_3_10
  use m_data_parameters, only :     &
    mpiint, iintegers, ireals, irealLUT, ireal_dp,     &
    one, zero, i1, default_str_len, &
    init_mpi_data_parameters
  use m_optprop_parameters, only : stddev_atol
  use m_boxmc_geometry, only : setup_default_unit_cube_geometry
  use m_geometric_coeffs, only : dir2dir3_geometric_coeffs
  use m_helper_functions, only : spherical_2_cartesian, cstr, toStr

  use pfunit_mod
  implicit none

  real(ireal_dp) :: bg(3), phi,theta,dx,dy,dz
  real(ireals) :: S(10),T(3), S_target(10), T_target(3)
  real(ireals) :: S_tol(10),T_tol(3)
  real(ireal_dp), allocatable :: vertices(:)

  type(t_boxmc_3_10) :: bmc_3_10

  integer(mpiint) :: myid,mpierr,numnodes,comm
  character(len=120) :: msg

  real(ireal_dp),parameter :: sigma = 3 ! normal test range for coefficients

  real(ireal_dp),parameter :: atol=1e-3, rtol=1e-2

  logical, parameter :: ldebug = .False.
contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()


    call init_mpi_data_parameters(comm)

    call bmc_3_10%init(comm)

    call mpi_comm_size(comm, numnodes, mpierr)
    if(myid.eq.0) print *,numnodes,'Testing Box-MonteCarlo model with tolerances atol/rtol :: ',atol,rtol

    phi   =  0
    theta =  0

    dx = 100
    dy = dx
    dz = 50
    call setup_default_unit_cube_geometry(dx, dy, dz, vertices)

    S_target = zero
    T_target = zero
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    if(myid.eq.0) print *,'Finishing boxmc tests module'
  end subroutine teardown


  !@test(npes =[1])
  subroutine test_boxmc_dir2dir3_geometric_coeffs_vs_online_monte_carlo(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), allocatable :: verts(:), verts_dtd(:)
    real( ireals), parameter :: dx=1, dy= dx, dz=1 * dx
    real(irealLUT) :: v(9), v_mc(9)
    real(ireals) :: sundir(3)
    integer(iintegers) :: itheta, iphi

    bg = [real(ireal_dp) :: - log(5e-1)/dz, 0e-0/dz, 1./2 ]
    S_target = zero
    iphi=30
    itheta=5

    ! 1. Test: Einzelbox
    ! ein Winkel: Wahrheit -> trgt(9)
    ! geometric coeffs -> assertequal
    ! weitere Winkel

    ! 2. Test: Symmetrie
    ! 2 Boxen unterschiedlich verschoben (in x/y geschert)
    ! Winkel drehen -> selbe Werte
    ! Spiegeln -> selbe Werte

    ! 3. Test: Analytische Werte in Schleife rechtwinklige Box
    ! schräg einstrahlen (45°) -> x = y

    ! 4. Test: Gaußhügel starten -> abkack Werte -> Special Case

    do iphi=180,360,30
      do itheta=10,50,20
        phi = real(iphi, ireals)
        theta = real(itheta, ireals)

        call setup_default_unit_cube_geometry(dx, dy, dz, verts)
        verts_dtd = verts
        verts_dtd([9,12,21,24]) = verts_dtd([9,12,21,24]) + dz / 4

        do src = 1,3
          call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
          v(src:3**2:3) = real(T, irealLUT)
        enddo

        sundir = spherical_2_cartesian(real(phi, ireals), real(theta, ireals))

        print *, cstr('regular not corrected', 'red')
        print *, 'src z', v(1:9:3)
        print *, 'src x', v(2:9:3)
        print *, 'src y', v(3:9:3)

        call dir2dir3_geometric_coeffs(verts_dtd, sundir * [-one, -one, one], bg, v)

        print *, cstr('regular corrected', 'blue')
        print *, 'src z', v(1:9:3)
        print *, 'src x', v(2:9:3)
        print *, 'src y', v(3:9:3)

        do src = 1,3
          call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts_dtd, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
          v_mc(src:3**2:3) = real(T, irealLUT)
        enddo

        print *, cstr('montecarlo distorted', 'green')
        print *, 'src z', v_mc(1:9:3)
        print *, 'src x', v_mc(2:9:3)
        print *, 'src y', v_mc(3:9:3)

        @assertEqual(v_mc, v, max(maxval(v_mc)*0.05_irealLUT, 1e-6_irealLUT), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
      enddo
    enddo

  end subroutine

  @test(npes =[1])
  subroutine test_geometric_coeffs_distorted_box_no_scatter(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src, i
    real(ireals), allocatable :: verts(:)
    real(irealLUT) :: v(9), v_symmetry(9), v_mc(9), v_mc_225_40(20,9)
    real(ireals) :: sundir(3), sundir_symmetry(3), verts_dtd(24), verts_dtd_symmetry(24)
    real(ireals), parameter :: dx=1, dy=dx, dz=dx

    !bg = [real(ireal_dp) :: 1e-6, 1, 0.85]
    bg = [real(ireal_dp) :: 1e-6, 0, 0.85]

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)

    S_target = zero

    phi = 225._ireals
    theta = 40._ireals


    sundir = spherical_2_cartesian(phi, theta) * [-one, -one, one]
    sundir_symmetry = spherical_2_cartesian(phi - 180, theta) * [-one, -one, one]
    print *, 'sundir', sundir
    print *, 'sundir_symmetry', sundir_symmetry

    v_mc_225_40(1,:) = [ real(irealLUT) :: &
      0.579161942, 0.880491078, 0.880592525, 0.210434467, 0.00000000, 0.119407229, 0.210403204, 0.119508661, 0.00000000 &
      ]
    v_mc_225_40(2,:) = [ real(irealLUT) :: &
      0.434624285, 0.829745471, 0.829779685, 0.282772988, 0.00000000, 0.170220003, 0.282602191, 0.170254186, 0.00000000 &
      ]
    v_mc_225_40(3,:) = [ real(irealLUT) :: &
      0.363792837, 0.801579118, 0.801497698, 0.318097949, 0.00000000, 0.198501915, 0.318108648, 0.198420525, 0.00000000 &
      ]
    v_mc_225_40(4,:) = [ real(irealLUT) :: &
      0.321538657, 0.783636987, 0.783575475, 0.339426547,  0.00000000, 0.216424137, 0.339034200, 0.216362625, 0.00000000 &
      ]
    v_mc_225_40(5,:) = [ real(irealLUT) :: &
      0.294336140, 0.771124899, 0.771123469, 0.352795750,  0.00000000, 0.228876129, 0.352867484, 0.228874698,  0.00000000 &
      ]
    v_mc_225_40(6,:) = [ real(irealLUT) :: &
      0.274791181, 0.762104392, 0.761982083, 0.362527847,  0.00000000, 0.238017485, 0.362680346, 0.237895191,  0.00000000 &
      ]
    v_mc_225_40(7,:) = [ real(irealLUT) :: &
      0.260479659, 0.755124331, 0.755310655, 0.369733512,  0.00000000, 0.244688913, 0.369786203, 0.244875193,  0.00000000 &
      ]
    v_mc_225_40(8,:) = [ real(irealLUT) :: &
      0.249499917, 0.749804854, 0.749615192, 0.375236362,  0.00000000, 0.250384361, 0.375263065, 0.250194699,  0.00000000 &
      ]
    v_mc_225_40(9,:) = [ real(irealLUT) :: &
      0.240674764, 0.745009124, 0.745452464, 0.379770517,  0.00000000, 0.254547089, 0.379554063, 0.254990399,  0.00000000 &
      ]
    v_mc_225_40(10,:) = [ real(irealLUT) :: &
      0.233765468, 0.741617262, 0.741670489, 0.383070976,  0.00000000, 0.258329034, 0.383162886, 0.258382291,  0.00000000 &
      ]
    v_mc_225_40(11,:) = [ real(irealLUT) :: &
      0.227727160, 0.738582730, 0.738496304, 0.386111259,  0.00000000, 0.261503220, 0.386160910, 0.261416823,  0.00000000 &
      ]
    v_mc_225_40(12,:) = [ real(irealLUT) :: &
      0.222842917, 0.736097336, 0.736135900, 0.388618767, 0.00000000, 0.263863593, 0.388537675, 0.263902217, 0.00000000 &
      ]
    v_mc_225_40(13,:) =  [ real(irealLUT) :: &
      0.218528301, 0.733748436, 0.733729243, 0.390776843, 0.00000000, 0.266270280, 0.390694201, 0.266251087, 0.00000000 &
      ]
    v_mc_225_40(14,:) = [ real(irealLUT) :: &
      0.214885697, 0.731690168, 0.731787622, 0.392676234, 0.00000000, 0.268211931, 0.392437398, 0.268309325, 0.00000000 &
      ]
    v_mc_225_40(15,:) = [ real(irealLUT) :: &
      0.211768627, 0.730218410, 0.730048716, 0.394240022, 0.00000000, 0.269950777, 0.393990666, 0.269781083, 0.00000000 &
      ]
    v_mc_225_40(16,:) = [ real(irealLUT) :: &
      0.208832026, 0.728453875, 0.728371441, 0.395480126, 0.00000000, 0.271628052, 0.395687193, 0.271545678, 0.00000000 &
      ]
    v_mc_225_40(17,:) = [ real(irealLUT) :: &
      0.206350222, 0.727174580, 0.727063775, 0.396706939, 0.00000000, 0.272935718, 0.396942168, 0.272824913, 0.00000000 &
      ]
    v_mc_225_40(18,:) = [ real(irealLUT) :: &
      0.204160988, 0.725816369, 0.725947320, 0.397873521,  0.00000000, 0.274052173, 0.397964835, 0.274183154, 0.00000000 &
      ]
    v_mc_225_40(19,:) = [real(irealLUT) :: &
      0.202239454, 0.724798799, 0.724844992, 0.398938626,  0.00000000, 0.275154531, 0.398821265, 0.275200695, 0.00000000 &
      ]
    v_mc_225_40(20,:) = [ real(irealLUT) :: &
      0.200536996, 0.723780394, 0.723834753, 0.399654120,  0.00000000, 0.276164711, 0.399808198, 0.276219100, 0.00000000 &
      ]

    do i = 1, 20
      verts_dtd = verts
      verts_dtd( 3) = verts( 3) + 2 * dz / i
      verts_dtd( 6) = verts( 6) + dz / i
      verts_dtd( 9) = verts( 9) + dz / (2 * i)
      verts_dtd(12) = verts( 6) + verts_dtd(12) - verts_dtd( 9)
      verts_dtd(15:24:3) = verts_dtd(3:12:3) + dz

      verts_dtd_symmetry = verts_dtd
      verts_dtd_symmetry( 3) = verts_dtd(12)
      verts_dtd_symmetry(15) = verts_dtd(24)
      verts_dtd_symmetry( 6) = verts_dtd( 9)
      verts_dtd_symmetry(18) = verts_dtd(21)
      verts_dtd_symmetry( 9) = verts_dtd( 6)
      verts_dtd_symmetry(21) = verts_dtd(18)
      verts_dtd_symmetry(12) = verts_dtd( 3)
      verts_dtd_symmetry(24) = verts_dtd(15)

      if (ldebug) then
        print *, cstr('verts distorted z', 'red')
        print *, verts_dtd(3:24:3)
        print *, cstr('verts symmetry z', 'red')
        print *, verts_dtd_symmetry(3:24:3)
      endif


      call dir2dir3_geometric_coeffs(verts_dtd, sundir, bg, v)
      call dir2dir3_geometric_coeffs(verts_dtd_symmetry, sundir_symmetry, bg, v_symmetry)

      if (ldebug) then
        print *, cstr('dtd', 'blue')
        print *, 'src z', v(1:9:3)
        print *, 'src x', v(2:9:3)
        print *, 'src y', v(3:9:3)
      endif

      !do src = 1,3
      !  call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts_dtd, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      !  v_mc(src:9:3) = real(T, irealLUT)
      !enddo
      v_mc = v_mc_225_40(i,:)

      if (ldebug) then
        print *, cstr('montecarlo distorted', 'green')
        !print *, v_mc
        print *, 'src z', v_mc(1:9:3)
        print *, 'src x', v_mc(2:9:3)
        print *, 'src y', v_mc(3:9:3)

        print *, cstr('dtd_symmetry', 'yellow')
        !print *, v_symmetry
        print *, 'src z', v_symmetry(1:9:3)
        print *, 'src x', v_symmetry(2:9:3)
        print *, 'src y', v_symmetry(3:9:3)
      endif

      @assertEqual(v_mc, v, max(maxval(v_mc)*0.05_irealLUT, 1e-6_irealLUT), 'failed for i='//toStr(i)//'; phi='//toStr(phi)//'; theta='//toStr(theta))
      @assertEqual(v_mc, v_symmetry, max(maxval(v_mc)*0.05_irealLUT, 1e-6_irealLUT), 'symmetry case failed for i='//toStr(i)//'; phi='//toStr(phi-180._ireals)//'; theta='//toStr(theta))

      if (ldebug) then
        print *, '___________________________________________________________________________________________'
        print *, '___________________________________________________________________________________________'
      endif
    enddo

  end subroutine

  @test(npes = [1])
  subroutine test_geometric_coeffs_regular_box_sundir_up_down(this)
    class(MpiTestMethod), intent(inout) :: this
    real(ireals), allocatable :: verts(:)
    real(irealLUT) :: v(9), v_ref(9), v_mc(9)
    real(ireals) :: sundir(3), verts_dtd(24)
    real(ireals), parameter :: dx=1, dy=dx, dz=dx
    integer(iintegers) :: src

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)

    bg = [real(ireal_dp) :: 1e-15, 0, 0.85]

    verts_dtd = verts
    verts_dtd([3,6,15,18]) = verts_dtd([3,6,15,18]) + dz

    if (ldebug) then
      print *, cstr('verts', 'red')
      print *, 'a', verts_dtd( 1: 3)
      print *, 'b', verts_dtd( 4: 6)
      print *, 'c', verts_dtd( 7: 9)
      print *, 'd', verts_dtd(10:12)
      print *, 'e', verts_dtd(13:15)
      print *, 'f', verts_dtd(16:18)
      print *, 'g', verts_dtd(19:21)
      print *, 'h', verts_dtd(22:24)
    endif

    if (ldebug) print *, cstr('CASE 1', 'red')

    phi = 180
    theta = 90

    sundir = spherical_2_cartesian(phi, theta) * [-one, -one, one]

    if (ldebug) then
      print *, cstr('sundir', 'yellow')
      print *, sundir
    endif

    call dir2dir3_geometric_coeffs(verts_dtd, sundir, bg, v)
    v_ref = [real(irealLUT) :: 0, 0.5, 1, 0, 0, 0, 1, 0.5, 0]

    if (ldebug) then
      print *, cstr('v_gomtrc', 'blue')
      !print *, v
      print *, 'src z', v(1:9:3)
      print *, 'src x', v(2:9:3)
      print *, 'src y', v(3:9:3)
      print *, cstr('v_ref', 'green')
      !print *, v_ref
      print *, 'src z', v_ref(1:9:3)
      print *, 'src x', v_ref(2:9:3)
      print *, 'src y', v_ref(3:9:3)
    endif

    @assertEqual(v_ref, v, max(maxval(v_ref)*0.05_irealLUT, 1e-6_irealLUT), 'failed for case 1'//'; phi='//toStr(phi)//'; theta='//toStr(theta))

    if (ldebug) print *, cstr('CASE 2', 'red')

    phi = 0
    theta = 90

    sundir = spherical_2_cartesian(phi, theta) * [-one, -one, one]

    if (ldebug) then
      print *, cstr('sundir', 'yellow')
      print *, sundir
    endif

    call dir2dir3_geometric_coeffs(verts_dtd, sundir, bg, v)
    v_ref = [real(irealLUT) :: 0, 0, 0, 0, 0.5, 0, 1, 0.5, 1]

    if (ldebug) then
      print *, cstr('v_gomtrc', 'blue')
      !print *, v
      print *, 'src z', v(1:9:3)
      print *, 'src x', v(2:9:3)
      print *, 'src y', v(3:9:3)
      print *, cstr('v_ref', 'green')
      !print *, v_ref
      print *, 'src z', v_ref(1:9:3)
      print *, 'src x', v_ref(2:9:3)
      print *, 'src y', v_ref(3:9:3)
    endif

    @assertEqual(v_ref, v, max(maxval(v_ref)*0.05_irealLUT, 1e-6_irealLUT), 'failed for case 2'//'; phi='//toStr(phi)//'; theta='//toStr(theta))

    ! CASE 3

    phi = 45
    theta = 90

    sundir = spherical_2_cartesian(phi, theta) * [-one, -one, one]

    if (ldebug) then
      print *, cstr('sundir', 'yellow')
      print *, sundir
    endif

    !do src = 1,3
    !  call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts_dtd, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    !  v_mc(src:9:3) = real(T, irealLUT)
    !enddo

    call dir2dir3_geometric_coeffs(verts_dtd, sundir, bg, v)
    v_ref = [real(irealLUT) :: 0.5, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5]

    if (ldebug) then
    !  print *, cstr('v_mc', 'yellow')
    !  print *, 'src z', v_mc(1:9:3)
    !  print *, 'src x', v_mc(2:9:3)
    !  print *, 'src y', v_mc(3:9:3)
      print *, cstr('v_gomtrc', 'blue')
      !print *, v
      print *, 'src z', v(1:9:3)
      print *, 'src x', v(2:9:3)
      print *, 'src y', v(3:9:3)
      print *, cstr('v_ref', 'green')
      !print *, v_ref
      print *, 'src z', v_ref(1:9:3)
      print *, 'src x', v_ref(2:9:3)
      print *, 'src y', v_ref(3:9:3)
    endif

    @assertEqual(v_ref, v, max(maxval(v_ref)*0.05_irealLUT, 1e-6_irealLUT), 'failed for case 1'//'; phi='//toStr(phi)//'; theta='//toStr(theta))

  end subroutine

  !@test(npes = [1])
  subroutine test_geometric_coeffs_distorted_box_no_abso(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src, i
    real(ireals), allocatable :: verts(:)
    real(irealLUT) :: v(9), v_mc(9), s_mc_dst(30), s_mc_reg(30), s_mc_reg_225_40(3,30), s_mc_dst_225_40(3,30)
    real(ireals) :: sundir(3), verts_dtd(24)
    real(ireals), parameter :: dx=1, dy=dx, dz=dx

    !bg = [real(ireal_dp) :: 1e-6, 1, 0.85]
    bg = [real(ireal_dp) :: 0, 1, 0.85]

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)

    S_target = zero

    phi = 225._ireals
    theta = 40._ireals

    sundir = spherical_2_cartesian(phi, theta) * [-one, -one, one]

    s_mc_dst_225_40(1,:) = [ real(irealLUT) :: &
      8.92252196e-03 , 1.74059323e-03 , 9.69395507e-04 , 0.212046579    , 0.170165360    , 0.170212224      , &
      4.39815409e-02 , 1.60874843e-04 , 1.65660270e-02 , 4.61620418e-03 , 9.30564199e-03 , 1.99979171e-03   , &
      4.86775208e-03 , 6.17795682e-04 , 1.74403004e-03 , 1.04445626e-03 , 3.07400990e-03 , 1.04887481e-03   , &
      4.10553925e-02 , 1.59284249e-02 , 6.19455614e-06 , 4.88061132e-03 , 2.86520529e-03 , 1.00863120e-02   , &
      4.42215893e-03 , 1.50971371e-03 , 2.58491491e-04 , 6.93272275e-04 , 8.93732358e-04 , 3.24867712e-03 &
      ]
    s_mc_dst_225_40(2,:) =  [ real(irealLUT) :: &
      8.35524127e-03 , 1.96852721e-03 , 1.51717372e-03 , 0.197576657    , 0.205781668    , 0.206131473      , &
      7.62645602e-02 , 6.24178385e-04 , 3.25699784e-02 , 6.28670771e-03 , 1.15256449e-02 , 2.78950064e-03   , &
      7.89808948e-03 , 1.86777406e-03 , 4.32133488e-03 , 1.51318708e-03 , 4.14777594e-03 , 1.68501784e-03   , &
      7.30539635e-02 , 3.15257721e-02 , 1.33042151e-04 , 6.68100920e-03 , 3.62757617e-03 , 1.22405626e-02   , &
      7.95461889e-03 , 4.10446431e-03 , 1.36837480e-03 , 1.16909470e-03 , 1.53356825e-03 , 4.22939751e-03 &
      ]
    s_mc_dst_225_40(3,:) = [ real(irealLUT) :: &
      8.53906758e-03 , 2.34087859e-03 , 1.96651346e-03 , 0.178379506    , 0.217701852    , 0.217965201       , &
      9.45552438e-02 , 1.33445719e-03 , 4.33570817e-02 , 6.79989019e-03 , 1.25034377e-02 , 3.11340019e-03    , &
      8.91354121e-03 , 2.64713424e-03 , 5.88732492e-03 , 1.73592858e-03 , 4.62617166e-03 , 1.98503397e-03    , &
      9.18332189e-02 , 4.25207727e-02 , 6.01716514e-04 , 7.05446163e-03 , 3.88485100e-03 , 1.29956184e-02    , &
      9.22008976e-03 , 5.74483350e-03 , 2.30736798e-03 , 1.47456466e-03 , 1.86840037e-03 , 4.68546804e-03  &
      ]

    s_mc_reg_225_40(1,:) = [ real(irealLUT) :: &
      1.17830206e-02 , 4.29459894e-03 , 4.22810391e-03 , 0.106779635    , 0.227237672    , 0.227556005      , &
      0.145700470    , 7.57512357e-03 , 9.04987380e-02 , 7.17611006e-03 , 1.45932129e-02 , 4.01723990e-03   , &
      9.34450421e-03 , 4.48263204e-03 , 9.67616774e-03 , 2.26292550e-03 , 5.86380530e-03 , 2.76703737e-03   , &
      0.145690143    , 9.06547606e-02 , 7.63494754e-03 , 7.25133158e-03 , 4.05412307e-03 , 1.46636609e-02   , &
      9.31688491e-03 , 9.71033145e-03 , 4.48918715e-03 , 2.24047154e-03 , 2.70990422e-03 , 5.80736576e-03 &
      ]
    s_mc_reg_225_40(2,:) = [ real(irealLUT) :: &
      1.17989006e-02 , 4.12377995e-03 , 4.21515619e-03 , 0.106731184    , 0.227561533    , 0.227508828      , &
      0.145324066    , 7.72332028e-03 , 9.07657593e-02 , 7.27583608e-03 , 1.46190245e-02 , 4.06450778e-03   , &
      9.34136380e-03 , 4.48409142e-03 , 9.73774213e-03 , 2.23024306e-03 , 5.85877663e-03 , 2.71193800e-03   , &
      0.145519182    , 9.05537307e-02 , 7.63165951e-03 , 7.19224941e-03 , 4.08802368e-03 , 1.47001557e-02   , &
      9.42092761e-03 , 9.77022015e-03 , 4.41311952e-03 , 2.23001977e-03 , 2.74304394e-03 , 5.83068002e-03 &
      ]
    s_mc_reg_225_40(3,:) = [ real(irealLUT) :: &
      1.17882425e-02 , 4.26216517e-03 , 4.22528014e-03 , 0.106809281    , 0.227495804    , 0.227499247      , &
      0.145535097    , 7.62372883e-03 , 9.07295942e-02 , 7.18160812e-03 , 1.46157164e-02 , 4.07528132e-03   , &
      9.38039087e-03 , 4.47310507e-03 , 9.72650107e-03 , 2.23162957e-03 , 5.84326964e-03 , 2.71494873e-03   , &
      0.145678475    , 9.08066705e-02 , 7.64530385e-03 , 7.23642809e-03 , 4.15736903e-03 , 1.45084914e-02   , &
      9.31425206e-03 , 9.70079750e-03 , 4.50785458e-03 , 2.21454003e-03 , 2.68698810e-03 , 5.86110959e-03 &
      ]

    do i = 1, 3
      verts_dtd = verts
      verts_dtd( 3) = verts( 3) + 2 * dz / i
      verts_dtd( 6) = verts( 6) + dz / i
      verts_dtd( 9) = verts( 9) + dz / (2 * i)
      verts_dtd(12) = verts( 6) + verts_dtd(12) - verts_dtd( 9)

      verts_dtd(15:24:3) = verts_dtd(3:12:3) + dz

      call dir2dir3_geometric_coeffs(verts_dtd, sundir, bg, v)

      !print *, cstr('regular corrected', 'blue')
      !print *, v
      !print *, 'src z', v(1:9:3)
      !print *, 'src x', v(2:9:3)
      !print *, 'src y', v(3:9:3)

      !do src = 1,3
      !  call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      !  s_mc_reg(src:30:3) = real(S, irealLUT)
      !enddo

      !do src = 1,3
      !  call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts_dtd, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      !  v_mc(src:9:3) = real(T, irealLUT)
      !  s_mc_dst(src:30:3) = real(S, irealLUT)
      !enddo

      s_mc_dst = s_mc_dst_225_40(i,:)
      s_mc_reg = s_mc_reg_225_40(i,:)

      print *, cstr('montecarlo distorted', 'green')
      print *, v_mc
      !print *, 'src z', v_mc(1:9:3)
      !print *, 'src x', v_mc(2:9:3)
      !print *, 'src y', v_mc(3:9:3)

      print *, cstr('montecarlo dir2diff', 'red')
      print *, cstr('dst', 'blue')
      print *, 'src z', s_mc_dst(1:30:3), ' : ', sum(s_mc_dst(1:30:3))
      print *, 'src x', s_mc_dst(2:30:3), ' : ', sum(s_mc_dst(2:30:3))
      print *, 'src y', s_mc_dst(3:30:3), ' : ', sum(s_mc_dst(3:30:3))
      print *, cstr('reg', 'green')
      print *, 'src z', s_mc_reg(1:30:3), ' : ', sum(s_mc_reg(1:30:3))
      print *, 'src x', s_mc_reg(2:30:3), ' : ', sum(s_mc_reg(2:30:3))
      print *, 'src y', s_mc_reg(3:30:3), ' : ', sum(s_mc_reg(3:30:3))
      print *, cstr('dst / reg', 'yellow')
      print *, 'src z', s_mc_dst(1:30:3) / s_mc_reg(1:30:3)
      print *, 'src x', s_mc_dst(2:30:3) / s_mc_reg(2:30:3)
      print *, 'src y', s_mc_dst(3:30:3) / s_mc_reg(3:30:3)

      !@assertEqual(v_mc, v, max(maxval(v_mc)*0.05_irealLUT, 1e-6_irealLUT), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
    enddo

  end subroutine

  !@test(npes =[1])
  subroutine test_geometric_coeffs_distorted_box_symmetry_no_abso(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src, i
    real(ireals), allocatable :: verts(:)
    real(irealLUT) :: v(9), s_mc_reg(30), s_mc_dst(30), s_mc_dst_225_40(3,30), s_mc_reg_225_40(3,30)
    real(ireals) :: sundir(3), verts_dtd(24)
    real(ireals), parameter :: dx=1, dy=dx, dz=dx

    s_mc_dst_225_40(1,:) = [ real(irealLUT) :: &
      8.92252196e-03 , 1.74059323e-03 , 9.69395507e-04 , 0.212046579    , 0.170165360    , 0.170212224      , &
      4.39815409e-02 , 1.60874843e-04 , 1.65660270e-02 , 4.61620418e-03 , 9.30564199e-03 , 1.99979171e-03   , &
      4.86775208e-03 , 6.17795682e-04 , 1.74403004e-03 , 1.04445626e-03 , 3.07400990e-03 , 1.04887481e-03   , &
      4.10553925e-02 , 1.59284249e-02 , 6.19455614e-06 , 4.88061132e-03 , 2.86520529e-03 , 1.00863120e-02   , &
      4.42215893e-03 , 1.50971371e-03 , 2.58491491e-04 , 6.93272275e-04 , 8.93732358e-04 , 3.24867712e-03 &
      ]
    s_mc_dst_225_40(2,:) =  [ real(irealLUT) :: &
      8.35524127e-03 , 1.96852721e-03 , 1.51717372e-03 , 0.197576657    , 0.205781668    , 0.206131473      , &
      7.62645602e-02 , 6.24178385e-04 , 3.25699784e-02 , 6.28670771e-03 , 1.15256449e-02 , 2.78950064e-03   , &
      7.89808948e-03 , 1.86777406e-03 , 4.32133488e-03 , 1.51318708e-03 , 4.14777594e-03 , 1.68501784e-03   , &
      7.30539635e-02 , 3.15257721e-02 , 1.33042151e-04 , 6.68100920e-03 , 3.62757617e-03 , 1.22405626e-02   , &
      7.95461889e-03 , 4.10446431e-03 , 1.36837480e-03 , 1.16909470e-03 , 1.53356825e-03 , 4.22939751e-03 &
      ]
    s_mc_dst_225_40(3,:) = [ real(irealLUT) :: &
      8.53906758e-03 , 2.34087859e-03 , 1.96651346e-03 , 0.178379506    , 0.217701852    , 0.217965201       , &
      9.45552438e-02 , 1.33445719e-03 , 4.33570817e-02 , 6.79989019e-03 , 1.25034377e-02 , 3.11340019e-03    , &
      8.91354121e-03 , 2.64713424e-03 , 5.88732492e-03 , 1.73592858e-03 , 4.62617166e-03 , 1.98503397e-03    , &
      9.18332189e-02 , 4.25207727e-02 , 6.01716514e-04 , 7.05446163e-03 , 3.88485100e-03 , 1.29956184e-02    , &
      9.22008976e-03 , 5.74483350e-03 , 2.30736798e-03 , 1.47456466e-03 , 1.86840037e-03 , 4.68546804e-03  &
      ]

    s_mc_reg_225_40(1,:) = [ real(irealLUT) :: &
      1.17830206e-02 , 4.29459894e-03 , 4.22810391e-03 , 0.106779635    , 0.227237672    , 0.227556005      , &
      0.145700470    , 7.57512357e-03 , 9.04987380e-02 , 7.17611006e-03 , 1.45932129e-02 , 4.01723990e-03   , &
      9.34450421e-03 , 4.48263204e-03 , 9.67616774e-03 , 2.26292550e-03 , 5.86380530e-03 , 2.76703737e-03   , &
      0.145690143    , 9.06547606e-02 , 7.63494754e-03 , 7.25133158e-03 , 4.05412307e-03 , 1.46636609e-02   , &
      9.31688491e-03 , 9.71033145e-03 , 4.48918715e-03 , 2.24047154e-03 , 2.70990422e-03 , 5.80736576e-03 &
      ]
    s_mc_reg_225_40(2,:) = [ real(irealLUT) :: &
      1.17989006e-02 , 4.12377995e-03 , 4.21515619e-03 , 0.106731184    , 0.227561533    , 0.227508828      , &
      0.145324066    , 7.72332028e-03 , 9.07657593e-02 , 7.27583608e-03 , 1.46190245e-02 , 4.06450778e-03   , &
      9.34136380e-03 , 4.48409142e-03 , 9.73774213e-03 , 2.23024306e-03 , 5.85877663e-03 , 2.71193800e-03   , &
      0.145519182    , 9.05537307e-02 , 7.63165951e-03 , 7.19224941e-03 , 4.08802368e-03 , 1.47001557e-02   , &
      9.42092761e-03 , 9.77022015e-03 , 4.41311952e-03 , 2.23001977e-03 , 2.74304394e-03 , 5.83068002e-03 &
      ]
    s_mc_reg_225_40(3,:) = [ real(irealLUT) :: &
      1.17882425e-02 , 4.26216517e-03 , 4.22528014e-03 , 0.106809281    , 0.227495804    , 0.227499247      , &
      0.145535097    , 7.62372883e-03 , 9.07295942e-02 , 7.18160812e-03 , 1.46157164e-02 , 4.07528132e-03   , &
      9.38039087e-03 , 4.47310507e-03 , 9.72650107e-03 , 2.23162957e-03 , 5.84326964e-03 , 2.71494873e-03   , &
      0.145678475    , 9.08066705e-02 , 7.64530385e-03 , 7.23642809e-03 , 4.15736903e-03 , 1.45084914e-02   , &
      9.31425206e-03 , 9.70079750e-03 , 4.50785458e-03 , 2.21454003e-03 , 2.68698810e-03 , 5.86110959e-03 &
      ]

    bg = [real(ireal_dp) :: 0, 1, 0.85]

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)

    S_target = zero

    phi = 45._ireals ! 225 -> 45
    theta = 40._ireals

    sundir = spherical_2_cartesian(phi, theta) * [-one, -one, one]
    print*, 'sundir', sundir

    do i = 1, 20
      verts_dtd = verts
      ! 6 <-> 9 && 3 <-> 12
      verts_dtd(12) = verts(12) + 2 * dz / i
      verts_dtd( 9) = verts( 9) + dz / i
      verts_dtd( 6) = verts( 6) + dz / (2 * i)
      verts_dtd( 3) = verts( 6) + verts_dtd( 9) - verts_dtd(12)

      verts_dtd(15:24:3) = verts_dtd(3:12:3) + dz

      call dir2dir3_geometric_coeffs(verts_dtd, sundir, bg, v)

      print *, cstr('regular corrected', 'blue')
      print *, v
      !print *, 'src z', v(1:9:3)
      !print *, 'src x', v(2:9:3)
      !print *, 'src y', v(3:9:3)

      !do src = 1,3
      !  call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      !  s_mc_reg(src:30:3) = real(T, irealLUT)
      !enddo

      !do src = 1,3
      !  call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts_dtd, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
      !  s_mc_dst(src:30:3) = real(T, irealLUT)
      !enddo

      s_mc_reg = s_mc_reg_225_40(i,:)
      s_mc_dst = s_mc_dst_225_40(i,:)

      !@assertEqual(v_mc, v, max(maxval(v_mc)*0.05_irealLUT, 1e-6_irealLUT), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
    enddo
  end subroutine

  !@test(npes =[1])
  subroutine test_geometric_coeffs_distorted_box_new(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src
    real(ireals), allocatable :: verts(:)
    real(irealLUT) :: v(9), v_mc(9), v_mc_dir2diff_reg(30), v_mc_dir2diff_dst(30)
    real(ireals) :: sundir(3), verts_dtd(24)
    real( ireals), parameter :: dx=1, dy=dx, dz=dx

    !bg = [real(ireal_dp) :: 1e-6, 1, 0.85]
    bg = [real(ireal_dp) :: 1e-6, 0, 0.85]

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)

    verts([9,12,21,24]) = verts([9,12,21,24]) + dz / 2

    S_target = zero

    phi = 225._ireals
    theta = 40._ireals

    sundir = spherical_2_cartesian(phi, theta)

    call dir2dir3_geometric_coeffs(verts_dtd, sundir, bg, v)

    print *, cstr('regular corrected', 'blue')
    print *, 'src z', v(1:9:3)
    print *, 'src x', v(2:9:3)
    print *, 'src y', v(3:9:3)

        do src = 1,3
          call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
          v_mc_dir2diff_reg(src:30:3) = real(S, irealLUT)
        enddo

        do src = 1,3
          call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts_dtd, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
          v_mc(src:3**2:3) = real(T, irealLUT)
          v_mc_dir2diff_dst(src:30:3) = real(S, irealLUT)
        enddo

        print *, cstr('montecarlo distorted', 'green')
        print *, 'src z', v_mc(1:9:3)
        print *, 'src x', v_mc(2:9:3)
        print *, 'src y', v_mc(3:9:3)

        print *, cstr('dir2diff', 'red')
        print *, 'srz z', v_mc_dir2diff_reg(1:30:3)
        print *, 'srz x', v_mc_dir2diff_reg(2:30:3)
        print *, 'srz y', v_mc_dir2diff_reg(3:30:3)

        print *, cstr('dir2diff', 'red')
        print *, 'srz z', v_mc_dir2diff_dst(1:30:3)
        print *, 'srz x', v_mc_dir2diff_dst(2:30:3)
        print *, 'srz y', v_mc_dir2diff_dst(3:30:3)

        @assertEqual(v_mc, v, max(maxval(v_mc)*0.05_irealLUT, 1e-6_irealLUT), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
!        @assertEqual(v_mc([2,3,4]), v([2,3,4]), max(maxval(v_mc([2,3,4]))*0.05_irealLUT, 1e-6_irealLUT), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
    !  enddo
    !enddo
    end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check(S_target,T_target, S,T, msg)
    real(ireals),intent(in),dimension(:) :: S_target,T_target, S,T

    character(len=*),optional :: msg
    character(default_str_len) :: local_msgS, local_msgT

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

      print*,'---------------------'
      write(*, FMT='( " diffuse ::  :: ",10(es12.5) )' ) S
      write(*, FMT='( " target  ::  :: ",10(es12.5) )' ) S_target
      write(*, FMT='( " diff    ::  :: ",10(es12.5) )' ) S_target-S
      print*,''
      write(*, FMT='( " direct  ::  :: ", 8(es12.5) )' ) T
      write(*, FMT='( " target  ::  :: ", 8(es12.5) )' ) T_target
      write(*, FMT='( " diff    ::  :: ", 8(es12.5) )' ) T_target-T
      print*,'---------------------'
      print*,''

      @assertEqual(S_target, S, real(atol*sigma, ireals), local_msgS )
      @assertLessThanOrEqual   (zero, S)
      @assertGreaterThanOrEqual(one , S)

      @assertEqual(T_target, T, real(atol*sigma, ireals), local_msgT )
      @assertLessThanOrEqual   (zero, T)
      @assertGreaterThanOrEqual(one , T)
    endif
  end subroutine

end module
