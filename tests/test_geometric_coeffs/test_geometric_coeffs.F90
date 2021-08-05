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

  real(ireals) :: phi,theta,dx,dy,dz,bg(3)
  real(ireals) :: S(10),T(3), S_target(10), T_target(3)
  real(ireals) :: S_tol(10),T_tol(3)
  real(ireals), allocatable :: vertices(:)

  type(t_boxmc_3_10) :: bmc_3_10

  integer(mpiint) :: myid,mpierr,numnodes,comm
  character(len=120) :: msg

  real(ireals),parameter :: sigma = 3 ! normal test range for coefficients

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

  @test(npes=[1])
  subroutine test_geometric_coeffs_lambert_beer(this)
    class (MpitestMethod), intent(inout) :: this
    real(ireals), allocatable :: verts(:)
    real(ireals) :: v(9), v_trgt(9), sundir(3), c_ext
    real(ireals), parameter :: dx=1, dy=dx, dz=dx

    bg = [real(ireals) :: one, zero, 0.85]
    c_ext = bg(1) + bg(2)
    call setup_default_unit_cube_geometry(dx, dy, dz, verts)
    phi = 180
    theta = 0
    sundir = spherical_2_cartesian(phi, theta) * [-one, -one, one]

    if (ldebug) print *, 'sundir', sundir
    call dir2dir3_geometric_coeffs(verts, sundir, c_ext, v)
    v_trgt = [real(ireals) :: &
      exp(-c_ext*dz), (one-exp(-c_ext*dz))/(c_ext*dz), (one-exp(-c_ext*dz))/(c_ext*dz), &
      zero, zero, zero, &
      zero, zero, zero &
      ]

    if (ldebug) then
      print *, cstr('v_gomtrc', 'blue')
      print *, 'src z', v(1:9:3)
      print *, 'src x', v(2:9:3)
      print *, 'src y', v(3:9:3)
      print *, cstr('v_trgt', 'green')
      print *, cstr('v_gomtrc', 'blue')
      print *, 'src z', v_trgt(1:9:3)
      print *, 'src x', v_trgt(2:9:3)
      print *, 'src y', v_trgt(3:9:3)
    endif


    @assertEqual(v_trgt, v, max(maxval(v_trgt)*0.05_ireals, 1e-6_ireals), 'failed for phi='//toStr(phi)//'; theta='//toStr(theta))
  end subroutine


  @test(npes=[1])
  subroutine test_geometric_coeffs_distorted_box_no_scatter(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers) :: src, i
    real(ireals), allocatable :: verts(:)
    real(ireals) :: v(9), v_symmetry(9), v_mc(9), v_mc_225_40(20,9)
    real(ireals) :: sundir(3), sundir_symmetry(3), verts_dtd(24), verts_dtd_symmetry(24), c_ext
    real(ireals), parameter :: dx=1, dy=dx, dz=dx

    bg = [real(ireals) :: tiny(bg), 0, 0.85]
    c_ext = bg(1) + bg(2)

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)

    S_target = zero

    phi = 225._ireals
    theta = 40._ireals

    sundir = spherical_2_cartesian(phi, theta) * [-one, -one, one]
    sundir_symmetry = spherical_2_cartesian(phi - 180, theta) * [-one, -one, one]

    if (ldebug) then
      print *, 'sundir', sundir
      print *, 'sundir_symmetry', sundir_symmetry
    endif

    v_mc_225_40(1,:) = [ real(ireals) :: &
      0.579161942, 0.880491078, 0.880592525, 0.210434467, 0.00000000, 0.119407229, 0.210403204, 0.119508661, 0.00000000 &
      ]
    v_mc_225_40(2,:) = [ real(ireals) :: &
      0.434624285, 0.829745471, 0.829779685, 0.282772988, 0.00000000, 0.170220003, 0.282602191, 0.170254186, 0.00000000 &
      ]
    v_mc_225_40(3,:) = [ real(ireals) :: &
      0.363792837, 0.801579118, 0.801497698, 0.318097949, 0.00000000, 0.198501915, 0.318108648, 0.198420525, 0.00000000 &
      ]
    v_mc_225_40(4,:) = [ real(ireals) :: &
      0.321538657, 0.783636987, 0.783575475, 0.339426547,  0.00000000, 0.216424137, 0.339034200, 0.216362625, 0.00000000 &
      ]
    v_mc_225_40(5,:) = [ real(ireals) :: &
      0.294336140, 0.771124899, 0.771123469, 0.352795750,  0.00000000, 0.228876129, 0.352867484, 0.228874698,  0.00000000 &
      ]
    v_mc_225_40(6,:) = [ real(ireals) :: &
      0.274791181, 0.762104392, 0.761982083, 0.362527847,  0.00000000, 0.238017485, 0.362680346, 0.237895191,  0.00000000 &
      ]
    v_mc_225_40(7,:) = [ real(ireals) :: &
      0.260479659, 0.755124331, 0.755310655, 0.369733512,  0.00000000, 0.244688913, 0.369786203, 0.244875193,  0.00000000 &
      ]
    v_mc_225_40(8,:) = [ real(ireals) :: &
      0.249499917, 0.749804854, 0.749615192, 0.375236362,  0.00000000, 0.250384361, 0.375263065, 0.250194699,  0.00000000 &
      ]
    v_mc_225_40(9,:) = [ real(ireals) :: &
      0.240674764, 0.745009124, 0.745452464, 0.379770517,  0.00000000, 0.254547089, 0.379554063, 0.254990399,  0.00000000 &
      ]
    v_mc_225_40(10,:) = [ real(ireals) :: &
      0.233765468, 0.741617262, 0.741670489, 0.383070976,  0.00000000, 0.258329034, 0.383162886, 0.258382291,  0.00000000 &
      ]
    v_mc_225_40(11,:) = [ real(ireals) :: &
      0.227727160, 0.738582730, 0.738496304, 0.386111259,  0.00000000, 0.261503220, 0.386160910, 0.261416823,  0.00000000 &
      ]
    v_mc_225_40(12,:) = [ real(ireals) :: &
      0.222842917, 0.736097336, 0.736135900, 0.388618767, 0.00000000, 0.263863593, 0.388537675, 0.263902217, 0.00000000 &
      ]
    v_mc_225_40(13,:) =  [ real(ireals) :: &
      0.218528301, 0.733748436, 0.733729243, 0.390776843, 0.00000000, 0.266270280, 0.390694201, 0.266251087, 0.00000000 &
      ]
    v_mc_225_40(14,:) = [ real(ireals) :: &
      0.214885697, 0.731690168, 0.731787622, 0.392676234, 0.00000000, 0.268211931, 0.392437398, 0.268309325, 0.00000000 &
      ]
    v_mc_225_40(15,:) = [ real(ireals) :: &
      0.211768627, 0.730218410, 0.730048716, 0.394240022, 0.00000000, 0.269950777, 0.393990666, 0.269781083, 0.00000000 &
      ]
    v_mc_225_40(16,:) = [ real(ireals) :: &
      0.208832026, 0.728453875, 0.728371441, 0.395480126, 0.00000000, 0.271628052, 0.395687193, 0.271545678, 0.00000000 &
      ]
    v_mc_225_40(17,:) = [ real(ireals) :: &
      0.206350222, 0.727174580, 0.727063775, 0.396706939, 0.00000000, 0.272935718, 0.396942168, 0.272824913, 0.00000000 &
      ]
    v_mc_225_40(18,:) = [ real(ireals) :: &
      0.204160988, 0.725816369, 0.725947320, 0.397873521,  0.00000000, 0.274052173, 0.397964835, 0.274183154, 0.00000000 &
      ]
    v_mc_225_40(19,:) = [real(ireals) :: &
      0.202239454, 0.724798799, 0.724844992, 0.398938626,  0.00000000, 0.275154531, 0.398821265, 0.275200695, 0.00000000 &
      ]
    v_mc_225_40(20,:) = [ real(ireals) :: &
      0.200536996, 0.723780394, 0.723834753, 0.399654120,  0.00000000, 0.276164711, 0.399808198, 0.276219100, 0.00000000 &
      ]

    do i = 1, 20
      verts_dtd = verts
      verts_dtd( 3) = verts( 3) + 2 * dz / real(i, ireals)
      verts_dtd( 6) = verts( 6) + dz / real(i, ireals)
      verts_dtd( 9) = verts( 9) + dz / (2 * real(i, ireals))
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


      call dir2dir3_geometric_coeffs(verts_dtd, sundir, c_ext, v)
      call dir2dir3_geometric_coeffs(verts_dtd_symmetry, sundir_symmetry, c_ext, v_symmetry)

      if (ldebug) then
        print *, cstr('dtd', 'blue')
        print *, 'src z', v(1:9:3)
        print *, 'src x', v(2:9:3)
        print *, 'src y', v(3:9:3)
      endif

      v_mc = v_mc_225_40(i,:)
      do src = 1,3
        !call bmc_3_10%get_coeff( &
        !  & comm, &
        !  & real(bg,ireal_dp), &
        !  & src, &
        !  & .True., &
        !  & real(phi,ireal_dp), &
        !  & real(theta, ireal_dp), &
        !  & real(verts_dtd, ireal_dp), &
        !  & S, T, S_tol, T_tol, &
        !  & inp_atol=atol, inp_rtol=rtol)
        !v_mc(src:9:3) = real(T, ireals)
      enddo

      if (ldebug) then
        print *, cstr('montecarlo distorted', 'green')
        print *, v_mc
        print *, 'src z', v_mc(1:9:3)
        print *, 'src x', v_mc(2:9:3)
        print *, 'src y', v_mc(3:9:3)

        print *, cstr('dtd_symmetry', 'yellow')
        print *, v_symmetry
        print *, 'src z', v_symmetry(1:9:3)
        print *, 'src x', v_symmetry(2:9:3)
        print *, 'src y', v_symmetry(3:9:3)
      endif

      @assertEqual(v_mc, v, max(maxval(v_mc)*0.05_ireals, 1e-6_ireals), 'failed for i='//toStr(i)//'; phi='//toStr(phi)//'; theta='//toStr(theta))
      @assertEqual(v_mc, v_symmetry, max(maxval(v_mc)*0.05_ireals, 1e-6_ireals), 'symmetry case failed for i='//toStr(i)//'; phi='//toStr(phi-180._ireals)//'; theta='//toStr(theta))
    enddo

  end subroutine

  @test(npes = [1])
  subroutine test_geometric_coeffs_regular_box_sundir_up_down(this)
    class(MpiTestMethod), intent(inout) :: this
    real(ireals), allocatable :: verts(:)
    real(ireals) :: v(9), v_mc(9)
    real(ireals) :: sundir(3), verts_dtd(24), c_ext
    real(ireals), parameter :: dx=1, dy=dx, dz=dx
    !integer(iintegers) :: src

    call setup_default_unit_cube_geometry(dx, dy, dz, verts)

    bg = [real(ireals) :: tiny(bg), 0, 0.85]
    c_ext = bg(1) + bg(2)

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

    call dir2dir3_geometric_coeffs(verts_dtd, sundir, c_ext, v)
    v_mc = [real(ireals) :: 0, 0.5, 1, 0, 0, 0, 1, 0.5, 0]

    if (ldebug) then
      print *, cstr('v_gomtrc', 'blue')
      print *, 'src z', v(1:9:3)
      print *, 'src x', v(2:9:3)
      print *, 'src y', v(3:9:3)
      print *, cstr('v_mc', 'green')
      print *, 'src z', v_mc(1:9:3)
      print *, 'src x', v_mc(2:9:3)
      print *, 'src y', v_mc(3:9:3)
    endif

    @assertEqual(v_mc, v, max(maxval(v_mc)*0.05_ireals, 1e-6_ireals), 'failed for case 1'//'; phi='//toStr(phi)//'; theta='//toStr(theta))

    if (ldebug) print *, cstr('CASE 2', 'red')

    phi = 0
    theta = 90

    sundir = spherical_2_cartesian(phi, theta) * [-one, -one, one]

    if (ldebug) then
      print *, cstr('sundir', 'yellow')
      print *, sundir
    endif

    call dir2dir3_geometric_coeffs(verts_dtd, sundir, c_ext, v)
    v_mc = [real(ireals) :: 0, 0, 0, 0, 0.5, 0, 1, 0.5, 1]

    if (ldebug) then
      print *, cstr('v_gomtrc', 'blue')
      print *, 'src z', v(1:9:3)
      print *, 'src x', v(2:9:3)
      print *, 'src y', v(3:9:3)
      print *, cstr('v_mc', 'green')
      print *, 'src z', v_mc(1:9:3)
      print *, 'src x', v_mc(2:9:3)
      print *, 'src y', v_mc(3:9:3)
    endif

    @assertEqual(v_mc, v, max(maxval(v_mc)*0.05_ireals, 1e-6_ireals), 'failed for case 2'//'; phi='//toStr(phi)//'; theta='//toStr(theta))

    ! CASE 3

    !phi = 45
    !theta = 90

    !sundir = spherical_2_cartesian(phi, theta) * [-one, -one, one]

    !if (ldebug) then
    !  print *, cstr('sundir', 'yellow')
    !  print *, sundir
    !endif

    !v_mc = [real(ireals) :: 0.5, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5]
    !do src = 1,3
    !  call bmc_3_10%get_coeff(comm,bg,src,.True.,phi,theta,real(verts_dtd, ireal_dp),S,T,S_tol,T_tol, inp_atol=atol, inp_rtol=rtol)
    !  v_mc(src:9:3) = real(T, ireals)
    !enddo

    !call dir2dir3_geometric_coeffs(verts_dtd, sundir, c_ext, v)

    !if (ldebug) then
    !  print *, cstr('v_gomtrc', 'blue')
    !  print *, 'src z', v(1:9:3)
    !  print *, 'src x', v(2:9:3)
    !  print *, 'src y', v(3:9:3)
    !  print *, cstr('v_mc', 'green')
    !  print *, 'src z', v_mc(1:9:3)
    !  print *, 'src x', v_mc(2:9:3)
    !  print *, 'src y', v_mc(3:9:3)
    !endif

    !@assertEqual(v_mc, v, max(maxval(v_mc)*0.05_ireals, 1e-6_ireals), 'failed for case 3'//'; phi='//toStr(phi)//'; theta='//toStr(theta))
  end subroutine
end module
