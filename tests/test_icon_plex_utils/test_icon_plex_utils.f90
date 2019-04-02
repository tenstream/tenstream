module test_icon_plex_utils
  use m_data_parameters, only : &
    mpiint, ireals, iintegers,  &
    one, zero, default_str_len, &
    i0, i1, i2, i3, i4,         &
    init_mpi_data_parameters

  use m_helper_functions, only : CHKERR, deg2rad, itoa

  use m_icon_plex_utils, only : date_to_julian_day, get_sun_vector

  use pfunit_mod
  implicit none

  integer(mpiint) :: comm, myid, numnodes, ierr

contains

  @before
  subroutine setup(this)
      class (MpiTestMethod), intent(inout) :: this
      comm     = this%getMpiCommunicator()
      numnodes = this%getNumProcesses()
      myid     = this%getProcessRank()

      call init_mpi_data_parameters(comm)
  end subroutine setup

  @after
  subroutine teardown(this)
      class (MpiTestMethod), intent(inout) :: this
      if(myid.eq.0) print *,'Finishing icon_plex_utils tests module'
  end subroutine teardown

  @test(npes =[1])
  subroutine test_gregorian_date_to_julian_day(this)
      class (MpiTestMethod), intent(inout) :: this

      @assertEqual(2451545.0_ireals, date_to_julian_day(2000_iintegers, i1, 1.5_ireals), 'Greenwich Noon')

      @assertEqual(2446113.75_ireals, date_to_julian_day(1985_iintegers, i2, 17.25_ireals))
      @assertEqual(0._ireals, date_to_julian_day(-4712_iintegers, i1, 1.5_ireals))
      @assertEqual(2458365.29749_ireals, date_to_julian_day(2018_iintegers, 9_iintegers, 3._ireals+real(19*3600+8*60+23, ireals)/86400._ireals), 1e-4_ireals)

  end subroutine

  @test(npes =[1])
  subroutine test_sundirection_computation(this)
      class (MpiTestMethod), intent(inout) :: this
      integer(iintegers) :: year, month
      real(ireals) :: day

      real(ireals) :: sundir(3)
      integer :: hour, iday
      real(ireals) :: max_z_component
      max_z_component = sin(deg2rad(23.439_ireals))

      year = 2018
      month = 1

      do hour = 6,18
        day = 5 + hour/24._ireals
        sundir = get_sun_vector(year, month, day)

        @assertGreaterThan(sundir(1), zero)
        @assertEqual(one, norm2(sundir), sqrt(epsilon(sundir)))
      enddo

      do hour = 1,10
        day = 5 + hour/24._ireals
        sundir = get_sun_vector(year, month, day)

        @assertGreaterThan(sundir(2), zero, sqrt(epsilon(sundir)))
        @assertEqual(one, norm2(sundir), sqrt(epsilon(sundir)))
      enddo

      do hour = 13,23
        day = 5 + hour/24._ireals
        sundir = get_sun_vector(year, month, day)

        @assertLessThan(sundir(2), zero, sqrt(epsilon(sundir)))
        @assertEqual(one, norm2(sundir), sqrt(epsilon(sundir)))
      enddo

      do iday = 1, 31
        do hour = 1, 24
          do month = 4, 8
            sundir = get_sun_vector(year, month, iday+(hour/24._ireals))
            @assertGreaterThan(sundir(3), zero)
          enddo
        enddo

        do month = 1, 2
          sundir = get_sun_vector(year, month, iday+(hour/24._ireals))
          @assertLessThan(sundir(3), zero)
        enddo
        do month = 10, 12
          sundir = get_sun_vector(year, month, iday+(hour/24._ireals))
          @assertLessThan(sundir(3), zero)
        enddo
      enddo

      do year=2000, 2010, 1
        do iday = 1, 31
          do month = 1, 12
            do hour = 1, 24
              sundir = get_sun_vector(year, month, iday+(hour/24._ireals))
              @assertLessThan(abs(sundir(3)), max_z_component, 'bad z_component '//itoa(year)//':'//itoa(month)//':'//itoa(iday))
            enddo
          enddo
        enddo
      enddo

      sundir = get_sun_vector(2018_iintegers, 9_iintegers, 3._ireals+real(19*3600+8*60+23, ireals)/86400._ireals)
      @assertEqual([-0.29155451713720204_ireals, -0.94795805644649656_ireals, 0.12795111080046895_ireals], sundir, sqrt(epsilon(sundir)))
  end subroutine
end module
