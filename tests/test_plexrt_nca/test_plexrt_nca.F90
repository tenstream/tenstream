module test_plexrt_nca
    use m_data_parameters, only: iintegers, ireals, default_str_len
    use m_helper_functions, only: triangle_area_by_edgelengths
    use m_tenstream_interpolation, only: interp_1d, interp_2d, &
      interp_vec_simplex_nd

    use m_plexrt_nca, only: plexrt_nca_init, plexrt_nca

    use pfunit_mod


    implicit none

contains
  @test(npes =[1])
  subroutine test_single_box_nca(this)
    class (MpiTestMethod), intent(inout) :: this
    real(ireals), parameter :: dx1=1e3, dx2=dx1, dx3=dx1, dz=1e3
    real(ireals) :: atop, abot, a1, a2, a3, vol, hr
    real(ireals) :: base_info(7), side_info(3*5)
    real(ireals), parameter :: kabs=1e-5

    atop = triangle_area_by_edgelengths(dx1, dx2, dx3)
    abot = atop
    a1   = dx1 * dz
    a2   = dx2 * dz
    a3   = dx3 * dz

    vol = atop * dz

    base_info = [        &
            kabs,        & !kabs    
            kabs,        & !kabs_top
            0._ireals,   & !Ldn_top 
            1._ireals,   & !Btop    
            kabs,        & !kabs_bot
            0._ireals,   & !Lup_bot 
            1._ireals    & !Bbot    
            ]
    side_info(1:5) = [   &
            kabs,        & !kabs
            0._ireals,   & !Ldn_top
            0._ireals,   & !Lup_top
            0._ireals,   & !Ldn_bot
            0._ireals    & !Lup_bot
            ]

    side_info( 6:10) = side_info(1:5)
    side_info(11:15) = side_info(1:5)

    call plexrt_nca_init()

    call plexrt_nca (dx1, dx2, dx3, &
            dz, atop, abot, a1, a2, a3, vol, &
            base_info, side_info, hr)

    print *,'hr:', hr
    @assertEqual(-9.9e-5, hr, epsilon(hr)*10, 'NCA Heating Test 1 false')

  end subroutine

end module
