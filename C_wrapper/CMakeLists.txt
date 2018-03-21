include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(ENABLE_RRTM)

  add_library( f2c_tenstream f2c_tenstream.F90 )
  target_link_libraries( f2c_tenstream pprts pprts_base tenstr_rrtmg tenstr_rrtm_lw tenstr_rrtm_sw tenstream )
  set_target_properties( f2c_tenstream PROPERTIES PUBLIC_HEADER "f2c_solver_ids.h;f2c_tenstream.h")

  install(TARGETS f2c_tenstream
    DESTINATION lib/
    LIBRARY DESTINATION lib/
    PUBLIC_HEADER DESTINATION include/)

  file(COPY ./
    DESTINATION include/
    FILES_MATCHING PATTERN "*.h" )

  add_executable( "Cwrapper_ex1" "ex1.c" )
  target_link_libraries( "Cwrapper_ex1" f2c_tenstream tenstr_rrtmg tenstr_rrtm_lw tenstr_rrtm_sw tenstream )

  add_executable( "Cwrapper_ex2" "ex2.c" )
  target_link_libraries( "Cwrapper_ex2" f2c_tenstream tenstr_rrtmg tenstr_rrtm_lw tenstr_rrtm_sw tenstream )
endif()
