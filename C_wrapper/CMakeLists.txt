include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(ENABLE_RRTM)

  add_library( f2c_pprts f2c_pprts.F90 )
  target_link_libraries( f2c_pprts pprts plexrt pprts_base )
  set_target_properties( f2c_pprts PROPERTIES PUBLIC_HEADER "f2c_solver_ids.h;f2c_pprts.h")

  add_library( f2c_pprts_rrtm f2c_pprts_rrtm.F90 )
  target_link_libraries( f2c_pprts_rrtm pprts_rrtmg pprts pprts_base )

  file(COPY ./
    DESTINATION ${CMAKE_BINARY_DIR}/include/
    FILES_MATCHING PATTERN "*.h" )

  add_executable( "Cwrapper_ex1" "ex1.c" )
  target_link_libraries( "Cwrapper_ex1" f2c_pprts_rrtm )

  add_executable( "Cwrapper_ex2" "ex2.c" )
  target_link_libraries( "Cwrapper_ex2" f2c_pprts )
endif()
