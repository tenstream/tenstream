include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(ENABLE_RRTM)

  add_library( f2c_pprts f2c_pprts.F90 )
  target_link_libraries( f2c_pprts pprts plexrt pprts_base ${EXTLIB})
  set_target_properties( f2c_pprts PROPERTIES PUBLIC_HEADER "f2c_solver_ids.h;f2c_pprts.h")

  add_library( f2c_pprts_rrtm f2c_pprts_rrtm.F90 )
  target_link_libraries( f2c_pprts_rrtm pprts_rrtmg pprts plexrt pprts_base ${EXTLIB})

  file(COPY ./
    DESTINATION ${CMAKE_BINARY_DIR}/include/
    FILES_MATCHING PATTERN "*.h" )

  add_executable( "cwrapper_opp" "opp.c" )
  target_link_libraries( "cwrapper_opp" ${TENSTR_SINGLE_LIB} ${EXTLIB})

  add_executable( "cwrapper_pprts" "pprts.c" )
  target_link_libraries( "cwrapper_pprts" ${TENSTR_SINGLE_LIB} ${EXTLIB})

  add_executable( "cwrapper_pprts_rrtmg" "pprts_rrtmg.c" )
  target_link_libraries( "cwrapper_pprts_rrtmg" ${TENSTR_SINGLE_LIB} ${EXTLIB})

endif()
