include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(ENABLE_RRTM)

  add_library( f2c_pprts_objs OBJECT f2c_pprts.F90 )
  add_dependencies( f2c_pprts_objs pprts pprts_base )
  if(ENABLE_PLEXRT)
    add_dependencies( f2c_pprts_objs plexrt )
  endif()

  add_library( f2c_pprts $<TARGET_OBJECTS:f2c_pprts_objs> )
  target_link_libraries( f2c_pprts pprts pprts_base ${EXTLIB})
  if(ENABLE_PLEXRT)
    target_link_libraries( f2c_pprts plexrt )
  endif()
  set_target_properties( f2c_pprts PROPERTIES PUBLIC_HEADER "f2c_solver_ids.h;f2c_pprts.h")

  add_library( f2c_pprts_rrtm_objs OBJECT f2c_pprts_rrtm.F90 )
  add_dependencies( f2c_pprts_rrtm_objs pprts_rrtmg pprts pprts_base )
  if(ENABLE_PLEXRT)
    add_dependencies( f2c_pprts_rrtm_objs plexrt )
  endif()

  add_library( f2c_pprts_rrtm $<TARGET_OBJECTS:f2c_pprts_rrtm_objs> )
  target_link_libraries( f2c_pprts_rrtm pprts_rrtmg pprts pprts_base ${EXTLIB})
  if(ENABLE_PLEXRT)
    target_link_libraries( f2c_pprts_rrtm plexrt )
  endif()

  file(COPY ./
    DESTINATION ${CMAKE_BINARY_DIR}/include/
    FILES_MATCHING PATTERN "*.h" )

  add_executable( "cwrapper_opp" "opp.c" )
  target_link_libraries( "cwrapper_opp" ${TENSTR_SINGLE_LIB} ${EXTLIB})

  if(WITH_PETSC)
    add_executable( "cwrapper_pprts" "pprts.c" )
    target_link_libraries( "cwrapper_pprts" ${TENSTR_SINGLE_LIB} ${EXTLIB})

    add_executable( "cwrapper_pprts_rrtmg" "pprts_rrtmg.c" )
    target_link_libraries( "cwrapper_pprts_rrtmg" ${TENSTR_SINGLE_LIB} ${EXTLIB})
  endif()

endif()
