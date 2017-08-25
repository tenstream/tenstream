include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(ENABLE_RRTM)

  add_library( c_tenstream c_tenstream.F90 )
  target_link_libraries( c_tenstream tenstr_rrtmg tenstr_rrtm_lw tenstr_rrtm_sw tenstream )

  add_executable( "Cwrapper_ex1" "ex1.c" )
  target_link_libraries( "Cwrapper_ex1" c_tenstream tenstr_rrtmg tenstr_rrtm_lw tenstr_rrtm_sw tenstream )

endif()
