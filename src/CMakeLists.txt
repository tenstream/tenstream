
#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${INCLUDE_DIRS})

set( general_sources
  "helper_functions.f90"
  "helper_functions_dp.f90"
  "netcdfio.f90"
  "data_parameters.f90"
  "eddington.f90"
  "interpolation.f90"
  "mersenne.f90"
  "ranlux.for"
  "twostream.f90"
  "schwarzschild.f90"
  "tenstream.f90"
  "tenstream_options.f90"
  "f2c_tenstream.f90"
  "nca.f90"
  )

set( optprop_sources
  "boxmc.f90"
  "optprop_parameters.f90"
  "optprop_LUT.f90"
  "optprop_ANN.f90"
  "optprop.f90"
  )

add_library( tenstream STATIC ${general_sources} ${optprop_sources} )
target_link_libraries( tenstream ${EXTLIB} ${USERLIB} ${LIBS})

file( GLOB streams createLUT_*.f90 )
foreach ( stream ${streams} )
  GET_FILENAME_COMPONENT(basename "${stream}" NAME_WE )

  add_executable("${basename}" "${basename}.f90" )
  target_link_libraries( "${basename}" tenstream )

  install(TARGETS ${basename} RUNTIME DESTINATION bin)
endforeach ()


add_executable( "Cwrapper_test" "f2c_tenstream.c" )
target_link_libraries( "Cwrapper_test" tenstream)

add_executable( "fitLUT" "fitLUT.f90" )
target_link_libraries( "fitLUT" tenstream )

add_executable( "test_eddington" "test_eddington.f90" )
target_link_libraries( "test_eddington" tenstream )

install(TARGETS tenstream DESTINATION lib)

