#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

set( general_sources 
  "helper_functions.f90" 
  "arrayio.f90"  
  "data_parameters.f90"  
  "eddington.f90"  
  "gridtransform.f90"  
  "interpolation.f90"  
  "kato_data.f90"  
  "mersenne.f90"  
  "twostream.f90"
  "tenstream.f90"
  "tenstream_options.f90"
  "f2c_tenstream.f90"
  )

set( optprop_sources
  "boxmc.f90"
  "optprop_parameters.f90"
  "optprop_LUT.f90"
  "optprop_ANN.f90"
  "optprop.f90"
  )

add_library( tenstream STATIC ${general_sources} ${optprop_sources} )
TARGET_LINK_LIBRARIES( tenstream "${EXTLIB}" )

#target_include_directories (tenstream PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

file( GLOB streams createLUT_*.f90 )
foreach ( stream ${streams} )
  GET_FILENAME_COMPONENT(basename "${stream}" NAME_WE )

  add_executable("${basename}" "${basename}.f90" )                                                                                                                                                                   
  target_link_libraries( "${basename}" tenstream )                                                                                                                                                                  
endforeach ()
                                                                                                                                                                                                                         
add_executable( "petsc_solver" "petsc_tsv3.F90" )
target_link_libraries( "petsc_solver" tenstream )

add_executable( "test_boxmc" "test_boxmc.f90" )
target_link_libraries( "test_boxmc" tenstream )

add_executable( "Cwrapper_test" "f2c_tenstream.c" )
target_link_libraries( "Cwrapper_test" tenstream )

add_executable( "fitLUT" "fitLUT.f90" )
target_link_libraries( "fitLUT" tenstream )


install(TARGETS tenstream DESTINATION lib)

