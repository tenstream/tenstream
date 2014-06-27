message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

set( general_sources 
  "helper_functions.f90" 
  "arrayio.f90"  
  "data_parameters.f90"  
  "eddington.f90"  
  "gridtransform.f90"  
  "interpolation.f90"  
  "kato_data.f90"  
  "mersenne.f95"  
  "twostream.f90"
  )

set( optprop_sources
  "boxmc.f90"
  "optprop_parameters.f90"
  "optprop_LUT.f90"
  "optprop.f90"
  )

file( GLOB streams createLUT_*.f90 )
foreach ( stream ${streams} )
  GET_FILENAME_COMPONENT(basename "${stream}" NAME_WE )

  add_executable("${basename}" ${optprop_sources} ${general_sources} "${basename}.f90" )                                                                                                                                                                   
  target_link_libraries( "${basename}" "${EXTLIB}" )                                                                                                                                                                  
endforeach ()
                                                                                                                                                                                                                         
# petsc_solver
list (APPEND solver_sources 
  "${general_sources}"
  "${optprop_sources}"
  "petsc_tsv3.F90" 
  )

add_executable( "petsc_solver" ${solver_sources} )
target_link_libraries( "petsc_solver" "${EXTLIB}" )

add_executable( "test_boxmc" ${optprop_sources} ${general_sources} "test_boxmc.f90" )
target_link_libraries( "test_boxmc" "${EXTLIB}" )

#
### createLUT 
#file( GLOB streams optprop_dir*_diff* )
#foreach ( stream ${streams} )
#    GET_FILENAME_COMPONENT(basename "${stream}" NAME ) 
#
#    unset( optprop_sources )
#    file (GLOB_RECURSE optprop_sources "${stream}/createLUT.f90" "boxmc_oop*.f90" "${stream}/optprop*.?90")
#    list( APPEND optprop_sources "helper_functions.f90" "arrayio.f90" "mersenne.f95" "interpolation.f90" "data_parameters.f90" )
#
#    message( "Building CreateLUT for ${basename} with ${optprop_sources} " )
#
#    add_executable( "createLUT_${basename}" ${optprop_sources}  )
#    target_link_libraries( "createLUT_${basename}" "${EXTLIB}" )
#
#endforeach ()

## test_boxmc
#foreach ( stream ${streams} )
#    GET_FILENAME_COMPONENT(basename "${stream}" NAME ) 
#
#    unset( optprop_sources )
#    file (GLOB_RECURSE optprop_sources "${stream}/optprop*.?90" "${stream}/boxmc*.f90" "${stream}/test_boxmc.f90" "boxmc*.f90" )
#    list( APPEND optprop_sources "helper_functions.f90" "arrayio.f90" "mersenne.f95" "interpolation.f90" "data_parameters.f90" "eddington.f90" )
#
#    message( "Building test_boxmc for ${basename} with ${optprop_sources} " )
#
#    add_executable( "test_boxmc_${basename}" ${optprop_sources}  )
#    target_link_libraries( "test_boxmc_${basename}" "${EXTLIB}" )
#
#endforeach ()
