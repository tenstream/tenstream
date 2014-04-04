message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

# test_cmake
list( APPEND test_cmake_SRCS test_cmake.f90 )
add_executable( "${APP_NAME}" "${test_cmake_SRCS}" )
target_link_libraries( "${APP_NAME}" "${MPI_Fortran_LIBRARIES}" )


# petsc_solver
file (GLOB_RECURSE optprop_sources "optprop_dir*_diff*/optprop*.?90" "optprop_dir*_diff*/boxmc*.f90" )

list (APPEND solver_sources "helper_functions.f90" "arrayio.f90"  "data_parameters.f90"  "eddington.f90"  "gridtransform.f90"  "interpolation.f90"  "kato_data.f90"  "mersenne.f95"  "petsc_tsv3.F90" )

add_executable( "petsc_solver" ${solver_sources} ${optprop_sources} )
target_link_libraries( "petsc_solver" "${EXTLIB}" )

# createLUT 
file( GLOB streams optprop_dir*_diff* )
foreach ( stream ${streams} )
    GET_FILENAME_COMPONENT(basename "${stream}" NAME ) 

    unset( optprop_sources )
    file (GLOB_RECURSE optprop_sources "${stream}/optprop*.?90" "${stream}/boxmc*.f90" "${stream}/createLUT.f90" )
    list( APPEND optprop_sources "helper_functions.f90" "arrayio.f90" "mersenne.f95" "interpolation.f90" "data_parameters.f90" )

    message( "Building CreateLUT for ${basename} with ${optprop_sources} " )

    add_executable( "createLUT_${basename}" ${optprop_sources}  )
    target_link_libraries( "createLUT_${basename}" "${EXTLIB}" )

endforeach ()

# test_boxmc
foreach ( stream ${streams} )
    GET_FILENAME_COMPONENT(basename "${stream}" NAME ) 

    unset( optprop_sources )
    file (GLOB_RECURSE optprop_sources "${stream}/optprop*.?90" "${stream}/boxmc*.f90" "${stream}/test_boxmc.f90" )
    list( APPEND optprop_sources "helper_functions.f90" "arrayio.f90" "mersenne.f95" "interpolation.f90" "data_parameters.f90" "eddington.f90" )

    message( "Building test_boxmc for ${basename} with ${optprop_sources} " )

    add_executable( "test_boxmc_${basename}" ${optprop_sources}  )
    target_link_libraries( "test_boxmc_${basename}" "${EXTLIB}" )

endforeach ()
