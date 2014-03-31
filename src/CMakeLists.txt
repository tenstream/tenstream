# test_cmake
list( APPEND test_cmake_SRCS test_cmake.f90 )
add_executable( "${APP_NAME}" "${test_cmake_SRCS}" )
target_link_libraries( "${APP_NAME}" "${MPI_Fortran_LIBRARIES}" )


# petsc_solver
file (GLOB_RECURSE optprop_sources "optprop_dir8_diff10/optprop*.?90" "optprop_dir8_diff10/boxmc*.f90" )

list (APPEND solver_sources "arrayio.f90"  "data_parameters.f90"  "eddington.f90"  "gridtransform.f90"  "interpolation.f90"  "kato_data.f90"  "mersenne.f95"  "petsc_tsv3.F90" )

add_executable( "petsc_solver" ${solver_sources} ${optprop_sources} )
target_link_libraries( "petsc_solver" "${MPI_Fortran_LIBRARIES}" "${HDF5_LIBRARIES}" "${PETSC_LIBRARIES}" )

# createLUT 
file( GLOB streams optprop_dir*_diff* )
foreach ( stream ${streams} )
    GET_FILENAME_COMPONENT(basename "${stream}" NAME ) 
    message( "Building CreateLUT for ${basename}" )
    file (GLOB_RECURSE optprop_sources "${stream}/optprop*.?90" "${stream}/boxmc*.f90" "createLUT.f90" )
    list( APPEND optprop_sources "arrayio.f90" "mersenne.f95" "interpolation.f90" )
    add_executable( "createLUT_${basename}" ${optprop_sources}  )
    target_link_libraries( "createLUT_${basename}" "${MPI_Fortran_LIBRARIES}" "${HDF5_LIBRARIES}" "${PETSC_LIBRARIES}" )
endforeach ()

#file (GLOB_RECURSE optprop_sources "optprop_dir8_diff10/optprop*.?90" "optprop_dir8_diff10/boxmc*.f90" )
#add_executable( "createLUT" ${solver_sources} ${optprop_sources} )
