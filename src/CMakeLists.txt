
#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${INCLUDE_DIRS})

if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  SET_SOURCE_FILES_PROPERTIES(ranlux.F90 PROPERTIES COMPILE_FLAGS -nowarn)
else()
  SET_SOURCE_FILES_PROPERTIES(ranlux.F90 PROPERTIES COMPILE_FLAGS -Wno-error)
endif()
SET_SOURCE_FILES_PROPERTIES(c_syscall_constants.c PROPERTIES COMPILE_FLAGS "-D_BSD_SOURCE")


include(CheckSymbolExists)
check_symbol_exists(nc_set_log_level "netcdf.h" HAVE_NC_SET_LOG_LEVEL)
if (${HAVE_NC_SET_LOG_LEVEL})
  SET_SOURCE_FILES_PROPERTIES(netcdfio.F90 PROPERTIES COMPILE_DEFINITIONS __HAVE_NC_SET_LOG_LEVEL__=TRUE)
else(${HAVE_NC_SET_LOG_LEVEL})
  SET_SOURCE_FILES_PROPERTIES(netcdfio.F90 PROPERTIES COMPILE_DEFINITIONS __HAVE_NC_SET_LOG_LEVEL__=FALSE)
  message(STATUS "Could not find C function for 'nc_set_log_level'
     -- this is however needed if you want to show hdf5 error messages through the netcdf functions.
     If you need to debug netcdf/IO, recompile your netcdf-c lib with --enable-logging")
endif(${HAVE_NC_SET_LOG_LEVEL})


set(preprocessed_files "")
file( GLOB preprocessing_files *.fypp )
foreach ( infile ${preprocessing_files} )
  GET_FILENAME_COMPONENT(basename "${infile}" NAME_WE )

  SET(outfile "${CMAKE_CURRENT_BINARY_DIR}/${basename}.F90")
  message(STATUS "Preprocessing <${FYPP_EXE}> ${infile} -> ${outfile}")

  ADD_CUSTOM_COMMAND(
    OUTPUT "${outfile}"
    COMMAND ${Python_EXECUTABLE} ${FYPP_EXE} "${infile}" "${outfile}"
    MAIN_DEPENDENCY "${infile}"
    DEPENDS ${FYPP_EXE}
    VERBATIM
    )

  list(APPEND preprocessed_files "${outfile}")
endforeach ()


set( general_sources
  ${preprocessed_files}
  "adaptive_spectral_integration.F90"
  "buildings.F90"
  "c_syscall_constants.c"
  "c_syscall_wrappers.F90"
  "data_parameters.F90"
  "eddington.F90"
  "interpolation.F90"
  "mcdmda.F90"
  "mmap.F90"
  "nca.F90"
  "petsc_helpers.F90"
  "pprts_base.F90"
  "ranlux.F90"
  "schwarzschild.F90"
  "search.F90"
  "tenstream_options.F90"
  "twostream.F90"
  "xdmf_export.F90"
  )

set( optprop_sources
  "boxmc_geometry.F90"
  "boxmc.F90"
  "optprop_parameters.F90"
  "LUT_param_phi.F90"
  "optprop_LUT.F90"
  "optprop_ANN.F90"
  "optprop_base.F90"
  "optprop.F90"
  "geometric_coeffs.F90"
  )

set( solver_sources
  "pprts.F90"
  "pprts_shell.F90"
  "pprts_explicit.F90"
  "pprts_external_solvers.F90"
  )

if(BUILD_FORNADO)
  SET_SOURCE_FILES_PROPERTIES(optprop_ANN.F90
    PROPERTIES COMPILE_DEFINITIONS HAVE_FORNADO)
endif()

include(CheckSymbolExists)
check_symbol_exists(OMPI_MPI_H "mpi.h" HAVE_OMPI)
if(HAVE_OMPI)
  SET_SOURCE_FILES_PROPERTIES(petsc_helpers.F90
    PROPERTIES COMPILE_DEFINITIONS HAVE_OMPI)
endif()

# pprts_base libs
add_library( pprts_base_objs OBJECT ${general_sources} ${optprop_sources} )
if(BUILD_FORNADO)
  add_dependencies(pprts_base_objs ${FORNADO_LIBS})
endif()
add_library( pprts_base $<TARGET_OBJECTS:pprts_base_objs> )
target_link_libraries( pprts_base ${FORNADO_LIBS} ${EXTLIB})

# pprts libs
add_library( pprts_objs OBJECT ${solver_sources} )
add_dependencies( pprts_objs pprts_base plexrt )
add_library( pprts $<TARGET_OBJECTS:pprts_objs> )
target_link_libraries( pprts    pprts_base plexrt   ${EXTLIB})

# CreateLUT binaries
file( GLOB streams createLUT*.[fF]90 )
foreach ( stream ${streams} )
  GET_FILENAME_COMPONENT(basename "${stream}" NAME_WE )
  GET_FILENAME_COMPONENT(fnameext "${stream}" EXT )

  add_executable("${basename}" "${basename}${fnameext}" )
  target_link_libraries( "${basename}" pprts_base ${EXTLIB})
  install(TARGETS ${basename} DESTINATION bin/)
endforeach ()
