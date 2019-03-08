
#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${INCLUDE_DIRS})

string(COMPARE EQUAL ${CMAKE_Fortran_COMPILER_ID} "Intel" _is_Intel)
if(${_is_Intel})
  SET_SOURCE_FILES_PROPERTIES(ranlux.for PROPERTIES COMPILE_FLAGS -nowarn)
else()
  SET_SOURCE_FILES_PROPERTIES(ranlux.for PROPERTIES COMPILE_FLAGS -Wno-error)
endif(${_is_Intel})
SET_SOURCE_FILES_PROPERTIES(c_syscall_constants.c PROPERTIES COMPILE_FLAGS "-D_BSD_SOURCE")


check_function_exists(nc_set_log_level HAVE_NC_SET_LOG_LEVEL)
if (${HAVE_NC_SET_LOG_LEVEL})
  SET_SOURCE_FILES_PROPERTIES(netcdfio.f90 PROPERTIES COMPILE_DEFINITIONS __HAVE_NC_SET_LOG_LEVEL__=TRUE)
else(${HAVE_NC_SET_LOG_LEVEL})
  SET_SOURCE_FILES_PROPERTIES(netcdfio.f90 PROPERTIES COMPILE_DEFINITIONS __HAVE_NC_SET_LOG_LEVEL__=FALSE)
  message(STATUS "Could not find C function for 'nc_set_log_level'
     -- this is however needed if you want to show hdf5 error messages through the netcdf functions.
     If you need to debug file reading or writing,
     please recompile your netcdf-c lib with --enable-logging")
endif(${HAVE_NC_SET_LOG_LEVEL})


set( general_sources
  "helper_functions.f90"
  "helper_functions_dp.f90"
  "c_syscall_wrappers.f90"
  "c_syscall_constants.c"
  "netcdfio.f90"
  "data_parameters.f90"
  "eddington.f90"
  "interpolation.f90"
  "ranlux.for" "kiss.F90"
  "mcrts_dmda.F90"
  "twostream.f90"
  "schwarzschild.f90"
  "tenstream_options.f90"
  "petsc_helpers.f90"
  "nca.f90"
  "mmap.f90"
  "pprts_base.F90"
  "adaptive_spectral_integration.f90"
  )

set( optprop_sources
  "boxmc_geometry.f90"
  "boxmc.f90"
  "optprop_parameters.f90"
  "LUT_param_phi.F90"
  "optprop_LUT.f90"
  "optprop_ANN.f90"
  "optprop.f90"
  )


add_library( pprts_base STATIC ${general_sources} ${optprop_sources} )
target_link_libraries( pprts_base ${EXTLIB} ${USERLIB} ${LIBS})

add_library( pprts STATIC "pprts.f90")
target_link_libraries( pprts pprts_base ${EXTLIB} ${USERLIB} ${LIBS})

file( GLOB streams createLUT*.[fF]90 )
foreach ( stream ${streams} )
  GET_FILENAME_COMPONENT(basename "${stream}" NAME_WE )
  GET_FILENAME_COMPONENT(fnameext "${stream}" EXT )

  add_executable("${basename}" "${basename}${fnameext}" )
  target_link_libraries( "${basename}" pprts_base)

  install(TARGETS ${basename} RUNTIME DESTINATION bin)
endforeach ()
