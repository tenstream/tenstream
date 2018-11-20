
#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${INCLUDE_DIRS})

SET_SOURCE_FILES_PROPERTIES(ranlux.for PROPERTIES COMPILE_FLAGS -Wno-error)
SET_SOURCE_FILES_PROPERTIES(pprts.f90 PROPERTIES COMPILE_FLAGS -Wno-unused-function)
SET_SOURCE_FILES_PROPERTIES(tenstream.f90 PROPERTIES COMPILE_FLAGS -Wno-error)
SET_SOURCE_FILES_PROPERTIES(c_syscall_constants.c PROPERTIES COMPILE_FLAGS "-D_BSD_SOURCE")

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
  "optprop_LUT.f90"
  "optprop_ANN.f90"
  "optprop.f90"
  )


add_library( pprts_base STATIC ${general_sources} ${optprop_sources} )
target_link_libraries( pprts_base ${EXTLIB} ${USERLIB} ${LIBS})

add_library( pprts STATIC "pprts.f90")
target_link_libraries( pprts pprts_base ${EXTLIB} ${USERLIB} ${LIBS})

file( GLOB streams createLUT_*.f90 )
foreach ( stream ${streams} )
  GET_FILENAME_COMPONENT(basename "${stream}" NAME_WE )

  add_executable("${basename}" "${basename}.f90" )
  target_link_libraries( "${basename}" pprts_base)

  install(TARGETS ${basename} RUNTIME DESTINATION bin)
endforeach ()
