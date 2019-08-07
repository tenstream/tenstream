
#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${INCLUDE_DIRS})

set( disort_sources
  disort.f
  ERRPACK.f
  BDREF.f
  DISOBRDF.f
  RDI1MACH.f
  DISOTESTAUX.f
  LINPAK.f
  tenstream_disort.F90
  )

# disort lib
add_library( tenstr_disort ${disort_sources} )
target_link_libraries( tenstr_disort pprts_base ${EXTLIB} ${USERLIB} ${LIBS})


add_executable( ex_disort_tenstr ex_disort_tenstr.F90 )
target_link_libraries( ex_disort_tenstr tenstr_disort ${EXTLIB} ${USERLIB} ${LIBS})


# Official DISORT tests
set( disort_test_sources
  disotest.f90
  DISOTESTAUX.f
  )

add_executable( disort_tests ${disort_test_sources} )
target_link_libraries( disort_tests tenstr_disort ${EXTLIB} ${USERLIB} ${LIBS})

add_test(disort_tests
  ${MPIEXEC} ${CTEST_MPIRUN_FLAGS} ${MPIEXEC_NUMPROC_FLAG} 1 ${PROJECT_BINARY_DIR}/bin/disort_tests)

# Ignore Warnings for disort files
string(COMPARE EQUAL ${CMAKE_Fortran_COMPILER_ID} "Intel" _is_Intel)
if(${_is_Intel})
  SET_SOURCE_FILES_PROPERTIES(${disort_sources} ${disort_test_sources}
    PROPERTIES COMPILE_FLAGS -nowarn)
else()
  SET_SOURCE_FILES_PROPERTIES(${disort_sources} ${disort_test_sources}
    PROPERTIES COMPILE_FLAGS "-Wno-error")
  string(REGEX REPLACE "-ffpe-trap=[a-z]+" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  string(REGEX REPLACE "-ffpe-trap=[a-z]+" "" CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}")
  string(REGEX REPLACE "-ffpe-trap=[a-z]+" "" CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}")
endif(${_is_Intel})
