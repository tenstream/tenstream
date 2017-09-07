get_filename_component(CDIR ${CMAKE_CURRENT_LIST_DIR} NAME)
set(TESTDIR "${PROJECT_BINARY_DIR}/generated/${CDIR}" )

message( STATUS "Creating Cmake Test wrapper for directory :: ${TESTDIR} :: CDIR ${CDIR}" )

file(MAKE_DIRECTORY  ${TESTDIR} )
include_directories( ${TESTDIR} )


file(WRITE ${TESTDIR}/testSuites.inc "")
set(_test_sources)

file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.[fF]90")
message(STATUS "Found following test files: ${files}")
foreach(_file ${files})
  get_filename_component (_test_NAME ${_file} NAME   )
  get_filename_component (_test_WE   ${_file} NAME_WE)
  get_filename_component (_test_EXT  ${_file} EXT    )

  message(STATUS "FNAMES are: ${_test_NAME} :: ${_test_WE} :: ${_test_EXT}")
  if (DEFINED PFUNIT_INSTALL)
    set(test_dependency ${_file})
  else()
    set(test_dependency pfunit ${_file})
  endif()
  add_custom_command(
    OUTPUT ${TESTDIR}/${_test_NAME}
    COMMAND python ${PFUNIT_DIR}/bin/pFUnitParser.py ${_file} ${TESTDIR}/${_test_NAME}
    DEPENDS ${test_dependency}
    )
  set(_test_sources ${_test_sources} ${TESTDIR}/${_test_NAME})
  file(APPEND ${TESTDIR}/testSuites.inc "ADD_TEST_SUITE(${_test_WE}_suite)\n")
endforeach()

set_source_files_properties(${PFUNIT_DIR}/include/driver.F90 PROPERTIES GENERATED 1)

message(STATUS "Source files determined as :: ${_test_sources}")

add_executable(
  pfunit_${CDIR}
  ${PFUNIT_DIR}/include/driver.F90
  ${_test_sources}
  )
target_link_libraries(
  pfunit_${CDIR}
  ${PFUNIT_DIR}/lib/libpfunit.a
  tenstream
  )

add_test(pfunit_${CDIR} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 8 ${PROJECT_BINARY_DIR}/bin/pfunit_${CDIR})

