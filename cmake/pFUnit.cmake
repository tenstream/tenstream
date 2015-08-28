enable_testing()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fprofile-arcs -ftest-coverage")
set(PFUNIT_INSTALL $ENV{PFUNIT_INSTALL})

if(DEFINED PFUNIT_INSTALL)
    message(STATUS "Manual setup of variable PFUNIT_INSTALL : ${PFUNIT_INSTALL}")
    set(PFUNIT_DIR ${PFUNIT_INSTALL})
else(DEFINED PFUNIT_INSTALL)
    message(ERROR "##### PFUNIT_INSTALL environment variable not set")
endif(DEFINED PFUNIT_INSTALL)

include_directories(${PFUNIT_DIR}/source)

include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_BINARY_DIR}/generated
    )

file(WRITE ${PROJECT_BINARY_DIR}/generated/testSuites.inc "")
set(_test_sources)
file(GLOB files "${PROJECT_SOURCE_DIR}/tests/*.F90")
foreach(_file ${files})
    get_filename_component (_test ${_file} NAME_WE)
    set(test_dependency ${PROJECT_SOURCE_DIR}/tests/${_test}.F90)
    add_custom_command(
        OUTPUT ${PROJECT_BINARY_DIR}/generated/${_test}.F90
        COMMAND python ${PFUNIT_DIR}/bin/pFUnitParser.py ${PROJECT_SOURCE_DIR}/tests/${_test}.F90 ${PROJECT_BINARY_DIR}/generated/${_test}.F90
        DEPENDS ${test_dependency}
        )
    set(_test_sources ${_test_sources} ${PROJECT_BINARY_DIR}/generated/${_test}.F90)
    file(APPEND ${PROJECT_BINARY_DIR}/generated/testSuites.inc "ADD_TEST_SUITE(${_test}_suite)\n")
endforeach()

set_source_files_properties(${PFUNIT_DIR}/include/driver.F90 PROPERTIES GENERATED 1)

file(GLOB tenstream_source "${PROJECT_SOURCE_DIR}/src/tenstream.[fF]90")

add_executable(
    pftest_alltests
    ${PFUNIT_DIR}/include/driver.F90
    ${_test_sources}
    )
target_link_libraries(
    pftest_alltests
    ${PFUNIT_DIR}/source/libpfunit.a
    tenstream
    )

add_test(pftest_alltests ${PROJECT_BINARY_DIR}/pftest_alltests)
