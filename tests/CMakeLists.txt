message(STATUS "Installing pfUnit  --  Installation taken from https://github.com/bast/pfunit-demo/blob/master/CMakeLists.txt")

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

enable_testing()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DUSE_MPI")
set(PFUNIT_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

string(COMPARE EQUAL ${CMAKE_Fortran_COMPILER_ID} "Intel" _is_Intel)
if(${_is_Intel})
# Find out if the compiler is Intel and if so add compile flag for parameter tests
  set(EXTRA_TESTDRIVER_COMPILE_FLAGS "-warn noerrors -warn nounused")
  set(EXTRA_TEST_COMPILE_FLAGS "-warn noerrors -warn nounused")
  set(PFUNIT_Fortran_FLAGS "${PFUNIT_Fortran_FLAGS} -assume realloc_lhs -check nobounds -nowarn")
else()
  set(EXTRA_TESTDRIVER_COMPILE_FLAGS "-Wno-error -w")
  set(EXTRA_TEST_COMPILE_FLAGS "-Wno-error -Wno-unused-dummy-argument -ffree-line-length-none")
  set(PFUNIT_Fortran_FLAGS "${PFUNIT_Fortran_FLAGS} -fprofile-arcs -ftest-coverage")
endif(${_is_Intel})


if(DEFINED PFUNIT_INSTALL)
    message(STATUS "Manual setup of variable PFUNIT_INSTALL : ${PFUNIT_INSTALL}")
    set(PFUNIT_DIR ${PFUNIT_INSTALL})
else()
    add_custom_target(
        git_update
        COMMAND git submodule init
        COMMAND git submodule update
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )

    include(ExternalProject)

    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/generated)

    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/pfunit
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_Fortran_FLAGS=${PFUNIT_Fortran_FLAGS}
        -DCMAKE_Fortran_FLAGS_RELEASE=
        -DMPI_USE_MPIEXEC=YES
        -DMPI=YES
        )
    ExternalProject_Add(pfunit
        DOWNLOAD_COMMAND git submodule update
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/pfunit
        BINARY_DIR ${PROJECT_BINARY_DIR}/external/pfunit-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/pfunit-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/pfunit-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${ExternalProjectCMakeArgs}
        )
    include_directories(${PROJECT_BINARY_DIR}/external/pfunit/mod)
    add_dependencies(pfunit git_update)
    set(PFUNIT_DIR ${PROJECT_BINARY_DIR}/external/pfunit)
endif()

include_directories(${PFUNIT_DIR}/mod)

macro(subdirlist result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child})
        list(APPEND dirlist ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

subdirlist(subdirs ${CMAKE_CURRENT_SOURCE_DIR} )

foreach(subdir ${subdirs})
  add_subdirectory(${subdir})
endforeach()
