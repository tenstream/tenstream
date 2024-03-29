if(ENABLE_PFUNIT)
  message(STATUS "Installing pfUnit  --  Installation taken from https://github.com/bast/pfunit-demo/blob/master/CMakeLists.txt")
  enable_testing()

  set(PFUNIT_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${PFUNIT_Fortran_FLAGS} -DUSE_MPI")

  string(COMPARE EQUAL ${CMAKE_Fortran_COMPILER_ID} "Intel" _is_Intel)
  if(${_is_Intel})
    # Find out if the compiler is Intel and if so add compile flag for parameter tests
    set(EXTRA_TESTDRIVER_COMPILE_FLAGS
      "${PFUNIT_Fortran_FLAGS} -warn noerrors -warn nounused"
      CACHE INTERNAL "FLAGS used to compile pfunit driver files")
    set(EXTRA_TEST_COMPILE_FLAGS
      "${PFUNIT_Fortran_FLAGS} -warn noerrors -warn nounused"
      CACHE INTERNAL "FLAGS used to compile files after they have been parsed by pfunit")
    set(PFUNIT_Fortran_FLAGS "${PFUNIT_Fortran_FLAGS} -assume realloc_lhs -check nobounds -nowarn")
  else()
    set(EXTRA_TESTDRIVER_COMPILE_FLAGS "${PFUNIT_Fortran_FLAGS} -Wno-error -w"
      CACHE INTERNAL "FLAGS used to compile pfunit driver files")
    set(EXTRA_TEST_COMPILE_FLAGS
      "${PFUNIT_Fortran_FLAGS} \
      -Wno-error \
      -Wno-unused-dummy-argument \
      -ffree-line-length-none"
      CACHE INTERNAL
      "FLAGS used to compile files after they have been parsed by pfunit")
    set(PFUNIT_Fortran_FLAGS
      "${PFUNIT_Fortran_FLAGS} \
      -Wno-compare-reals \
      -Wno-conversion \
      -Wno-intrinsic-shadow \
      -Wno-maybe-uninitialized \
      -Wno-unused-dummy-argument \
      -Wno-unused-function \
      -Wno-unused-parameter \
      -Wno-unused-variable \
    ")
  endif(${_is_Intel})

  message(STATUS "Using PFUNIT_Fortran_FLAGS: ${PFUNIT_Fortran_FLAGS}")

  if(DEFINED PFUNIT_INSTALL)
    message(STATUS "Manual setup of variable PFUNIT_INSTALL : ${PFUNIT_INSTALL}")
    set(PFUNIT_DIR ${PFUNIT_INSTALL} CACHE INTERNAL "PFUNIT install directory")
  else()
    add_custom_target(
      git_update
      COMMAND git submodule init
      COMMAND git submodule update --init --recursive
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      )

    include(ExternalProject)

    file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/generated)

    set(ExternalProjectCMakeArgs
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external/pfunit
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_Fortran_FLAGS=${PFUNIT_Fortran_FLAGS}
      -DCMAKE_Fortran_FLAGS_RELEASE=${PFUNIT_Fortran_FLAGS}
      -DCMAKE_Fortran_FLAGS_DEBUG=${PFUNIT_Fortran_FLAGS}
      -DPYTHON_EXECUTABLE=${Python_EXECUTABLE}
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
    #include_directories(${PROJECT_BINARY_DIR}/external/pfunit/mod)
    add_dependencies(pfunit git_update)
    set(PFUNIT_DIR ${PROJECT_BINARY_DIR}/external/pfunit CACHE INTERNAL "PFUNIT install directory")
    set(PFUNIT_PARSER ${PFUNIT_DIR}/bin/pFUnitParser.py CACHE INTERNAL "PFUNIT path to parser")
    #set(PFUNIT_DIR ${PROJECT_BINARY_DIR}/external/pfunit/PFUNIT-4.2 CACHE INTERNAL "PFUNIT install directory")
    #set(PFUNIT_PARSER ${PFUNIT_DIR}/bin/funitproc CACHE INTERNAL "PFUNIT path to parser")
  endif()
endif(ENABLE_PFUNIT)


message(STATUS "Installing fypp - Fortran Preprocessor")
#ExternalProject_Add(fypp
#  GIT_REPOSITORY https://github.com/aradi/fypp.git
#  SOURCE_DIR "${PROJECT_BINARY_DIR}/external/fypp"
#  CONFIGURE_COMMAND ""
#  BUILD_COMMAND ""
#  INSTALL_COMMAND ""
#)
set(FYPP_EXE ${PROJECT_SOURCE_DIR}/external/fypp CACHE STRING "Python Script Location for Fortran PreProcessor")


# Build RayLi Raytracer
option(BUILD_RAYLI "Build RayLi - Raytracer" OFF)
set(BUILD_RAYLI_GIT_URL "https://gitlab.com/d70-t/rayli.git"
  CACHE STRING
  "Git url which is used when downloading RayLi (e.g. to git@gitlab.com:d70-t/rayli.git)"
  )

set(BUILD_RAYLI_GIT_TAG origin/master
  CACHE STRING
  "Git branch which is used when downloading RayLi (e.g. to origin/master)"
  )

if(BUILD_RAYLI)

  message(STATUS "Building RayLi - Raytracer (${BUILD_RAYLI_GIT_TAG})")
  cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
  include(FetchContent)
  set(RAYLI_DIR ${CMAKE_BINARY_DIR}/external/rayli CACHE INTERNAL "RAYLI install directory")

  FetchContent_Declare(rayli
    GIT_REPOSITORY ${BUILD_RAYLI_GIT_URL}
    SOURCE_DIR "${CMAKE_BINARY_DIR}/external/rayli-src"
    PREFIX ${RAYLI_DIR}
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${RAYLI_DIR}
      -DRAYLI_BUILD_EXAMPLES=OFF
    GIT_TAG ${BUILD_RAYLI_GIT_TAG}
    )
  FetchContent_GetProperties(rayli)
  if(NOT rayli_POPULATED)
    FetchContent_Populate(rayli)
    add_subdirectory(${rayli_SOURCE_DIR} ${rayli_BINARY_DIR})
  endif()

  set(RAYLI_LIBS rayli_c_wrapper rayli_common spdlog fmt CACHE INTERNAL "RAYLI Libraries")
  message(STATUS " Rayli Libraries => ${RAYLI_LIBS}")
  set(EXTLIB "${EXTLIB};stdc++;pthread;m" PARENT_SCOPE)
  message(STATUS "EXTLIB after user Rayli: ${EXTLIB}")
else(BUILD_RAYLI) # else, check if we have a custom RAYLI_DIR set as environment variable

  set(RAYLI_DIR "$ENV{RAYLI_DIR}")
  if(RAYLI_DIR)
    message(STATUS "Configuring with RayliDIR: ${RAYLI_DIR}")

    include_directories(${RAYLI_DIR}/include)
    find_library(RAYLI_LIB1 rayli_c_wrapper ${RAYLI_DIR} ${RAYLI_DIR}/lib )
    find_library(RAYLI_LIB2 rayli_common    ${RAYLI_DIR} ${RAYLI_DIR}/lib )

    set(RAYLI_LIBS "${RAYLI_LIB1};${RAYLI_LIB2}" CACHE INTERNAL "RAYLI Libraries")
    list(APPEND RAYLI_LIBS stdc++ pthread m)
    message(STATUS " Rayli Libraries => ${RAYLI_LIBS}")
  endif(RAYLI_DIR)

endif(BUILD_RAYLI)


option(BUILD_FORNADO "Build Fornado - Fortran tensorflow inference module" OFF)
if(BUILD_FORNADO)
  set(BUILD_FORNADO_GIT_URL "git@gitlab.com:LeonhardScheck/fornado.git"
    CACHE STRING
    "Git url which is used when downloading Fornado (e.g. https://gitlab.com/LeonhardScheck/fornado/)"
    )
  set(BUILD_FORNADO_GIT_TAG origin/FJ/refactoring_to_use_modules
    CACHE STRING
    "Git branch which is used when downloading Fornado (e.g. to origin/master)"
    )

  cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
  include(FetchContent)
  set(FORNADO_DIR ${CMAKE_BINARY_DIR}/external/fornado CACHE INTERNAL "FORNADO install directory")
  set(ENABLE_OPENMP OFF CACHE INTERNAL "Disable fornado openmp")

  FetchContent_Declare(fornado
    GIT_REPOSITORY ${BUILD_FORNADO_GIT_URL}
    SOURCE_DIR "${CMAKE_BINARY_DIR}/external/fornado-src"
    PREFIX ${FORNADO_DIR}
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${FORNADO_DIR} -DENABLE_PYTHON=OFF -DENABLE_OPENMP=OFF"
    GIT_TAG ${BUILD_FORNADO_GIT_TAG}
    )
  FetchContent_GetProperties(fornado)
  if(NOT fornado_POPULATED)
    FetchContent_Populate(fornado)
    add_subdirectory(${fornado_SOURCE_DIR} ${fornado_BINARY_DIR})
  endif()

  set(FORNADO_LIBS fornado CACHE INTERNAL "FORNADO Libraries")
  message(STATUS " Fornado Libraries => ${FORNADO_LIBS}")
endif(BUILD_FORNADO)
