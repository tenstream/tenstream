
message(STATUS "Installing pfUnit  --  Installation taken from https://github.com/bast/pfunit-demo/blob/master/CMakeLists.txt")

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
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
set(BUILD_RAYLI_GIT_URL "git@gitlab.lrz.de:ri96kit/rayli.git"
  CACHE STRING
  "Git url which is used when downloading RayLi (e.g. to https://gitlab.lrz.de/ri96kit/rayli.git)"
  )

set(BUILD_RAYLI_GIT_TAG origin/feature/city_like_mystic
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
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${RAYLI_DIR}
    GIT_TAG ${BUILD_RAYLI_GIT_TAG}
    )
  FetchContent_GetProperties(rayli)
  if(NOT rayli_POPULATED)
    FetchContent_Populate(rayli)
    add_subdirectory(${rayli_SOURCE_DIR} ${rayli_BINARY_DIR})
  endif()

  set(RAYLI_LIBS rayli_c_wrapper rayli_common CACHE INTERNAL "RAYLI Libraries")
  message(STATUS " Rayli Libraries => ${RAYLI_LIBS}")

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
