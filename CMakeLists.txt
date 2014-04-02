cmake_minimum_required (VERSION 2.6.3)
project (tenstream)
enable_language(Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set (CMAKE_Fortran_FLAGS "-fbacktrace -finit-real=nan -ffree-line-length-none -W -Wall -Wuninitialized " )
set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -Ofast -ffast-math -march=native -mtune=native")
set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -fcheck=all -fbacktrace -O0 -g -pedantic -Wsurprising -g -pg -fbounds-check -ffpe-trap=invalid,zero,overflow")

#set(CMAKE_VERBOSE_MAKEFILE ON)
MESSAGE("Using PETSC $ENV{PETSC_DIR} / $ENV{PETSC_ARCH} ")

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/petsc_cmake_modules")
list (APPEND CMAKE_MODULE_PATH "/usr/users/jakub/ticalib/ompi1.7.3/hdf5/hdf5-latest/build/install/share/cmake/hdf5")
message( "Searching for cmake Files in ${CMAKE_MODULE_PATH}" )

#zlib
set( ZLIB_DIR "/usr/users/jakub/ticalib/ompi1.7.3/hdf5/hdf5-1.8.12/build/ZLIB-prefix/src/ZLIB-build/zlib-config.cmake")

#find_package(HDF5 REQUIRED)
FIND_PACKAGE (HDF5 NAMES hdf5 REQUIRED)
# FIND_PACKAGE (HDF5) # Find non-cmake built HDF5
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
SET (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})

# Petsc
find_package(PETSc REQUIRED)
add_definitions (${PETSC_DEFINITIONS})

# Set MPI Stuff
find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})

include_directories (${PETSC_INCLUDES} ${ITAPS_INCLUDES} ${HDF5_INCLUDE_DIRS})
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Output directory for Dohp archives")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Output directory for Dohp libraries")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin CACHE PATH "Output directory for Dohp executables")
set (CMAKE_INSTALL_RPATH "lib")

list (APPEND EXTLIB "${MPI_Fortran_LIBRARIES}" )
list (APPEND EXTLIB "${HDF5_LIBRARIES}" )
list (APPEND EXTLIB "${PETSC_LIBRARIES}" )

add_subdirectory (src)

