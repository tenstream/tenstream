cmake_minimum_required (VERSION 2.8.10)
project (tenstream)
enable_language(Fortran)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set (CMAKE_Fortran_FLAGS "-cpp -ffree-line-length-none -W -Wall -Wuninitialized" )
set (CMAKE_Fortran_FLAGS_RELEASE "-fno-backtrace -ffast-math -O3 -funroll-loops ")#-march=native -mtune=native")
set (CMAKE_Fortran_FLAGS_DEBUG   "-finit-real=nan -fcheck=all -fbacktrace -O0 -g -pedantic -Wsurprising -g -pg -fbounds-check ")#-ffpe-trap=invalid,zero,overflow,denormal")

#set(CMAKE_VERBOSE_MAKEFILE ON)
#MESSAGE("Using PETSC $ENV{PETSC_DIR} / $ENV{PETSC_ARCH} ")

list (APPEND CMAKE_MODULE_PATH "$ENV{PETSC_DIR}/cmake-modules")
#list (APPEND CMAKE_MODULE_PATH "$ENV{HDF5_DIR}/share/cmake")
#message( "Searching for cmake Files in ${CMAKE_MODULE_PATH}" )

#zlib
#set( ZLIB_DIR "$ENV{PETSC_DIR}/../hdf5/hdf5-latest/build/ZLIB-prefix/src/ZLIB-build/zlib-config.cmake")

set( HDF5_ROOT_DIR_HINT "$ENV{HDF5_DIR}/share/cmake/hdf5")
FIND_PACKAGE (HDF5 NAMES hdf5 REQUIRED)
#INCLUDE_DIRECTORIES( $ENV{HDF5_DIR}/include )
#SET (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
#message(" HDF5 Includes are:${HDF5_INCLUDE_DIR}" )
#set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_Fortran_LINK_FLAGS} )
set(CMAKE_INCLUDE_SYSTEM_FLAG_Fortran "-I")

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

#message(" External Libraries are: ${EXTLIB}")

add_subdirectory (src)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${PROJECT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)
