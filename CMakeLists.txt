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

set(CMAKE_VERBOSE_MAKEFILE ON)
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

set(APP_NAME test_cmake)
message("Starting cmake for ${APP_NAME} ...")

add_subdirectory (src)

# Normally PETSc is built with MPI, if not, use CC=mpicc, etc
#find_package (PETSc REQUIRED)
#find_package (ITAPS COMPONENTS MESH GEOM REL)
#find_package (GDAL)
#find_package (HDF5)
#find_library (MATH_LIB m)
#if (ITAPS_REL_FOUND)
#  set (dHAVE_ITAPS_REL TRUE)
#endif ()
#if (ITAPS_GEOM_FOUND)
#  set (dHAVE_ITAPS_GEOM TRUE)
#endif ()
#if (GDAL_FOUND)
#  set (dHAVE_GDAL TRUE)
#endif ()
#if (NOT HDF5_INCLUDE_DIRS) # The CMake API changed sometime between 1.8.0 and 1.8.4
#  set (HDF5_INCLUDE_DIRS "${HDF5_INCLUDE_DIR}")
#endif ()
#find_path (MEMCHECK_HEADER valgrind/memcheck.h)
#if (MEMCHECK_HEADER)
#  set (dUSE_VALGRIND TRUE)
#endif ()
#include (CheckCSourceCompiles) # The name is misleading, this also tries to link
#check_c_source_compiles ("
##define PragmaQuote(a) _Pragma(#a)
##define PragmaGCC(a) PragmaQuote(GCC a)
#int main(int argc,char *argv[]) {
#PragmaGCC(diagnostic ignored \"-Wconversion\")
#char c = (int)argv[0][0] + argv[argc-1][0];
#return c;
#}" dHAVE_PRAGMA_GCC)
#
##find_package (Doxygen)
#
#option (dUSE_DEBUG "Compile Dohp with extra debugging" ON)
#option (dUSE_LOG_FINEGRAIN "Log fine grained operations for profiling (slower)" OFF)
#option (dUSE_PARALLEL_HDF5 "Write HDF5 files in parallel" OFF)
#option (Dohp_BUILD_TESTS "Build tests" ON)
#if (Dohp_BUILD_TESTS)
#  option (Dohp_TESTS_MALLOC_DUMP "Run all tests with -malloc_dump" ON)
#  mark_as_advanced (Dohp_TESTS_MALLOC_DUMP)
#endif ()
#option (Dohp_BUILD_TOOLS "Build tools" ON)
#option (Dohp_BUILD_SANDBOX "Build programs in sandbox" OFF)
#option (Dohp_PEDANTIC_WARNINGS "Compile with pedantic warnings" ON)
#option (Dohp_WERROR "Treat warnings as errors (-Werror)" OFF)
#option (BUILD_SHARED_LIBS "Build shared Dohp libraries" ON)
#
#add_definitions (-std=c99)
#
#if (Dohp_BUILD_TESTS)
#  enable_testing ()
#  set (Dohp_DATA_DIR "${Dohp_SOURCE_DIR}/data")
#endif (Dohp_BUILD_TESTS)
#
#if (Dohp_PEDANTIC_WARNINGS)
#  set (DEFAULT_PEDANTIC_FLAGS "-pedantic -Wall -Wextra -Wundef -Wshadow -Wpointer-arith -Wbad-function-cast -Wcast-align -Wwrite-strings -Wconversion -Wlogical-op -Wsign-compare -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Winline -Wno-long-long -Wmissing-format-attribute -Wmissing-noreturn -Wpacked -Wdisabled-optimization -Wmultichar -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wendif-labels -Wdeclaration-after-statement -Wold-style-definition -Winvalid-pch -Wmissing-field-initializers -Wvariadic-macros -Wunsafe-loop-optimizations -Wvolatile-register-var -Wstrict-aliasing -funit-at-a-time -Wno-sign-conversion")
#  #set (DEFAULT_PEDANTIC_FLAGS "-Wunreachable-code -Wfloat-equal -Wc++-compat")
#  #set (DEFAULT_PEDANTIC_FLAGS "-pedantic -Wall -Wextra -Winline -Wshadow -Wconversion -Wlogical-op -Wmissing-prototypes -Wvla")
#  #set (DEFAULT_PEDANTIC_FLAGS "${DEFAULT_PEDANTIC_FLAGS} -Wno-sign-conversion -Wwrite-strings -Wstrict-aliasing -Wcast-align -fstrict-aliasing")
#  #set (DEFAULT_PEDANTIC_FLAGS "${DEFAULT_PEDANTIC_FLAGS} -Wdisabled-optimization -funit-at-a-time")
#  #set (DEFAULT_PEDANTIC_FLAGS "${DEFAULT_PEDANTIC_FLAGS} -Wpadded")
#  set (Dohp_PEDANTIC_FLAGS ${DEFAULT_PEDANTIC_FLAGS} CACHE STRING "Compiler flags to enable pedantic warnings")
#  add_definitions (${Dohp_PEDANTIC_FLAGS})
#endif ()
#if (Dohp_WERROR)
#  add_definitions (-Werror)
#endif ()
#
#configure_file (${Dohp_SOURCE_DIR}/dohpconfig.h.in
#  ${Dohp_BINARY_DIR}/include/dohpconfig.h @ONLY)
#
#add_custom_target (etags
#  etags `find include src sandbox -name *.[ch]`
#  WORKING_DIRECTORY ${Dohp_SOURCE_DIR})
#
## Essential: include our directories first otherwise we can get internal headers from some installed path
#include_directories (${PETSC_INCLUDES} ${ITAPS_INCLUDES} ${HDF5_INCLUDE_DIRS})
#add_definitions (${PETSC_DEFINITIONS})
#
#set (Dohp_ITAPS_LIBRARIES)
#if (dHAVE_ITAPS_REL)
#  list (APPEND Dohp_ITAPS_LIBRARIES "${ITAPS_REL_LIBRARIES}")
#endif ()
#if (dHAVE_ITAPS_GEOM)
#  list (APPEND Dohp_ITAPS_LIBRARIES "${ITAPS_GEOM_LIBRARIES}")
#endif ()
#list (APPEND Dohp_ITAPS_LIBRARIES "${ITAPS_MESH_LIBRARIES}")
#set (Dohp_DEPENDENT_LIBRARIES "${Dohp_ITAPS_LIBRARIES}" "${PETSC_LIBRARIES}" "${HDF5_LIBRARIES}" "${MATH_LIB}")
#set (Dohp_LIBRARIES dohp ${Dohp_DEPENDENT_LIBRARIES})
#
#mark_as_advanced (CMAKE_ARCHIVE_OUTPUT_DIRECTORY CMAKE_LIBRARY_OUTPUT_DIRECTORY CMAKE_RUNTIME_OUTPUT_DIRECTORY)
#
#set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#
## Convenience for testing
#function (dohp_link_executable name source)
#  add_executable ("${name}" "${source}" ${ARGN})
#  target_link_libraries ("${name}" dohp)
#endfunction ()
#function (dohp_add_test name np)
#  if (Dohp_TESTS_MALLOC_DUMP)
#    set (extra_args "-malloc_dump")
#  endif ()
#  add_test ("${name}" "${Dohp_TestWithReference}" "${PETSC_MPIEXEC}" "${np}" "${CMAKE_CURRENT_SOURCE_DIR}/refout/${name}.refout" ${ARGN} ${extra_args})
#endfunction ()
#set (Dohp_TestWithReference "${Dohp_SOURCE_DIR}/TestWithReference.sh")
#
#add_subdirectory (include)
#add_subdirectory (src)
#
#if (Dohp_BUILD_SANDBOX)
#  add_subdirectory (${Dohp_SOURCE_DIR}/sandbox)
#endif (Dohp_BUILD_SANDBOX)
