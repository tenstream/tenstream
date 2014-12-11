set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/config)
include(CMakeForceCompiler) #Necessary to change between MPI/parallel/profiler compilers, without having to do a make clean
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/petsc_cmake_modules/")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/hdf5/")

cmake_minimum_required (VERSION 2.8.1)

#Load system specific settings
if (NOT SYST)
  set (SYST default CACHE STRING
      "Choose the location: mpipc thunder fedora archlinux supermuc jugene" FORCE )
endif ()
include (${SYST} OPTIONAL)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

########################### Tenstream Project: ###################################
project (tenstream)
enable_language(Fortran)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

#set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_C_FLAGS               ${USER_C_FLAGS} )
set(CMAKE_Fortran_FLAGS         ${USER_Fortran_FLAGS} )
set(CMAKE_Fortran_FLAGS_RELEASE ${USER_Fortran_FLAGS_RELEASE} )
set(CMAKE_Fortran_FLAGS_DEBUG   ${USER_Fortran_FLAGS_DEBUG})

set(CMAKE_INCLUDE_SYSTEM_FLAG_Fortran "-I")

MESSAGE(STATUS "Build Type " ${CMAKE_BUILD_TYPE})
if (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE})
else()
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG})
endif()

# Add HDF5
set(HDF5_ROOT_DIR_HINT ${USER_HDF5_ROOT_DIR_HINT})
find_package(HDF5 REQUIRED)
message(STATUS "Found HDF5            :: ${HDF5_FOUND} ")
message(STATUS "HDF5 Includes         :: ${HDF5_INCLUDE_DIRS} ")
message(STATUS "HDF5 Libraries        :: ${HDF5_LIBRARIES} ")

include_directories (${HDF5_INCLUDE_DIRS})
list (APPEND EXTLIB "${HDF5_LIBRARIES}" )

# Add Lapack
find_package( LAPACK REQUIRED )
list (APPEND EXTLIB "${LAPACK_LIBRARIES}")

# Petsc
MESSAGE(STATUS "Using PETSC $ENV{PETSC_DIR} / $ENV{PETSC_ARCH} ")
find_package(PETSc REQUIRED)
message(STATUS "Found Petsc:    ${PETSC_FOUND}" )
message(STATUS "Petsc Includes: ${PETSC_INCLUDES}" )
message(STATUS "Petsc Library : ${PETSC_LIBRARIES}")
message(STATUS "Petsc compiler: ${PETSC_COMPILER}" )
message(STATUS "Petsc VERSION : ${PETSC_VERSION}" )
message(STATUS "Petsc definiti: ${PETSC_DEFINITIONS}")

add_definitions (${PETSC_DEFINITIONS})
include_directories (${PETSC_INCLUDES})
list (APPEND EXTLIB "${PETSC_LIBRARIES}" )

# Set MPI Stuff
find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})
list (APPEND EXTLIB "${MPI_Fortran_LIBRARIES}" )

include_directories(${ITAPS_INCLUDES})


# Set Output Directories
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Output directory for Dohp archives")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Output directory for Dohp libraries")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin CACHE PATH "Output directory for Dohp executables")

# Add Source Directories
add_subdirectory (src)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${PROJECT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)
