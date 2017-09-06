list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/config)
include(CMakeForceCompiler) #Necessary to change between MPI/parallel/profiler compilers, without having to do a make clean
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/petsc_cmake_modules/")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

cmake_minimum_required (VERSION 2.8.10)

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

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.")
endif()

########################### Tenstream Project: ###################################
project (tenstream C Fortran)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

# For python interfaces need -fPIC:
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_C_FLAGS               ${USER_C_FLAGS} )
set(CMAKE_Fortran_FLAGS         ${USER_Fortran_FLAGS} )
set(CMAKE_Fortran_FLAGS_RELEASE ${USER_Fortran_FLAGS_RELEASE} )
set(CMAKE_Fortran_FLAGS_DEBUG   ${USER_Fortran_FLAGS_DEBUG})

#set(CMAKE_INCLUDE_SYSTEM_FLAG_Fortran "-I")

MESSAGE(STATUS "Build Type " ${CMAKE_BUILD_TYPE})
MESSAGE(STATUS "Compiler ID "  ${CMAKE_Fortran_COMPILER_ID})
if (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE})
else()
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG})
endif()

# HDF5
#find_package( HDF5 ) # use this for non-cmake install of HDF5
#FIND_PACKAGE (HDF5 NAMES hdf5 ) # or this if HDF5 was built with cmake support
#message(STATUS "Found HDF5    : ${HDF5_FOUND}" )
#if(${HDF5_FOUND})
#  message(STATUS "HDF5  Includes: ${HDF5_INCLUDE_DIR}" )
#  message(STATUS "HDF5 F90      : ${HDF5_BUILD_FORTRAN}" )
#  message(STATUS "HDF5 Library  : ${HDF5_LIBRARIES}")
#  include_directories(${HDF5_INCLUDE_DIR})
#  include_directories(${HDF5_FORTRAN_INCLUDE_DIR})
#  list (APPEND EXTLIB ${HDF5_LIBRARIES})
#else(${HDF5_FOUND})
#  message(WARNING "CMake could not find your HDF5 install, you may need to add the hdf5 libraries in config/<loc>.cmake by hand." )
#endif(${HDF5_FOUND})

# NetCDF
set (NETCDF_F90 "YES")
find_package( NetCDF )
if(${NETCDF_FOUND})
  message(STATUS "Found NetCDF:    ${NETCDF_FOUND}" )
  message(STATUS "NetCDF  Includes: ${NETCDF_INCLUDES}" )
  message(STATUS "NetCDFF Includes: ${NETCDF_INCLUDES_F90}" )
  message(STATUS "NetCDF  Library : ${NETCDF_LIBRARIES}")
  include_directories(${NETCDF_INCLUDES})
  include_directories(${NETCDF_INCLUDES_F90})
  list (APPEND EXTLIB ${NETCDF_LIBRARIES})
else(${NETCDF_FOUND})
  message(WARNING "CMake could not find your NETCDF install, you may need to add the NETCDF libraries in config/<loc>.cmake by hand." )
endif(${NETCDF_FOUND})

# Add Lapack
find_package( LAPACK )
if(${LAPACK_FOUND})
  list (APPEND EXTLIB "${LAPACK_LIBRARIES}")
endif(${LAPACK_FOUND})

# Petsc
MESSAGE(STATUS "Using PETSC $ENV{PETSC_DIR} / $ENV{PETSC_ARCH} ")
find_package( PETSc )
message(STATUS "Found Petsc:    ${PETSC_FOUND}" )
message(STATUS "Petsc Includes: ${PETSC_INCLUDES}" )
message(STATUS "Petsc Library : ${PETSC_LIBRARIES}")
message(STATUS "Petsc compiler: ${PETSC_COMPILER}" )
message(STATUS "Petsc VERSION : ${PETSC_VERSION}" )
message(STATUS "Petsc definiti: ${PETSC_DEFINITIONS}")

list(APPEND CMAKE_MODULE_PATH "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/petsc/conf/")
include(PETScConfig OPTIONAL)

if(PETSC_HDF5_LIB)
  message(STATUS "Petsc HDF5:     ${PETSC_HDF5_LIB}")
  add_definitions(-D_PETSC_HAVE_HDF5)
endif(PETSC_HDF5_LIB)

add_definitions (${PETSC_DEFINITIONS})
include_directories (${PETSC_INCLUDES})
list (APPEND EXTLIB "${PETSC_LIBRARIES}" )

## Set MPI Stuff
option(AUTODETECT_MPI "Let CMAKE try to find MPI variables" OFF)
if(AUTODETECT_MPI)
  find_package( MPI REQUIRED )
  if(${MPI_Fortran_FOUND})
    add_definitions(${MPI_Fortran_COMPILE_FLAGS})
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    link_directories(${MPI_Fortran_LIBRARIES})
    list (APPEND EXTLIB "${MPI_Fortran_LIBRARIES}" )
    include_directories(${ITAPS_INCLUDES})
    message(STATUS "+++++++ MPI FOUND")
  else(${MPI_Fortran_FOUND})
    message(WARNING "!!!!!!! MPI MISSING")
  endif(${MPI_Fortran_FOUND})
else(AUTODETECT_MPI)
    set(MPIEXEC_NUMPROC_FLAG -np) # This is just a wild guess but needed to run tests succesfully
endif(AUTODETECT_MPI)

# Set Output Directories
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Output directory for Dohp archives")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib CACHE PATH "Output directory for Dohp libraries")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin CACHE PATH "Output directory for Dohp executables")

# Install Fortran *.mod files into include dir
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include )
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

# Set default install dir to local path
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "default install path" FORCE )
endif()

# Find out if the compiler is XLF and if so give this to the C preprocessor
string(COMPARE EQUAL ${CMAKE_Fortran_COMPILER_ID} "XL" _is_XLF)
if(${_is_XLF})
  add_definitions(-D_XLF)
  # and add the directory to the include path, as IBM xlf compiler will otherwise not find it. This is already an open bug in cmake and was supposed to be fixed.... grml
  include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
endif(${_is_XLF})

# Add Source Directories
add_subdirectory (src)
add_subdirectory (plexrt)

# Documentation targets
add_custom_target(doc DEPENDS doxygen)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${PROJECT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
  add_custom_target(doxygen
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
else(FOXYGEN_FOUND)
  message(STATUS "Could not find Doxygen to create source code documentation")
endif(DOXYGEN_FOUND)

# Add pfUnit for testing
option(ENABLE_PFUNIT "Enable pfUnit testing Framework" OFF)
if(ENABLE_PFUNIT)
	enable_testing()
	add_subdirectory(tests)
	add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --verbose)
endif(ENABLE_PFUNIT)

option(ENABLE_RRTM "Enable use of RRTM" OFF)
if(ENABLE_RRTM)
	add_subdirectory(src_rrtmg)
endif(ENABLE_RRTM)

option(ENABLE_EXAMPLES "Enable compilation of EXAMPLES" OFF)
if(ENABLE_EXAMPLES)
	add_subdirectory(examples)
endif(ENABLE_EXAMPLES)

find_program(CTAGS_PATH ctags)
if(CTAGS_PATH)
    message(STATUS "Found ctags: ${CTAGS_PATH}")
    set_source_files_properties( tags PROPERTIES GENERATED true)
    add_custom_target(tags ALL
        COMMAND ctags -R ${CMAKE_CURRENT_SOURCE_DIR}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Creating ctags" VERBATIM
        )
endif(CTAGS_PATH)

option(ENABLE_PYTHON "Enable creation of Python Interfaces via f2py" OFF)
if(ENABLE_PYTHON)
  find_program(F2PY_PATH f2py REQUIRED)
  message(STATUS "Found f2py: ${F2PY_PATH}")
  add_custom_target(python_interface ALL)
	add_subdirectory(python)
endif(ENABLE_PYTHON)

option(ENABLE_C_WRAPPER "Enable creation of C Wrapper " ON)
if(ENABLE_C_WRAPPER)
	add_subdirectory(C_wrapper)
endif(ENABLE_C_WRAPPER)

GET_DIRECTORY_PROPERTY(output INCLUDE_DIRECTORIES)
MESSAGE(WARNING ${output})
