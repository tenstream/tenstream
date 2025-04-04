# Add our custom cmake find scripts
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

message(STATUS "CMake Version ${CMAKE_VERSION}")
cmake_minimum_required (VERSION 3.6.0...3.26.5)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR
    " In-source builds not allowed. Please make a new directory (e.g. called build) and run CMake from there.\n \n"
    " CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}\n"
    " cannot be\n"
    " CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}\n \n"
    " However, before that, you have to remove the already created files,\n  e.g. with `rm -rf CMakeCache.txt CMakeFiles/`")
endif()

#Load system specific settings
if ($ENV{FC} MATCHES "gfortran|mpif90")
  set (SYST gcc CACHE STRING "Location specific settings from config/<SYST>.cmake files")
elseif ($ENV{FC} MATCHES "ifort|mpiifort")
  set (SYST intel CACHE STRING "Location specific settings from config/<SYST>.cmake files")
elseif ($ENV{FC} MATCHES "ifx|mpiifx")
  set (SYST intel CACHE STRING "Location specific settings from config/<SYST>.cmake files")
endif()

if(SYST)
  list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/config)
  include ("${SYST}" OPTIONAL RESULT_VARIABLE FOUND_SYST)
  if(NOT FOUND_SYST)
    message(FATAL_ERROR "Could not find SYST config file for :: SYST=`${SYST}` :: ${FOUND_SYST}")
  else()
    message(STATUS "Using System config file: ${FOUND_SYST}")
  endif(NOT FOUND_SYST)
endif(SYST)

########################### Tenstream Project: ###################################
project (tenstream
  LANGUAGES C Fortran CXX
  HOMEPAGE_URL https://gitlab.com/jakubfabian/tenstream.git
  DESCRIPTION  "Library to solve Radiative Transer for participating media"
  VERSION 1.0 )
enable_testing()

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CMAKE_C_STANDARD 99)

# For python interfaces need -fPIC:
set(PythonCompat_FIND_VERSION 3)
find_package(PythonCompat REQUIRED)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10)
  set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -std=legacy")
endif()

# Add IPO support if compiler allows it -- default is off because it A) doesnt work or B) compiles super slow
option(ENABLE_IPO "enable link time optimizations, also called INTERPROCEDURAL_OPTIMIZATIONs" OFF)
if(ENABLE_IPO)
  cmake_minimum_required (VERSION 3.9.6)
  cmake_policy(SET CMP0069 NEW)
  set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT ipo_supported OUTPUT ipo_errmsg)
  if(ipo_supported)
    message(STATUS "IPO Enabled")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(WARNING "IPO could not be enabled ${ipo_errmsg}")
  endif()
endif()

option(ENABLE_GCOV "enable code coverage flags" OFF)
if(ENABLE_GCOV)
  set(GCC_COVERAGE_COMPILE_FLAGS "-g -coverage -fprofile-arcs -ftest-coverage")
  set(GCC_COVERAGE_LINK_FLAGS    "-coverage -lgcov")
  message(INFO " Enabling Gcov flags: ${GCC_COVERAGE_COMPILE_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
  set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
  set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
endif()

set(CMAKE_C_FLAGS               "${USER_C_FLAGS} ${CMAKE_C_FLAGS} " )
set(CMAKE_Fortran_FLAGS         "${USER_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS} " )
set(CMAKE_Fortran_FLAGS_RELEASE "${USER_Fortran_FLAGS_RELEASE} ${CMAKE_Fortran_FLAGS_RELEASE} " )
set(CMAKE_Fortran_FLAGS_DEBUG   "${USER_Fortran_FLAGS_DEBUG} ${CMAKE_Fortran_FLAGS_DEBUG} ")

# Add PGO Options
set(ENABLE_PGO "OFF" CACHE STRING "Use PGO in build, first compile with PGO set to `generate`, run make check, and then set to `use`")
set_property(CACHE ENABLE_PGO PROPERTY STRINGS OFF GENERATE USE)

if(NOT ENABLE_PGO STREQUAL "OFF")
  enable_testing()
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(PGO_GEN_FLAGS "-fprofile-generate=${CMAKE_BINARY_DIR}/pgo")
    set(PGO_USE_FLAGS "-fprofile-use=${CMAKE_BINARY_DIR}/pgo")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(PGO_GEN_FLAGS "-prof-gen -prof-dir=${CMAKE_BINARY_DIR}/pgo")
    set(PGO_USE_FLAGS "-prof-use -prof-dir=${CMAKE_BINARY_DIR}/pgo")
  else()
    message(FATAL_ERROR " Dont know the PGO command line options for ${CMAKE_Fortran_COMPILER_ID}, please tell me")
  endif()

  if(ENABLE_PGO STREQUAL "GENERATE")
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pgo)
    set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS}       ${PGO_GEN_FLAGS} ")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${PGO_GEN_FLAGS} ")
  elseif(ENABLE_PGO STREQUAL "USE")
    enable_testing()
    if(NOT EXISTS ${CMAKE_BINARY_DIR}/pgo)
      message(FATAL_ERROR " Could not find profiling data in ${CMAKE_BINARY_DIR}/pgo, please configure first with -DENABLE_PGO=GENERATE,"
        "then run examples and only then configure with -DENABLE_PGO=USE")
    endif()
    set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS}       ${PGO_USE_FLAGS} ")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${PGO_USE_FLAGS} ")
  endif()
endif()
add_custom_target(pgobuild
  COMMAND ${CMAKE_COMMAND} ${CMAKE_BINARY_DIR} -DENABLE_PGO=GENERATE
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --
  COMMAND ${CMAKE_CTEST_COMMAND} -E 'boxmc|LUT' -R 'pprts|plexrt'
  COMMAND ${CMAKE_COMMAND} ${CMAKE_BINARY_DIR} -DENABLE_PGO=USE
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --
  )

MESSAGE(STATUS "Build Type " ${CMAKE_BUILD_TYPE})
MESSAGE(STATUS "Compiler ID "  ${CMAKE_Fortran_COMPILER_ID})
if (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE})
  add_definitions(-D__RELEASE_BUILD__)
elseif(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  MESSAGE(STATUS "Compiler Flags " ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG})
else()
  message(FATAL " -- CMAKE_BUILD_TYPE should be either RELEASE or DEBUG")
endif()

## Set MPI Stuff
option(AUTODETECT_MPI "Let CMAKE try to find MPI variables" OFF)
if(AUTODETECT_MPI)
  find_package( MPI REQUIRED )

  add_definitions(${MPI_C_COMPILE_FLAGS})
  include_directories(${MPI_C_INCLUDE_PATH})
  list (APPEND EXTLIB "${MPI_C_LIBRARIES}" )

  add_definitions(${MPI_CXX_COMPILE_FLAGS})
  include_directories(${MPI_CXX_INCLUDE_PATH})
  list (APPEND EXTLIB "${MPI_CXX_LIBRARIES}" )

  if(${MPI_Fortran_FOUND})
    add_definitions(${MPI_Fortran_COMPILE_FLAGS})
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    link_directories(${MPI_Fortran_LIBRARIES})
    list (APPEND EXTLIB "${MPI_Fortran_LIBRARIES}" )
  else(${MPI_Fortran_FOUND})
    message(FATAL "!!!!!!! Fotran MPI MISSING")
  endif(${MPI_Fortran_FOUND})

  set( CMAKE_EXE_LINKER_FLAGS
    ${CMAKE_EXE_LINKER_FLAGS}
    ${MPI_CXX_LINK_FLAGS}
    ${MPI_Fortran_LIBRARIES}
    ${MPI_C_LIBRARIES}
    )
  string (REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

  list (APPEND CMAKE_INSTALL_RPATH
    ${MPI_C_LIBRARY_DIRS}
    ${MPI_CXX_LIBRARY_DIRS}
    ${MPI_Fortran_LIBRARY_DIRS})

  include_directories(${ITAPS_INCLUDES})
else(AUTODETECT_MPI)
endif(AUTODETECT_MPI)
set(MPIEXEC mpirun CACHE STRING "Command to run MPI programs, e.g. used to run tests")
if(DEFINED USER_MPIEXEC)
  set(MPIEXEC ${USER_MPIEXEC})
endif()
set(MPIEXEC_NUMPROC_FLAG -n CACHE STRING "options to set number of MPI tasks/nodes")
MESSAGE(STATUS "Using MPIEXEC: ${MPIEXEC}")

# Petsc
if(NOT PETSC_INSTALL)
  if(NOT DEFINED ENV{PETSC_DIR})
    message(VERBOSE "Environment variable $PETSC_DIR is not set.
                     Either export $PETSC_DIR and $PETSC_ARCH or supply -DPETSC_INSTALL=<path_to_petsc>
                     Ignore this warning in case that you are using pkg-config to find PETSc")
  else()
    message(STATUS "Using environment vars to find PETSc: PETSC_DIR=$ENV{PETSC_DIR} PETSC_ARCH=$ENV{PETSC_ARCH}")
    set(PETSC_INSTALL "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}"
      CACHE STRING "Path to PETSc Installation, usually $PETSC_DIR/$PETSC_ARCH")
    MESSAGE(STATUS "Using PETSC ${PETSC_INSTALL} ")
    if(NOT EXISTS "${PETSC_INSTALL}" )
      message(FATAL "Directory ${PETSC_INSTALL} does not exist!")
    endif()
  endif(NOT DEFINED ENV{PETSC_DIR})
endif()
if(PETSC_INSTALL)
  set(CMAKE_PREFIX_PATH ${PETSC_INSTALL} ${CMAKE_PREFIX_PATH})
  message(STATUS "Added ${PETSC_INSTALL} to CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
endif(PETSC_INSTALL)

include(FindPkgConfig)
pkg_search_module(PETSC REQUIRED PETSc)
if(EXISTS ${PETSC_LIBDIR})
  set (PETSC_LDFLAGS "-L${PETSC_LIBDIR} ${PETSC_LDFLAGS}")
endif()
message(STATUS "PETSC_INCLUDEDIR     : ${PETSC_INCLUDEDIR}    ")
message(STATUS "PETSC_INCLUDE_DIRS   : ${PETSC_INCLUDE_DIRS}  ")
message(STATUS "PETSC_LIBDIR         : ${PETSC_LIBDIR}        ")
message(STATUS "PETSC_LINK_LIBRARIES : ${PETSC_LINK_LIBRARIES}")
message(STATUS "PETSC_LDFLAGS        : ${PETSC_LDFLAGS}       ")
list (APPEND EXTLIB ${PETSC_LINK_LIBRARIES} ${PETSC_LDFLAGS})
list (APPEND CMAKE_INSTALL_RPATH ${PETSC_LIBDIR})
include_directories(${PETSC_INCLUDEDIR} ${PETSC_INCLUDE_DIRS})

if(PETSC_STATIC_LIBS)
  message(STATUS "PETSC_STATIC_LIBS         : ${PETSC_STATIC_LIBS}")
  message(STATUS "PETSC_STATIC_LIBRARIES    : ${PETSC_STATIC_LIBRARIES}")
  message(STATUS "PETSC_STATIC_LDFLAGS      : ${PETSC_STATIC_LDFLAGS}")
  list (APPEND EXTLIB ${PETSC_STATIC_LDFLAGS} ${PETSC_STATIC_LIBS} ${PETSC_STATIC_LIBRARIES})
endif()
if(PETSC_VERSION VERSION_GREATER_EQUAL "3.16.99")
  add_definitions(-DPETSC_DMPlexDistributeSetDefault)
endif()

# NetCDF
set (NETCDF_FOUND FALSE)
set (NETCDF_F90 "YES")
find_package( NetCDF QUIET )
if(${NETCDF_FOUND})
  message(STATUS "Found NetCDF:     ${NETCDF_FOUND}" )
  message(STATUS "NetCDF  Includes: ${NETCDF_INCLUDES}" )
  message(STATUS "NetCDFF Includes: ${NETCDF_INCLUDES_F90}" )
  message(STATUS "NetCDF  Library : ${NETCDF_LIBRARIES}")
  include_directories(${NETCDF_INCLUDES})
  include_directories(${NETCDF_INCLUDES_F90})
  list (APPEND EXTLIB ${NETCDF_LIBRARIES})
  list (APPEND CMAKE_INSTALL_RPATH ${NETCDF_LIBRARY_DIRS})
endif(${NETCDF_FOUND})

# if netcdf not found through find_package, try pkg-config
if(NOT ${NETCDF_FOUND})
  message(STATUS "Looking for NETCDF via pkg-config modules: netcdf netcdf-fortran")
  include(FindPkgConfig)
  pkg_check_modules(NETCDF REQUIRED netcdf netcdf-fortran)
  if(${NETCDF_FOUND})
    message(STATUS "NetCDF  Includes: ${NETCDF_netcdf-fortran_INCLUDEDIR} ${NETCDF_netcdf_INCLUDEDIR}" )
    message(STATUS "NetCDF  LDFLAGS : ${NETCDF_LDFLAGS} ${NETCDF_LDFLAGS_OTHER} ${NETCDF_STATIC_LDFLAGS} ${NETCDF_STATIC_LDFLAGS_OTHER}")
    message(STATUS "NetCDF  LibDir  : ${NETCDF_LIBRARY_DIRS}")

    include_directories(${NETCDF_netcdf-fortran_INCLUDEDIR}  ${NETCDF_netcdf_INCLUDEDIR})
    list (APPEND EXTLIB
      ${NETCDF_LDFLAGS}
      ${NETCDF_LDFLAGS_OTHER}
      ${NETCDF_STATIC_LDFLAGS}
      ${NETCDF_STATIC_LDFLAGS_OTHER})
    list (APPEND CMAKE_INSTALL_RPATH ${NETCDF_LIBRARY_DIRS})

    # check if netcdf supports parallel calls
    file(STRINGS "${NETCDF_netcdf_INCLUDEDIR}/netcdf_meta.h" _netcdf_lines REGEX "#define[ \t]+NC_HAS_PARALLEL[ \t]")
    string(REGEX REPLACE ".*NC_HAS_PARALLEL[ \t]*([0-1]+).*" "\\1" _netcdf_has_parallel "${_netcdf_lines}")
    if (_netcdf_has_parallel)
      add_definitions(-D__HAVE_NC_PARALLEL__)
      message(STATUS "Netcdf supports parallel ${_netcdf_has_parallel}")
    else()
      message(STATUS "Netcdf does NOT supports parallel ${_netcdf_has_parallel}")
    endif()


  endif(${NETCDF_FOUND})
endif(NOT ${NETCDF_FOUND})

if(NOT ${NETCDF_FOUND})
  message(FATAL_ERROR "CMake could not find your NETCDF install, you may need to add the NETCDF libraries in config/<loc>.cmake by hand." )
endif(NOT ${NETCDF_FOUND})

# HDF5
find_package( HDF5 REQUIRED ) # use this for non-cmake install of HDF5
message(STATUS "Found HDF5    : ${HDF5_FOUND}" )
if(${HDF5_FOUND})
  message(STATUS "HDF5 Includes: ${HDF5_INCLUDE_DIR}" )
  message(STATUS "HDF5 Library : ${HDF5_LIBRARIES}")
  include_directories(${HDF5_INCLUDE_DIR})
  link_directories("${HDF5_BUILD_DIR}/lib")
  list (APPEND ${EXTLIB} ${HDF5_C_LIBRARY} ${HDF5_LIBRARIES})
  list (APPEND CMAKE_INSTALL_RPATH ${HDF5_BUILD_DIR}/lib ${HDF5_LIBRARY_DIRS})
else(${HDF5_FOUND})
  message(WARNING "CMake could not find your HDF5 install, you may need to add the hdf5 libraries in config/<loc>.cmake by hand." )
endif(${HDF5_FOUND})


# Add Lapack
option(AUTODETECT_LAPACK "Let CMAKE try to find LAPACK" ON)
if(AUTODETECT_LAPACK)
  find_package( LAPACK )
  if(${LAPACK_FOUND})
    list (APPEND EXTLIB "${LAPACK_LIBRARIES}")
  endif(${LAPACK_FOUND})
endif(AUTODETECT_LAPACK)

## Add User Libs
if(DEFINED USER_LIBS)
    list (APPEND EXTLIB "${USER_LIBS}" )
    message(STATUS "EXTLIB after user LIBS: ${USER_LIBS}")
endif()

set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

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
install(DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ DESTINATION lib)

# set path to '/share', i.e. common data files
set (TENSTREAM_SHARE_DIR ${CMAKE_BINARY_DIR}/share CACHE STRING
  "Directory where we find tenstream related data files"
  FORCE)
set (TENSTREAM_SHARE_DIR_FULL ${TENSTREAM_SHARE_DIR})
install(DIRECTORY ${TENSTREAM_SHARE_DIR} DESTINATION .)

# split share dir string into multiple to adhere to fortran char line limits
string(LENGTH "${TENSTREAM_SHARE_DIR_FULL}" TENSTREAM_SHARE_DIR_LEN)
set(TENSTREAM_SHARE_DIR_CNT 0)
set(TENSTREAM_SHARE_DIR_SUBLEN 100)
while(TENSTREAM_SHARE_DIR_LEN GREATER 0)
  MATH(EXPR TENSTREAM_SHARE_DIR_CNT "${TENSTREAM_SHARE_DIR_CNT}+1")
  string(SUBSTRING "${TENSTREAM_SHARE_DIR_FULL}" 0 ${TENSTREAM_SHARE_DIR_SUBLEN} SUBSTR)
  add_definitions(-DTENSTREAM_SHARE_DIR${TENSTREAM_SHARE_DIR_CNT}=\"${SUBSTR}\")
  MESSAGE(STATUS "TENSTREAM_SHARE_DIR${TENSTREAM_SHARE_DIR_CNT}= ${SUBSTR}")
  if(TENSTREAM_SHARE_DIR_LEN GREATER ${TENSTREAM_SHARE_DIR_SUBLEN})
    string(SUBSTRING "${TENSTREAM_SHARE_DIR_FULL}" ${TENSTREAM_SHARE_DIR_SUBLEN} -1 TENSTREAM_SHARE_DIR_FULL)
  else()
    string(SUBSTRING "${TENSTREAM_SHARE_DIR_FULL}" ${TENSTREAM_SHARE_DIR_LEN} -1 TENSTREAM_SHARE_DIR_FULL)
  endif()
  string(LENGTH "${TENSTREAM_SHARE_DIR_FULL}" TENSTREAM_SHARE_DIR_LEN)
endwhile()


# Find out if the compiler is XLF and if so give this to the C preprocessor
string(COMPARE EQUAL ${CMAKE_Fortran_COMPILER_ID} "XL" _is_XLF)
if(${_is_XLF})
  add_definitions(-D_XLF)
  # and add the directory to the include path, as IBM xlf compiler will otherwise not find it. This is already an open bug in cmake and was supposed to be fixed.... grml
  include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
endif(${_is_XLF})

# Set preprocessor Flag if Compiler is out of GNU Family
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  add_definitions(-D_GNU)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  add_definitions(-D_INTEL)
endif()

option(ENABLE_C_WRAPPER "Enable creation of C Wrapper " ON)
option(ENABLE_ECCKD "Enable ECCKD" ON)
option(ENABLE_EXAMPLES "Enable compilation of EXAMPLES" ON)
option(ENABLE_PFUNIT "Enable pfUnit testing Framework" ON)
option(ENABLE_PLEXRT "Enable PLEXRT, Radiative transfer on wedges" ON)
option(ENABLE_REPWVL "Enable use of REPWVL" ON)
option(ENABLE_RRTM "Enable use of RRTM" ON)
option(ENABLE_SPECINT "Enable use of SPECINT" ON)

option(ENABLE_ALL "Enable use of all additional modules " OFF)
if(ENABLE_ALL)
  set(ENABLE_C_WRAPPER ON CACHE BOOL "Enable C_WRAPPER" FORCE)
  set(ENABLE_ECCKD  ON CACHE BOOL "Enable ECCKD" FORCE)
  set(ENABLE_EXAMPLES ON CACHE BOOL "Enable EXAMPLES" FORCE)
  set(ENABLE_PFUNIT ON CACHE BOOL "Enable PFUNIT" FORCE)
  set(ENABLE_PLEXRT ON CACHE BOOL "Enable PLEXRT" FORCE)
  set(ENABLE_REPWVL ON CACHE BOOL "Enable REPWVL" FORCE)
  set(ENABLE_RRTM   ON CACHE BOOL "Enable RRTM" FORCE)
  set(ENABLE_SPECINT ON CACHE BOOL "Enable SPECINT" FORCE)
endif(ENABLE_ALL)


# Add Source Directories
add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(disort)
add_subdirectory(rayli_wrapper)
add_subdirectory(misc)

if(ENABLE_PLEXRT)
  add_subdirectory (plexrt)
endif(ENABLE_PLEXRT)

set(TENSTR_LIBS
  pprts_base
  pprts
  tenstr_disort
  f2c_rayli
  )
if(ENABLE_C_WRAPPER)
  list(APPEND TENSTR_LIBS f2c_pprts f2c_pprts_rrtm)
endif()
if(ENABLE_PLEXRT)
  list(APPEND TENSTR_LIBS plexrt )
  if(ENABLE_RRTM)
    list(APPEND TENSTR_LIBS plexrt_rrtmg )
  endif()
endif()
if(ENABLE_RRTM)
  list(APPEND TENSTR_LIBS pprts_rrtm_sw pprts_rrtm_lw pprts_rrtmg)
endif()
if(ENABLE_REPWVL)
  list(APPEND TENSTR_LIBS repwvl_pprts)
endif()
if(ENABLE_ECCKD)
  list(APPEND TENSTR_LIBS ecckd_pprts)
endif()
if(ENABLE_SPECINT)
  list(APPEND TENSTR_LIBS specint)
endif()
if(BUILD_RAYLI)
  list(APPEND TENSTR_LIBS ${RAYLI_LIBS})
endif()
if(BUILD_FORNADO)
  list(APPEND TENSTR_LIBS ${FORNADO_LIBS})
endif()

# Generate Master Library which contains all symbols
if(NOT BUILD_SHARED_LIBS)
  function(combine_archives output_archive list_of_input_archives)
    set(obj_folder ${CMAKE_BINARY_DIR}/objs/)
    file(MAKE_DIRECTORY ${obj_folder})
    set(FULL_OUTPUT_PATH ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/lib${output_archive}.a)
    message(STATUS "Generating one lib to rule them all: ${FULL_OUTPUT_PATH}")

    set(output_archive_dummy_file ${CMAKE_BINARY_DIR}/${output_archive}.dummy.F90)
    add_custom_command(OUTPUT ${output_archive_dummy_file}
      COMMAND touch ${output_archive_dummy_file}
      DEPENDS ${list_of_input_archives})

    add_library(${output_archive} STATIC ${output_archive_dummy_file})

    FOREACH(in_archive ${list_of_input_archives})
      message(STATUS "Adding static lib: ${in_archive}")
      add_custom_command(TARGET ${output_archive}
        PRE_BUILD ${output_archive}
        WORKING_DIRECTORY ${obj_folder}
        COMMAND ${CMAKE_AR} x ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/lib${in_archive}.a)
    ENDFOREACH()

    add_custom_command(TARGET ${output_archive}
      POST_BUILD
      WORKING_DIRECTORY ${obj_folder}
      COMMAND ${CMAKE_AR} -qc ${FULL_OUTPUT_PATH} *.o)
  endfunction(combine_archives)


  SET(TENSTR_SINGLE_LIB "tenstream")
	combine_archives(${TENSTR_SINGLE_LIB} "${TENSTR_LIBS}")
endif()

if(ENABLE_RRTM)
	add_subdirectory(rrtmg)
endif(ENABLE_RRTM)

if(ENABLE_REPWVL)
	add_subdirectory(repwvl)
endif(ENABLE_REPWVL)

if(ENABLE_ECCKD)
  add_subdirectory(ecckd)
endif(ENABLE_ECCKD)

if(ENABLE_SPECINT)
  add_subdirectory(specint)
endif(ENABLE_SPECINT)

if(ENABLE_EXAMPLES)
	add_subdirectory(examples)
endif(ENABLE_EXAMPLES)

# Add pfUnit for testing
if(ENABLE_PFUNIT)
  add_subdirectory(tests)
  add_custom_target(check    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -E 'boxmc|LUT')
  add_custom_target(check_verbose COMMAND ${CMAKE_CTEST_COMMAND} -V -E 'boxmc|LUT')
  add_custom_target(checkLUT COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -R 'LUT')
  set(CTEST_MPIRUN_FLAGS "" CACHE STRING "Additional mpiexec flags to run tests with")
endif(ENABLE_PFUNIT)

if(ENABLE_C_WRAPPER)
	add_subdirectory(c_wrapper)
endif(ENABLE_C_WRAPPER)


find_program(CTAGS_PATH ctags)
if(CTAGS_PATH)
    message(STATUS "Found ctags: ${CTAGS_PATH}")

    execute_process(COMMAND ctags --help
        COMMAND grep "langmap"
        OUTPUT_VARIABLE CTAGS_LANGMAP)
    if(CTAGS_LANGMAP STREQUAL "")
      set(CTAGS_LANGMAP "")
    else()
      set(CTAGS_LANGMAP "--langmap=Fortran:.fypp.f.for.ftn.f77.f90.f95.F.FOR.FTN.F77.F90.F95")
    endif()
    message(STATUS "ctags supports langmap?: ${CTAGS_LANGMAP}")

    set_source_files_properties( tags PROPERTIES GENERATED true)
    add_custom_target(tags ALL
      COMMAND ctags -R
      ${CTAGS_LANGMAP}
      ${CMAKE_CURRENT_BINARY_DIR}/external
      ${CMAKE_CURRENT_SOURCE_DIR}/c_wrapper
      ${CMAKE_CURRENT_SOURCE_DIR}/disort
      ${CMAKE_CURRENT_SOURCE_DIR}/examples
      ${CMAKE_CURRENT_SOURCE_DIR}/plexrt
      ${CMAKE_CURRENT_SOURCE_DIR}/python
      ${CMAKE_CURRENT_SOURCE_DIR}/rayli_wrapper
      ${CMAKE_CURRENT_SOURCE_DIR}/rrtmg
      ${CMAKE_CURRENT_SOURCE_DIR}/repwvl
      ${CMAKE_CURRENT_SOURCE_DIR}/specint
      ${CMAKE_CURRENT_SOURCE_DIR}/src
      ${CMAKE_CURRENT_SOURCE_DIR}/tests
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
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

GET_DIRECTORY_PROPERTY(output INCLUDE_DIRECTORIES)

add_custom_target(download_common_luts
  COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/misc/download_LUT.sh _diffuse_10.tau31.w020.aspect_zx23.g6.ds1000.nc.Sdiff.mmap4 ${TENSTREAM_SHARE_DIR}
  COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/misc/download_LUT.sh _direct_3_10.tau31.w020.aspect_zx23.g6.phi19.theta19.ds1000.nc.Sdir.mmap ${TENSTREAM_SHARE_DIR}
  COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/misc/download_LUT.sh _direct_3_10.tau31.w020.aspect_zx23.g6.phi19.theta19.ds1000.nc.Tdir.mmap ${TENSTREAM_SHARE_DIR}
  COMMENT "Download most common lookuptables" VERBATIM
  )
add_custom_target(download_all_luts
  COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/misc/download_LUT.sh mmap4 ${TENSTREAM_SHARE_DIR}
  COMMENT "Download all lookuptables" VERBATIM
  )
add_custom_target(download_all_nc_luts
  COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/misc/download_LUT.sh nc ${TENSTREAM_SHARE_DIR}
  COMMENT "Download all lookuptables, including netcdf info" VERBATIM
  )

# Like a distclean
add_custom_target(distclean_clean
  COMMAND make clean
  COMMENT "Deleting compiler and linker files" VERBATIM
  )
add_custom_target(distclean
  DEPENDS distclean_clean
  COMMAND rm -rf
  ${CMAKE_BINARY_DIR}/CMakeCache.txt
  ${CMAKE_BINARY_DIR}/*.cmake
  ${CMAKE_BINARY_DIR}/CMakeFiles
  ${CMAKE_BINARY_DIR}/*/CMakeFiles
  COMMENT "Deleting CMake files" VERBATIM
  )

# Documentation targets
add_custom_target(doc_doxygen DEPENDS doxygen)

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


# list_all_cmake_vars
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#	message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

# create pkg-config files
string(REPLACE ";" " " EXTLIB_SPLIT "${EXTLIB}")
#set(pc_req_public "PETSc netcdf-fortran lapack") # would be nice but we dont want to fail if a lib is not in the PKG_CONFIG_PATH
#set(pc_req_private "")
CONFIGURE_FILE("pkgconfig/tenstream.pc.in" "lib/pkgconfig/tenstream.pc" @ONLY)
CONFIGURE_FILE("pkgconfig/tenstream.pc.in" "lib/pkgconfig/TenStream.pc" @ONLY)
install(DIRECTORY ${CMAKE_BINARY_DIR}/lib/pkgconfig DESTINATION lib)
