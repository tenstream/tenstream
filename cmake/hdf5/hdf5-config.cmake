#-----------------------------------------------------------------------------
# HDF5 Config file for compiling against hdf5 install directory
#-----------------------------------------------------------------------------
GET_FILENAME_COMPONENT (SELF_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${SELF_DIR}" PATH)
GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
if (NOT WIN32)
  GET_FILENAME_COMPONENT(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
endif (NOT WIN32)

#-----------------------------------------------------------------------------
# User Options
#-----------------------------------------------------------------------------
set (HDF5_ENABLE_PARALLEL ON)
set (HDF5_BUILD_FORTRAN   ON)
set (HDF5_ENABLE_F2003    ON)
set (HDF5_BUILD_CPP_LIB   OFF)
set (HDF5_BUILD_TOOLS     ON)
set (HDF5_BUILD_HL_LIB    ON)
set (HDF5_ENABLE_Z_LIB_SUPPORT ON)
set (HDF5_ENABLE_SZIP_SUPPORT  ON)
set (HDF5_ENABLE_SZIP_ENCODING ON)
set (HDF5_BUILD_SHARED_LIBS    ON)
set (HDF5_PACKAGE_EXTLIBS      ON)

#-----------------------------------------------------------------------------
# Dependencies
#-----------------------------------------------------------------------------
IF(HDF5_ENABLE_PARALLEL)
  SET(HDF5_MPI_C_INCLUDE_PATH "/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/include")
  SET(HDF5_MPI_C_LIBRARIES    "/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/lib64/libmpi.so")
ENDIF(HDF5_ENABLE_PARALLEL)

#-----------------------------------------------------------------------------
# Directories
#-----------------------------------------------------------------------------
set (HDF5_INCLUDE_DIR "${_IMPORT_PREFIX}/include" "${HDF5_MPI_C_INCLUDE_PATH}" )

if (HDF5_BUILD_FORTRAN)
  set (HDF5_INCLUDE_DIR_FORTRAN "${_IMPORT_PREFIX}/include/fortran" )
endif (HDF5_BUILD_FORTRAN)
  
if (HDF5_BUILD_CPP_LIB)
  set (HDF5_INCLUDE_DIR_CPP "${_IMPORT_PREFIX}/include/cpp" )
endif (HDF5_BUILD_CPP_LIB)

if (HDF5_BUILD_HL_LIB)
  set (HDF5_INCLUDE_DIR_HL "${_IMPORT_PREFIX}/include/hl" )
endif (HDF5_BUILD_HL_LIB)

if (HDF5_BUILD_HL_LIB AND HDF5_BUILD_CPP_LIB)
  set (HDF5_INCLUDE_DIR_HL_CPP "${_IMPORT_PREFIX}/include/hl/cpp" )
endif (HDF5_BUILD_HL_LIB AND HDF5_BUILD_CPP_LIB)

if (HDF5_BUILD_TOOLS)
  set (HDF5_INCLUDE_DIR_TOOLS "${_IMPORT_PREFIX}/include" )
  set (HDF5_TOOLS_DIR "${_IMPORT_PREFIX}/bin" )
endif (HDF5_BUILD_TOOLS)

#-----------------------------------------------------------------------------
# Version Strings
#-----------------------------------------------------------------------------
set (HDF5_VERSION_STRING 1.8.13)
set (HDF5_VERSION_MAJOR  1.8)
set (HDF5_VERSION_MINOR  13)

#-----------------------------------------------------------------------------
# Don't include targets if this file is being picked up by another
# project which has already built hdf5 as a subproject
#-----------------------------------------------------------------------------
if (NOT TARGET "hdf5")
  if (HDF5_ENABLE_Z_LIB_SUPPORT AND HDF5_PACKAGE_EXTLIBS AND NOT TARGET "zlib")
    include (${SELF_DIR}/../ZLIB/zlib-targets.cmake)
  endif (HDF5_ENABLE_Z_LIB_SUPPORT AND HDF5_PACKAGE_EXTLIBS AND NOT TARGET "zlib")
  if (HDF5_ENABLE_SZIP_SUPPORT AND HDF5_PACKAGE_EXTLIBS AND NOT TARGET "szip")
    include (${SELF_DIR}/../SZIP/szip-targets.cmake)
  endif (HDF5_ENABLE_SZIP_SUPPORT AND HDF5_PACKAGE_EXTLIBS AND NOT TARGET "szip")
  include (${SELF_DIR}/hdf5-targets.cmake)
  set (HDF5_LIBRARIES "hdf5;hdf5_f90cstub;hdf5_fortran;hdf5_hl_f90cstub;hdf5_hl_fortran;hdf5_tools;hdf5_hl")
endif (NOT TARGET "hdf5")

