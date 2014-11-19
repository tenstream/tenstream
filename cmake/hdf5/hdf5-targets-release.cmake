#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "hdf5" for configuration "Release"
set_property(TARGET hdf5 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5 PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "m;dl;/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/lib64/libmpi.so;zlib;szip;dl;/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/lib64/libmpi.so"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhdf5.so.1.8.13"
  IMPORTED_SONAME_RELEASE "libhdf5.so.8.0.2"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5 )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5 "${_IMPORT_PREFIX}/lib/libhdf5.so.1.8.13" )

# Import target "hdf5_f90cstub" for configuration "Release"
set_property(TARGET hdf5_f90cstub APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5_f90cstub PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "hdf5;m;dl;/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/lib64/libmpi.so;zlib;szip;/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/lib64/libmpi_usempi.so;/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/lib64/libmpi_mpifh.so;/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/lib64/libmpi.so"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhdf5_f90cstub.so.1.8.13"
  IMPORTED_SONAME_RELEASE "libhdf5_f90cstub.so.8.0.2"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5_f90cstub )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5_f90cstub "${_IMPORT_PREFIX}/lib/libhdf5_f90cstub.so.1.8.13" )

# Import target "hdf5_fortran" for configuration "Release"
set_property(TARGET hdf5_fortran APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5_fortran PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "hdf5_f90cstub;hdf5;/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/lib64/libmpi_usempi.so;/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/lib64/libmpi_mpifh.so;/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/lib64/libmpi.so"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhdf5_fortran.so.1.8.13"
  IMPORTED_SONAME_RELEASE "libhdf5_fortran.so.8.0.2"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5_fortran )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5_fortran "${_IMPORT_PREFIX}/lib/libhdf5_fortran.so.1.8.13" )

# Import target "hdf5_hl_f90cstub" for configuration "Release"
set_property(TARGET hdf5_hl_f90cstub APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5_hl_f90cstub PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "hdf5_f90cstub;hdf5_hl"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhdf5_hl_f90cstub.so.1.8.13"
  IMPORTED_SONAME_RELEASE "libhdf5_hl_f90cstub.so.8.0.2"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5_hl_f90cstub )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5_hl_f90cstub "${_IMPORT_PREFIX}/lib/libhdf5_hl_f90cstub.so.1.8.13" )

# Import target "hdf5_hl_fortran" for configuration "Release"
set_property(TARGET hdf5_hl_fortran APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5_hl_fortran PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "hdf5_hl_f90cstub;hdf5_fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhdf5_hl_fortran.so.1.8.13"
  IMPORTED_SONAME_RELEASE "libhdf5_hl_fortran.so.8.0.2"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5_hl_fortran )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5_hl_fortran "${_IMPORT_PREFIX}/lib/libhdf5_hl_fortran.so.1.8.13" )

# Import target "hdf5_tools" for configuration "Release"
set_property(TARGET hdf5_tools APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5_tools PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "hdf5"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhdf5_tools.so.1.8.13"
  IMPORTED_SONAME_RELEASE "libhdf5_tools.so.8.0.2"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5_tools )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5_tools "${_IMPORT_PREFIX}/lib/libhdf5_tools.so.1.8.13" )

# Import target "hdf5_hl" for configuration "Release"
set_property(TARGET hdf5_hl APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(hdf5_hl PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "hdf5"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libhdf5_hl.so.1.8.13"
  IMPORTED_SONAME_RELEASE "libhdf5_hl.so.8.0.2"
  )

list(APPEND _IMPORT_CHECK_TARGETS hdf5_hl )
list(APPEND _IMPORT_CHECK_FILES_FOR_hdf5_hl "${_IMPORT_PREFIX}/lib/libhdf5_hl.so.1.8.13" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
