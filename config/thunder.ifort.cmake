# Thunder
# module switch intel intel/15.0.1
# module switch openmpi openmpi/1.6.5-static-intel15

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_Fortran_COMPILER "mpif90")

set(USER_C_FLAGS "-std=c99")
set(USER_Fortran_FLAGS "-cpp -traceback -extend_source -g ")
set(USER_Fortran_FLAGS_RELEASE " -O3 -no-prec-div -xAVX -fp-model source ")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")


set(NETCDF_DIR      "/sw/squeeze-x64/netcdf-4.2-static/")
set(NETCDF_DIR_F90  "/sw/squeeze-x64/netcdf_fortran-4.2-static-intel15/")
#set(HDF_ROOT        "/sw/squeeze-x64/hdf5-1.8.8-static/") # we should set this if hdf5 was made with cmake... otherwise set libraries by hand:

set(HDF5_LIB_1         "/sw/squeeze-x64/hdf5-1.8.8-static/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/sw/squeeze-x64/hdf5-1.8.8-static/lib/libhdf5.a")
set(SZIP_LIB           "/sw/squeeze-x64/szip-2.1-static/lib/libsz.a")
set(LIBS ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB})
