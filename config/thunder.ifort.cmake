# Thunder
# module switch intel intel/15.0.1
# module switch openmpi openmpi/1.6.5-static-intel15

set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_Fortran_COMPILER "mpif90")

set(USER_C_FLAGS "-std=c99")
set(USER_Fortran_FLAGS "-cpp -traceback -r8 -ftz -extend_source")
set(USER_Fortran_FLAGS_RELEASE "-O3 -no-prec-div -xAVX -fp-model source")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

set(NETCDF_INCLUDE_DIR "/scratch/mpi/mpiaes/m300362/libs/netcdf-fortran-intel15/include/")
set(NETCDF_LIB_1       "/scratch/mpi/mpiaes/m300362/libs/netcdf-fortran-intel15/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/scratch/mpi/mpiaes/m300362/libs/netcdf-intel15/lib/libnetcdf.a")

set(HDF5_INCLUDE_DIRS       "/scratch/mpi/mpiaes/m300362/libs/hdf5-intel15/include")
list(APPEND HDF5_LIBRARIES  "/scratch/mpi/mpiaes/m300362/libs/hdf5-intel15/lib/libhdf5hl_fortran.a")
list(APPEND HDF5_LIBRARIES  "/scratch/mpi/mpiaes/m300362/libs/hdf5-intel15/lib/libhdf5_fortran.a")

set(SZIP_LIB           "/sw/squeeze-x64/szip-latest-static/lib/libsz.a")

set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${SZIP_LIB} ${HDF5_LIBRARIES} m z curl jpeg)

