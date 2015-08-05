# Thunder
# For further installation instruction please look into the wiki documentation:
# <https://github.com/tenstream/tenstream/wiki/Install-instructions-for-thunder-at-ZMAW,-Hamburg>

# For IntelMPI
#set(CMAKE_C_COMPILER "mpiicc")
#set(CMAKE_CXX_COMPILER "mpiicpc")
#set(CMAKE_Fortran_COMPILER "mpiifort")

# For OpenMPI
set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")
set(CMAKE_Fortran_COMPILER "mpif90")

set(USER_C_FLAGS "-std=c99")
set(USER_Fortran_FLAGS "-cpp -traceback -extend-source -g -mkl ")
set(USER_Fortran_FLAGS_RELEASE " -O3 -no-prec-div -xAVX -fp-model source -fno-omit-frame-pointer")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

# From bashrc:
#    export NETCDFCROOT=/sw/squeeze-x64/netcdf-4.3.3.1-static/
#    export NETCDFFROOT=/sw/squeeze-x64/netcdf_fortran-4.2-static-intel15/
#    export HDF5ROOT=/sw/squeeze-x64/hdf5-latest-static/
#    export SZIPROOT=/sw/squeeze-x64/szip-latest-static/

set(NETCDF_DIR      "$ENV{NETCDFCROOT}")
set(NETCDF_DIR_F90  "$ENV{NETCDFFROOT}")

set(HDF5_LIB_1         "$ENV{HDF5ROOT}/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "$ENV{HDF5ROOT}/lib/libhdf5.a")
set(SZIP_LIB           "$ENV{SZIPROOT}/lib/libsz.a")
set(LIBS ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB})
