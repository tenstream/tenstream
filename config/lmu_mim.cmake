# Default GCC
set(CMAKE_Fortran_COMPILER "/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/bin/mpif90")
set(Fortran_COMPILER_WRAPPER "/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/bin/mpif90")

set(USER_C_FLAGS               "-cpp -W -Wall -Wuninitialized --std=c99") 
set(USER_Fortran_FLAGS         "-cpp -ffree-line-length-none -W -Wall -Wuninitialized -g") 
set(USER_Fortran_FLAGS_RELEASE "-fno-backtrace -fno-range-check -O3") 
set(USER_Fortran_FLAGS_DEBUG   "-fbacktrace -finit-real=nan -W -Wall -Wuninitialized -g -pg -fcheck=all -fbounds-check -pedantic -Wsurprising")

set(NETCDF_INCLUDE_DIR "/home/opt/cosmo_tica_lib//ompi1.8.1/netcdf-fortran-4.2/install/include/")
set(NETCDF_LIB_1       "/home/opt/cosmo_tica_lib//ompi1.8.1/netcdf-fortran-4.2/install/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/home/opt/cosmo_tica_lib//ompi1.8.1/netcdf-4.3.0/install/lib64/libnetcdf.a")

set(HDF5_INCLUDE_DIRS       "/home/opt/cosmo_tica_lib//ompi1.8.1/hdf5/HDF5-1.8.13-Linux/HDF_Group/HDF5/1.8.13/include")
list(APPEND HDF5_LIBRARIES  "/home/opt/cosmo_tica_lib//ompi1.8.1/hdf5/HDF5-1.8.13-Linux/HDF_Group/HDF5/1.8.13/lib/libhdf5_hl_fortran.so")
list(APPEND HDF5_LIBRARIES  "/home/opt/cosmo_tica_lib//ompi1.8.1/hdf5/HDF5-1.8.13-Linux/HDF_Group/HDF5/1.8.13/lib/libhdf5_fortran.so")

set(SZIP_LIB           "/home/opt/cosmo_tica_lib//ompi1.8.1/hdf5/HDF5-1.8.13-Linux/HDF_Group/HDF5/1.8.13/lib/libszip.so")

set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${SZIP_LIB} ${HDF5_LIBRARIES} m z curl)

