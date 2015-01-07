# Default GCC
set(CMAKE_Fortran_COMPILER "/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/bin/mpif90")
set(Fortran_COMPILER_WRAPPER "/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/bin/mpif90")

set(USER_C_FLAGS               "-cpp -W -Wall -Wuninitialized --std=c99") 
set(USER_Fortran_FLAGS         "-cpp -ffree-line-length-none -W -Wall -Wuninitialized -g") 
set(USER_Fortran_FLAGS_RELEASE "-fno-backtrace -fno-range-check -O3") 
set(USER_Fortran_FLAGS_DEBUG   "-fbacktrace -finit-real=nan -W -Wall -Wuninitialized -g -pg -fcheck=all -fbounds-check -pedantic -Wsurprising")

set(NETCDF_INCLUDE_DIR "/home/opt/cosmo_tica_lib//ompi1.8.1/netcdf-fortran-4.2/install/include/")
include_directories (${NETCDF_INCLUDE_DIR} )
set(NETCDF_INCLUDE_DIR "/home/opt/cosmo_tica_lib//ompi1.8.1/netcdf-4.3.0/install/include/")
include_directories (${NETCDF_INCLUDE_DIR} )

set(NETCDF_LIB_1       "/home/opt/cosmo_tica_lib//ompi1.8.1/netcdf-fortran-4.2/install/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/home/opt/cosmo_tica_lib//ompi1.8.1/netcdf-4.3.0/install/lib64/libnetcdf.a")

set(USERLIB ${NETCDF_LIB_1} ${NETCDF_LIB_2} m z curl)

