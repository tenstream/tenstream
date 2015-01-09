# Default GCC
#
# Config script for the linux environment at the Meteorological Institute Munich 08.Jan.2015
#
# There may be a maintained(Fabian Jakub) library structure at /home/opt/cosmo_tica_lib/
# if so, you might get away with copying the following into your .bashrc:
# 
#> TICAP=/home/opt/cosmo_tica_lib/
#> 
#> OMPI=1.8.1
#> export OMPIDIR=$TICAP/ompi$OMPI/openmpi-$OMPI/install/
#> 
#> if [ -n "$OMPIDIR" ]
#>     then # echo "setting OpenMPI environment to: $OMPIDIR"
#>             export LD_RUN_PATH=$OMPIDIR/lib64/:$LD_RUN_PATH
#>             export PATH=$OMPIDIR/bin/:$PATH
#>             export LD_LIBRARY_PATH=$OMPIDIR/lib64/:$LD_LIBRARY_PATH
#>             export MANPATH=$OMPIDIR/share/man:$MANPATH
#>     fi      
#>             
#> 
#> #PETSC 
#> export PETSC_DIR=$TICAP/ompi$OMPI/petsc-maint
#> export PETSc_DIR=$PETSC_DIR
#> export PETSC_ARCH=fast-single
#> 
#> # NETCDF
#> export NETCDFFROOT=$MPI/netcdf-fortran-4.2/install
#> export NETCDFROOT=$MPI/netcdf-4.3.0/install
#> export PATH=$NETCDFFROOT/bin:$NETCDFROOT/bin:$PATH
#> export LD_LIBRARY_PATH=$NETCDFFROOT/lib64:$NETCDFROOT/lib64:$LD_LIBRARY_PATH
#
#
#
# and use this config file with cmake <tenstream_root_dir> -DSYST:STRING=lmu_mim
#

set(CMAKE_Fortran_COMPILER "/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/bin/mpif90")
set(Fortran_COMPILER_WRAPPER "/home/opt/cosmo_tica_lib/ompi1.8.1/openmpi-1.8.1/install/bin/mpif90")

set(USER_C_FLAGS               "-cpp -W -std=c99") 
set(USER_Fortran_FLAGS         "-cpp -ffree-line-length-none -g") 
set(USER_Fortran_FLAGS_RELEASE "-fno-backtrace -fno-range-check -O3") 
set(USER_Fortran_FLAGS_DEBUG   "-fbacktrace -finit-real=nan -W -Wall -Wuninitialized -g -pg -fcheck=all -fbounds-check -pedantic -Wsurprising")

set(NETCDF_DIR      "/home/opt/cosmo_tica_lib//ompi1.8.1/netcdf-c-4.3.3-rc2/build/install/")
set(NETCDF_DIR_F90  "/home/opt/cosmo_tica_lib//ompi1.8.1/netcdf-fortran-4.4.1/build/install/")
set(HDF_ROOT        "/home/opt/cosmo_tica_lib//ompi1.8.1/hdf5/HDF5-1.8.13-Linux/HDF_Group/HDF5/1.8.13/")
