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
#> export MPI=$TICAP/ompi$OMPI/
#> export MPIDIR=$MPI/openmpi-$OMPI/install/
#> 
#> if [ -n "$MPIDIR" ]
#> then 
#>   #echo "setting OpenMPI environment to: $OMPIDIR"
#>   export PATH=$MPIDIR/bin/:$PATH
#>   export LD_RUN_PATH=$MPIDIR/lib64/:$LD_RUN_PATH
#>   export LD_LIBRARY_PATH=$MPIDIR/lib64/:$LD_LIBRARY_PATH
#>   export MANPATH=$MPIDIR/share/man:$MANPATH
#> fi
#>             
#> # PETSC
#> export PETSC_DIR=$MPI/petsc
#> export PETSC_ARCH=fast_single
#> 
#> # NETCDF
#> export NETCDFFROOT=$MPI/netcdf-latest/
#> export NETCDFCROOT=$MPI/netcdf-latest/
#> export NETCDFC_DIR=$NETCDFCROOT
#> export NETCDFF_DIR=$NETCDFFROOT
#> export PATH=$NETCDFFROOT/bin:$NETCDFCROOT/bin:$PATH
#> export LD_LIBRARY_PATH=$NETCDFFROOT/lib64:$NETCDFCROOT/lib64:$LD_LIBRARY_PATH
#> export LD_RUN_PATH=$NETCDFFROOT/lib64:$NETCDFCROOT/lib64:$LD_RUN_PATH
#> 
#> # HDF5
#> export HDF5_DIR=$MPI/hdf5/HDF5-1.8.13-Linux/HDF_Group/HDF5/1.8.13/
#> export HDF5_ROOT=$HDF5_DIR
#> export PATH=$HDF5_DIR/bin:$PATH
#> export LD_LIBRARY_PATH=$HDF5_DIR/lib:$LD_LIBRARY_PATH
#> export LD_RUN_PATH=$HDF5_DIR/lib:$LD_RUN_PATH

# and use this config file with `cmake <tenstream_root_dir> -DSYST:STRING=lmu_mim`

set(CMAKE_C_COMPILER   "mpiicc")
set(CMAKE_Fortran_COMPILER   "mpiifort")
set(Fortran_COMPILER_WRAPPER "mpiifort")

set(USER_C_FLAGS       "-cpp -W -std=c99 ") 
set(USER_Fortran_FLAGS "-cpp -traceback -extend_source -g -sox -no-wrap-margin ")
set(USER_Fortran_FLAGS_RELEASE "-O3 -fp-model source ")
set(USER_Fortran_FLAGS_DEBUG "-warn -fpe0 -O2 -g -check all -check nopointers -check noarg_temp_created ")

set(NETCDF_DIR      "$ENV{NETCDF}")
set(NETCDF_DIR_F90  "$ENV{NETCDF}")
#set(HDF5_DIR        "$ENV{HDF5ROOT}")
