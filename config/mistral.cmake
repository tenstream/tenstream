# Mistral supercomputer at DKRZ Hamburg
# modules loaded:

# module load intel/18.0.4 openmpi/2.0.2p2_hpcx-intel14 cmake
#
# export OMPI_MCA_pml=cm         # sets the point-to-point management layer
# export OMPI_MCA_mtl=mxm        # sets the matching transport layer (MPI-2 one-sided comm.)
# export MXM_RDMA_PORTS=mlx5_0:1
# export MXM_LOG_LEVEL=ERROR
# export MXM_HANDLE_ERRORS=bt
# export UCX_HANDLE_ERRORS=bt
#
# # enable HCOLL based collectives
# export OMPI_MCA_coll=^fca              # disable FCA for collective MPI routines
# export OMPI_MCA_coll_hcoll_enable=1    # enable HCOLL for collective MPI routines
# export OMPI_MCA_coll_hcoll_priority=95
# export OMPI_MCA_coll_hcoll_np=8        # use HCOLL for all communications with more than 8 tasks
# export HCOLL_MAIN_IB=mlx5_0:1
# export HCOLL_ENABLE_MCAST=1
# export HCOLL_ENABLE_MCAST_ALL=1
#
# # disable specific HCOLL functions (strongly depends on the application)
# export HCOLL_ML_DISABLE_BARRIER=1
# export HCOLL_ML_DISABLE_IBARRIER=1
# export HCOLL_ML_DISABLE_BCAST=1
# export HCOLL_ML_DISABLE_REDUCE=1
#
# export MPIMODE=openmpi-bull
# export MYmpif90=mpiifort
# export MYmpicc=mpiicc
# export MYmpicxx=mpiicpc
#
# export NETCDFCROOT=$PETSC_DIR/$PETSC_ARCH
# export NETCDFFROOT=$PETSC_DIR/$PETSC_ARCH
# export HDF5ROOT=$PETSC_DIR/$PETSC_ARCH

# Petsc is installed with:
# export PETSC_ARCH=fast-$MPIMODE
# NETCDFFILE="netcdf-c-4.6.2.tar.gz"
# if [ ! -e $NETCDFFILE ]
# then
#   wget ftp://ftp.unidata.ucar.edu/pub/netcdf/$NETCDFFILE
# fi
# ./configure                          \
#  --with-make-np=10                   \
#  --with-cc=$(which  $MYmpicc)        \
#  --with-fc=$(which  $MYmpif90)       \
#  --with-cxx=$(which $MYmpicxx)       \
#  --with-shared-libraries=1           \
#  --with-blas-lapack-dir=$MKL         \
#  --with-fortran                      \
#  --with-fortran-interfaces           \
#  --with-cmake=$(which cmake)         \
#  --with-precision=single             \
#  --download-hdf5                     \
#  --download-netcdf=$NETCDFFILE       \
#  --download-metis                    \
#  --download-parmetis                 \
#  --download-ptscotch                 \
#  --with-debugging=0                  \
#  --with-valgrind-dir=/usr            \
#  COPTFLAGS='-fPIC -O3 -xCORE-AVX2 -mkl' \
#  FOPTFLAGS='-fPIC -O3 -xCORE-AVX2 -mkl' \
#  && make all
#
# And additionally install netcdf-fortran into the PETSC_ARCH dir
#
# NETCDFF90FILE="netcdf-fortran-4.4.4.tar.gz"
# if [ ! -e $NETCDFF90FILE ]
# then
#   wget ftp://ftp.unidata.ucar.edu/pub/netcdf/$NETCDFF90FILE
# fi
# cd $PETSC_DIR/$PETSC_ARCH/externalpackages/
# tar -xzf $PETSC_DIR/$NETCDFF90FILE
# cd $(basename $NETCDFF90FILE .tar.gz)
# export LD_LIBRARY_PATH=${PETSC_DIR}/${PETSC_ARCH}/lib:${LD_LIBRARY_PATH}
# export LIBRARY_PATH=${PETSC_DIR}/${PETSC_ARCH}/lib:${LIBRARY_PATH}
# export LD_RUN_PATH=${PETSC_DIR}/${PETSC_ARCH}/lib:${LD_RUN_PATH}
# export PATH=${PETSC_DIR}/${PETSC_ARCH}/bin/:${PATH}
# autoconf
# CC=$(which $MYmpicc) FC=$(which $MYmpif90) F90=$(which $MYmpif90) CPPFLAGS=-I$PETSC_DIR/$PETSC_ARCH/include LDFLAGS=-L$PETSC_DIR/$PETSC_ARCH/lib \
# ./configure --prefix=$PETSC_DIR/$PETSC_ARCH && make -j install
#

set(NETCDF_DIR      "$ENV{NETCDFCROOT}")
set(NETCDF_DIR_F90  "$ENV{NETCDFFROOT}")
set(HDF5_DIR        "$ENV{HDF5ROOT}")

set(CMAKE_C_COMPILER       "$ENV{MYmpicc}") # see definition of compiler vars above in comments
set(CMAKE_CXX_COMPILER     "$ENV{MYmpicxx}")
set(CMAKE_Fortran_COMPILER "$ENV{MYmpif90}")

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

set(USER_C_FLAGS "-Wall -std=c99")
set(USER_Fortran_FLAGS "-warn all -cpp -extend_source -g -sox -no-wrap-margin -fp-model source -mkl -mkl=sequential")
set(USER_Fortran_FLAGS_RELEASE "-march=native -O3 -ftz -pc64 -xCORE-AVX2 -fp-model fast=2 -no-prec-div -no-prec-sqrt -fast-transcendentals")
set(USER_Fortran_FLAGS_DEBUG "-traceback -warn error -ftrapuv -fpe0 -O2 -g -check all -check nopointers -check noarg_temp_created")
