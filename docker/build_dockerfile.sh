#!/bin/bash

DOCKERBASEFILE=$1
BUILD_TYPE=$2
WORKDIR=$3
CC=$4
FC=$5
CXX=$6
DOCKER_TARGET=$7

cat $DOCKERBASEFILE > $DOCKER_TARGET

CURRENT_PETSC_HASH=$(git ls-remote https://gitlab.com/petsc/petsc | head -n 1 | awk '{print $1}')

export PETSC_DIR=$WORKDIR/petsc
export PETSC_ARCH=$BUILD_TYPE
echo "Using Work Dir :: $WORKDIR"
echo "Installing PETSc :: $PETSC_DIR // $PETSC_ARCH // $CC // $FC // $CXX"

#PETSC_OPT="--with-cc=$CC --with-fc=$FC --with-cxx=$CXX \
#  --with-fortran --with-fortran-interfaces \
#  --with-valgrind --download-hdf5 --download-zlib --download-openmpi"

PETSC_OPT="\
  --with-mpiexec=/bin/true \
  --with-fortran --with-fortran-interfaces \
  --with-valgrind --download-hdf5 --download-zlib"

if grep openmpi $DOCKERBASEFILE; then
  echo "Found MPI installed"
  PETSC_OPT="$PETSC_OPT --with-cc=$CC --with-fc=$FC --with-cxx=$CXX"
else
  echo "Cannot find mpirun ... will petsc install openmpi"
  PETSC_OPT="$PETSC_OPT --download-openmpi"
fi

CMAKE_BUILD_TYPE="RELEASE"
[[ $PETSC_ARCH = *"DEBUG"* ]] && CMAKE_BUILD_TYPE="DEBUG"
[[ $PETSC_ARCH = *"DEBUG"* ]] && PETSC_OPT="$PETSC_OPT --with-debugging=1"
[[ $PETSC_ARCH = *"int64"* ]] && PETSC_OPT="$PETSC_OPT --with-64-bit-indices=1"
[[ $PETSC_ARCH = *"single"* ]] && PETSC_OPT="$PETSC_OPT --with-precision=single"

echo "PETSc_Configure_Options: $PETSC_OPT"

NCFC=netcdf-c-4.7.0
NCFF=netcdf-fortran-4.4.5
export NETCDF_DIR=$PETSC_DIR/$PETSC_ARCH

cat >> $DOCKER_TARGET << EOF
RUN echo "export PETSC_DIR=$PETSC_DIR" >> $WORKDIR/.profile && \
    echo "export PETSC_ARCH=$PETSC_ARCH" >> $WORKDIR/.profile && \
    echo "export PATH=$PETSC_DIR/$PETSC_ARCH/bin:\$PATH" >> $WORKDIR/.profile && \
    echo "export NETCDF_DIR=$NETCDF_DIR" >> $WORKDIR/.profile && \
    echo "export CMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE" >> $WORKDIR/.profile && \
    cat $WORKDIR/.profile && \
  \
  cd $WORKDIR && . $WORKDIR/.profile && \
  git clone --depth=1 https://gitlab.com/petsc/petsc -b master \$PETSC_DIR && \
  cd \$PETSC_DIR && git checkout $CURRENT_PETSC_HASH && \
  ./configure $PETSC_OPT || (cat configure.log; false) && make && \
  \
  cd $WORKDIR && . $WORKDIR/.profile && \
  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/${NCFC}.tar.gz && \
  tar -xzf ${NCFC}.tar.gz && cd $NCFC && \
  CC=$CC CPPFLAGS=-I\$PETSC_DIR/\$PETSC_ARCH/include LDFLAGS=-L\$PETSC_DIR/\$PETSC_ARCH/lib ./configure --prefix=\$NETCDF_DIR && \
  make -j install && \
  \
  cd $WORKDIR && . $WORKDIR/.profile && \
  wget ftp://ftp.unidata.ucar.edu/pub/netcdf/${NCFF}.tar.gz && \
  tar -xzf ${NCFF}.tar.gz && cd $NCFF && \
  CC=$CC FC=$FC F70=$FC CPPFLAGS=-I\$PETSC_DIR/\$PETSC_ARCH/include LDFLAGS=-L\$PETSC_DIR/\$PETSC_ARCH/lib ./configure --prefix=\$NETCDF_DIR && \
  make -j install
EOF
