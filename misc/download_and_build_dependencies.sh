#!/bin/bash
set -euo pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT_ROOT="$SCRIPTDIR/../"

echo ""
echo "** Downloading Petsc, and netcdf libs and install them in ${PROJECT_ROOT}/external"
echo ""

if [ -z ${FC:-} ]; then
  echo " Need to define FC, e.g. set with FC=mpif90"
  exit 1
fi
if [ -z ${CC:-} ]; then
  echo " Need to define CC, e.g. set with CC=mpicc"
  exit 2
fi
if [ -z ${CXX:-} ]; then
  echo " Need to define CXX, e.g. set with CXX=mpicxx"
  exit 3
fi

printf "Using Compilers:\n\
  C-Compiler:   ${CC}\n\
  F-Compiler:   ${FC}\n\
  C++ Compiler: ${CXX}\n"

PETSC_URL=https://gitlab.com/petsc/petsc.git
PETSC_BRANCH=master

PETSC_DIR=$PROJECT_ROOT/external/petsc
if [ -e "$PETSC_DIR" ]
then
  echo "Using PETSC_DIR: $PETSC_DIR"
else
  git clone $PETSC_URL -b $PETSC_BRANCH $PETSC_DIR
fi

PETSC_ARCH=default
PETSC_OPTIONS="\
  --with-debugging=0 \
  --download-hdf5 \
  --download-szlib \
  --download-zlib \
  "

echo "Using Petsc Options: $PETSC_OPTIONS"
cd $PETSC_DIR
./configure CC=$CC FC=$FC F90=$FC CXX=$CXX $PETSC_OPTIONS
make


# Download and install netcdf libs:
function download_file() {
  URL=$1
  DST=$2
  if [ -e $DST ]; then
    echo "Skipping Download '$URL' to '$DST' because it already exists"
    return
  fi
  echo "Download '$URL' to '$DST'"
  WGET_BIN=$(command -v wget || true)
  if [ -z $WGET_BIN ]; then
    echo "Could not find wget but I need it to download a file"
    echo "Either download it yourself or make sure we have wget available"
    exit 1
  fi
  wget $URL -O $DST
  }

function install_netcdf() {
  FILE=$1
  PREFIX=$2
  SRCDIR=$PREFIX/externalpackages
  mkdir -p $SRCDIR
  cd $SRCDIR

  if [ -e "$(basename $FILE .tar.gz)" ]; then
    echo "Skipping netcdf install for $FILE because src directory already present"
    return
  fi

  URL="ftp://ftp.unidata.ucar.edu/pub/netcdf/$FILE"
  download_file "$URL" "$SRCDIR/$FILE"
  tar -xzf $FILE
  cd $(basename $FILE .tar.gz)
  export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH:-}
  export PATH=${PREFIX}/bin:${PATH:-}
  CC=$CC FC=$FC F90=$FC CXX=$CXX CPPFLAGS=-I$PREFIX/include LDFLAGS=-L$PREFIX/lib ./configure --prefix=$PREFIX && make -j install
  echo "Installed NetCDF lib $FILE into $PREFIX -- CC $CC FC $FC CXX $CXX"
}

install_netcdf "netcdf-c-4.7.4.tar.gz"       "$PETSC_DIR/$PETSC_ARCH/"
install_netcdf "netcdf-fortran-4.5.3.tar.gz" "$PETSC_DIR/$PETSC_ARCH/"
install_netcdf "netcdf-cxx4-4.3.1.tar.gz"    "$PETSC_DIR/$PETSC_ARCH/"

printf "\n** Make sure to export PETSC_DIR and PETSC_ARCH before cmake'ing TenStream, i.e. set \n\
  \n\
  export PETSC_DIR=$PETSC_DIR\n\
  export PETSC_ARCH=$PETSC_ARCH\n\
  \n"
