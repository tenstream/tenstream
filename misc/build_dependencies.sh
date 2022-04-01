#!/bin/bash
set -euo pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT_ROOT="$(readlink -f $SCRIPTDIR/../)"

export PETSC_DIR=${PETSC_DIR:-${1:-$PROJECT_ROOT/external/petsc}}
export PETSC_ARCH=${PETSC_ARCH:-${2:-default}}
PETSC_PRECISION=${3:-single}
PETSC_DEBUGGING=${4:-0}
PETSC_64_INTEGERS=${5:-0}
PETSC_OPTS=${@:6}

if [[ "x$PETSC_ARCH" == x"debug"* ]]; then
  PETSC_DEBUGGING=1
  PETSC_PRECISION=double
fi

if [[ -z ${PETSC_OPTS} ]]; then
  PETSC_OPTS="--download-hdf5 --download-szlib --download-zlib"
fi

echo ""
echo "** Downloading Petsc, and netcdf libs and install them"

if [[ -z ${FC:-} ]]; then
  echo " Need to define FC, e.g. set with FC=mpif90"
  exit 1
fi
if [[ -z ${CC:-} ]]; then
  echo " Need to define CC, e.g. set with CC=mpicc"
  exit 2
fi
if [[ -z ${CXX:-} ]]; then
  echo " Need to define CXX, e.g. set with CXX=mpicxx"
  exit 3
fi

printf "Using:\n\
  PETSC_DIR:          $PETSC_DIR \n\
  PETSC_ARCH          $PETSC_ARCH \n\
  PETSC_PRECISION     $PETSC_PRECISION \n\
  PETSC_DEBUGGING     $PETSC_DEBUGGING \n\
  PETSC_64_INTEGERS   $PETSC_64_INTEGERS \n\
  PETSC_OPTS          $PETSC_OPTS \n\

  C-Compiler:   ${CC}\n\
  F-Compiler:   ${FC}\n\
  C++ Compiler: ${CXX}\n\
  \n"

PETSC_URL=https://gitlab.com/petsc/petsc.git
PETSC_BRANCH=main

if [ -e "$PETSC_DIR" ]
then
  echo "Using PETSC_DIR: $PETSC_DIR"
else
  git clone $PETSC_URL -b $PETSC_BRANCH $PETSC_DIR
fi

PETSC_OPTIONS="\
  --with-cc=${CC} \
  --with-cxx=${CXX} \
  --with-fc=${FC} \
  --COPTFLAGS=\"-g\" \
  --CXXOPTFLAGS=\"-g\" \
  --FOPTFLAGS=\"-g\" \
  --with-debugging=$PETSC_DEBUGGING \
  --with-precision=$PETSC_PRECISION \
  --with-64-bit-indices=$PETSC_64_INTEGERS \
  $PETSC_OPTS \
  "

cd $PETSC_DIR
COPT="$PETSC_OPTIONS PETSC_DIR=$(pwd)"
echo "Running configure with: $COPT"
./configure $COPT
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
  CURL_BIN=$(command -v curl || true)
  if [ -z $CURL_BIN ]; then
    echo "Could not find curl but I need it to download a file"
    echo "Either download it yourself or make sure we have curl available"
    exit 1
  fi
  curl -L --progress-bar $URL --output $DST
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

  URL="https://github.com/Unidata/$FILE"
  download_file "$URL" "$SRCDIR/$(basename $FILE)"
  tar -xzf $(basename $FILE)
  cd $(basename $FILE .tar.gz)
  export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH:-}
  export PATH=${PREFIX}/bin:${PATH:-}
  CC=$CC FC=$FC F90=$FC CXX=$CXX CPPFLAGS=-I$PREFIX/include LDFLAGS=-L$PREFIX/lib ./configure --prefix=$PREFIX --disable-dap && make -j install
  echo "Installed NetCDF lib $FILE into $PREFIX -- CC $CC FC $FC CXX $CXX"
}

install_netcdf "netcdf-c/archive/refs/tags/v4.8.1.tar.gz"       "$PETSC_DIR/$PETSC_ARCH/"
install_netcdf "netcdf-fortran/archive/refs/tags/v4.5.4.tar.gz" "$PETSC_DIR/$PETSC_ARCH/"

printf "\n** Make sure to export PETSC_DIR and PETSC_ARCH before cmake'ing TenStream, i.e. set \n\
  \n\
  export PETSC_DIR=$PETSC_DIR\n\
  export PETSC_ARCH=$PETSC_ARCH\n\
  \n"
