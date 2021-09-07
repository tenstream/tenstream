#!/bin/bash
set -eu -o pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT_ROOT=$(readlink -f "$SCRIPTDIR/../")

PYTHON=$(which "python3")
if [[ ! -e $PYTHON ]]; then
  echo "Could not find python3! ($PYTHON)"
  exit 1
else
  echo "Using python3: $PYTHON"
fi

BIN=$PROJECT_ROOT/external/fprettify/fprettify.py
if [[ ! -e $BIN ]]; then
  git clone git@github.com:pseewald/fprettify.git $(dirname $BIN)
fi
if [[ ! -e $BIN ]]; then
  echo "Could not find fprettify at $BIN"
  exit 2
else
  echo "Using fprettify: $BIN"
fi

cd $PROJECT_ROOT

FILES="\
  src/*.fypp \
  src/*.F90 \
  src/*.inc \
  rrtmg/rrtmg/*.F90 \
  plexrt/*.F90 \
  c_wrapper/*.F90 \
  examples/**/*.F90 \
  "

OPT="\
  --indent 2 \
  --line-length 132 \
  --whitespace 3 \
  --strict-indent \
  --enable-decl \
  --case 1 1 1 1 \
  "

NPROC=$(nproc)
count=0

for f in $FILES
do
  (
  echo "prettify => $f"
  $PYTHON $BIN $OPT $f
  ) &
  count=$(( $count + 1 ))
  if [ $count -eq $NPROC ]; then
    count=0
    wait
  fi
done
wait
echo done
