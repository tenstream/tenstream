#!/bin/bash
set -eu -o pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
PROJECT_ROOT="$SCRIPTDIR/../"
DST=$PROJECT_ROOT/misc/LUT

MASK="*${1:-LUT}*"
echo "Looking for Tables with mask: $MASK and saving files to $DST"

mkdir -p $DST
cd $DST

WGET_BIN=$(command -v wget || true)
if [ -z $WGET_BIN ]; then
  echo "Could not find wget but I need it to download files"
  echo "Make sure we have wget available"
  exit 1
fi

BASEURL="https://www.meteo.physik.uni-muenchen.de/~Fabian.Jakub/TenstreamLUT/"
BIN="$WGET_BIN -q --show-progress --progress=bar:force:noscroll -c -r --no-parent -nH --cut-dirs=3 -A "

echo "Downloading LUT's for: ${MASK}"
$BIN $MASK $BASEURL

echo ""
echo "For convenience you may add the location of the luts to your global petscrc, i.e. run:"
echo ""
echo "  echo -lut_basename $DST/LUT >> $HOME/.petscrc"
echo ""
