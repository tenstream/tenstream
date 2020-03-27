#!/bin/bash
set -eu -o pipefail

MASK="*${1:-'LUT'}*"

[ $(which wget) ] || echo "Could not find wget"

BASEURL="https://www.meteo.physik.uni-muenchen.de/~Fabian.Jakub/TenstreamLUT/"
BIN="$(which wget) -c -r --no-parent -nH --cut-dirs=3 -A"

echo "Downloading LUT's for: ${MASK}"
$BIN $MASK $BASEURL
