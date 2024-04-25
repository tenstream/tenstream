#!/bin/bash
set -euo pipefail

. $HOME/.profile
py local

BASE=bmcdb_ex_ey
rm -f $BASE.nc $BASE.annotated.nc

./client.py -d --Ndir 3 -o $BASE.nc
./annotate.py -v -s -r -tf .05 $BASE.nc $BASE.annotated.nc
