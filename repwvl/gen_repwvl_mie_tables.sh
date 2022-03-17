#!/bin/bash
set -euo pipefail

SOLARDB=${1:-missing 1st argument}
THERMALDB=${2:-missing 2nd argument}

if [ ! -e "$SOLARDB"   ]; then echo "Solar database file (arg 1) does not exist: $SOLARDB"; exit 1; fi
if [ ! -e "$THERMALDB" ]; then echo "Thermal database file (arg 2) does not exist: $THERMALDB"; exit 1; fi

WVLVAR=ChosenWvls

ALL_WVL=
for DB in $SOLARDB $THERMALDB; do
  WVL=$(python -c "import xarray as xr; D=xr.open_dataset('$DB'); print(' '.join(map(str, D['$WVLVAR'].data)))")
  ALL_WVL="$WVL $ALL_WVL"
done

echo "Combined wavelengths are $ALL_WVL"

REFF_WC=$(seq -s " " 2 .1 40)
REFF_IC=$(seq -s " " 5 1 140)

WC=mie.wc.table.nc
if [ ! -e $WC ]; then
${MPIEXEC:-} python gen_mie_tables.py $WC -p water \
  --wvls $ALL_WVL \
  --reffs $REFF_WC \
  -v -m $WORK/lib/libRadtran/bin/mie
else
  echo "Skipping computation of $WC because file already exists"
fi

IC=mie.ic.table.nc
if [ ! -e $IC ]; then
${MPIEXEC:-} python gen_mie_tables.py $IC -p ice \
  --wvls $ALL_WVL \
  --reffs $REFF_IC \
  -v -m $WORK/lib/libRadtran/bin/mie
else
  echo "Skipping computation of $IC because file already exists"
fi

echo done
