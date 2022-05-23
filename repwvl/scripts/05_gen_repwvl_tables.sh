set -euo pipefail

. config.sh

REFF_WC=$(seq -s " " 2 .5 40)

if [ ! -e $MIE_REPWVL ]; then
  READWVL_PY="import xarray as xr; import numpy as np; import sys;\
    wvl = xr.open_dataset(sys.argv[1]).wvl.data; \
    print(' '.join(map(str, wvl))); \
    "
  MIE_WVL_SW=$(python -c "$READWVL_PY" $REPWVL_SW)
  MIE_WVL_LW=$(python -c "$READWVL_PY" $REPWVL_LW)

  ${MPIEXEC:-} python $SCRIPTS/gen_mie_tables.py \
    $MIE_REPWVL \
    -p water \
    --wvls $MIE_WVL_SW $MIE_WVL_LW \
    --reffs $REFF_WC \
    --wdir $WORK \
    -v \
    -m $LIBRADTRAN_ROOT/bin/mie
else
  echo "Skipping computation of $MIE_REPWVL because file already exists"
fi

if [ ! -e $FU_REPWVL ]; then
  python $SCRIPTS/gen_fu_ice_tables.py \
    $LIBRADTRAN_ROOT \
    $FU_REPWVL \
    --verbose \
    --repwvl_file $REPWVL_SW

  python $SCRIPTS/gen_fu_ice_tables.py \
    $LIBRADTRAN_ROOT \
    $FU_REPWVL \
    --verbose \
    --repwvl_file $REPWVL_LW
else
  echo "Skipping computation of $FU_REPWVL because file already exists"
fi
