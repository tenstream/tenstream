set -euo pipefail

. config.sh

MIE_WVL=$(python -c "import numpy as np; \
WVLLOW, WVLHIGH = 250, 200e3; \
wvl2wvnr = lambda wvl: 1e7 / wvl; \
WVNHI, WVNLO = [ wvl2wvnr(_) for _ in (WVLLOW, WVLHIGH) ]; \
wvl = wvl2wvnr(np.linspace(WVNLO, WVNHI, 30));
print(' '.join(map(str, wvl))); \
")

REFF_WC=$(seq -s " " 2 .1 40)

if [ ! -e $MIE_GENERAL ]; then
  cd $WORK
${MPIEXEC:-} python $SCRIPTS/gen_mie_tables.py \
  $MIE_GENERAL \
  -p water \
  --wvls $MIE_WVL \
  --reffs $REFF_WC \
  -v \
  -m $LIBRADTRAN_ROOT/bin/mie
else
  echo "Skipping computation of $MIE_GENERAL because file already exists"
fi

if [ ! -e $FU_GENERAL ]; then
  python $SCRIPTS/gen_fu_ice_tables.py $LIBRADTRAN_ROOT -v $FU_GENERAL
else
  echo "Skipping computation of $FU_GENERAL because file already exists"
fi
