set -euo pipefail

BIN=ex_uclales_cld_file
CLDDIR=$SCRATCH/compare_specint_uclales/
CLD=$CLDDIR/phillip1_0d_25m_srfc.strided_4.nc; OUTPREFIX="$WORK/strided4"
CLD=$CLDDIR/phillip1_0d_25m_srfc.strided_8.nc; OUTPREFIX="$WORK/strided8"
CLD=$CLDDIR/phillip1_0d_25m_srfc.merged.nc;    OUTPREFIX="$WORK/full"

NP=64
export MPIEXEC=${MPIEXEC:-srun -n $NP -N$(($NP / 32)) -c1 --mem=180G --time=48:00:00 --exclusive --cpu-bind=cores,verbose --mpi=pmi2}

mkdir -p $OUTPREFIX

TSTART=800
TEND=900

BASEOPT="\
  -specint repwvl \
  -cld $CLD \
  -dx 25 \
  -dy 25 \
  -twostr_ratio 3 \
  -log_view \
  -tod_offset .25 \
  -pprts_view_suninfo \
  "
