#!/bin/bash
set -euo pipefail

BIN=ex_pprts_buildings

make -j -C $TENSTREAMROOT $BIN

BASEOPT="\
  -dtau 0 \
  -Ag=0.25 \
  -thermal no \
  -phi 225 \
  -theta 40 \
  -Nz 21 \
  -Nx 40 \
  -Ny 40 \
  -BNz 20 \
  -BNx 10 \
  -BNy 10 \
  -dx 1 \
  -dy 1 \
  -dz 1 \
  -BAg .5\
  -solar_dir_ksp_atol 1e-60 \
  -solar_dir_ksp_rtol 1e-7 \
  -solar_diff_ksp_atol 1e-60 \
  -solar_diff_ksp_rtol 1e-7 \
  "

PREFIX="buildings_benchmark_"
function run_ex() {
  NP=$1
  NC=$2
  OUT=${PREFIX}$3
  OPTS=$4

  MPIEXEC="srun --mpi=pmix -n $NP -c $NC"
  SNAP=rayli_snaphots.nc
  if [ ! -e $OUT ]; then
    rm -f $SNAP
    $MPIEXEC $TENSTREAMROOT/bin/$BIN $BASEOPT \
      -out $OUT \
      -pprts_xdmf $OUT.srfc \
      -pprts_xdmf_buildings $OUT.bldgs \
      $OPTS || true
    if [ -e $SNAP ] ; then
      mv $SNAP snaps_$OUT
    fi
  else
    echo "Skipping $OUT"
  fi
}

NP=1; NC=1
run_ex $NP $NC "tenstr_3_10.nc" "-solver 3_10 ${1:-}"
run_ex $NP $NC "tenstr_3_24.nc" "-solver 3_24 ${1:-}"
run_ex $NP $NC "tenstr_8_10.nc" "-solver 8_10 ${1:-}"
run_ex $NP $NC "tenstr_8_16.nc" "-solver 8_16 ${1:-}"
run_ex $NP $NC "twostr.nc" "-twostr_only ${1:-}"

NP=1; NC=64
ROPT="\
  -pprts_use_rayli \
  -rayli_cyclic_bc \
  -rayli_diff_flx_origin 0,0,-inf \
  -pprts_rayli_photons 1e8 \
  -rayli_nthreads $(($NP * $NC)) \
  ${1:-}\
  "
run_ex $NP $NC "rayli.nc" "$ROPT"


cat > benchmark_buildings.py <<EOF
from pylab import *
import glob
import xml.etree.ElementTree as ET

rmse = lambda a,b: np.sqrt(np.mean((a-b)**2))
rrmse = lambda a,b: rmse(a,b) / np.mean(b)

idents = ['rayli', 'tenstr_3_10', 'tenstr_3_24', 'tenstr_8_10', 'tenstr_8_16', 'twostr']
modes = ['bldgs', 'srfc']
r={}
for ident in idents:
    for mode in modes:
        fname = glob.glob(f'buildings_benchmark_{ident}.nc.{mode}*xmf')[0]
        for a in ET.parse(fname).iter('Attribute'):
            k = '.'.join((ident,mode,a.attrib['Name']))
            d = np.array(a[0].text.split()).astype(float)
            r[k] = d


vars = ['bldgs.edir', 'bldgs.incoming', 'bldgs.outgoing', 'srfc.edn']
for var in vars:
    print('-')
    for ident in idents:
        err = rrmse(r[f'{ident}.{var}'], r[f'rayli.{var}'])
        print(f'{var:20s} {ident:20s}: {err*100:12.1f}')

EOF
