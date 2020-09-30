#!/bin/bash
set -euo pipefail

BIN=ex_pprts_rrtm_buildings

make -j -C $TENSTREAMROOT $BIN

ATM=$TENSTREAMROOT/generated/test_plexrt_rrtmg_lw_sw/atm.dat

BASEOPT="-atm_filename $ATM -thermal no -phi -55 -theta 40 -Nz 12 -BAg .25"


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
    $MPIEXEC $TENSTREAMROOT/bin/$BIN $BASEOPT -out $OUT $OPTS || true
    if [ -e $SNAP ] ; then
      mv $SNAP snaps_$OUT
    fi
  fi
}

NP=1; NC=1
run_ex $NP $NC "tenstr_3_10.nc" "-solver 3_10 ${1:-}"
run_ex $NP $NC "tenstr_8_10.nc" "-solver 8_10 ${1:-}"
run_ex $NP $NC "twostr.nc" "-twostr_only ${1:-}"

NP=1; NC=128
run_ex $NP $NC "rayli.nc" "-rayli_opts -pprts_rayli_photons 1e6 -rayli_nthreads $(($NP * $NC))  ${1:-}"

cat > compare_buildings_benchmark.py <<EOF
import xarray as xr
import numpy as np

import glob
files = sorted(glob.glob("$PREFIX*.nc"))

ftruth = "${PREFIX}rayli.nc"

Dtruth = xr.open_dataset(ftruth)

dir = lambda D: D['rank0_buildings_edir'].data
inc = lambda D: D['rank0_buildings_edir'].data + D['rank0_buildings_incoming'].data
out = lambda D: D['rank0_buildings_outgoing'].data

rmse = lambda a,b: np.sqrt(np.mean((a-b)**2))
rrmse = lambda a,b: rmse(a,b)/np.mean(b) * 100

print("")
print("|{:40s} | {:>10s} | {:>10s} | {:>10s}".format("Mean abs. deviation [W/m2]", "Edir", "Eincoming", "Eoutgoing"))
print("|---|---:|---:|---:")
for f in files:
  fn = f.replace('$PREFIX','').replace('.nc','')
  with xr.open_dataset(f) as D:
    err = [
      np.mean(np.abs(dir(D) - dir(Dtruth))),
      np.mean(np.abs(inc(D) - inc(Dtruth))),
      np.mean(np.abs(out(D) - out(Dtruth))),
    ]
    print("|{:40s} | {:10.2f} | {:10.2f} | {:10.2f}".format(fn, err[0], err[1], err[2]))
print("|<img width=300/>|||")

print("")
print("|{:40s} | {:>10s} | {:>10s} | {:>10s}".format("RMSE [W/m2]", "Edir", "Eincoming", "Eoutgoing"))
print("|---|---:|---:|---:")
for f in files:
  fn = f.replace('$PREFIX','').replace('.nc','')
  with xr.open_dataset(f) as D:
    err = [
      rmse(dir(D), dir(Dtruth)),
      rmse(inc(D), inc(Dtruth)),
      rmse(out(D), out(Dtruth)),
    ]
    print("|{:40s} | {:10.2f} | {:10.2f} | {:10.2f}".format(fn, err[0], err[1], err[2]))
print("|<img width=300/>|||")

print("")
print("|{:40s} | {:>10s} | {:>10s} | {:>10s}".format("Rel. RMSE [%]", "Edir", "Eincoming", "Eoutgoing"))
print("|---|---:|---:|---:")
for f in files:
  fn = f.replace('$PREFIX','').replace('.nc','')
  with xr.open_dataset(f) as D:
    err = [
      rrmse(dir(D), dir(Dtruth)),
      rrmse(inc(D), inc(Dtruth)),
      rrmse(out(D), out(Dtruth)),
    ]
    print("|{:40s} | {:10.2f} | {:10.2f} | {:10.2f}".format(fn, err[0], err[1], err[2]))
print("|<img width=300/>|||")


print("")
D10 = xr.open_dataset('buildings_benchmark_tenstr_3_10.nc')
D1  = xr.open_dataset('buildings_benchmark_twostr.nc')
for a,b,c in zip(dir(Dtruth), dir(D10), dir(D1)):
  print("dir", a,b,c)
EOF

cat > imshow_buildings_benchmark.py <<EOF
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
D=xr.open_dataset('snaps_buildings_benchmark_rayli.nc')
im = np.mean([ v.data for v in  D.values()] , axis=0)

plt.imshow (im, origin='lower', cmap='Spectral_r')
plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
plt.savefig('pyramid_snap.pdf', bbox_inches='tight')
EOF
