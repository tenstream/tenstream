#!/bin/bash

set -euo pipefail

BIN=ex_uclales_cld_file
CLD=/project/meteo/public/jakub/phillip1_0d_25m_srfc/phillip1_0d_25m_srfc.merged.nc
CLD=/project/meteo/public/jakub/phillip1_0d_25m_srfc/phillip1_0d_25m_srfc.strided_8.nc

make -j -C $TENSTREAMROOT $BIN

NP=8

MPIEXEC=${MPIEXEC:-srun -n $NP -N1 -c1 --mem=120G --time=08:00:00 --exclusive --cpu-bind=cores,verbose}

TSTART=800
TEND=900
TINC=1

DT=10

cat > plot_timings.py << EOF
import re
import matplotlib.pyplot as plt
import sys
import os

outfile = sys.argv[1]
files = sys.argv[2:]

if os.path.exists(outfile):
    raise Exception(f"outfile {outfile} already exists")

res = []

for f in sorted(files):
    print(f"Reading {f}")
    with open(f) as fh:
        lines = fh.read()

    tsw = tlw = None
    try:
        tsw = re.search(": repwvl_pprts_repwvl_solar: (.*?) ", lines).groups()[0]
        tlw = re.search(": repwvl_pprts_repwvl_thermal: (.*?) ", lines).groups()[0]
    except:
        pass

    try:
        tsw = re.search(": pprts_rrtmg_solar: (.*?) ", lines).groups()[0]
        tlw = re.search(": pprts_rrtmg_thermal: (.*?) ", lines).groups()[0]
    except:
        pass


    label=".".join(os.path.basename(f).split(".")[3:5])
    print(f"{label} {tsw} {tlw}")

    res.append(dict(l=label, sw=float(tsw), lw=float(tlw)))

print(res)

plt.clf()
w=.2
for i, v in enumerate(res):
    print(i,v)
    plt.bar(i - w, v["sw"], width=w, align="edge")
    plt.text(i - w, v["sw"], f'sw {int(v["sw"])}s', rotation=45, fontsize='x-small')
    plt.bar(i, v["lw"], width=w, align="edge")
    plt.text(i, v["lw"], f'lw {int(v["lw"])}s', rotation=45, fontsize='x-small')
    plt.bar(i+w, v["sw"] + v["lw"], width=w, align="edge")
    plt.text(i+w, v["sw"] + v["lw"], f'total {int(v["sw"] + v["lw"])}s', rotation=45, fontsize='x-small')

labels = [_['l'] for _ in res]
plt.gca().set_xticks(range(len(labels)))
plt.gca().set_xticklabels(labels)
plt.xticks(rotation=20)

print(f"Plotting to {outfile}")
plt.savefig(outfile, bbox_inches='tight')
EOF

cat > plot_profiles.py << EOF
#!/usr/bin/env python3
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys

fname = sys.argv[1]
outdir = sys.argv[2]

D = xr.open_dataset(fname)
timesteps = np.unique([int(k.split(".")[1]) for k in list(filter(lambda k: "zlev" in k, D.keys()))])
E = xr.Dataset(
    data_vars=dict(
        zlay=xr.concat([D[f"zlay.{t}"] for t in timesteps], dim='time'),
        zlev=xr.concat([D[f"zlev.{t}"] for t in timesteps], dim='time'),
        edir=xr.concat([D[f"edir.{t}"] for t in timesteps], dim='time'),
        edn =xr.concat([D[f"edn.{t}"]  for t in timesteps], dim='time'),
        eup =xr.concat([D[f"eup.{t}"]  for t in timesteps], dim='time'),
        abso=xr.concat([D[f"abso.{t}"] for t in timesteps], dim='time'),
    ),
    coords=dict(
        time = timesteps,
    )
)

for t in E.time.data:
  Et = E.sel(time=t).mean(dim=("nx", "ny"))

  sub = plt.subplot(141)
  plt.plot(Et.edir, Et.zlev, label=f't{t}')
  plt.xlabel('edir [W/m2]')

  sub = plt.subplot(142)
  plt.plot(Et.edn, Et.zlev, label=f't{t}')
  plt.xlabel('edn [W/m2]')

  sub = plt.subplot(143)
  plt.plot(Et.eup, Et.zlev, label=f't{t}')
  plt.xlabel('eup [W/m2]')

  sub = plt.subplot(144)
  plt.plot(Et.abso, Et.zlay, label=f't{t}')
  plt.xlabel('abso [W/m3]')

for ax in plt.gcf().get_axes():
     ax.legend()
     ax.set_ylim(0, 1e4)
     ax.label_outer()

plt.ylabel('height [m]')
plt.ylim(0,1e4)

outfile = f"{outdir}/profiles.pdf"
print(f"plotting to file: {outfile}")
plt.savefig(outfile, bbox_inches='tight')
EOF

for specint in repwvl rrtmg
do
for solver in 2str 3_10 3_10_incomplete
do
  ID=$(basename $CLD .nc).np${NP}.$specint.$solver
  OUT=$SCRATCH/compare_specint_uclales/$ID.nc
  OUT=
  LOG=$SCRATCH/compare_specint_uclales/$ID.log
  PLOTS=$SCRATCH/compare_specint_uclales/$ID.plots/

  OPT=""
  if [[ $solver == *"incomplete"* ]]; then
    OPT="$OPT \
      -solver 3_10 \
      -solar_dir_explicit \
      -Xsolar_dir_ksp_view \
      -solar_dir_ksp_max_it 1 \
      -solar_dir_ksp_skip_residual \
      -solar_dir_ksp_ignore_max_it $(($TSTART * $DT -1)) \
      -solar_diff_explicit \
      -Xsolar_diff_ksp_view \
      -solar_diff_ksp_max_it 1 \
      -solar_diff_ksp_skip_residual \
      -solar_diff_ksp_ignore_max_it $(($TSTART * $DT -1)) \
      -thermal_diff_explicit \
      -Xthermal_diff_ksp_view \
      -thermal_diff_ksp_max_it 1 \
      -thermal_diff_ksp_skip_residual \
      -thermal_diff_ksp_ignore_max_it $(($TSTART * $DT -1)) \
      -accept_incomplete_solve \
      -absorption_by_coeff_divergence \
    "
  fi

  if [ ! -z "$OUT" ]; then
    mkdir -p $(dirname $OUT)
  fi

  if [ ! -e $LOG ]; then
  ${MPIEXEC} \
    $TENSTREAMROOT/bin/$BIN \
    -cld $CLD \
    -dx 25 -dy 25 \
    -specint $specint \
    -solver $solver \
    -tstart $TSTART \
    -tend $TEND \
    -tinc $TINC \
    -log_view \
    ${OUT:+'-out'} $OUT \
    $OPT \
    &> >(tee $LOG)
  else
    echo "Skipping $ID because logfile $LOG already exists"
  fi

  if [ ! -z "$OUT" ]; then
    mkdir -p $PLOTS
    python plot_profiles.py $OUT $PLOTS
  fi
done
done

python plot_timings.py $SCRATCH/compare_specint_uclales/timings.pdf $SCRATCH/compare_specint_uclales/$(basename $CLD .nc).np${NP}.*.log
