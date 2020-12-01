#!/bin/bash

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

bin="ex_pprts_hill"

make -j $bin || exit 1
baseopt="\
  -dx 100 \
  -dy 100 \
  -atm_filename $SCRIPTDIR/atm.dat \
  -thermal no \
  -lwc .1 \
  -cld_width 5 \
  -cld_bot 700 \
  -cld_top 600 \
  -Ny 64 \
  -Nz 50 \
  -hill_dP 100 \
  -hill_shape 10 \
  -theta0 40 -phi0 180 \
  -log_view \
  -rrtmg_bands 1,1 \
  $1"

rayli_opt="\
  -pprts_use_rayli \
  -pprts_rayli_photons 1e7 \
  -rayli_cyclic_bc \
  -show_rayli_dm3d hdf5:dm.h5 \
  "

snap_opt="\
-rayli_snapshot \
-rayli_snap_photons 1e7 \
-rayli_snap_Nx 1024 \
-rayli_snap_Ny 1024 \
-visit_image_zoom 5.5 \
-visit_parallel_scale 15e3 \
-visit_focus 150,3200,2500
-visit_view_normal 0.6212888997511092,0.3246068765558223,0.7131833416020943 \
-visit_view_up -0.6270958235924525,-0.3397617048659849,0.7009370955652606 \
  "

cat > plot_snapshot.py << EOF
import xarray as xr
from pylab import *
from matplotlib.colors import LogNorm
from contextlib import closing
with closing(xr.open_dataset('rayli_snaphots.nc')) as D:
  im = np.sum([ D[k].data for k in D.keys() ], axis=0)
im[im==0] = np.NaN
plt.imshow(im, origin='lower', cmap=cm.viridis, norm=LogNorm(vmin=1e-2))
plt.xlabel('px')
plt.ylabel('px')
cbar = plt.colorbar()
cbar.set_label('radiance [arbitrary units]')
plt.tight_layout()
plt.savefig('rayli_snaphots.pdf', bbox_inches='tight')
EOF

cat > plot_cross_section.py << EOF
import argparse
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('inp', type=str)
parser.add_argument('out', type=str)
args = parser.parse_args()

print(f"Plotting irradiance cross_sections for {args.inp} -> {args.out}")

D = xr.open_dataset(args.inp)

def get_grid_coords(var):
  Ny, Nz = D[var].mean(axis=1).shape
  z      = D['hhl'].mean(axis=1).data[:,::-1][:,:Nz][:,::-1]
  y = np.linspace(-(D.Ny/2-.5)*D.dy, (D.Ny/2-.5)*D.dy, D.Ny)
  y = np.tile(y, Nz).reshape(Nz, Ny).T
  return y,z

y, z = get_grid_coords('edir')

plt.figure(figsize=(7,8))
yrng = (0,5e3)

plt.suptitle(args.inp)

plt.subplot(411)
plt.pcolormesh(y,z, D['edir'].mean(axis=1))
cbar = plt.colorbar(); cbar.set_label('edir [W/m^2]'); plt.ylim(*yrng)

plt.subplot(412)
plt.pcolormesh(y,z, D['edn'].mean(axis=1))
cbar = plt.colorbar(); cbar.set_label('edn [W/m^2]'); plt.ylim(*yrng)

plt.subplot(413)
plt.pcolormesh(y,z, D['eup'].mean(axis=1))
cbar = plt.colorbar(); cbar.set_label('eup [W/m^2]'); plt.ylim(*yrng)

y, z = get_grid_coords('abso')

plt.subplot(414)
plt.pcolormesh(y,z, D['abso'].mean(axis=1))
cbar = plt.colorbar(); cbar.set_label('abso [W/m^3]'); plt.ylim(*yrng)

plt.savefig(args.out, bbox_inches='tight')
EOF

cat > plot_distorted.py << EOF
import glob
from pylab import *
import xarray as xr
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('basename', help='files prefix which should be read')
parser.add_argument('-o', '--out', help='plot output path')
args = parser.parse_args()

print("Plotting distorted surface irradiance")

basename = args.basename
out = 'plot_srfc_'+basename+'.pdf' if args.out is None else args.out

files = glob.glob(f'{basename}*.nc')

D3 = xr.open_dataset(f'{basename}rayli_ac.nc')


bias = lambda a,b: (mean(a)/mean(b) - 1.) * 100
rmse = lambda a,b: sqrt(mean((a-b)**2))
rrmse = lambda a,b: rmse(a,b)/mean(b) * 100

def plot(E, title):
  plt.title(title)
  for f in sorted(files):
      print(f)
      with xr.open_dataset(f) as D:
          kw = {
                  'label': f+'({:.2f}, {:.2f})'.format(rrmse(E(D),E(D3)), bias(E(D),E(D3))),
                  'linestyle': '--' if '3_10' in f else '-',
                  'alpha': .7 if '_ac' in f else 1,
                  }
          plt.plot(E(D), **kw)
  plt.legend(loc='best')


#E = lambda D: D['edir'][:,0,-1].data + D['edn'][:,0,-1].data
E = lambda D: D['edir'][:,0,-1].data

plt.figure(figsize=(12,9))

plt.subplot(4,1,1)
plot(lambda D: D['edir'][:,0,-1].data, 'edir')

plt.subplot(4,1,2)
plot(lambda D: D['edn' ][:,0,-1].data, 'edn')

plt.subplot(4,1,3)
plot(lambda D: D['eup' ][:,0,-1].data, 'eup')

plt.subplot(4,1,4)
plot(lambda D: D['abso' ][:,0,-1].data, 'abso')

print("saving to", out)
plt.savefig(out)
EOF

function runex {
  EXEC=$1
  OUT=$2
  OPT=$3

  if [ -e $OUT ]; then
    echo "Skipping $OUT"
  else
    ($EXEC bin/$bin $baseopt -out $OUT $OPT) &
  fi
}

mpiexec="srun --time=08:00:00 -n 40 -N 1-4 -p met-ws,cluster --mpi=pmix"
rayexec="srun --time=08:00:00 -n 2 -N 2 -c 10 -p met-ws,cluster --mpi=pmix"

for SZA in 0 20 40 60
do
  runex "$mpiexec" "res_${SZA}_1d.nc                  "  "-theta0 $SZA  -twostr_only"
  runex "$rayexec" "res_${SZA}_rayli_ac.nc            "  "-theta0 $SZA  -pprts_atm_correction $rayli_opt"
  runex "$mpiexec" "res_${SZA}_3_10str.nc             "  "-theta0 $SZA "
  runex "$mpiexec" "res_${SZA}_3_10str_distorted_ac.nc"  "-theta0 $SZA  -pprts_atm_correction -bmc_online"
  wait

  python plot_distorted.py res_${SZA}_
  for f in res_${SZA}_*.nc; do
    python plot_cross_section.py $f plot_cross_$f.pdf
  done
done
