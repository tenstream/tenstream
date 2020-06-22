make -j ex_pprts_hill || exit 1
baseopt="\
  -atm_filename $HOME/tenstream/examples/pprts_hill/afglus_100m.dat \
  -dx 100 \
  -dy 100 \
  -thermal no \
  -lwc .1 \
  -cld_width 5 \
  -cld_bot 700 \
  -cld_top 650 \
  -Ny 64 \
  -hill_dP 150 \
  -hill_shape 10 \
  -theta0 40 -phi0 0 \
  $1"

rayli_opt="\
  -pprts_use_rayli \
  -rayli_diff_flx_origin 0,0,-inf \
  -rayli_cylic_bc \
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

bin=bin/ex_pprts_hill
out="out_pprts_hill.nc"
rm -f $out

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

[ ! -e res_rayli.nc ] && $bin $baseopt $rayli_opt $snap_opt && mv $out res_rayli.nc

[ ! -e res_1d.nc ] && $bin $baseopt -twostr_only -pprts_slope_correction && mv $out res_1d.nc
[ ! -e res_10str.nc ] && $bin $baseopt && mv $out res_10str.nc
[ ! -e res_10str_topo.nc ] && $bin $baseopt -topography && mv $out res_10str_topo.nc
