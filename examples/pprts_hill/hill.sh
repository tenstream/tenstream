rm out_pprts_hill.nc
make -j -C $HOME/tenstream/build ex_pprts_hill && \
  $HOME/tenstream/build/bin/ex_pprts_hill \
  -atm_filename $HOME/tenstream/examples/pprts_hill/afglus_100m.dat \
  -twostr_only \
  -pprts_slope_correction \
  -lwc .1 \
  -Ny 64 \
  -hill_dP 50 \
  -hill_shape 5 \
  -theta0 40 -phi0 0 \
  $1

cat > plot_hill_srfc_edir.py << EOF
from pylab import *
import xarray as xr
D=xr.open_dataset("out_pprts_hill.nc")
plot(D["0edir"][:,1,-1])
#ylim(750,820)
waitforbuttonpress()
EOF
