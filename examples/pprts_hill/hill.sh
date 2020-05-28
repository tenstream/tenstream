make -j -C $HOME/tenstream/build ex_pprts_hill || exit 1
baseopt="\
  -atm_filename $HOME/tenstream/examples/pprts_hill/afglus_100m.dat \
  -thermal no \
  -lwc 1 \
  -Ny 64 \
  -hill_dP 200 \
  -hill_shape 5 \
  -theta0 40 -phi0 0 \
  $1"

bin=$HOME/tenstream/build/bin/ex_pprts_hill
out="out_pprts_hill.nc"
rm $out
$bin $baseopt -twostr_only -pprts_slope_correction && mv $out res_1d.nc
$bin $baseopt && mv $out res_10str.nc
$bin $baseopt -topography && mv $out res_10str_topo.nc

