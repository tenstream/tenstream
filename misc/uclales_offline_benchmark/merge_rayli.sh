#!/bin/bash
. config.sh

ITER=${1:-}
IDENT=rayli${ITER}

# pack results
for SPEC in sw lw
do
  python -c "\
import xarray as xr; import glob; \
files = glob.glob('${OUTPREFIX}/$IDENT/$IDENT.*.$SPEC.nc'); \
print(f'Merging {files}'); \
Ds = [ xr.open_dataset(f) for f in glob.glob('${OUTPREFIX}/$IDENT/$IDENT.*.$SPEC.nc') ]; \
xr.merge(Ds).to_netcdf('${OUTPREFIX}/$IDENT/$IDENT.$SPEC.nc'); \
"
done
