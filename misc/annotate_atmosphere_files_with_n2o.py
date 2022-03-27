#!/usr/bin/env python3

def main(garand_file, atm_file, out_file):
    import numpy as np
    import xarray as xr
    import os
    import sys
    from scipy.interpolate import interp1d

    if not os.path.exists(garand_file):
        raise FileNotFoundError(f"garand_file does not exist: {garand_file}")

    if not os.path.exists(atm_file):
        raise FileNotFoundError(f"atm_file does not exist: {atm_file}")

    atm = np.loadtxt(atm_file)
    if atm.shape[1] != 9:
        raise Exception(f"atm_file does not have 9 columns: shape {atm.shape}")

    G = garand = xr.open_dataset(garand_file)
    p = G.Atmosphere.isel(numCols=0).mean("numAtm").data
    n2o = G.Atmosphere.isel(numCols=6).mean("numAtm").data
    fn2o = interp1d(p * 1e-2, n2o, bounds_error=False, fill_value="extrapolate")

    new_atm = np.c_[ atm,  fn2o(atm[:, 1]) * atm[:, 3]]

    header = []
    with open(atm_file) as fh:
        for line in fh:
            if "#" in line:
                header.append(line)

    new_header = []
    for h in header:
        h = h.replace("\n", "")
        if "no2" in h:
            h += "     n2o(cm-3)\n"
        new_header.append(h)

    if out_file is None:
        fh = sys.stdout
    else:
        fh = open(out_file, 'w')

    for h in new_header:
        fh.writelines(h + '\n')
    for k in range(new_atm.shape[0]):
        fh.writelines(" ".join(map( lambda f: f"{f:14.8g}", new_atm[k])) + '\n')

    if fh is not sys.stdout:
        fh.close()


import argparse

parser = argparse.ArgumentParser()
parser.add_argument('garand_file', help='garand input file (garand_profiles.nc)')
parser.add_argument('atm_file', help='atmosphere file')
parser.add_argument('--out', '-o', help='output atmosphere file')
parser.add_argument('--inplace', '-i', action='store_true', default=False, help='change input file inplace')
args = parser.parse_args()

if args.inplace:
    args.out = args.atm_file

main(args.garand_file, args.atm_file, args.out)
