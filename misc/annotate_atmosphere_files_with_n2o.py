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

    header = []
    with open(atm_file) as fh:
        for line in fh:
            if "#" in line:
                header.append(line)

    have_n2o = np.any(["n2o" in _ for _ in header])
    have_ch4 = np.any(["ch4" in _ for _ in header])

    if (not have_n2o) and (not have_ch4):
        if atm.shape[1] != 9:
            raise Exception(f"atm_file does not have 9 columns: shape {atm.shape}")

    if (have_n2o) and (not have_ch4):
        if atm.shape[1] != 10:
            raise Exception(f"atm_file does not have 10 columns: shape {atm.shape}")

    if (have_n2o) and (have_ch4) and atm.shape[1] == 11:
        print('seems we have all the data we want... returning')

    # Garand profiles (9 columns)
    # p, T, z, H2O, O3, CO2, N2O, CO, CH4
    G = garand = xr.open_dataset(garand_file)
    p = G.Atmosphere.isel(numCols=0).mean("numAtm").data

    if not have_n2o:
        n2o = G.Atmosphere.isel(numCols=6).mean("numAtm").data
        fn2o = interp1d(p * 1e-2, n2o, bounds_error=False, fill_value="extrapolate")
        atm = np.c_[ atm,  fn2o(atm[:, 1]) * atm[:, 3]]

    if not have_ch4:
        ch4 = G.Atmosphere.isel(numCols=8).mean("numAtm").data
        fch4 = interp1d(p * 1e-2, ch4, bounds_error=False, fill_value="extrapolate")
        atm = np.c_[ atm,  fch4(atm[:, 1]) * atm[:, 3]]


    header = []
    with open(atm_file) as fh:
        for line in fh:
            if "#" in line:
                header.append(line)

    new_header = []
    for h in header:
        h = h.replace("\n", "")
        if (not have_n2o) and ("no2" in h):
            h += "     n2o(cm-3)\n"
        if (not have_ch4) and ("no2" in h):
            h += "     ch4(cm-3)\n"
        new_header.append(h)

    if out_file is None:
        fh = sys.stdout
    else:
        fh = open(out_file, 'w')

    for h in new_header:
        fh.writelines(h + '\n')
    for k in range(atm.shape[0]):
        fh.writelines(" ".join(map( lambda f: f"{f:14.8g}", atm[k])) + '\n')

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
