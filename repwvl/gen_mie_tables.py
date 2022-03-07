#!/usr/bin/env python3
import glob
import subprocess
import shutil
import os
import tempfile
from jinja2 import Template
import numpy as np
import xarray as xr

def run_mie_calc(wdir, wvls, phase='water', reff=10., mie_bin=None, verbose=False):
    if verbose:
        print('workdir', wdir)
    os.chdir(wdir)

    mie_template = Template("""
mie_program MIEV0
refrac             {{phase}} # water or ice
distribution GAMMA {{alpha}} # water 7, ice 1
r_eff              {{reff}}
wavelength         {{wvl_grid_file}}
temperature        {{temperature}}
basename           {{basename}}
verbose
output_user netcdf
""")

    if phase == 'water':
        alpha = 7
    elif phaes == 'ice':
        alpha = 1
    else:
        raise ValueError(f"Dont know which alpha value should be used for phase = {phase}")

    mie_input = {
            'phase': phase,
            'reff':  reff,
            'alpha': alpha,
            'wvl_grid_file': 'wvl_grid_file.tmp',
            'temperature': 280,
            'basename': f'mie.{phase}.reff{reff}'
            }
    np.savetxt(mie_input['wvl_grid_file'], list(zip(range(1,len(wvls)+1), wvls)))
    inp = mie_template.render(mie_input)
    if verbose:
        print(inp)
    with open("mie.inp", "w") as fh:
        fh.write(inp)

    cmd = [mie_bin, "-f", "mie.inp"]
    if verbose:
        print(f"Running: {cmd}")
    subout = subprocess.run(cmd, capture_output=True)

    outfile = glob.glob(mie_input['basename']+'*.cdf')[0]
    with xr.open_dataset(outfile) as D:
        R = xr.Dataset(
                data_vars=dict(
                    ext=(["wvl", "reff"], D.ext.data, D.ext.attrs),
                    ssa=(["wvl", "reff"], D.ssa.data, D.ssa.attrs),
                    g  =(["wvl", "reff"], D.gg.data, D.gg.attrs),
                    ),
                coords=dict(
                    reff=(["reff",], D.reff.data, D.reff.attrs),
                    wvl=(["wvl",], D.wavelen.data, D.wavelen.attrs),
                    ),
                )
    return R

def run_exp(wvls, reffs, mie_bin, verbose=False):
    ret = []
    for reff in reffs:
        with tempfile.TemporaryDirectory() as tmpdir:
            ret.append( run_mie_calc(tmpdir, wvls, reff=reff, mie_bin=mie_bin, verbose=verbose) )
    return xr.merge(ret)

def _main():
    import argparse

    parser = argparse.ArgumentParser(description="Run the libRadtran mie tool to create input data for representative wavelengths")
    parser.add_argument('output', type=str, help='output path')
    parser.add_argument('--wvls', '-w', type=float, nargs='+', help='wavelengths to compute [nm]')
    parser.add_argument('--reffs', '-r', type=float, nargs='+', help='reffs to compute [mu]')
    parser.add_argument('--mietool', '-m', type=str, default=None, help='path to the libRadtran mietool')
    parser.add_argument('--verbose', '-v', action='store_true', default=False, help='verbose output')
    args = parser.parse_args()

    args.output = os.path.abspath(args.output)

    if args.mietool is None:
        args.mietool = shutil.which("mie")

    if args.mietool is None:
        raise ValueError(f"cant find mietool, please provide the path to the libRadtran mietool")

    if args.verbose:
        print(args)

    ret = run_exp(args.wvls, args.reffs, mie_bin=args.mietool, verbose=args.verbose)

    if args.verbose:
        print(ret)

    if args.verbose:
        print(f"Writing output to {args.output}")
    ret.to_netcdf(args.output)


if __name__ == '__main__':
    _main()
