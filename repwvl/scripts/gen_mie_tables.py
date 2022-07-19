#!/usr/bin/env python3
from itertools import repeat
from jinja2 import Template
import glob
import multiprocessing
import numpy as np
import os
import shutil
import subprocess
import sys
import tempfile
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
        temperature = 280
    elif phase == 'ice':
        alpha = 1
        temperature = 260
    else:
        raise ValueError(f"Dont know which alpha value should be used for phase = {phase}")

    mie_input = {
            'phase': phase,
            'reff':  reff,
            'alpha': alpha,
            'wvl_grid_file': 'wvl_grid_file.tmp',
            'temperature': temperature,
            'basename': f'mie.{phase}.reff{reff}'
            }
    np.savetxt(mie_input['wvl_grid_file'], list(zip(range(1,len(wvls)+1), wvls)))
    inp = mie_template.render(mie_input)
    if verbose:
        print(f"Wavelenghts {wvls}")
        print(inp)

    with open("mie.inp", "w") as fh:
        fh.write(inp)

    cmd = [mie_bin, "-f", "mie.inp"]

    if verbose:
        print(f"Running: {cmd} {mie_input}")

    subout = subprocess.run(cmd, capture_output=True)

    if verbose:
        print(f"Finished: {cmd} {mie_input}")

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
    if verbose:
        print(f"Written output for {mie_input} to {outfile}")
    return R

def pool_fct(reff, wvls, phase, mie_bin, verbose, wdir='./'):
    fct_one_dir = lambda directory: run_mie_calc(directory, wvls, phase=phase, reff=reff, mie_bin=mie_bin, verbose=verbose)
    if wdir == 'tmp':
        with tempfile.TemporaryDirectory() as tmpdir:
            return fct_one_dir(tmpdir)
    else:
        dirname = os.path.abspath(f"{wdir}/mie_calcs/{phase}/{reff}/{np.sum(wvls)}/")
        if not os.path.exists(dirname):
            os.makedirs(dirname, exist_ok=True)

        if not os.path.exists(f"{dirname}/result.nc"):
            R = fct_one_dir(dirname)
            R.to_netcdf(f"{dirname}/result.nc")

        with xr.open_dataset(f"{dirname}/result.nc") as R:
            return R.load()

def parse_omp_envvar(env_value):
    return int(env_value.strip().split(",")[0])

def estimate_cpus():
    limit = float("inf")
    if "OMP_THREAD_LIMIT" in os.environ:
        limit = parse_omp_envvar(os.environ["OMP_THREAD_LIMIT"])

    if "OMP_NUM_THREADS" in os.environ:
        cpus = parse_omp_envvar(os.environ["OMP_NUM_THREADS"])
    else:
        try:
            cpus = len(os.sched_getaffinity(os.getpid()))
        except (AttributeError, OSError) as e:
            cpus = os.cpu_count()

    return min(cpus, limit)

def run_exp(wvls, reffs, phase, mie_bin, verbose=False, wdir='./'):
    pool = multiprocessing.Pool(processes=estimate_cpus())
    ret = list( pool.starmap(pool_fct, zip(reffs, repeat(wvls), repeat(phase), repeat(mie_bin), repeat(verbose), repeat(wdir))) )
    return xr.merge(ret).sortby("reff").sortby("wvl")

def _main():
    import argparse

    parser = argparse.ArgumentParser(description="Run the libRadtran mie tool to create input data for representative wavelengths")
    parser.add_argument('output', type=str, help='output path')
    parser.add_argument('--phase', '-p', default='water', choices=['water', 'ice'], help='aggregate phase (default: %(default)s)')
    parser.add_argument('--wvls', '-w', type=float, nargs='+', help='wavelengths to compute [nm]')
    parser.add_argument('--reffs', '-r', type=float, nargs='+', help='reffs to compute [mu]')
    parser.add_argument('--mietool', '-m', type=str, default=None, help='path to the libRadtran mietool')
    parser.add_argument('--verbose', '-v', action='store_true', default=False, help='verbose output')
    parser.add_argument('--wdir', '-C', type=str, default="./", help='work directory where to compute mie calculations, if "tmp" then temporary dirs are used (no restarts though)')
    args = parser.parse_args()

    args.output = os.path.abspath(args.output)

    if args.mietool is None:
        args.mietool = shutil.which("mie")

    if args.mietool is None:
        raise ValueError(f"cant find mietool, please provide the path to the libRadtran mietool")

    if args.verbose:
        print(args)

    ret = run_exp(args.wvls, args.reffs, args.phase, mie_bin=args.mietool, verbose=args.verbose, wdir=args.wdir)

    if args.verbose:
        print(ret)

    if args.verbose:
        print(f"Writing output to {args.output}")
    ret.to_netcdf(args.output)


if __name__ == '__main__':
    _main()
