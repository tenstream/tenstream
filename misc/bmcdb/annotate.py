#!/usr/bin/env python3

import numpy as np
import os
import xarray as xr

from eddington import eddington_ec

def load_data(ncfile, f_input, f_feature):
    with xr.load_dataset(ncfile) as D:
        return f_input(D), f_feature(D)

def annotate(inpfile, var='C', lsolar=False, renormalization=True, shuffle=None, test_fraction=0, verbose=False, **kwargs):
    if verbose: print(f" ... loading data for {var}")
    inp, out = load_data(inpfile,
            lambda D: D.inp,
            lambda D: D[[var,f'{var}tol']],
            )
    dim_N, dim_src, dim_dst = out.dims

    theta = inp.sel(inpvars='theta') * lsolar
    phi   = inp.sel(inpvars='phi') * lsolar

    t, r, rdir, sdir, tdir = eddington_ec(
            dtau = (inp.sel(inpvars='kabs') + inp.sel(inpvars='ksca')) * inp.sel(inpvars='dz'),
            w0 = inp.sel(inpvars='ksca') / (inp.sel(inpvars='kabs') + inp.sel(inpvars='ksca')),
            g = inp.sel(inpvars='g'),
            mu0 = np.cos(np.deg2rad(theta)),
            ).T

    X = np.vstack([
        inp.sel(inpvars='kabs'),
        inp.sel(inpvars='ksca'),
        -np.expm1(-inp.sel(inpvars='kabs') * inp.sel(inpvars='dz')),
        -np.expm1(-inp.sel(inpvars='ksca') * inp.sel(inpvars='dz')),
        inp.sel(inpvars='g'),
        inp.sel(inpvars='dz'),
        phi,
        theta,
        inp.sel(inpvars='tanx'),
        inp.sel(inpvars='tany'),
        t,
        r,
        rdir,
        sdir,
        tdir,
        ]).T

    X_descr = [
        'kabs',
        'ksca',
        '1-exp(tau_abs)',
        '1-exp(tau_sca)',
        'g',
        'dz',
        'phi',
        'theta',
        'tanx',
        'tany',
        'eddington_t',
        'eddington_r',
        'eddington_rdir',
        'eddington_sdir',
        'eddington_tdir',
        ]

    coeffs = out[var]
    norm   = np.maximum(1e-40, out[var].sum(dim_dst))
    stddev = out[f"{var}tol"]
    # never trust any coefficient more than the stddev exit condition of the montecarlo solver
    #stddev = np.maximum(1e-4, stddev)
    stddev = np.maximum(1e-8, stddev)
    stddev = np.where(coeffs <= 0, 1e-7, stddev)

    if renormalization:
        if verbose: print(" ... renormalization")
        coeffs = coeffs / norm

    Y = np.hstack([
        coeffs.data.reshape(inp.N.size,-1),
        norm.data,
        stddev.reshape(inp.N.size,-1),
        ])

    if shuffle is not None:
        print(f"Shuffling data")
        X, Y = X[shuffle], Y[shuffle]

    split = int(inp.N.size * test_fraction)
    sl_test = slice(None, split)
    sl_train  = slice(split, None)

    D = xr.Dataset(
        data_vars={
            f"{var}_X"      : (("N", f"inpvars"), X[sl_train]),
            f"{var}_Y"      : (("N", f"{var}_outvars"), Y[sl_train]),
        },
        coords={
            f"inpvars" : X_descr,
        },
        attrs={
            f"Nsrc_{var}"    : out.dims[dim_src],
            f"Ndst_{var}"    : out.dims[dim_dst],
            f"Ncoeff_{var}"  : coeffs.data.reshape(inp.N.size,-1).shape[-1],
            f"Nnorm_{var}"   : norm.data.shape[-1],
            f"Nstddev_{var}" : stddev.reshape(inp.N.size,-1).shape[-1],
            },
    )

    if test_fraction > 0:
        if verbose: print(f" ... splitting apart a test fraction {test_fraction}")
        D.attrs['test_fraction'] = test_fraction
        D[f"{var}_X_test"] = (("Ntest", f"inpvars"), X[sl_test])
        D[f"{var}_Y_test"] = (("Ntest", f"{var}_outvars"), Y[sl_test])

    return D

def _main():
    import argparse
    parser = argparse.ArgumentParser(description='Annotate bmc database output with additional features')
    parser.add_argument('input', type=str, help='bmcdatabse dump')
    parser.add_argument('out', type=str, help='output filename')
    parser.add_argument('-tf', '--test_fraction', type=float, default=.05, help='fraction to set apart a test dataset which should never be used to train or validate a model but rather to check against before emitting a model')
    parser.add_argument('-v', '--verbose', action="store_true", help='more verbose output')
    parser.add_argument('-s', '--shuffle', action="store_true", help='shuffle training part of dataset')
    parser.add_argument('-r', '--renormalization', action="store_true", help='apply renormalization to the coefficients')
    args = parser.parse_args()

    if args.verbose:
        print(f"Reading from < {args.input}\nWriting to   > {args.out}")

    if os.path.exists(args.out):
        raise ValueError(f"Output file already exists... will not overwrite!")

    if args.shuffle: # shuffling with permutation
        if args.verbose: print(" ... shuffling")
        with xr.open_dataset(args.input) as D:
            perm = np.arange(D.N.size)
        np.random.shuffle(perm)
        args.shuffle = perm
        #D = D.isel(N=perm)

    D_C = annotate(args.input, 'C', **vars(args))
    D_S = annotate(args.input, 'S', lsolar=True, **vars(args))
    D_T = annotate(args.input, 'T', lsolar=True, **vars(args))
    D = xr.merge([D_C, D_S, D_T], combine_attrs='no_conflicts')


    if args.verbose: print(D)

    if args.verbose: print(f" ... writing data to netcdf: {args.out}")
    D.to_netcdf(args.out)


if __name__ == "__main__":
    _main()
