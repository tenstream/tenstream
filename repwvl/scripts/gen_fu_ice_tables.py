#!/usr/bin/env python3

import os
import xarray as xr
import numpy as np
import datetime


def read_fu_data(libRadtran_dir, **kwargs):
    fu96_asy = np.loadtxt(f"{libRadtran_dir}/data/ic/fu96/fu96.asy")
    fu96_del = np.loadtxt(f"{libRadtran_dir}/data/ic/fu96/fu96.del")
    fu96_ext = np.loadtxt(f"{libRadtran_dir}/data/ic/fu96/fu96.ext")
    fu96_ssa = np.loadtxt(f"{libRadtran_dir}/data/ic/fu96/fu96.ssa")
    fu96_wvl = np.hstack((fu96_asy[0,0], fu96_asy[:,1]))


    fu98_abs = np.loadtxt(f"{libRadtran_dir}/data/ic/fu98/fu98.abs")
    fu98_asy = np.loadtxt(f"{libRadtran_dir}/data/ic/fu98/fu98.asy")
    fu98_ext = np.loadtxt(f"{libRadtran_dir}/data/ic/fu98/fu98.ext")


    R = xr.Dataset(
                data_vars=dict(
                    fu96_asy=(['fu96_bands', 'c_asy'], fu96_asy[:,2:],
                                dict(
                                    units='1',
                                    long_name='Asymmetry factor g  =  c0  +  c1 * D  +  c2 * D**2  +c3 * D**3, eq. 3.9c',
                                    )),
                    fu96_ext=(['fu96_bands', 'c_ext'], fu96_ext[:,2:],
                                dict(
                                    units='1',
                                    long_name='Extinction coefficient b = IWC * (a0 + a1/D), eq. 3.9a',
                                    )),
                    fu96_ssa=(['fu96_bands', 'c_ssa'], fu96_ssa[:,2:],
                                dict(
                                    units='1',
                                    long_name='Coalbedo  1-w  =  b0  +  b1 * D  +  b2 * D**2  +b3 * D**3, eq. 3.9b',
                                    )),
                    #fu96_del=(['fu96_bands', 'c_del'], fu96_del[:,2:],
                    #            dict(
                    #                units='1',
                    #                long_name='delta-fraction  f  =  d0  +  d1 * D  +  d2 * D**2  +d3 * D**3, eq. 3.9d',
                    #                )),
                    fu98_abs=(['fu98_wvl', 'cc_abs'], fu98_abs[:,1:],
                                dict(
                                    units='1',
                                    long_name='Absorption coefficient b = IWC/D * (b0 + b1*D + b2*D**2 + b3*D**3)',
                                    )),
                    fu98_asy=(['fu98_wvl', 'cc_asy'], fu98_asy[:,1:],
                                dict(
                                    units='1',
                                    long_name='Asymmetry parameter g = c0 +  c1 * D  + c2 * D**2  +  c3 * D**3',
                                    )),
                    fu98_ext=(['fu98_wvl', 'cc_ext'], fu98_ext[:,1:],
                                dict(
                                    units='1',
                                    long_name='Extinction coefficient b = IWC * (a0 + a1/D + a2/D**2)',
                                    )),
                    ),
                coords=dict(
                    fu96_wvl=(["fu96_wvl",], fu96_wvl[1:], dict(units='um', long_name='upper bound of wavelength band')),
                    fu98_wvl=(["fu98_wvl",], fu98_abs[:,0], dict(units='um', long_name='wavelength limits')),
                    ),
                attrs=dict(
                    created = datetime.datetime.now().isoformat(),
                    libRadtran_data_sourcedir = os.path.realpath(libRadtran_dir),
                    )
            )
    return R


def gen_repwvl_fu_data(repwvl_file, fu):
    repwvl = xr.open_dataset(repwvl_file)
    if repwvl.wvl.data[0] < 1e3:
        print("First wavelength of repwvl file is smaller than 1mu, assuming this is a solar file")

        bands = []
        for w in repwvl.wvl.data * 1e-3:
            bands.append(np.argmin(fu.fu96_wvl.data < w) -1)

        repwvl_fu = fu.isel(fu96_bands=bands)

        R = xr.Dataset(
                data_vars=dict(
                    fu96_ext = repwvl_fu.fu96_ext,
                    fu96_ssa = repwvl_fu.fu96_ssa,
                    fu96_asy = repwvl_fu.fu96_asy,
                    #fu96_del = repwvl_fu.fu96_del,
                    ),
                coords=dict(
                    fu96_wvl = ("fu96_wvl", repwvl.wvl.data*1e-3, dict(units='um', long_name='wavelength from repwvl file')),
                    ),
                attrs=dict(
                    created_96 = datetime.datetime.now().isoformat(),
                    repwvl_file_96 = os.path.realpath(repwvl_file),
                    )
                )

    else:
        print("First wavelength of repwvl file is larger than 1mu, assuming this is a thermal file")

        repwvl_fu = fu.sel(fu98_wvl=repwvl.wvl.data*1e-3, method='nearest')

        R = xr.Dataset(
                data_vars=dict(
                    fu98_ext = repwvl_fu.fu98_ext,
                    fu98_abs = repwvl_fu.fu98_abs,
                    fu98_asy = repwvl_fu.fu98_asy,
                    ),
                coords=dict(
                    fu98_wvl = repwvl_fu.fu98_wvl,
                    ),
                attrs=dict(
                    created_98 = datetime.datetime.now().isoformat(),
                    repwvl_file_98 = os.path.realpath(repwvl_file),
                    )
                )
    return R


def _main():
    import argparse

    parser = argparse.ArgumentParser(description="creates a data file for fu ice optical properties computations")
    parser.add_argument('libRadtran_dir', type=str, help='path to a libRadtran installation, used to read the FU tables from')
    parser.add_argument('output', type=str, help='output path')
    parser.add_argument('--repwvl_file', type=str, help='path to a repwvl file that is used for its wavelength grid')
    parser.add_argument('--verbose', '-v', action='store_true', default=False, help='verbose output')
    args = parser.parse_args()

    args.output = os.path.abspath(args.output)

    if args.verbose:
        print(args)

    fu = read_fu_data(**vars(args))

    if args.repwvl_file is None:
        if args.verbose:
            print("fu_data:", fu)
            print(f"Writing output to {args.output}")
        fu.to_netcdf(args.output)
    else:
        repwvl_fu = gen_repwvl_fu_data(args.repwvl_file, fu)
        print("repwvl_fu data:", repwvl_fu)

        if os.path.exists(args.output):
            print(f"Merging with already present output file: {args.output}")
            with xr.open_dataset(args.output) as old:
                repwvl_fu = xr.merge((repwvl_fu, old.load()), compat='no_conflicts', combine_attrs='override')
            print("merged repwvl_fu data:", repwvl_fu)

        if args.verbose:
            print(f"Writing output to {args.output}")
        repwvl_fu.to_netcdf(args.output)

if __name__ == '__main__':
    _main()
