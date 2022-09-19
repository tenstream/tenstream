import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import glob

f2ts = lambda f: int(f.split('.')[1])

rmse = lambda a, b: np.sqrt(np.mean((a - b) ** 2))
rel_rmse = lambda a, b: rmse(a, b) / np.abs(np.mean(b)) * 100
bias = lambda a, b: np.mean(a - b)
rel_bias = lambda a, b: bias(a, b) / np.mean(b) * 100
mae = lambda a, b: np.mean(np.abs(a - b))
rel_mae = lambda a, b: mae(a, b) / np.abs(np.mean(b)) * 100

def main(prefix, spec='sw', tmin=800, tmax=900, verbose=False):
    D2str      = xr.open_dataset(f"{prefix}/2str/2str.{spec}.nc")
    D2str_rrtmg= xr.open_dataset(f"{prefix}/2str_rrtmg/2str_rrtmg.{spec}.nc")
    D3_10      = xr.open_dataset(f"{prefix}/3_10/3_10.{spec}.nc")
    #D3_10lt    = xr.open_dataset(f"{prefix}/3_10_lowtol/3_10_lowtol.{spec}.nc")
    D3_10i1_2  = xr.open_dataset(f"{prefix}/3_10_incomplete_1_2/3_10_incomplete_1_2.{spec}.nc")
    D3_10i1_4  = xr.open_dataset(f"{prefix}/3_10_incomplete_1_4/3_10_incomplete_1_4.{spec}.nc")
    D3_10i1_8  = xr.open_dataset(f"{prefix}/3_10_incomplete_1_8/3_10_incomplete_1_8.{spec}.nc")
    D3_10i6_2  = xr.open_dataset(f"{prefix}/3_10_incomplete_6_2/3_10_incomplete_6_2.{spec}.nc")
    D3_10i6_4  = xr.open_dataset(f"{prefix}/3_10_incomplete_6_4/3_10_incomplete_6_4.{spec}.nc")
    D3_10i6_8  = xr.open_dataset(f"{prefix}/3_10_incomplete_6_8/3_10_incomplete_6_8.{spec}.nc")
    D3_10i6_12 = xr.open_dataset(f"{prefix}/3_10_incomplete_6_12/3_10_incomplete_6_12.{spec}.nc")
    Drayli0    = xr.open_dataset(f"{prefix}/rayli0/rayli0.{spec}.nc")
    Drayli1    = xr.open_dataset(f"{prefix}/rayli1/rayli1.{spec}.nc")

    def plot_rmse(varfunc, errfunc, D, label, nstep=1, **kwargs):
        errors = []
        times = []

        tselect = lambda t: tmin + ((t - tmin) // nstep ) * nstep

        for t in range(tmin, tmax+1):
            tsel = tselect(t)
            ttruth = t

            try:
                ytruth = varfunc(Drayli0, ttruth)
            except KeyError:
                print(f"Skipping because ttruth is not there {ttruth}")
                continue

            try:
                y = varfunc(D, tsel)
            except KeyError:
                print(f"{label} Skipping because tsel is not there {tsel}")
                continue


            err = errfunc(y.data, ytruth.data)
            print(f"{label} {tsel} {ttruth} {err}")

            times.append(float(ttruth))
            errors.append(float(err))
        plt.plot(times, errors, '-', label=f"{label} ({np.mean(errors):.0f}%)", **kwargs)
        return times, errors

    def get_abso(X, t, maxz=4.5e3):
        zmask = (X[f"zlay.{t}"] < maxz).data
        a = X[f"abso.{t}"].isel(zlay=zmask)
        return a

    rayli_plot_kwargs = dict(marker='.', linestyle=':', color='.1', linewidth=.8)

    def plot_sims(varfunc, errfunc):
        plot_rmse(varfunc, errfunc, Drayli1    , 'mc'                , nstep=1 , **rayli_plot_kwargs)
        plot_rmse(varfunc, errfunc, D2str      , '2str'              , nstep=1)
        plot_rmse(varfunc, errfunc, D2str_rrtmg, '2str_rrtmg'        , nstep=1)
        plot_rmse(varfunc, errfunc, D2str      , '2str dt=1min'      , nstep=6 , linestyle='--')
        plot_rmse(varfunc, errfunc, D3_10      , '3_10'              , nstep=1)
        #plot_rmse(varfunc, errfunc, D3_10lt    , '3_10_lt'           , nstep=1)
        plot_rmse(varfunc, errfunc, D3_10      , '3_10 dt=1min'      , nstep=6 , linestyle='--')
        plot_rmse(varfunc, errfunc, D3_10i1_2  , '3_10i it=2'        , nstep=1 , linestyle=':')
        plot_rmse(varfunc, errfunc, D3_10i1_4  , '3_10i it=4'        , nstep=1 , linestyle=':')
        plot_rmse(varfunc, errfunc, D3_10i1_8  , '3_10i it=8'        , nstep=1 , linestyle=':')
        plot_rmse(varfunc, errfunc, D3_10i6_2  , '3_10i it=2  dt=1min', nstep=6 , linestyle='-.')
        plot_rmse(varfunc, errfunc, D3_10i6_4  , '3_10i it=4  dt=1min', nstep=6 , linestyle='-.')
        plot_rmse(varfunc, errfunc, D3_10i6_8  , '3_10i it=8  dt=1min', nstep=6 , linestyle='-.')
        plot_rmse(varfunc, errfunc, D3_10i6_12 , '3_10i it=12 dt=1min', nstep=6 , linestyle='-.')

    plt.figure(1)
    plt.clf()
    plot_sims(get_abso, rel_rmse)
    plt.legend()
    plt.xlabel('timestep [dt=10s]')
    plt.ylabel('absorption rel. RMSE [%]')
    plt.savefig(f'{prefix}/{spec}_Abso_rmse.pdf', bbox_inches='tight')


    plt.figure(2)
    plt.clf()
    plot_sims(get_abso, rel_bias)
    plt.legend()
    plt.xlabel('timestep [dt=10s]')
    plt.ylabel('absorption rel. bias [%]')
    plt.savefig(f'{prefix}/{spec}_Abso_bias.pdf', bbox_inches='tight')

    plt.figure(3)
    plt.clf()
    plot_sims(get_abso, rel_mae)
    plt.legend()
    plt.xlabel('timestep [dt=10s]')
    plt.ylabel('absorption rel. mae [%]')
    plt.savefig(f'{prefix}/{spec}_Abso_mae.pdf', bbox_inches='tight')


    def get_Enet(X, t):
       if f'edir.{t}' in X:
           return X[f'edir.{t}'].isel(zlev=-1) + X[f'edn.{t}'].isel(zlev=-1) - X[f'eup.{t}'].isel(zlev=-1)
       else:
           return X[f'edn.{t}'].isel(zlev=-1) - X[f'eup.{t}'].isel(zlev=-1)

    plt.figure(4)
    plt.clf()
    plot_sims(get_Enet, rel_rmse)
    plt.legend()
    plt.xlabel('timestep [dt=10s]')
    plt.ylabel('Enet @ surface rel. RMSE [%]')
    plt.savefig(f'{prefix}/{spec}_Enet_rmse.pdf', bbox_inches='tight')



    plt.figure(5)
    plt.clf()
    plot_sims(get_Enet, rel_bias)
    plt.legend()
    plt.xlabel('timestep [dt=10s]')
    plt.ylabel('Enet @ surface rel. bias [%]')
    plt.savefig(f'{prefix}/{spec}_Enet_bias.pdf', bbox_inches='tight')

    plt.figure(6)
    plt.clf()
    plot_sims(get_Enet, rel_mae)
    plt.legend()
    plt.xlabel('timestep [dt=10s]')
    plt.ylabel('Enet @ surface rel. MAE [%]')
    plt.savefig(f'{prefix}/{spec}_Enet_mae.pdf', bbox_inches='tight')



def _main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('prefix', help='output folder prefix, e.g. strided4 or strided8')
    parser.add_argument('--tmin', default=800, type=int, help='start timestep')
    parser.add_argument('--tmax', default=900, type=int, help='start timestep')

    parser.add_argument('--verbose', '-v', action='store_true', default=False, help='verbose output')
    parser.add_argument('--skip_thermal', '-t', action='store_true', default=False, help='dont run thermal plots')
    parser.add_argument('--skip_solar', '-s', action='store_true', default=False, help='dont run solar plots')

    args = parser.parse_args()

    if not args.skip_thermal:
        main(args.prefix, 'lw', tmin=args.tmin, tmax=args.tmax, verbose=args.verbose)

    if not args.skip_solar:
        main(args.prefix, 'sw', tmin=args.tmin, tmax=args.tmax, verbose=args.verbose)

if __name__ == '__main__':
    _main()
