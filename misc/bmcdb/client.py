#!/usr/bin/env python3
import peewee as pg
import numpy as np
import os

np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)

DB = pg.PostgresqlDatabase("bmcdb", user="bmc", password="bmcdb", host="kiwi.meteo.physik.uni-muenchen.de", port=5432)


class Entry(pg.Model):
    inp = pg.BlobField()
    C = pg.BlobField()
    S = pg.BlobField()
    T = pg.BlobField()
    Ctol = pg.BlobField()
    Stol = pg.BlobField()
    Ttol = pg.BlobField()

    class Meta:
        database = DB


DB.create_tables(
    [
        Entry,
    ]
)
DB.close()

def select_mc_func(Ndir, Ndif, use_verts=False):
    from py_boxmc import m_py_boxmc as P
    if use_verts:
        if (Ndir==3) and (Ndif==10):
            return P.py_get_coeff_3_10_verts
    else:
        if (Ndir==3) and (Ndif==10):
            return P.py_get_coeff_3_10
        if (Ndir==8) and (Ndif==10):
            return P.py_get_coeff_8_10

    raise ValueError(f"No suitable solver configured for {Ndir=} and {Ndif=} and {use_verts=}")

def run(kabs, ksca, g, dz, phi, theta, dx=1.0, atol=1e-4, rtol=1e-2, tanx=0, tany=0, Ndir=8, Ndif=10):
    from mpi4py import MPI

    comm = MPI.COMM_WORLD.py2f()

    inp = np.array([kabs, ksca, g, dz, phi, theta, tanx, tany])
    print(f"Computing {inp}", flush=True)

    C = np.empty((Ndif, Ndif))
    S = np.empty((Ndir, Ndif))
    T = np.empty((Ndir, Ndir))
    Ctol = np.empty((Ndif, Ndif))
    Stol = np.empty((Ndir, Ndif))
    Ttol = np.empty((Ndir, Ndir))

    verts = np.zeros((8,3))
    verts[0] = [0,0,0]
    verts[1] = [1,0,0+tanx]
    verts[2] = [0,1,0+tany]
    verts[3] = [1,1,0+tanx+tany]
    verts[4] = [0,0,dz]
    verts[5] = [1,0,dz+tanx]
    verts[6] = [0,1,dz+tany]
    verts[7] = [1,1,dz+tanx+tany]

    mc_func = select_mc_func(Ndir, Ndif, use_verts=True)

    if True:
        for src in range(Ndir):
            _S, _T, _S_tol, _T_tol = mc_func(comm, dx, dx, dz, verts.flatten(), [kabs, ksca, g], src + 1, True, phi, theta, atol, rtol)
            S[src, :] = _S
            Stol[src, :] = _S_tol
            T[src, :] = _T
            Ttol[src, :] = _T_tol
    else:
        S   [:] = -1
        Stol[:] = -1
        T   [:] = -1
        Ttol[:] = -1

    for src in range(Ndif):
        _S, _T, _S_tol, _T_tol = mc_func(
            comm, dx, dx, dz, verts.flatten(), [kabs, ksca, g], src + 1, False, -1.0, -1.0, inp_atol=atol, inp_rtol=rtol
        )
        C[src, :] = _S
        Ctol[src, :] = _S_tol

    ret = dict(inp=inp, C=C, S=S, T=T, Ctol=Ctol, Stol=Stol, Ttol=Ttol)

    from pprint import pprint
    pprint(ret)

    return ret



def rlog(s=1e-10, e=1e2):
    lb = np.log10(s)
    ub = np.log10(e)
    sample = np.random.random() * (ub - lb) + lb
    return 10**sample


def sample(i=None, **kwargs):
    print(f"Running sample {i}", flush=True)
    param = dict(
        kabs=rlog(),
        ksca=rlog() * np.random.randint(2), # half the time use ksca=0 which is important to sample for w_0 being 0, particularly in the thermal spectral range
        g=(np.random.random()**2) * np.random.randint(2),
        dz=np.random.random() * (5.0 - 0.01) + 0.01,
        phi=np.random.random() * 90,
        theta=np.random.random() * 90,
        tanx=(np.random.random()**2 * np.random.choice([-1, 1])) * np.tan(np.deg2rad(45)),
        tany=(np.random.random()**2 * np.random.choice([-1, 1])) * np.tan(np.deg2rad(45)),
    )
    param.update(kwargs)

    return run(**param)


def gen_dataset(Ndir=8, Ndif=10):
    import xarray as xr
    from playhouse.shortcuts import model_to_dict

    N = Entry.select().count()
    first = Entry.select().first()
    for k, v in model_to_dict(first).items():
        try:
            print(k, 'shape', np.frombuffer(v).shape)
        except:
            print(k, v)


    D = xr.Dataset(
        data_vars=dict(
            inp=(("N", "inpvars"),              np.empty((N, len(np.frombuffer(first.inp))), dtype=np.float32)),
            C=(("N", "dif_src", "dif_dst"),     np.empty((N, Ndif, Ndif), dtype=np.float32)),
            Ctol=(("N", "dif_src", "dif_dst"),  np.empty((N, Ndif, Ndif), dtype=np.float32)),
            S=(("N", "dir_src", "dif_dst"),     np.empty((N, Ndir, Ndif), dtype=np.float32)),
            Stol=(("N", "dir_src", "dif_dst"),  np.empty((N, Ndir, Ndif), dtype=np.float32)),
            T=(("N", "dir_src", "dir_dst"),     np.empty((N, Ndir, Ndir), dtype=np.float32)),
            Ttol=(("N", "dir_src", "dir_dst"),  np.empty((N, Ndir, Ndir), dtype=np.float32)),
        ),
        coords=dict(
            inpvars=(("inpvars",), ["kabs", "ksca", "g", "dz", "phi", "theta", "tanx", "tany"]),
        ),
        attrs=dict(),
    )


    from alive_progress import alive_bar
    with alive_bar(N) as bar:
        for i, e in enumerate(Entry.select().iterator()):
            if i > N-1:
                break

            D["inp"][i, :] = np.frombuffer(e.inp)
            D["C"][i, :, :] = np.frombuffer(e.C).reshape((Ndif, Ndif))
            D["Ctol"][i, :, :] = np.frombuffer(e.Ctol).reshape((Ndif, Ndif))
            D["S"][i, :, :] = np.frombuffer(e.S).reshape((Ndir, Ndif))
            D["Stol"][i, :, :] = np.frombuffer(e.Stol).reshape((Ndir, Ndif))
            D["T"][i, :, :] = np.frombuffer(e.T).reshape((Ndir, Ndir))
            D["Ttol"][i, :, :] = np.frombuffer(e.Ttol).reshape((Ndir, Ndir))
            bar()
    return D


def _main():
    import argparse
    import glob

    parser = argparse.ArgumentParser(description="Compute random samples of bmc coeffs and store them in a database")
    parser.add_argument("--iter", "-N", type=int, default=100, help="number of bmc calls to run")
    parser.add_argument("--nproc", "-j", type=int, default=None, help="number of processes to use")
    parser.add_argument("--parallel", "-p", default=False, action="store_true", help="if running in parallel")
    parser.add_argument("--dry", "-d", default=False, action="store_true", help="dry run, dont run compute loop")
    parser.add_argument("--out", "-o", default=None, type=str, help="output filename to dump data from postgres to xarray file")
    parser.add_argument("--verbose", "-v", default=False, action="store_true", help="verbose output")
    parser.add_argument("--atol", "-a", default=1e-4, type=float, help="absolute tolerance for montecarlo")
    parser.add_argument("--rtol", "-r", default=1e-2, type=float, help="relative tolerance for montecarlo")
    parser.add_argument("--Ndir", default=8, type=int, help="Number of direct streams")
    parser.add_argument("--Ndif", default=10, type=int, help="Number of diffuse streams")

    args = parser.parse_args()

    if args.nproc is not None:
        args.parallel = True

    sample_args = dict(
            atol = args.atol,
            rtol = args.rtol,
            Ndir = args.Ndir,
            Ndif = args.Ndif,
            )

    if not args.dry:
        if args.parallel:
            from joblib import Parallel, delayed
            numproc = 2
            numproc = int(os.environ.get('SLURM_JOB_CPUS_PER_NODE', numproc))
            numproc = int(os.environ.get('OMP_NUM_THREADS', numproc))
            if args.nproc is not None:
                numproc = args.nproc
            print(f"Using {numproc} processes")
            parallel = Parallel(n_jobs=numproc, return_as="generator")
            output_generator = parallel(delayed(sample)(i, **sample_args) for i in range(args.iter))
            for ret in output_generator:
                print(f"Adding results for {ret}")
                Entry(**{k: v.tobytes() for k, v in ret.items()}).save()
                DB.close()
        else:
            i = 0
            while i < args.iter:
                ret = sample(**sample_args)
                Entry(**{k: v.tobytes() for k, v in ret.items()}).save()
                DB.close()
                i += 1

    print(f"Have {Entry.select().count()} entries in db")
    if args.verbose:
        for e in Entry.select().iterator():
            print(e, np.frombuffer(e.inp))  # , np.frombuffer(e.C))
    DB.close()

    if args.out:
        if args.verbose:
            print(f"Dumping Data to output file {args.out}")
        D = gen_dataset(args.Ndir, args.Ndif)
        print(D)

        if args.out.endswith('zarr'):
            import zarr
            compressor = zarr.Blosc(cname="zstd", clevel=3, shuffle=2)
            enc = {x: {"compressor": compressor} for x in D}
            D.to_zarr(args.out, encoding=enc)
        else:
            D.to_netcdf(args.out)


if __name__ == "__main__":
    _main()
