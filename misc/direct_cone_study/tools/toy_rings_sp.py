import sys, os, itertools, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toy
toy.RING_SPACING = sys.argv[3]
N = int(sys.argv[1]); out = open(sys.argv[2], 'w')
layouts = [(2, 2), (3, 2), (4, 2), (4, 4)]
methods = []
for tb, se in layouts:
    L = dict(reset="rings", nring_tb=tb, side_elev=se)
    ns = 2 * tb + 8 * se
    methods += [(f"r{tb}e{se}_n{ns}+c20", dict(L, cone_deg=20))]
out.write("tau_lay,sza,method,R,T,A,R_ex,T_ex,A_ex,alay_rms,alay_ex_mean\n")
for tl, sza in itertools.product((0.2, 1.0, 4.0), (0, 30, 60, 75)):
    t0 = time.time()
    base = dict(dx=200.0, sza=sza, albedo=0.1, N=N, g=0.9, tau=10 * tl, ncld=10)
    ex = toy.slab(reset="none", **base)
    for name, kw in methods:
        r = toy.slab(**base, **kw)
        e = (r['Alay'] - ex['Alay'])[5:15]
        out.write(f"{tl},{sza},{name},{r['R']:.5f},{r['T']:.5f},{r['A']:.5f},{ex['R']:.5f},{ex['T']:.5f},{ex['A']:.5f},"
                  f"{np.sqrt(np.mean(e**2)):.6f},{np.mean(ex['Alay'][5:15]):.6f}\n")
    out.flush(); print(tl, sza, f"{time.time()-t0:.1f}s", flush=True)
