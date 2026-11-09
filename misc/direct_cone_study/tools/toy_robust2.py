import sys, os, itertools, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toy
N = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
out = open(sys.argv[2] if len(sys.argv) > 2 else 'toy_robust2.csv', 'w')
methods = [("delta", dict(delta=True))]
methods += [(f"c{c}p0", dict(cone_deg=c)) for c in (10, 15, 20, 30, 40, 50)]
methods += [(f"d+c{c}", dict(delta=True, cone_deg=c)) for c in (10, 20, 30)]
out.write("g,tau_lay,sza,dx,method,R,T,A,Tdir,R_ex,T_ex,A_ex,alay_rms,alay_ex_mean\n")
for g, tl, sza, dx in itertools.product((0.85, 0.9), (0.2, 1.0, 4.0), (0, 30, 60, 75), (200.0, 1000.0)):
    t0 = time.time()
    base = dict(dx=dx, sza=sza, albedo=0.1, N=N, g=g, tau=10 * tl, ncld=10)
    ex = toy.slab(reset="none", **base)
    for name, kw in methods:
        r = toy.slab(**base, **kw)
        e = (r['Alay'] - ex['Alay'])[5:15]
        out.write(f"{g},{tl},{sza},{dx},{name},{r['R']:.5f},{r['T']:.5f},{r['A']:.5f},{r['Tdir']:.5f},"
                  f"{ex['R']:.5f},{ex['T']:.5f},{ex['A']:.5f},{np.sqrt(np.mean(e**2)):.6f},{np.mean(ex['Alay'][5:15]):.6f}\n")
    out.flush()
    print(g, tl, sza, dx, f"{time.time()-t0:.1f}s", flush=True)
