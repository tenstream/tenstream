import sys, os, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toy
N = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
out = open(sys.argv[2] if len(sys.argv) > 2 else 'toy_robust.csv', 'w')
methods = [("delta", dict(delta=True)), ("cone0", dict()),
           ("c30p0", dict(cone_deg=30)), ("c45p0", dict(cone_deg=45)), ("c60p0", dict(cone_deg=60)),
           ("c30p1", dict(cone_deg=30, phase1=True))]
out.write("g,tau_lay,sza,dx,method,R,T,A,Tdir,R_ex,T_ex,A_ex\n")
for g, tl, sza, dx in itertools.product((0.7, 0.85, 0.9), (0.2, 1.0, 4.0), (0, 30, 60, 75), (200.0, 1000.0)):
    t0 = time.time()
    base = dict(dx=dx, sza=sza, albedo=0.1, N=N, g=g, tau=10 * tl, ncld=10)
    ex = toy.slab(reset="none", **base)
    for name, kw in methods:
        r = toy.slab(**base, **kw)
        out.write(f"{g},{tl},{sza},{dx},{name},{r['R']:.5f},{r['T']:.5f},{r['A']:.5f},{r['Tdir']:.5f},{ex['R']:.5f},{ex['T']:.5f},{ex['A']:.5f}\n")
    out.flush()
    print(g, tl, sza, dx, f"{time.time()-t0:.1f}s", flush=True)
