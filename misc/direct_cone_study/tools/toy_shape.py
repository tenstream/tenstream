import sys, os, itertools, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toy
N = int(sys.argv[1]); out = open(sys.argv[2], 'w'); tl = float(sys.argv[3])
methods = [("3_10", {}), ("3_10+delta", dict(delta=True)), ("3_10+c20", dict(cone_deg=20)), ("3_10+d+c20", dict(delta=True, cone_deg=20)),
           ("r3_unsplit", dict(reset="rings", nring_tb=3, side_elev=-1))]
for p, q in itertools.product((1, 2, 3, 5), (0, 1, 2, 3)):
    if p == 1 and q == 0:
        continue
    L = dict(reset="shaped", p_tb=p, q_side=q)
    methods += [(f"p{p}q{q}", L), (f"p{p}q{q}+delta", dict(L, delta=True)), (f"p{p}q{q}+c20", dict(L, cone_deg=20))]
out.write("tau_lay,sza,method,R,T,A,R_ex,T_ex,A_ex,alay_rms,alay_ex_mean\n")
for sza in (0, 30, 60, 75):
    base = dict(dx=200.0, sza=sza, albedo=0.1, N=N, g=0.9, tau=10 * tl, ncld=10)
    ex = toy.slab(reset="none", **base)
    for name, kw in methods:
        r = toy.slab(**base, **kw)
        e = (r['Alay'] - ex['Alay'])[5:15]
        out.write(f"{tl},{sza},{name},{r['R']:.5f},{r['T']:.5f},{r['A']:.5f},{ex['R']:.5f},{ex['T']:.5f},{ex['A']:.5f},"
                  f"{np.sqrt(np.mean(e**2)):.6f},{np.mean(ex['Alay'][5:15]):.6f}\n")
    out.flush()
