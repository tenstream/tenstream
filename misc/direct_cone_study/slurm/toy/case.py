"""one case of the large toy study: python3 case.py <case_index> <outdir>"""
import sys, os, itertools, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toy
N = 40000
cases = list(itertools.product((0.85, 0.9), (0.1, 0.3, 1.0, 3.0), (0, 20, 40, 60, 75), (100.0, 200.0, 1000.0)))
methods = [("delta", dict(delta=True)), ("cone0", dict())]
methods += [(f"c{c}p0", dict(cone_deg=c)) for c in (15, 20, 25, 30, 35, 40, 50)]
methods += [(f"d+c{c}", dict(delta=True, cone_deg=c)) for c in (10, 20, 30)]
methods += [(f"c{c}p1", dict(cone_deg=c, phase1=True)) for c in (15, 20)]
methods += [("exact_seed2", dict(reset="none", seed=2))]
i = int(sys.argv[1]); out = sys.argv[2]
g, tl, sza, dx = cases[i]
base = dict(dx=dx, sza=sza, albedo=0.1, N=N, g=g, tau=10 * tl, ncld=10)
t0 = time.time()
ex = toy.slab(reset="none", **base)
res = dict(case=dict(g=g, tau_lay=tl, sza=sza, dx=dx), exact={k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in ex.items()}, methods={})
for m, kw in methods:
    kw2 = dict(base); kw2.update(kw)
    r = toy.slab(**kw2)
    res["methods"][m] = {k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in r.items()}
res["time"] = time.time() - t0
json.dump(res, open(os.path.join(out, f"case{i:03d}.json"), "w"))
