import sys, os, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toy
N = 40000
meths = [("delta", dict(delta=True)), ("c20p0", dict(cone_deg=20)), ("c30p0", dict(cone_deg=30)), ("d+c30", dict(delta=True, cone_deg=30)),
         ("exact_seed2", dict(reset="none", seed=2))]
res = {m: [] for m, _ in meths}
for g, tl, sza in itertools.product((0.85,), (0.2, 1.0, 4.0), (0, 30, 60, 75)):
    base = dict(dx=200.0, sza=sza, albedo=0.1, N=N, g=g, tau=10 * tl, ncld=10)
    ex = toy.slab(reset="none", **base)
    line = f"tau_lay {tl} sza {sza}:"
    for m, kw in meths:
        kw2 = dict(base); kw2.update(kw)
        r = toy.slab(**kw2)
        e = 100 * np.sqrt(np.mean(((r['Alay'] - ex['Alay'])[5:15])**2)) / np.mean(ex['Alay'][5:15])
        res[m].append(e)
        line += f" {m} {e:5.2f}"
    print(line, flush=True)
print("mean rel rms layer abs err [%]:", {m: round(float(np.mean(v)), 2) for m, v in res.items()})
