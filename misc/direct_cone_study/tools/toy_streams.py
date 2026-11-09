import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, toy
for sza in (0, 30, 60):
    base = dict(dx=200.0, sza=sza, albedo=0.1, N=20000)
    ex = toy.slab(reset="none", **base)
    def rep(n, r):
        e = 100*np.sqrt(np.mean(((r['Alay']-ex['Alay'])[5:15])**2))/np.mean(ex['Alay'][5:15])
        print(f"  {n:16s} dR {100*(r['R']-ex['R']):+5.1f} dT {100*(r['T']-ex['T']):+5.1f} AB {e:5.1f}", flush=True)
    print(f"SZA {sza}: exact R {100*ex['R']:.1f} T {100*ex['T']:.1f}")
    rep("3_10 delta", toy.slab(delta=True, **base))
    rep("3_10 cone25", toy.slab(cone_deg=25, **base))
    rep("3_30 nods", toy.slab(reset="3_30", **base))
    rep("3_30 delta", toy.slab(reset="3_30", delta=True, **base))
    rep("3_30 cone25", toy.slab(reset="3_30", cone_deg=25, **base))
    rep("suncap41 nods", toy.slab(reset="suncap", **base))
    rep("suncap25 nods", toy.slab(reset="suncap", suncap_deg=25, **base))
    rep("suncap41 delta", toy.slab(reset="suncap", delta=True, **base))
