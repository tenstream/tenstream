import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import toy
for sza in (0, 60):
    base = dict(dx=200.0, sza=sza, albedo=0.1, N=20000)
    ex = toy.slab(reset="none", **base)
    print(f"SZA {sza}: exact R {100*ex['R']:.1f} T {100*ex['T']:.1f}")
    rows = [("delta", dict(delta=True))]
    for c in (20, 30, 45, 60, 75):
        rows.append((f"cone{c} ph0", dict(cone_deg=c)))
        rows.append((f"cone{c} ph1", dict(cone_deg=c, phase1=True)))
        rows.append((f"cone{c} ph1all", dict(cone_deg=c, phase1=True, phase1_scattered_only=False)))
    for name, kw in rows:
        r = toy.slab(**base, **kw)
        print(f"  {name:14s} R {100*r['R']:5.1f} T {100*r['T']:5.1f} Tdir {100*r['Tdir']:5.1f} A {100*r['A']:5.2f} | dR {100*(r['R']-ex['R']):+5.1f} dT {100*(r['T']-ex['T']):+5.1f}")
