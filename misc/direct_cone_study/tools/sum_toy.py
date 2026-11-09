import sys, csv, collections
import numpy as np
rows = list(csv.DictReader(open(sys.argv[1])))
agg = collections.defaultdict(list)
for r in rows:
    f = {k: float(v) for k, v in r.items() if k != 'method'}
    agg[r['method']].append((100 * max(abs(f['R'] - f['R_ex']), abs(f['T'] - f['T_ex'])), 100 * f['alay_rms'] / f['alay_ex_mean'],
                             100 * (f['R'] - f['R_ex']), float(r['sza'])))
print(f"{'method':22s} {'RT':>5s} {'AB%':>6s} {'dR':>6s} | {'RT<=30':>6s} {'AB<=30':>6s} | {'RT>=60':>6s} {'AB>=60':>6s}")
for m, v in agg.items():
    a = np.array(v); lo = a[:, 3] <= 30; hi = ~lo
    print(f"{m:22s} {a[:,0].mean():5.2f} {a[:,1].mean():6.2f} {a[:,2].mean():+6.2f} | {a[lo,0].mean():6.2f} {a[lo,1].mean():6.2f} | {a[hi,0].mean():6.2f} {a[hi,1].mean():6.2f}")
