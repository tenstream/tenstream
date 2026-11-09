import json, glob, sys, collections
import numpy as np
D = sys.argv[1]
res = [json.load(open(f)) for f in sorted(glob.glob(D + '/case*.json'))]
meths = list(res[0]['methods'].keys())
def errs(r, m):
    ex, x = r['exact'], r['methods'][m]
    rt = 100 * max(abs(x['R'] - ex['R']), abs(x['T'] - ex['T']))
    a0 = np.array(ex['Alay'][5:15]); a = np.array(x['Alay'][5:15])
    ae = 100 * np.sqrt(np.mean((a - a0) ** 2)) / np.mean(a0)
    return rt, ae
E = {m: np.array([errs(r, m) for r in res]) for m in meths}
key = lambda r, k: r['case'][k]
print(f"{len(res)} cases;  metric: RT = max(|dR|,|dT|) [pts], AB = rel rms layer absorption err [%]")
print(f"{'method':12s} {'RT mean':>8s} {'RT p90':>7s} {'RT max':>7s} {'AB mean':>8s} {'AB p90':>7s} {'AB max':>7s}")
for m in meths:
    e = E[m]
    print(f"{m:12s} {e[:,0].mean():8.2f} {np.percentile(e[:,0],90):7.2f} {e[:,0].max():7.2f} {e[:,1].mean():8.2f} {np.percentile(e[:,1],90):7.2f} {e[:,1].max():7.2f}")
for k in ('sza', 'tau_lay', 'dx', 'g'):
    vals = sorted({key(r, k) for r in res})
    for col, lab in ((0, 'RT'), (1, 'AB')):
        print(f"--- {lab} mean by {k}")
        print(f"{k:>8s} " + ' '.join(f'{m[:7]:>7s}' for m in meths))
        for v in vals:
            sel = [i for i, r in enumerate(res) if key(r, k) == v]
            print(f"{v:8g} " + ' '.join(f"{E[m][sel, col].mean():7.2f}" for m in meths))
