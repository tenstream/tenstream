"""parse vertical profiles printed by ex_pprts_ex1 -v (rank 0, column x=1,y=1) and compare slab runs
usage: slab.py <ref.log> <run.log> ...
prints TOA reflectance, srfc transmittance (dir+diff), atmospheric absorptance and the max abs error of
the down/up profiles relative to the reference (normalized by the TOA incoming flux)
"""
import re
import sys
import os
import numpy as np


def parse(fname):
    prof = {"fdir": {}, "fdn": {}, "fup": {}}
    pat = re.compile(r"^\s*(\d+)\s+\x1b\[\d+m(fdir|fdn |fup ) x=1\x1b\[0m\s+([-\d.Ee+]+)")
    for line in open(fname, errors="replace"):
        m = pat.match(line)
        if m:
            k, var, val = int(m.group(1)), m.group(2).strip(), float(m.group(3))
            prof[var].setdefault(k, val)  # first rank only
    return {v: np.array([d[k] for k in sorted(d)]) for v, d in prof.items()}


def main(ref, runs):
    r0 = parse(ref)
    S = r0["fdir"][0]
    print(f"{'run':20s} {'R_toa':>7s} {'T_srfc':>7s} {'T_dir':>7s} {'A_atm':>7s} | {'maxerr_dn':>9s} {'maxerr_up':>9s}  [% of S0*mu]")
    for n in [ref] + runs:
        p = parse(n)
        dn = p["fdir"] + p["fdn"]
        R = p["fup"][0] / S
        T = dn[-1] / S
        A = (dn[0] - p["fup"][0] - (dn[-1] - p["fup"][-1])) / S
        e_dn = np.max(np.abs(dn - (r0["fdir"] + r0["fdn"]))) / S
        e_up = np.max(np.abs(p["fup"] - r0["fup"])) / S
        print(
            f"{os.path.basename(n):20s} {100*R:7.2f} {100*T:7.2f} {100*p['fdir'][-1]/S:7.2f} {100*A:7.2f} | {100*e_dn:9.2f} {100*e_up:9.2f}"
        )


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
