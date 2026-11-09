"""compare absorption profiles / TOA eup printed by ex_pprts_box_cld (v2) against averaged mcdmda runs
usage: abso3d.py <ref1.log,ref2.log,...> <run.log> ...
"""
import sys, os, re
import numpy as np


def parse(fn):
    toa, prof = None, []
    for line in open(fn, errors="replace"):
        if "TOA eup domain mean" in line:
            toa = [float(v) for v in line.split()[-2:]]
        m = re.match(r"\s*abso_profile\s+(\d+)\s+(\S+)\s+(\S+)", line)
        if m:
            prof.append([float(m.group(2)), float(m.group(3))])
    return np.array(toa), np.array(prof)


def avg(names):
    rs = [parse(n) for n in names.split(",")]
    return np.mean([r[0] for r in rs], axis=0), np.mean([r[1] for r in rs], axis=0), rs


ref_toa, ref_p, rs = avg(sys.argv[1])
cld = slice(5, 15)  # cloud layers 6..15
dz = 100.
print(f"reference: {len(rs)} runs; TOA eup dom/cld {ref_toa[0]:.2f}/{ref_toa[1]:.2f}; "
      f"cloud abs (column sum W/m2) dom/cld {dz*ref_p[cld,0].sum():.2f}/{dz*ref_p[cld,1].sum():.2f}")
if len(rs) > 1:
    sp = np.std([dz * r[1][cld, 1].sum() for r in rs])
    pe = np.mean([np.sqrt(np.mean((r[1][cld, 1] - ref_p[cld, 1])**2)) / np.mean(ref_p[cld, 1]) * 100 for r in rs])
    print(f"  run-to-run std of cloud column abs {sp:.3f} W/m2, single run rel rms profile deviation {pe:.2f}%")
print(f"{'run':16s} {'dTOAdom':>8s} {'dTOAcld':>8s} | {'dAbs_dom':>8s} {'dAbs_cld':>8s} | {'rmsP_dom%':>9s} {'rmsP_cld%':>9s} {'belowcld%':>9s}")
for n in sys.argv[2:]:
    toa, p = parse(n)
    e_dom = 100 * np.sqrt(np.mean((p[cld, 0] - ref_p[cld, 0])**2)) / np.mean(ref_p[cld, 0])
    e_cld = 100 * np.sqrt(np.mean((p[cld, 1] - ref_p[cld, 1])**2)) / np.mean(ref_p[cld, 1])
    below = slice(15, None)
    e_bel = 100 * (p[below, 0].sum() - ref_p[below, 0].sum()) / ref_p[below, 0].sum()
    print(f"{os.path.basename(n)[:16]:16s} {toa[0]-ref_toa[0]:+8.2f} {toa[1]-ref_toa[1]:+8.2f} | "
          f"{dz*(p[cld,0].sum()-ref_p[cld,0].sum()):+8.2f} {dz*(p[cld,1].sum()-ref_p[cld,1].sum()):+8.2f} | {e_dom:9.2f} {e_cld:9.2f} {e_bel:+9.2f}")
