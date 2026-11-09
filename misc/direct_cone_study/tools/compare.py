"""compare surface fluxes of tenstream runs against an mcdmda reference
usage: compare.py <ref.xmf> <run1.xmf> [<run2.xmf> ...]
the shadow mask is taken from the reference: edir < 0.5 * max(edir)
"""
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import xmf  # noqa: E402


def stats(r, ref, shadow):
    tot, t0 = r["edir"] + r["edn"], ref["edir"] + ref["edn"]
    return dict(
        edir=np.mean(r["edir"]),
        edn=np.mean(r["edn"]),
        eup=np.mean(r["eup"]),
        sh_edir=np.mean(r["edir"][shadow]),
        sh_edn=np.mean(r["edn"][shadow]),
        sh_tot=np.mean(tot[shadow]),
        lit_tot=np.mean(tot[~shadow]),
        rmse_tot=np.sqrt(np.mean((tot - t0) ** 2)),
        rmse_eup=np.sqrt(np.mean((r["eup"] - ref["eup"]) ** 2)),
    )


def read_avg(name):
    """read one file or average a comma separated list of files"""
    fs = name.split(",")
    rs = [xmf.read(f) for f in fs]
    return {k: np.mean([r[k] for r in rs], axis=0) for k in rs[0]}, rs


def main(ref_name, names, profile_x=None):
    ref, refs = read_avg(ref_name)
    if len(refs) > 1:
        sh = refs[0]["edir"] < 0.5 * np.nanmax(refs[0]["edir"])
        tots = [np.mean((r["edir"] + r["edn"])[sh]) for r in refs]
        print(f"reference: mean of {len(refs)} runs, run-to-run std of shadow total: {np.std(tots):.2f}")
    shadow = ref["edir"] < 0.5 * np.nanmax(ref["edir"])
    hdr = f"{'run':28s} {'edir':>7s} {'edn':>6s} {'eup':>6s} | {'sh_edir':>7s} {'sh_edn':>6s} {'sh_tot':>6s} {'lit_tot':>7s} | {'rmse_tot':>8s} {'rmse_eup':>8s}"
    print(f"shadow px {shadow.sum()} of {shadow.size}")
    print(hdr)
    for n in [ref_name] + names:
        s = stats(read_avg(n)[0], ref, shadow)
        print(
            f"{os.path.basename(n.split(',')[0]):28s} {s['edir']:7.2f} {s['edn']:6.2f} {s['eup']:6.2f} | "
            f"{s['sh_edir']:7.2f} {s['sh_edn']:6.2f} {s['sh_tot']:6.2f} {s['lit_tot']:7.2f} | {s['rmse_tot']:8.2f} {s['rmse_eup']:8.2f}"
        )
    if profile_x is not None:
        np.set_printoptions(linewidth=250, precision=0, suppress=True)
        print(f"total down along y at x={profile_x}")
        for n in [ref_name] + names:
            r = read_avg(n)[0]
            print(f"{os.path.basename(n)[:12]:12s}", (r["edir"] + r["edn"])[profile_x, :])


if __name__ == "__main__":
    px = None
    args = sys.argv[1:]
    if args and args[0].startswith("--x="):
        px = int(args[0][4:])
        args = args[1:]
    main(args[0], args[1:], px)
