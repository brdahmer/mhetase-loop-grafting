#!/usr/bin/env python3
"""
qregion.py -- regional native-contact retention from the output of perresq.py.

Reads <run>_qres.csv (per-residue q over time) and <run>_qres_summary.csv
(per-residue native contact count) and reports, for one residue range, the
mean q of each variant at each temperature and the difference of differences.
This is the contrast behind Figure 4B.

The two averages answer different questions:

  --weight residues   plain mean over residues: "the region's residues retain,
                      on average, X% of their native contacts"
  --weight contacts   contacts retained / contacts present, each residue
                      counting in proportion to its own n_a: "the region
                      retains X% of its native contacts"

perresq.py makes the same distinction for global Q; see the header of any
*_Q_global.csv. Chimeric residue numbers are mapped onto WT numbering (+9
above residue 169) so one range applies to both variants.

  python3 qregion.py --dir production/q_res --region 163 182
  python3 qregion.py --dir production/q_res --region 128 140 --weight contacts
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np


def chim_to_wt(resid):
    return resid if resid <= 169 else resid + 9


def read_run(d, run, is_chim, window):
    """(resids on WT numbering, mean q over the window, n_a)."""
    na = {}
    with open(os.path.join(d, f"{run}_qres_summary.csv")) as fh:
        for line in fh:
            if not line.startswith(("#", "resid")):
                f = line.split(",")
                na[int(f[0])] = float(f[3])

    with open(os.path.join(d, f"{run}_qres.csv")) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    resids = np.array([int(h) for h in lines[0].split(",")[1:]])

    total, n = np.zeros(len(resids)), 0
    for line in lines[1:]:
        f = line.split(",")
        if float(f[0]) >= window:
            total += [np.nan if v.strip() == "nan" else float(v) for v in f[1:]]
            n += 1
    if n == 0:
        raise SystemExit(f"{run}: no frames at or after {window} ns")

    counts = np.array([na[int(r)] for r in resids])
    if is_chim:
        resids = np.array([chim_to_wt(r) for r in resids])
    return resids, total / n, counts


def region_mean(resids, q, na, lo, hi, weighted):
    m = (resids >= lo) & (resids <= hi) & np.isfinite(q)
    if not m.any():
        return np.nan
    return float(np.sum(q[m] * na[m]) / np.sum(na[m])) if weighted \
        else float(np.mean(q[m]))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="qregion.py",
                                 description=__doc__.strip().split("\n\n")[1])
    ap.add_argument("--dir", required=True, help="directory of perresq.py output")
    ap.add_argument("--region", nargs=2, type=int, required=True,
                    metavar=("LO", "HI"), help="residue range, WT numbering")
    ap.add_argument("--window", type=float, default=150.0, metavar="NS",
                    help="average from this time on (default: 150)")
    ap.add_argument("--weight", choices=["residues", "contacts"],
                    default="residues", help="how to average (default: residues)")
    ap.add_argument("--wt-prefix", default="wt")
    ap.add_argument("--chim-prefix", default="graft")
    ap.add_argument("--temps", nargs=2, type=int, default=[303, 363],
                    metavar=("LOW", "HIGH"))
    ap.add_argument("--boot", type=int, default=10000,
                    help="bootstrap resamples (default: 10000)")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)

    lo, hi = a.region
    tlo, thi = a.temps
    weighted = a.weight == "contacts"

    means = {}
    for prefix, label in ((a.wt_prefix, "WT"), (a.chim_prefix, "CHIM")):
        for T in (tlo, thi):
            runs = sorted(os.path.basename(p)[:-len("_qres.csv")]
                          for p in glob.glob(
                              os.path.join(a.dir, f"{prefix}_{T}_*_qres.csv")))
            if not runs:
                raise SystemExit(f"no runs matched {prefix}_{T}_*_qres.csv "
                                 f"in {a.dir}")
            vals = [region_mean(*read_run(a.dir, r, label == "CHIM", a.window),
                                lo, hi, weighted) for r in runs]
            means[(label, T)] = np.array(vals)
            print(f"  {label:4s} {T} K  " +
                  "  ".join(f"{r.split('_')[-1]}={v:.4f}"
                            for r, v in zip(runs, vals)))

    m = {k: v.mean() for k, v in means.items()}
    print(f"\nregion {lo}-{hi}, {a.window:g}-end ns, weighted by "
          f"{'native contacts' if weighted else 'residues equally'}\n")
    print(f'{"":6s} {tlo:>9d} K {thi:>9d} K {"change":>10s}')
    for label in ("WT", "CHIM"):
        x, y = m[(label, tlo)], m[(label, thi)]
        print(f"{label:6s} {x:11.4f} {y:11.4f} {100*(y-x)/x:9.1f} %")

    rng = np.random.default_rng(a.seed)
    bs = np.empty(a.boot)
    for i in range(a.boot):
        s = {k: rng.choice(v, v.size, replace=True).mean()
             for k, v in means.items()}
        bs[i] = ((s[("CHIM", thi)] - s[("CHIM", tlo)])
                 - (s[("WT", thi)] - s[("WT", tlo)]))
    ci = np.percentile(bs, [2.5, 97.5])
    dWT = m[("WT", thi)] - m[("WT", tlo)]
    dCH = m[("CHIM", thi)] - m[("CHIM", tlo)]
    print(f"\nDelta WT   {dWT:+.4f}")
    print(f"Delta CHIM {dCH:+.4f}")
    print(f"DD         {dCH-dWT:+.4f}  95% CI [{ci[0]:+.3f}, {ci[1]:+.3f}]"
          f"   {'excludes' if ci[0] > 0 or ci[1] < 0 else 'includes'} zero")
    return 0


if __name__ == "__main__":
    sys.exit(main())
