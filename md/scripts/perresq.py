#!/usr/bin/env python3
"""
perresq.py -- per-residue fraction of native contacts, q_a(t), from an MD
trajectory, using the Best-Hummer-Eaton definition (PNAS 110:17874, 2013).

Global Q(t) is written alongside, both computed directly and rebuilt from the
per-residue terms as a check that the decomposition is exact.

Requires mdtraj, numpy and scipy.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

__version__ = "1.0"          # output format; unchanged by the 2026-11 tidy-up

DEF_BETA = 50.0     # nm^-1
DEF_LAMBDA = 1.8
DEF_CUTOFF = 0.45   # nm
DEF_SEQSEP = 3      # contacts need |res_i - res_j| > SEQSEP


def find_native_contacts(ref_xyz, res_of_atom, cutoff, seqsep):
    pairs = np.array(sorted(cKDTree(ref_xyz).query_pairs(cutoff)),
                     dtype=np.int64)
    if pairs.size == 0:
        return (np.empty(0, np.int64), np.empty(0, np.int64),
                np.empty(0, np.float64))
    ia, ja = pairs[:, 0], pairs[:, 1]
    keep = np.abs(res_of_atom[ia] - res_of_atom[ja]) > seqsep
    ia, ja = ia[keep], ja[keep]
    return ia, ja, np.linalg.norm(ref_xyz[ia] - ref_xyz[ja], axis=1)


def switching(d, r0, beta, lam):
    return 1.0 / (1.0 + np.exp(np.clip(beta * (d - lam * r0), -60.0, 60.0)))


def build_incidence(pi, pj, res_of_atom, n_res):
    """Accumulator mapping per-contact values onto both their residues."""
    rows = np.concatenate([np.arange(len(pi)), np.arange(len(pi))])
    cols = np.concatenate([res_of_atom[pi], res_of_atom[pj]])
    M = sparse.csr_matrix((np.ones(len(rows), np.float64), (rows, cols)),
                          shape=(len(pi), n_res))
    return lambda f: np.asarray((M.T @ f.T).T)


def compute(traj_files, top_file, out_prefix, selection="protein",
            cutoff=DEF_CUTOFF, seqsep=DEF_SEQSEP, beta=DEF_BETA,
            lam=DEF_LAMBDA, stride=1, chunk=100, ref_file=None, ref_frame=0,
            periodic=True, time_unit="ns", verbose=True):
    import mdtraj as md

    say = print if verbose else (lambda *a, **k: None)

    # ---- topology and atom selection --------------------------------------
    top = md.load_topology(top_file)
    sel_user = top.select(selection)
    if len(sel_user) == 0:
        raise SystemExit(f"selection '{selection}' matched no atoms in {top_file}")
    heavy_all = set(top.select_atom_indices("heavy").tolist())
    sel = np.array(sorted(set(sel_user.tolist()) & heavy_all), dtype=np.int64)
    if len(sel) == 0:
        raise SystemExit(f"selection '{selection}' contains no heavy atoms")

    residues = sorted({top.atom(int(i)).residue for i in sel},
                      key=lambda r: r.index)
    res_pos = {r.index: k for k, r in enumerate(residues)}
    n_res = len(residues)
    res_of_atom = np.array([res_pos[top.atom(int(i)).residue.index] for i in sel],
                           dtype=np.int64)
    resids = np.array([r.resSeq for r in residues], dtype=np.int64)
    resnames = [r.name for r in residues]
    chains = [r.chain.index for r in residues]

    say(f"  system      : {len(sel)} heavy atoms in {n_res} residues "
        f"(selection: {selection})")

    # ---- reference state --------------------------------------------------
    if ref_file:
        try:                       # pdb/gro/h5 carry their own topology
            ref = md.load(ref_file, atom_indices=sel)
        except (TypeError, ValueError):
            ref = md.load(ref_file, top=top_file, atom_indices=sel)
        say(f"  reference   : {os.path.basename(ref_file)} (frame 0)")
    else:
        ref = md.load_frame(traj_files[0], ref_frame, top=top_file,
                            atom_indices=sel)
        say(f"  reference   : frame {ref_frame} of "
            f"{os.path.basename(traj_files[0])}")
    ref_xyz = ref.xyz[0].astype(np.float64)

    # ---- native contact set ----------------------------------------------
    pi, pj, r0 = find_native_contacts(ref_xyz, res_of_atom, cutoff, seqsep)
    n_pairs = len(pi)
    if n_pairs == 0:
        raise SystemExit(
            "no native contacts found -- check --selection, --cutoff and that "
            "the reference structure is not garbled")
    n_a = np.zeros(n_res, np.float64)
    np.add.at(n_a, res_of_atom[pi], 1.0)
    np.add.at(n_a, res_of_atom[pj], 1.0)
    say(f"  native set  : |S| = {n_pairs} contacts  "
        f"(cutoff {cutoff*10:.1f} A, |i-j| > {seqsep})")
    say(f"                sum(n_a) = {int(n_a.sum())} = 2|S| "
        f"({'ok' if int(n_a.sum()) == 2 * n_pairs else 'MISMATCH'}); "
        f"{int((n_a == 0).sum())} residues have none (reported as NaN)")

    accumulate = build_incidence(pi, pj, res_of_atom, n_res)
    r0_row = r0[None, :]

    # ---- sweep the trajectory in chunks -----------------------------------
    num_blocks, q_direct_blocks, times = [], [], []
    n_frames = 0
    for tf in traj_files:
        for tr in md.iterload(tf, top=top_file, chunk=chunk, stride=stride,
                              atom_indices=sel):
            d = md.compute_distances(tr, np.column_stack([pi, pj]),
                                     periodic=periodic)
            f = switching(d.astype(np.float64), r0_row, beta, lam)
            q_direct_blocks.append(f.mean(axis=1))
            num_blocks.append(accumulate(f))
            times.append(np.asarray(tr.time, dtype=np.float64))
            n_frames += len(tr)
            if verbose:
                sys.stdout.write(f"\r  frames read : {n_frames}      ")
                sys.stdout.flush()
    if verbose:
        sys.stdout.write("\n")
    if n_frames == 0:
        raise SystemExit("no frames were read -- check the trajectory files")

    num = np.concatenate(num_blocks, axis=0)
    q_direct = np.concatenate(q_direct_blocks)
    t = np.concatenate(times)
    del num_blocks, q_direct_blocks

    # MDTraj reports time in ps; if the file carried none, the frame index is
    # used instead and the CSV header says so.
    have_time = np.any(np.diff(t) != 0) if len(t) > 1 else bool(t[0])
    if not have_time:
        t = np.arange(n_frames, dtype=np.float64) * stride
        t_label, t_out = "frame", t
    elif time_unit == "ns":
        t_label, t_out = "time_ns", t / 1000.0
    else:
        t_label, t_out = "time_ps", t

    # ---- q_a and reconstruction of global Q -------------------------------
    has = n_a > 0
    q_res = np.full((n_frames, n_res), np.nan)
    q_res[:, has] = num[:, has] / n_a[None, has]
    q_rebuilt = np.nansum(q_res[:, has] * n_a[None, has], axis=1) / n_a.sum()

    resid = float(np.max(np.abs(q_direct - q_rebuilt)))
    say(f"  Q(t)        : mean {q_direct.mean():.4f}  sd {q_direct.std():.4f}  "
        f"min {q_direct.min():.4f}")
    say(f"  consistency : max |Q_direct - Q_from_q_a| = {resid:.2e} "
        f"({'ok' if resid < 1e-6 else 'CHECK THIS'})")

    # ---- write ------------------------------------------------------------
    outdir = os.path.dirname(out_prefix)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    hdr_common = (
        f"# qres.py v{__version__} -- per-residue fraction of native contacts\n"
        f"# Best, Hummer & Eaton, PNAS 110:17874 (2013)\n"
        f"# topology={os.path.basename(top_file)} "
        f"traj={','.join(os.path.basename(x) for x in traj_files)}\n"
        f"# selection='{selection}' heavy_atoms={len(sel)} residues={n_res}\n"
        f"# cutoff={cutoff} nm  seqsep={seqsep}  beta={beta} nm^-1  "
        f"lambda={lam}  stride={stride}  periodic={periodic}\n"
        f"# reference={'file:' + os.path.basename(ref_file) if ref_file else 'frame ' + str(ref_frame)}"
        f"  n_native_contacts={n_pairs}\n")

    p_q = out_prefix + "_qres.csv"
    with open(p_q, "w") as fh:
        fh.write(hdr_common)
        fh.write("# columns after the first are residues, labelled by resSeq; "
                 "NaN = residue has no native contacts\n")
        fh.write(t_label + "," + ",".join(str(r) for r in resids) + "\n")
        for k in range(n_frames):
            row = ",".join("nan" if np.isnan(v) else f"{v:.5f}"
                           for v in q_res[k])
            fh.write(f"{t_out[k]:.4f},{row}\n")

    # statistics only over residues that actually form contacts
    q_mean = np.full(n_res, np.nan)
    q_sd = np.full(n_res, np.nan)
    q_min = np.full(n_res, np.nan)
    q_mean[has] = q_res[:, has].mean(axis=0)
    q_sd[has] = q_res[:, has].std(axis=0)
    q_min[has] = q_res[:, has].min(axis=0)
    p_s = out_prefix + "_qres_summary.csv"
    with open(p_s, "w") as fh:
        fh.write(hdr_common)
        fh.write("resid,resname,chain,n_native_contacts,q_mean,q_sd,q_min\n")
        for k in range(n_res):
            if n_a[k] == 0:
                fh.write(f"{resids[k]},{resnames[k]},{chains[k]},0,nan,nan,nan\n")
            else:
                fh.write(f"{resids[k]},{resnames[k]},{chains[k]},"
                         f"{int(n_a[k])},{q_mean[k]:.5f},{q_sd[k]:.5f},"
                         f"{q_min[k]:.5f}\n")

    p_g = out_prefix + "_Q_global.csv"
    with open(p_g, "w") as fh:
        fh.write(hdr_common)
        fh.write("# Q_direct   = mean of f over all native contacts\n"
                 "# Q_from_qres= sum(n_a * q_a) / sum(n_a), the contact-weighted\n"
                 "#              mean; equals Q_direct by construction. The\n"
                 "#              UNWEIGHTED mean over residues does not.\n")
        fh.write(f"{t_label},Q_direct,Q_from_qres,abs_diff\n")
        for k in range(n_frames):
            fh.write(f"{t_out[k]:.4f},{q_direct[k]:.6f},{q_rebuilt[k]:.6f},"
                     f"{abs(q_direct[k]-q_rebuilt[k]):.2e}\n")

    say(f"  written     : {os.path.basename(p_q)}, "
        f"{os.path.basename(p_s)}, {os.path.basename(p_g)}")

    return dict(resids=resids, q_res=q_res, q_direct=q_direct,
                q_rebuilt=q_rebuilt, n_a=n_a, t=t_out, t_label=t_label,
                n_pairs=n_pairs, q_mean=q_mean, q_sd=q_sd)


def selftest():
    """Check the decomposition on a case small enough to verify by hand."""
    contacts = [(1, 5, 0.97), (1, 5, 0.95), (1, 5, 0.96),
                (1, 5, 0.98), (2, 4, 0.10)]
    f = np.array([c[2] for c in contacts])
    n_res = 6
    pi = np.arange(0, 2 * len(contacts), 2)
    pj = pi + 1
    lut = np.zeros(2 * len(contacts), np.int64)
    for k, (a, b, _) in enumerate(contacts):
        lut[2 * k], lut[2 * k + 1] = a, b
    num = build_incidence(pi, pj, lut, n_res)(f[None, :])[0]
    den = np.zeros(n_res)
    np.add.at(den, lut[pi], 1.0)
    np.add.at(den, lut[pj], 1.0)
    has = den > 0
    q = np.full(n_res, np.nan)
    q[has] = num[has] / den[has]

    Q_direct = f.mean()
    Q_weighted = np.nansum(q[has] * den[has]) / den.sum()
    Q_plain = np.nanmean(q[has])

    print(f"  Q direct (mean of the {len(f)} f values)   = {Q_direct:.6f}")
    print(f"  Q rebuilt, contact-weighted             = {Q_weighted:.6f}"
          f"   diff {abs(Q_direct - Q_weighted):.1e}")
    print(f"  plain mean over residues (WRONG)        = {Q_plain:.6f}"
          f"   off by {Q_plain - Q_direct:+.4f}")
    print(f"  sum(n_a) = {int(den.sum())} = 2 x {len(f)} contacts   "
          f"{'ok' if int(den.sum()) == 2 * len(f) else 'MISMATCH'}")

    r0 = np.array([0.40])
    assert switching(np.array([0.40]), r0, DEF_BETA, DEF_LAMBDA)[0] > 0.98
    assert switching(np.array([2.00]), r0, DEF_BETA, DEF_LAMBDA)[0] < 1e-20
    assert abs(switching(np.array([DEF_LAMBDA * 0.40]), r0,
                         DEF_BETA, DEF_LAMBDA)[0] - 0.5) < 1e-12
    ok = abs(Q_direct - Q_weighted) < 1e-12
    print(f"  RESULT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="perresq.py", description=__doc__.strip().split("\n\n")[1],
        epilog="q depends only on internal distances, so the trajectory need "
               "not be fitted -- but molecules must be whole, e.g. "
               "gmx trjconv -pbc mol -ur compact.")
    ap.add_argument("-f", "--traj", nargs="+", metavar="FILE",
                    help="trajectory file(s): xtc, trr, dcd, nc, h5, ...")
    ap.add_argument("-s", "--top", metavar="FILE",
                    help="topology/structure: gro, pdb, prmtop, psf, ...")
    ap.add_argument("-o", "--out", default="qres", metavar="PREFIX")
    ap.add_argument("--selection", default="protein", metavar="STR",
                    help="MDTraj atom selection; heavy atoms are taken from it")
    ap.add_argument("--ref", metavar="FILE", default=None,
                    help="structure to use as the native state instead of a "
                         "trajectory frame")
    ap.add_argument("--ref-frame", type=int, default=0, metavar="N")
    ap.add_argument("--cutoff", type=float, default=DEF_CUTOFF, metavar="NM")
    ap.add_argument("--seqsep", type=int, default=DEF_SEQSEP, metavar="N",
                    help=f"require |i-j| > N residues (default: {DEF_SEQSEP})")
    ap.add_argument("--beta", type=float, default=DEF_BETA, metavar="X",
                    help=f"switching steepness, nm^-1 (default: {DEF_BETA})")
    ap.add_argument("--lam", "--lambda", dest="lam", type=float,
                    default=DEF_LAMBDA, metavar="X")
    ap.add_argument("--stride", type=int, default=1, metavar="N")
    ap.add_argument("--chunk", type=int, default=100, metavar="N",
                    help="frames held in memory at once (default: 100)")
    ap.add_argument("--nopbc", action="store_true",
                    help="do not apply the minimum image convention")
    ap.add_argument("--time-unit", choices=["ns", "ps"], default="ns")
    ap.add_argument("-q", "--quiet", action="store_true")
    ap.add_argument("--selftest", action="store_true",
                    help="verify the maths on a synthetic case and exit")
    ap.add_argument("--version", action="version",
                    version=f"perresq.py {__version__}")
    args = ap.parse_args(argv)

    if args.selftest:
        return selftest()
    if not args.traj or not args.top:
        ap.error("-f/--traj and -s/--top are required (or use --selftest)")
    for p in list(args.traj) + [args.top] + ([args.ref] if args.ref else []):
        if not os.path.exists(p):
            ap.error(f"file not found: {p}")

    compute(traj_files=args.traj, top_file=args.top, out_prefix=args.out,
            selection=args.selection, cutoff=args.cutoff, seqsep=args.seqsep,
            beta=args.beta, lam=args.lam, stride=args.stride, chunk=args.chunk,
            ref_file=args.ref, ref_frame=args.ref_frame,
            periodic=not args.nopbc, time_unit=args.time_unit,
            verbose=not args.quiet)
    return 0


if __name__ == "__main__":
    sys.exit(main())
