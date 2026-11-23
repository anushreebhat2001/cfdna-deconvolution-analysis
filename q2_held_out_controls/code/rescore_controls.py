#!/usr/bin/env python3
"""
Re-score existing deconvolution estimates against a subset of controls.

Reads the all_estimates.tsv that collect_all.py already writes, plus the
splits.tsv from make_splits.py, and recomputes the per-bin limit-of-detection
AUROC using only the controls a condition names as negatives.

No mixtures are built and nothing is deconvolved: this reuses estimates that
already exist. Its purpose is to measure how much AUROC moves purely because
the negative set is smaller, so that the setup (b) result can be read as
leakage rather than as a sample-size artifact.

Binning follows collect_all.py section 2 exactly: the 120 CRC virtual samples
(20 undiluted + 100 diluted) are sorted by est_TF and split into 10 bins of
12. That ordering does not depend on the controls, so bins are identical
across every condition and the comparison is like-for-like.

Usage
-----
  python3 rescore_controls.py \
      --estimates /scratch/mxf4642/results/summary/all_estimates.tsv \
      --splits    /scratch/$USER/mrd/data/splits.tsv \
      --out       /scratch/$USER/mrd/results/rescore

Needs numpy/pandas/sklearn, so run it inside the container.
"""

import argparse
import fnmatch
import os
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

N_BINS = 10


def auc(pos, neg):
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    return roc_auc_score(y, np.r_[pos, neg])


def load_conditions(splits_path, pattern, role):
    """Return {condition: [control sample names used as negatives]}."""
    sp = pd.read_csv(splits_path, sep="\t")
    for col in ("condition", "role", "sample"):
        if col not in sp.columns:
            sys.exit(f"{splits_path}: missing column {col!r}")

    sel = sp[(sp.role == role)
             & sp.condition.map(lambda c: fnmatch.fnmatch(c, pattern))]
    if sel.empty:
        have = sorted(sp.condition.unique())[:5]
        sys.exit(f"No conditions match {pattern!r} with role {role!r}. "
                 f"Conditions present (first few): {have}")
    return {c: sorted(g["sample"].tolist()) for c, g in sel.groupby("condition")}


def bin_table(crc):
    """Sort the CRC virtual samples by est_TF and split into N_BINS bins."""
    return np.array_split(crc.sort_values("est_TF"), N_BINS)


def _label(b):
    return f"{100 * b.est_TF.min():.1f}-{100 * b.est_TF.max():.1f}%"


def floor_bin(aurocs, bins, threshold):
    """Lowest bin (ascending est_TF) whose AUROC reaches the threshold.

    This is the permissive reading. AUROC is not monotonic in est_TF here --
    each bin holds 12 virtual samples from only 3-5 patients, so bin
    composition moves the number as much as tumor fraction does -- which
    means a bin can cross the threshold and later bins fall back below it.
    """
    for a, b in zip(aurocs, bins):
        if not np.isnan(a) and a >= threshold:
            return _label(b)
    return "not reached"


def floor_bin_sustained(aurocs, bins, threshold):
    """Lowest bin from which EVERY higher bin also reaches the threshold.

    The conservative reading. Where this disagrees with floor_bin, the
    difference is a direct measure of how fragile the floor is.
    """
    ok = [(not np.isnan(a)) and a >= threshold for a in aurocs]
    floor = None
    for i in range(len(bins) - 1, -1, -1):
        if ok[i]:
            floor = i
        else:
            break
    return _label(bins[floor]) if floor is not None else "not reached"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--estimates", required=True,
                    help="all_estimates.tsv written by collect_all.py")
    ap.add_argument("--splits", required=True, help="splits.tsv")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--pattern", default="a_rescore_*",
                    help="which conditions to use (fnmatch glob)")
    ap.add_argument("--role", default="score_heldout",
                    help="which splits.tsv role supplies the negatives")
    ap.add_argument("--threshold", type=float, default=0.95,
                    help="AUROC counted as detection, for the floor column")
    ap.add_argument("--method", default=None,
                    help="restrict to one method, e.g. celfieish-clip")
    ap.add_argument("--markers", default=None,
                    help="restrict to one marker set, e.g. max")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    est = pd.read_csv(args.estimates, sep="\t")
    conditions = load_conditions(args.splits, args.pattern, args.role)

    sizes = {len(v) for v in conditions.values()}
    print(f"{len(conditions)} conditions matching {args.pattern!r}, "
          f"negatives per condition: {sorted(sizes)}")

    if args.method:
        est = est[est.method == args.method]
    if args.markers:
        est = est[est.markers == args.markers]

    rows = []
    for (method, markers), g in est.groupby(["method", "markers"]):
        ctl_all = g[(g.group == "control") & (g.dilution == 1.0)]
        crc = g[g.group == "crc"]
        if len(ctl_all) != 22 or len(crc) != 120:
            print(f"  skipping {method} {markers}: "
                  f"{len(ctl_all)} controls, {len(crc)} CRC rows "
                  f"(need 22 and 120)")
            continue

        bins = bin_table(crc)
        full = [auc(b.epi.values, ctl_all.epi.values) for b in bins]

        print(f"\n=== {method} | {markers}")
        print(f"{'bin (est_TF)':>16}  {'all 22':>7}  "
              f"{'subset median':>13}  {'2.5-97.5 pct':>16}  "
              f"{'delta':>7}  patients")

        per_cond = {}
        for cond, members in conditions.items():
            sub = ctl_all[ctl_all["sample"].isin(members)]
            if len(sub) != len(members):
                missing = set(members) - set(ctl_all["sample"])
                sys.exit(f"{cond}: controls missing from estimates: "
                         f"{sorted(missing)}")
            per_cond[cond] = [auc(b.epi.values, sub.epi.values) for b in bins]

        arr = np.array([per_cond[c] for c in sorted(per_cond)])  # conds x bins

        for i, b in enumerate(bins):
            col = arr[:, i]
            med = np.nanmedian(col)
            lo, hi = np.nanpercentile(col, [2.5, 97.5])
            label = f"{100 * b.est_TF.min():.1f}-{100 * b.est_TF.max():.1f}%"
            print(f"{label:>16}  {full[i]:>7.3f}  {med:>13.3f}  "
                  f"{lo:>6.3f}-{hi:<6.3f}  {med - full[i]:>+7.3f}  "
                  f"{b.patient.nunique()}p")

            rows.append(dict(
                method=method, markers=markers, bin_index=i, bin_label=label,
                est_TF_min=b.est_TF.min(), est_TF_max=b.est_TF.max(),
                n_patients=b.patient.nunique(),
                auroc_all_controls=full[i], n_controls_all=len(ctl_all),
                auroc_subset_median=med, auroc_subset_p2_5=lo,
                auroc_subset_p97_5=hi, delta_median_minus_all=med - full[i],
                n_controls_subset=len(next(iter(conditions.values()))),
                n_conditions=len(conditions)))

        med_by_bin = np.nanmedian(arr, axis=0)
        f_first_all = floor_bin(full, bins, args.threshold)
        f_sust_all = floor_bin_sustained(full, bins, args.threshold)
        print(f"  floor at AUROC >= {args.threshold}")
        print(f"    first bin to cross:        all 22 -> {f_first_all};  "
              f"subset median -> {floor_bin(med_by_bin, bins, args.threshold)}")
        print(f"    sustained from:            all 22 -> {f_sust_all};  "
              f"subset median -> "
              f"{floor_bin_sustained(med_by_bin, bins, args.threshold)}")
        if f_first_all != f_sust_all:
            print(f"    ** the two definitions disagree: AUROC dips back "
                  f"below {args.threshold} above the first crossing, so this "
                  f"floor is fragile")

        # Per-condition floors, so the spread in the floor itself is visible.
        for name, fn in (("first-crossing", floor_bin),
                         ("sustained", floor_bin_sustained)):
            tally = pd.Series([fn(per_cond[c], bins, args.threshold)
                               for c in sorted(per_cond)]).value_counts()
            print(f"    per-condition {name} floors: {tally.to_dict()}")

    if not rows:
        sys.exit("No (method, markers) group had the expected 22 controls "
                 "and 120 CRC rows. Check --estimates.")

    out_tsv = os.path.join(args.out, "rescore_auroc.tsv")
    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)
    print(f"\nWrote {out_tsv}")
    print("\nThe spread across conditions is the sample-size artifact: it is "
          "what AUROC does when the negative set halves, with no leakage "
          "removed. Setup (b) has to move the floor by more than this to mean "
          "anything.")


if __name__ == "__main__":
    main()
