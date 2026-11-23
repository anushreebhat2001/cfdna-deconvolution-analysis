#!/usr/bin/env python3
"""
Score setup (b)/(c) conditions: the same mixtures, judged two ways.

For each condition, the mixtures built from that condition's control pool are
scored against:
    score_leaky    the controls that went into the pool   (full overlap)
    score_heldout  the disjoint controls                  (no overlap)

Because the mixtures, the pool size, the scoring-set size and the batch
balance are identical between the two arms, the paired difference isolates
leakage. Setup (a)'s numbers (all 22 controls, paper pool) are printed
alongside as a reference.

Only the mixture estimates are condition-specific. The undiluted CRC and
control estimates do not depend on the pool, so they are reused from the
existing all_estimates.tsv rather than recomputed.

Usage
-----
  python3 score_b.py \
      --base-estimates /scratch/mxf4642/results/summary/all_estimates.tsv \
      --est-dir   /scratch/$USER/mrd/results/b/est \
      --mixinfo   /scratch/$USER/mrd/results/b/mixinfo \
      --splits    /scratch/$USER/mrd/data/splits.tsv \
      --samples   /scratch/mxf4642/code/project/meta/samples.tsv \
      --atlas-dir /scratch/mxf4642/data/katsman_run \
      --out       /scratch/$USER/mrd/results/b/score

Run inside the container.
"""

import argparse
import fnmatch
import glob
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

N_BINS = 10

# Epithelial cell types, as grouped in collect_all.py.
EPI = {"Bladder-Ep", "Breast-Basal-Ep", "Breast-Luminal-Ep", "Colon-Ep",
       "Fallopian-Ep", "Gastric-Ep", "Head-Neck-Ep", "Kidney-Ep",
       "Lung-Ep-Alveo", "Lung-Ep-Bron", "Ovary-Ep", "Prostate-Ep",
       "Small-Int-Ep", "Thyroid-Ep", "Pancreas-Acinar", "Pancreas-Duct",
       "Epithelial"}

ATLAS = {"u250": "31CT_U250_Netanel_atlas_over_regions.txt",
         "max": "blood_epithelial_atlas_over_regions_max_markers.txt"}


def atlas_names(atlas_dir):
    out = {}
    for key, fn in ATLAS.items():
        path = os.path.join(atlas_dir, fn)
        if not os.path.exists(path):
            continue
        header = open(path).readline().rstrip("\n").split("\t")
        out[key] = [h[:-5] for h in header if h.endswith("_METH")]
    return out


def epi_of(path, names):
    v = np.loadtxt(path).ravel()
    total = v.sum()
    if total <= 0:
        return np.nan
    v = v / total
    return sum(x for n, x in zip(names, v) if n in EPI)


def auc(pos, neg):
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    y = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    return roc_auc_score(y, np.r_[pos, neg])


def _label(b):
    return f"{100 * b.est_TF.min():.1f}-{100 * b.est_TF.max():.1f}%"


def floor_first(aurocs, bins, thr):
    for a, b in zip(aurocs, bins):
        if not np.isnan(a) and a >= thr:
            return _label(b)
    return "not reached"


def floor_sustained(aurocs, bins, thr):
    ok = [(not np.isnan(a)) and a >= thr for a in aurocs]
    floor = None
    for i in range(len(bins) - 1, -1, -1):
        if ok[i]:
            floor = i
        else:
            break
    return _label(bins[floor]) if floor is not None else "not reached"


def load_mixinfo(mixinfo_dir):
    """(condition, mix_id) -> (sample, dilution)."""
    out = {}
    files = sorted(glob.glob(os.path.join(mixinfo_dir, "*.mixinfo.tsv")))
    if not files:
        sys.exit(f"No .mixinfo.tsv files in {mixinfo_dir}")
    for f in files:
        stem = os.path.basename(f)[:-len(".mixinfo.tsv")]
        cond = stem.split("__", 1)[0]
        df = pd.read_csv(f, sep="\t")
        for _, r in df.iterrows():
            out[(cond, r["mix_id"])] = (r["sample"], float(r["dilution"]))
    return out


def load_condition_estimates(est_dir, mixinfo, names_by_ms):
    """Parse this run's estimate files into rows."""
    rows = []
    for path in sorted(glob.glob(os.path.join(est_dir, "*.txt"))):
        base = os.path.basename(path)[:-4]
        try:
            method, ms, ident = base.split("@")
        except ValueError:
            continue
        if ms not in names_by_ms:
            continue
        # ident is "<condition>__<mix_id>"; conditions contain no "__".
        cond, mix_id = ident.split("__", 1)
        key = (cond, mix_id)
        if key not in mixinfo:
            continue
        sample, dilution = mixinfo[key]
        rows.append(dict(condition=cond, method=method, markers=ms,
                         mix_id=mix_id, sample=sample, dilution=dilution,
                         epi=epi_of(path, names_by_ms[ms])))
    if not rows:
        sys.exit(f"No usable estimate files parsed from {est_dir}")
    return pd.DataFrame(rows)


def roles_from_splits(splits_path):
    sp = pd.read_csv(splits_path, sep="\t")
    out = defaultdict(dict)
    for (cond, role), g in sp.groupby(["condition", "role"]):
        out[cond][role] = sorted(g["sample"].tolist())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-estimates", required=True)
    ap.add_argument("--est-dir", required=True)
    ap.add_argument("--mixinfo", required=True)
    ap.add_argument("--splits", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--atlas-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--conditions", default="*",
                    help="fnmatch pattern over condition names, e.g. "
                         "'b_split_*' or 'c_batch_*'. Setups (b) and (c) share "
                         "one estimate directory, so scoring them together "
                         "would pool unlike conditions into one median.")
    ap.add_argument("--threshold", type=float, default=0.95)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    names_by_ms = atlas_names(args.atlas_dir)
    meta = pd.read_csv(args.samples, sep="\t")
    base = pd.read_csv(args.base_estimates, sep="\t")
    roles = roles_from_splits(args.splits)

    mixinfo = load_mixinfo(args.mixinfo)
    mine = load_condition_estimates(args.est_dir, mixinfo, names_by_ms)
    all_conds = sorted(mine.condition.unique())
    conds = [c for c in all_conds if fnmatch.fnmatch(c, args.conditions)]
    if not conds:
        sys.exit(f"No conditions match {args.conditions!r}. "
                 f"Present: {', '.join(all_conds)}")
    skipped = [c for c in all_conds if c not in conds]
    print(f"{len(conds)} conditions with estimates: {', '.join(conds)}")
    if skipped:
        print(f"  ({len(skipped)} not matching {args.conditions!r} excluded)")
    mine = mine[mine.condition.isin(conds)]

    # Undiluted rows are pool-independent, so reuse them.
    undiluted = base[base.dilution == 1.0].copy()
    meta_cols = meta[["sample", "group", "patient", "batch", "ichorCNA_TF"]]

    summary_rows, percond_rows = [], []
    for (method, markers), g_mine in mine.groupby(["method", "markers"]):
        und = undiluted[(undiluted.method == method)
                        & (undiluted.markers == markers)]
        ctl_all = und[und.group == "control"]
        crc_und = und[und.group == "crc"]
        if len(ctl_all) != 22 or len(crc_und) != 20:
            print(f"  skipping {method} {markers}: baseline has "
                  f"{len(ctl_all)} controls, {len(crc_und)} undiluted CRC")
            continue

        print(f"\n=== {method} | {markers}")
        per_cond_leaky, per_cond_held, per_cond_ref = {}, {}, {}
        bins_ref = None

        for cond in conds:
            gm = g_mine[g_mine.condition == cond]
            if len(gm) != 100:
                print(f"  {cond}: {len(gm)} mixture estimates, expected 100 "
                      f"-- skipping")
                continue

            mix = gm.merge(meta_cols, on="sample", how="left")
            mix["est_TF"] = mix.ichorCNA_TF * mix.dilution
            crc = pd.concat([
                crc_und[["epi", "patient", "est_TF"]],
                mix[["epi", "patient", "est_TF"]],
            ], ignore_index=True)
            if len(crc) != 120:
                print(f"  {cond}: {len(crc)} CRC virtual samples, expected 120")
                continue

            bins = np.array_split(crc.sort_values("est_TF"), N_BINS)
            bins_ref = bins_ref or bins

            leaky = ctl_all[ctl_all["sample"].isin(
                roles[cond].get("score_leaky", []))]
            held = ctl_all[ctl_all["sample"].isin(
                roles[cond].get("score_heldout", []))]
            if len(held) == 0:
                print(f"  {cond}: no score_heldout controls in splits.tsv")
                continue

            per_cond_held[cond] = [auc(b.epi.values, held.epi.values)
                                   for b in bins]
            # Bin membership is set by est_TF (ichorCNA x dilution), which
            # does not depend on the pool, so bins are identical across
            # conditions and only the epi values differ. The all-22 arm must
            # still be recomputed per condition to be comparable.
            per_cond_ref[cond] = [auc(b.epi.values, ctl_all.epi.values)
                                  for b in bins]
            if len(leaky):
                per_cond_leaky[cond] = [auc(b.epi.values, leaky.epi.values)
                                        for b in bins]

        if not per_cond_held:
            continue

        ref = list(np.nanmedian(
            np.array([per_cond_ref[c] for c in sorted(per_cond_ref)]), axis=0))
        held_arr = np.array([per_cond_held[c] for c in sorted(per_cond_held)])
        has_leaky = bool(per_cond_leaky)
        leaky_arr = (np.array([per_cond_leaky[c] for c in sorted(per_cond_leaky)])
                     if has_leaky else None)

        # NOT setup (a): these are THIS condition's mixtures (built from an
        # 11-control pool) scored against all 22 controls, i.e. partial
        # overlap. The paper's own design is 20-control pool vs 22 scored --
        # get that from rescore_controls.py on the baseline estimates.
        hdr = f"{'bin (est_TF)':>16}  {'vs all 22':>9}"
        if has_leaky:
            hdr += f"  {'leaky':>7}"
        hdr += f"  {'held out':>8}"
        if has_leaky:
            hdr += f"  {'paired delta (2.5-97.5)':>26}"
        print(hdr + "  patients")

        for i, b in enumerate(bins_ref):
            held_med = np.nanmedian(held_arr[:, i])
            line = f"{_label(b):>16}  {ref[i]:>9.3f}"
            if has_leaky:
                line += f"  {np.nanmedian(leaky_arr[:, i]):>7.3f}"
            line += f"  {held_med:>8.3f}"
            if has_leaky:
                d = held_arr[:, i] - leaky_arr[:, i]
                lo, hi = np.nanpercentile(d, [2.5, 97.5])
                line += (f"  {np.nanmedian(d):>+8.3f} "
                         f"({lo:>+6.3f},{hi:>+6.3f})")
            print(line + f"  {b.patient.nunique()}p")

            row = dict(method=method, markers=markers, bin_index=i,
                       bin_label=_label(b), n_patients=b.patient.nunique(),
                       auroc_partial_overlap_all22=ref[i],
                       auroc_heldout_median=held_med,
                       n_conditions=held_arr.shape[0])
            if has_leaky:
                d = held_arr[:, i] - leaky_arr[:, i]
                row.update(auroc_leaky_median=np.nanmedian(leaky_arr[:, i]),
                           paired_delta_median=np.nanmedian(d),
                           paired_delta_p2_5=np.nanpercentile(d, 2.5),
                           paired_delta_p97_5=np.nanpercentile(d, 97.5))
            summary_rows.append(row)

        thr = args.threshold
        print(f"  floor at AUROC >= {thr}")
        print(f"    partial overlap (this pool, all 22 scored):  "
              f"first {floor_first(ref, bins_ref, thr)}"
              f" / sustained {floor_sustained(ref, bins_ref, thr)}")
        for name, arr in (("leaky", leaky_arr), ("held out", held_arr)):
            if arr is None:
                continue
            med = np.nanmedian(arr, axis=0)
            tally = pd.Series([floor_sustained(r, bins_ref, thr) for r in arr]
                              ).value_counts().to_dict()
            print(f"    {name:>9}, median:  "
                  f"first {floor_first(med, bins_ref, thr)}"
                  f" / sustained {floor_sustained(med, bins_ref, thr)}")
            print(f"    {name:>9}, per-condition sustained floors: {tally}")

        # ---- paired analysis -------------------------------------------
        # The two arms share mixtures within a condition, so compare them
        # condition by condition rather than comparing marginal medians.
        if has_leaky:
            order = sorted(per_cond_held)
            print("\n  PAIRED, within condition:")
            print(f"{'bin':>16}  {'worse':>5} {'better':>6} {'tied':>5}  "
                  f"{'median delta':>12}  {'Wilcoxon p':>10}")
            for i, b in enumerate(bins_ref):
                d = held_arr[:, i] - leaky_arr[:, i]
                d = d[~np.isnan(d)]
                worse = int((d < 0).sum())
                better = int((d > 0).sum())
                tied = int((d == 0).sum())
                p = np.nan
                if worse + better > 0:
                    try:
                        from scipy.stats import wilcoxon
                        p = wilcoxon(d, zero_method="zsplit").pvalue
                    except Exception:
                        pass
                pstr = "     n/a" if np.isnan(p) else f"{p:>10.4f}"
                print(f"{_label(b):>16}  {worse:>5} {better:>6} {tied:>5}  "
                      f"{np.nanmedian(d):>+12.3f}  {pstr}")

            # Floor shift, paired by condition.
            def fidx(aurocs):
                lab = floor_sustained(aurocs, bins_ref, thr)
                labels = [_label(x) for x in bins_ref]
                return labels.index(lab) if lab in labels else len(labels)

            shifts = [fidx(per_cond_held[c]) - fidx(per_cond_leaky[c])
                      for c in order if c in per_cond_leaky]
            shifts = np.array(shifts)
            print(f"\n  sustained-floor shift (held out minus leaky), "
                  f"in bins, over {len(shifts)} conditions:")
            print(f"    worse in {int((shifts > 0).sum())}, "
                  f"same in {int((shifts == 0).sum())}, "
                  f"better in {int((shifts < 0).sum())};  "
                  f"median shift {np.median(shifts):+.1f} bins")

            for ci, c in enumerate(order):
                if c not in per_cond_leaky:
                    continue
                for i, b in enumerate(bins_ref):
                    percond_rows.append(dict(
                        method=method, markers=markers, condition=c,
                        bin_index=i, bin_label=_label(b),
                        auroc_leaky=per_cond_leaky[c][i],
                        auroc_heldout=per_cond_held[c][i],
                        auroc_all22=per_cond_ref[c][i],
                        delta=per_cond_held[c][i] - per_cond_leaky[c][i]))

    if not summary_rows:
        sys.exit("Nothing scored. Check --est-dir and --base-estimates.")

    out_tsv = os.path.join(args.out, "setup_b_auroc.tsv")
    pd.DataFrame(summary_rows).to_csv(out_tsv, sep="\t", index=False)
    print(f"\nWrote {out_tsv}")
    if percond_rows:
        pc = os.path.join(args.out, "setup_b_per_condition.tsv")
        pd.DataFrame(percond_rows).to_csv(pc, sep="\t", index=False)
        print(f"Wrote {pc}")
    print("\nThe paired delta is held-out minus leaky on the SAME mixtures. "
          "Negative means removing the overlap makes detection worse, which "
          "is the leakage effect. Compare its size against the <=0.03 "
          "sample-size artifact measured by rescore_controls.py.")


if __name__ == "__main__":
    main()
