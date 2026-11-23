import collections
import re

import numpy as np
from scipy.stats import rankdata

S = "results/scores"
N_BOOT = 2000
BATCHES = ("all", "HU012", "HU083")
CHECK1 = ["A1_all_strings", "A2_2plus_strings", "A3_ladder_whole_vs_single",
          "A3_ladder_3bead_vs_single", "A3_ladder_2bead_vs_single",
          "B3_background_HU012", "B3_background_HU083", "D1_simulated"]
LADDER = [("whole", "A3_ladder_whole_vs_single"), ("3-bead pieces", "A3_ladder_3bead_vs_single"),
          ("2-bead pieces", "A3_ladder_2bead_vs_single")]
CHECK4 = CHECK1[:5] + ["A4_bucket2", "A4_bucket3", "A5_matched_bucket2", "A5_matched_bucket3",
                       "A5_matched_bucket4+"]


def corr(x, y):
    return np.corrcoef(x, y)[0, 1]


def check1(label, clip):
    series = collections.defaultdict(list)
    ctrl = collections.defaultdict(list)
    for line in list(open(f"{S}/{label}_{clip}_estimates.tsv"))[1:]:
        i, patient, batch, cancer, _, exp, a, b = line.rstrip("\n").split("\t")
        if cancer == "1":
            series[(patient, re.sub(r"__d[\d.]+$", "", i))].append((float(exp), float(a), float(b)))
        else:
            ctrl[batch].append((float(a), float(b)))
    by_pat = collections.defaultdict(list)
    for (patient, _), files in series.items():
        by_pat[patient].append(np.array(files))
    pats = list(by_pat.values())
    ctrl = {k: np.array(v) for k, v in ctrl.items()}

    def score(pick_p, pick_h):
        ser = [s for i in pick_p for s in pats[i]]
        exp, a, b = np.concatenate(ser).T
        h = {k: ctrl[k][pick_h[k]] for k in ctrl}
        h["all"] = np.concatenate(list(h.values()))
        v = [corr(rankdata(a), rankdata(exp)), corr(rankdata(b), rankdata(exp)),
             np.mean([corr(s[:, 1], s[:, 0]) for s in ser]), np.mean([corr(s[:, 2], s[:, 0]) for s in ser])]
        for k in BATCHES:
            base_a, base_b = np.median(h[k], axis=0)
            v += [-np.abs(a - base_a - exp).mean(), -np.abs(b - base_b - exp).mean()]
        return np.array(v)

    rng = np.random.default_rng(0)
    point = score(range(len(pats)), {k: slice(None) for k in ctrl})
    boot = np.array([score(rng.integers(0, len(pats), len(pats)),
                           {k: rng.integers(0, len(ctrl[k]), len(ctrl[k])) for k in ctrl})
                     for _ in range(N_BOOT)])
    names = ["rank correlation", "within-series corr"] + [f"error minus {k} baseline" for k in BATCHES]
    print(f"\n{label} | {clip} | cancer files {sum(len(s) for p in pats for s in p)}, patients {len(pats)}")
    print("measure\tA\tB\tgain\tci_lo\tci_hi\tmoved")
    for j, name in enumerate(names):
        a, b = point[2 * j], point[2 * j + 1]
        lo, hi = np.percentile(boot[:, 2 * j] - boot[:, 2 * j + 1], [2.5, 97.5])
        if name.startswith("error"):
            a, b = -a, -b
        moved = "yes" if lo > 0 or hi < 0 else "no"
        print(f"{name}\t{a:.4f}\t{b:.4f}\t{point[2 * j] - point[2 * j + 1]:+.4f}\t{lo:+.4f}\t{hi:+.4f}\t{moved}")


def summary(label):
    out, clip, head, part = {}, None, None, 0
    for line in open(f"{S}/{label}_summary.tsv"):
        f = line.rstrip("\n").split("\t")
        if line.startswith("==="):
            clip = "unclipped" if "ISH unclipped" in line else "clipped"
            out[clip], head, part = ([], []), None, 0
        elif line.startswith("healthy controls"):
            head, part = None, 1
        elif clip and len(f) > 3:
            if head is None:
                head = f
            else:
                out[clip][part].append(dict(zip(head, f)))
    return out


def main():
    print("CHECK 1: how much, without the false baseline")
    for label in CHECK1:
        for clip in ("unclipped", "clipped"):
            check1(label, clip)

    print("\n\nCHECK 2: detection by tumor level, background and controls from the same batch (B3)")
    print("batch\tclip\tsubset\tcancer_files\tauroc_A\tauroc_B\tauroc_gain\tci_lo\tci_hi\tmoved")
    for batch in ("HU012", "HU083"):
        tab = summary(f"B3_background_{batch}")
        for clip in ("unclipped", "clipped"):
            for r in tab[clip][0]:
                if r["controls_used"] == batch and r["subset"].startswith("tumor_"):
                    print(f"{batch}\t{clip}\t{r['subset']}\t{r['cancer_files']}\t{r['auroc_A']}\t{r['auroc_B']}\t"
                          f"{r['auroc_gain']}\t{r['auc_ci_lo']}\t{r['auc_ci_hi']}\t{r['auroc_moved']}")

    print("\n\nCHECK 3: the size ladder on the 20 undiluted samples (A3, 4+ bead strings)")
    print("clip\tcontrols_used\tversion\tauroc\tgain_vs_single\tci_lo\tci_hi\tmoved\terror\terror_gain_vs_single\terror_moved")
    for clip in ("unclipped", "clipped"):
        for batch in BATCHES:
            for name, label in LADDER:
                r = next(r for r in summary(label)[clip][0]
                         if r["subset"] == "dilution_1" and r["controls_used"] == batch)
                print(f"{clip}\t{batch}\t{name}\t{r['auroc_A']}\t{r['auroc_gain']}\t{r['auc_ci_lo']}\t"
                      f"{r['auc_ci_hi']}\t{r['auroc_moved']}\t{r['error_A']}\t{r['error_gain']}\t{r['error_moved']}")
            print(f"{clip}\t{batch}\tsingle beads\t{r['auroc_B']}\t\t\t\t\t{r['error_B']}\t\t")

    print("\n\nCHECK 4: healthy false signal, A against B, unclipped and clipped")
    print("comparison\tclip\tcontrols\tmedian_A\tmedian_B\tmean_A_minus_B\tci_lo\tci_hi")
    for label in CHECK4:
        tab = summary(label)
        for clip in ("unclipped", "clipped"):
            for r in tab[clip][1]:
                print(f"{label}\t{clip}\t{r['group']}\t{r['median_A']}\t{r['median_B']}\t"
                      f"{r['mean_A_minus_B']}\t{r['ci_lo']}\t{r['ci_hi']}")


if __name__ == "__main__":
    main()
