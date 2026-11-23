import collections
import contextlib
import io
import os
import re
import subprocess
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "common"))
sys.path.insert(0, os.path.join(HERE, "..", "step1"))
import score_pair
import step1_lookalikes as s1

N_BOOT = 2000
COLON = "results/checks/colon_ep"
PAIRS = {
    "A1_all_strings": ("results/step2/{id}/whole", "results/step2/{id}/cut"),
    "A2_2plus_strings": ("results/step2_2plus/{id}/whole", "results/step2_2plus/{id}/cut"),
    "A3_ladder_whole_vs_single": ("results/step4/{id}/whole", "results/step4/{id}/size1"),
    "A3_ladder_3bead_vs_single": ("results/step4/{id}/size3", "results/step4/{id}/size1"),
    "A3_ladder_2bead_vs_single": ("results/step4/{id}/size2", "results/step4/{id}/size1"),
    "A4_bucket2": ("results/step3/{id}/b2_whole", "results/step3/{id}/b2_cut"),
    "A4_bucket3": ("results/step3/{id}/b3_whole", "results/step3/{id}/b3_cut"),
    "D1_simulated": ("results/d1/{id}/whole", "results/d1/{id}/cut"),
}
SCORES = {"epithelial": "results/scores", "Colon-Ep": COLON}
CLIPS = ("unclipped", "clipped")
MAIN = ("A1_all_strings", "A3_ladder_whole_vs_single")


def rescore_colon():
    score_pair.EPI["u250"] = {"Colon-Ep"}
    for label, (a, b) in PAIRS.items():
        sys.argv = ["score_pair.py", "--meta", "meta/samples.tsv", "--markers", "u250", "--a", a, "--b", b,
                    "--label", label, "--out", COLON]
        with contextlib.redirect_stdout(io.StringIO()):
            score_pair.main()


def load(score, label, clip):
    rows = []
    for line in list(open(f"{SCORES[score]}/{label}_{clip}_estimates.tsv"))[1:]:
        i, patient, batch, cancer, dil, exp, a, b = line.rstrip("\n").split("\t")
        rows.append(dict(id=i, sample=re.sub(r"__d[\d.]+$", "", i), patient=patient, batch=batch,
                         cancer=cancer == "1", dil=float(dil), exp=float(exp), a=float(a), b=float(b)))
    return rows


def auroc(pos, neg):
    return float(score_pair.auroc(np.array(pos), np.array(neg)))


def e1():
    print("\n\nE1: batch alone (healthy HU012 against healthy HU083) next to cancer against controls, undiluted")
    print("score\tcomparison\tclip\thealthy_vs_healthy_A\thealthy_vs_healthy_B\tgain\tci_lo\tci_hi\t"
          "cancer_vs_controls_A\tcancer_vs_controls_B\tgain")
    rng = np.random.default_rng(0)
    for score in SCORES:
        for label in PAIRS:
            for clip in CLIPS:
                rows = load(score, label, clip)
                h = {b: np.array([(r["a"], r["b"]) for r in rows if r["batch"] == b]) for b in ("HU012", "HU083")}
                und = [r for r in rows if r["cancer"] and r["dil"] == 1]
                ctl = [r for r in rows if not r["cancer"]]
                ha, hb = (auroc(h["HU012"][:, k], h["HU083"][:, k]) for k in (0, 1))
                boot = []
                for _ in range(N_BOOT):
                    x, y = (h[b][rng.integers(0, len(h[b]), len(h[b]))] for b in ("HU012", "HU083"))
                    boot.append(auroc(x[:, 0], y[:, 0]) - auroc(x[:, 1], y[:, 1]))
                lo, hi = np.percentile(boot, [2.5, 97.5])
                ca, cb = (auroc([r[k] for r in und], [r[k] for r in ctl]) for k in "ab")
                print(f"{score}\t{label}\t{clip}\t{ha:.3f}\t{hb:.3f}\t{ha - hb:+.3f}\t{lo:+.3f}\t{hi:+.3f}\t"
                      f"{ca:.3f}\t{cb:.3f}\t{ca - cb:+.3f}")


def e2():
    print("\n\nE2: the two cancer samples with tumor share under 2%")
    print("score\tcomparison\tclip\tsample\ttumor_share\test_A\tcontrols_below_A\test_B\tcontrols_below_B")
    out = ["score\tcomparison\tclip\tcontrols_used\tauroc_A_all20\tauroc_B_all20\tgain\t"
           "auroc_A_without\tauroc_B_without\tgain"]
    for score in SCORES:
        for label in MAIN:
            for clip in CLIPS:
                rows = load(score, label, clip)
                und = [r for r in rows if r["cancer"] and r["dil"] == 1]
                ctl = [r for r in rows if not r["cancer"]]
                for r in und:
                    if r["exp"] < 0.02:
                        below = [sum(c[k] < r[k] for c in ctl) for k in "ab"]
                        print(f"{score}\t{label}\t{clip}\t{r['id']}\t{r['exp']:.4f}\t{r['a']:.4f}\t{below[0]} of 22\t"
                              f"{r['b']:.4f}\t{below[1]} of 22")
                for batch in ("all", "HU012", "HU083"):
                    c = [r for r in ctl if batch in ("all", r["batch"])]
                    v = []
                    for keep in (und, [r for r in und if r["exp"] >= 0.02]):
                        v += [auroc([r[k] for r in keep], [r[k] for r in c]) for k in "ab"]
                    out.append(f"{score}\t{label}\t{clip}\t{batch}\t{v[0]:.3f}\t{v[1]:.3f}\t{v[0] - v[1]:+.3f}\t"
                               f"{v[2]:.3f}\t{v[3]:.3f}\t{v[2] - v[3]:+.3f}")
    print("\nundiluted AUROC with all 20 samples and without those two")
    print("\n".join(out))


def e3():
    print("\n\nE3: one line per patient through the undiluted samples (estimate = intercept + slope x tumor share)")
    print("score\tcomparison\tclip\tpatient\tbatch\tsamples\ttumor_share_range\tslope_A\tintercept_A\tslope_B\tintercept_B")
    for score in SCORES:
        for label in MAIN:
            for clip in CLIPS:
                rows = load(score, label, clip)
                und = [r for r in rows if r["cancer"] and r["dil"] == 1]
                groups = collections.defaultdict(list)
                for r in und:
                    groups[(r["patient"], r["batch"])].append(r)
                groups[("all 20", "-")] = und
                for (patient, batch), v in sorted(groups.items()):
                    x = np.array([r["exp"] for r in v])
                    rng_ = f"{x.min():.3f}-{x.max():.3f}"
                    if len(v) < 2:
                        print(f"{score}\t{label}\t{clip}\t{patient}\t{batch}\t1\t{rng_}\tone sample, no line\t\t\t")
                        continue
                    (sa, ia), (sb, ib) = (np.polyfit(x, [r[k] for r in v], 1) for k in "ab")
                    print(f"{score}\t{label}\t{clip}\t{patient}\t{batch}\t{len(v)}\t{rng_}\t{sa:.2f}\t{ia:.3f}\t"
                          f"{sb:.2f}\t{ib:.3f}")
                for batch in ("all", "HU012", "HU083"):
                    c = [r for r in rows if not r["cancer"] and batch in ("all", r["batch"])]
                    print(f"{score}\t{label}\t{clip}\thealthy median\t{batch}\t{len(c)}\t\t\t"
                          f"{np.median([r['a'] for r in c]):.3f}\t\t{np.median([r['b'] for r in c]):.3f}")


def e4():
    print("\n\nE4: tracking of each sample's dilution series (6 levels)")
    print("A series whose estimate is the same at every level has no correlation; it is counted in flat_series and left out of the mean and median.")
    print("score\tcomparison\tclip\tversion\tseries_used\tmean_corr\tmedian_corr\tseries_with_corr_below_0.5\tflat_series\tmedian_slope")
    for score in SCORES:
        for label in ("A1_all_strings", "A3_ladder_whole_vs_single", "A3_ladder_3bead_vs_single",
                      "A3_ladder_2bead_vs_single"):
            for clip in CLIPS:
                series = collections.defaultdict(list)
                for r in load(score, label, clip):
                    if r["cancer"]:
                        series[r["sample"]].append(r)
                for k in ("a", "b") if label in MAIN[:1] + ("A3_ladder_whole_vs_single",) else ("a",):
                    for name, keep in (("all 20", 0.0), ("tumor share 5% up", 0.05)):
                        cor, slo = [], []
                        for v in series.values():
                            x = np.array([r["exp"] for r in v])
                            if x.max() < keep:
                                continue
                            y = np.array([r[k] for r in v])
                            cor.append(np.corrcoef(x, y)[0, 1])
                            slo.append(np.polyfit(x, y, 1)[0])
                        print(f"{score}\t{label}\t{clip}\t{'A' if k == 'a' else 'B'}\t{name} ({len(cor)})\t"
                              f"{np.nanmean(cor):.3f}\t{np.nanmedian(cor):.3f}\t{sum(c < 0.5 for c in cor)}\t"
                              f"{int(np.isnan(cor).sum())}\t{np.median(slo):.2f}")


def e5():
    print("\n\nE5: look-alike recheck, real healthy controls against simulated ones")
    os.makedirs("results/checks/step1_sim", exist_ok=True)
    starts, ends = s1.load_places("data/places/epithelial_u250.bed")
    meta = [line.rstrip("\n").split("\t") for line in list(open("meta/samples.tsv"))[1:]]
    print("batch\tdata\tbeads_in_strings\tstrings\tshare_all_off\tseen_over_chance\tbead_off_rate_in_these_strings")
    summary = []
    for batch in ("HU012", "HU083"):
        subprocess.run([sys.executable, "-I", os.path.join(HERE, "..", "step1", "step1_lookalikes.py"),
                        "--meta", "meta/samples.tsv", "--marker-dir", "results/d1/sim",
                        "--places", "data/places/epithelial_u250.bed", "--label", f"sim_u250_{batch}",
                        "--out", "results/checks/step1_sim", "--batch", batch], check=True, capture_output=True)
        tables = {"real": f"results/step1/u250_{batch}_pooled.tsv",
                  "simulated": f"results/checks/step1_sim/sim_u250_{batch}_pooled.tsv"}
        folders = {"real": "data/marker_pat", "simulated": "results/d1/sim"}
        rate = {}
        for data in ("real", "simulated"):
            n, off = collections.Counter(), collections.Counter()
            for f in meta:
                if f[1] == "control" and f[3] == batch:
                    for _, idx, o in s1.strings(f"{folders[data]}/{f[0]}.pat", starts, ends):
                        n[s1.bucket(len(idx))] += len(idx)
                        off[s1.bucket(len(idx))] += sum(o)
            for line in list(open(tables[data]))[1:]:
                t = line.split("\t")
                if t[0] == "healthy":
                    rate[(data, t[1])] = (off[t[1]] / n[t[1]], float(t[4]))
                    print(f"{batch}\t{data}\t{t[1]}\t{t[2]}\t{t[4]}\t{t[6]}\t{off[t[1]] / n[t[1]]:.4f}")
            rate[(data, "all")] = sum(off.values()) / sum(n.values())
        r, s = rate[("real", "4+")], rate[("simulated", "4+")]
        summary.append(f"{batch}\t{rate[('real', 'all')]:.4f}\t{rate[('simulated', 'all')]:.4f}\t{r[0]:.4f}\t{s[0]:.4f}\t"
                       f"{r[0] / s[0]:.2f}\t{(r[0] / s[0]) ** 4:.1f}\t{r[1] / s[1]:.1f}")
    print("\nbatch\toff_rate_all_beads_real\tsimulated\toff_rate_beads_on_4+_strings_real\tsimulated\treal_over_simulated\t"
          "fold_more_all_off_4_bead_strings_that_rate_alone_predicts\tfold_seen_(share_all_off_4+_real_over_simulated)")
    print("\n".join(summary))


def e6():
    print("\n\nE6: healthy baseline against the number of strings, 22 controls (A1, all strings)")
    est = {(score, clip): {r["id"]: r for r in load(score, lab, clip)}
           for score in SCORES for clip in CLIPS for lab in ("A1_all_strings",)}
    sim = {r["id"]: r for r in load("epithelial", "D1_simulated", "unclipped")}
    ctl = [r for r in est[("epithelial", "unclipped")].values() if not r["cancer"]]
    n = {r["id"]: int(open(f"results/step2/{r['id']}/whole.nlines").read()) for r in ctl}
    cols = [("epithelial unclipped A", lambda i: est[("epithelial", "unclipped")][i]["a"]),
            ("epithelial unclipped B", lambda i: est[("epithelial", "unclipped")][i]["b"]),
            ("epithelial clipped A", lambda i: est[("epithelial", "clipped")][i]["a"]),
            ("epithelial clipped B", lambda i: est[("epithelial", "clipped")][i]["b"]),
            ("Colon-Ep unclipped A", lambda i: est[("Colon-Ep", "unclipped")][i]["a"]),
            ("Colon-Ep unclipped B", lambda i: est[("Colon-Ep", "unclipped")][i]["b"]),
            ("simulated unclipped A", lambda i: sim[i]["a"]),
            ("simulated unclipped B", lambda i: sim[i]["b"])]
    print("control\tbatch\tstrings\t" + "\t".join(c[0] for c in cols))
    for r in sorted(ctl, key=lambda r: (r["batch"], n[r["id"]])):
        print(f"{r['id']}\t{r['batch']}\t{n[r['id']]}\t" + "\t".join(f"{f(r['id']):.4f}" for _, f in cols))
    print("\nrank correlation of the estimate with the number of strings")
    print("controls\t" + "\t".join(c[0] for c in cols))
    for batch in ("all", "HU012", "HU083"):
        ids = [r["id"] for r in ctl if batch in ("all", r["batch"])]
        print(f"{batch} ({len(ids)})\t" + "\t".join(
            f"{spearmanr([n[i] for i in ids], [f(i) for i in ids])[0]:+.2f}" for _, f in cols))
    dil = [int(open(f"results/step2/{d}/whole.nlines").read()) for d in os.listdir("results/step2") if "__d" in d]
    und = [int(open(f"results/step2/{r['id']}/whole.nlines").read())
           for r in est[("epithelial", "unclipped")].values() if r["cancer"] and r["dil"] == 1]
    for name, v in (("controls", list(n.values())), ("undiluted cancer samples", und), ("dilutions", dil)):
        print(f"strings per file, {name}: smallest {min(v)}, median {int(np.median(v))}, largest {max(v)}")


def main():
    print("EXPLORATORY CHECKS (chosen after seeing the planned results)")
    rescore_colon()
    for f in (e1, e2, e3, e4, e5, e6):
        f()


if __name__ == "__main__":
    main()
