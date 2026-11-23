import collections
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "common"))
import score_pair
import stage1_checks as s1

O = "results/factors"
META, BREAST, LOW, auroc = s1.META, s1.BREAST, s1.LOW, s1.auroc
CTRL = [m for m in META if not m["cancer"]]
CANC = [m for m in META if m["cancer"]]
LEVELS = ("50", "25", "12.5", "6.25", "3.125")
XS = np.array([float(p) / 100 for p in LEVELS])


def vec(path):
    v = np.loadtxt(path).ravel()
    return v / v.sum()


def pair(folder, i, clip="unclipped"):
    return vec(f"{folder}/{i}/whole.ish_{clip}.txt"), vec(f"{folder}/{i}/cut.ish_{clip}.txt")


def scores(folder_of, mask, ids, clip="unclipped"):
    out = {}
    for i in ids:
        w, c = pair(folder_of(i), i, clip)
        out[i] = np.array([w[mask].sum(), c[mask].sum()])
    return out


def detection(sc, rng, breast=True):
    neg = {k: np.array([sc[m["id"]] for m in CTRL if k in ("all", m["batch"])]) for k in ("all", "HU012", "HU083")}
    rows = []
    for gname, keep in (("all20", lambda m: True), ("low6", lambda m: m["tf"] < LOW)):
        by_pat = collections.defaultdict(list)
        for m in CANC:
            if keep(m):
                by_pat[m["patient"]].append(sc[m["id"]])
        for k in neg:
            rows.append(f"{gname}\t{k}\t" + s1.fmt(s1.boot_gain([np.array(v) for v in by_pat.values()], neg[k], rng)))
    if breast:
        for k in neg:
            rows.append(f"breast4\t{k}\t" + s1.fmt(s1.boot_gain([sc[i][None, :] for i in BREAST], neg[k], rng)))
    rows.append("healthy_HU012\tHU083\t" + s1.fmt(s1.boot_gain([v[None, :] for v in neg["HU012"]], neg["HU083"], rng)))
    return rows


def medians(sc, groups):
    return "\t".join(f"{np.median([sc[i][0] for i in g]):.3f}\t{np.median([sc[i][1] for i in g]):.3f}" for g in groups.values())


GROUPS = {"healthy HU012": [m["id"] for m in CTRL if m["batch"] == "HU012"],
          "healthy HU083": [m["id"] for m in CTRL if m["batch"] == "HU083"],
          "breast4": BREAST,
          "cancer low6": [m["id"] for m in CANC if m["tf"] < LOW],
          "cancer 20%+": [m["id"] for m in CANC if m["tf"] >= LOW]}
HEAD = "positives\tnegatives\tauroc_whole\tauroc_single\tgain\tci_lo\tci_hi\tmoved"


def t6():
    print("T6: pan-epithelial markers (step 7; 8 cell types, one of them megakaryocytes). The breast samples were not fitted on this atlas.")
    names = [h[:-5] for h in open(score_pair.ATLAS["max"]).readline().rstrip("\n").split("\t") if h.endswith("_METH")]
    ids = [m["id"] for m in META]
    groups = {k: v for k, v in GROUPS.items() if k != "breast4"}
    print("median share per group\ntype\t" + "\t".join(f"{g} whole\t{g} single" for g in groups))
    for j, n in enumerate(names):
        mask = np.arange(len(names)) == j
        print(f"{n}\t" + medians(scores(lambda i: "results/step7", mask, ids), groups))
    epi = np.array([n == "Epithelial" for n in names])
    for clip in ("unclipped", "clipped"):
        print(f"\nepithelial score, {clip}, undiluted samples\n" + HEAD)
        print("\n".join(detection(scores(lambda i: "results/step7", epi, ids, clip), np.random.default_rng(0), breast=False)))
    print("\nfor comparison, U250 epithelial score (16 types), unclipped\n" + HEAD)
    print("\n".join(detection(scores(lambda i: "results/step2", s1.MASK["epithelial16"], ids), np.random.default_rng(0), breast=False)))
    sc = scores(lambda i: "results/step7", epi, ids)
    print("\none line per patient through the undiluted samples (epithelial score = intercept + slope x tumor share), unclipped")
    print("patient\tsamples\twhole slope\twhole intercept\tsingle slope\tsingle intercept")
    by = collections.defaultdict(list)
    for m in CANC:
        by[m["patient"]].append((m["tf"], *sc[m["id"]]))
    for p, v in by.items():
        v = np.array(v)
        if len(v) >= 2:
            (sw, iw), (ss, is_) = np.polyfit(v[:, 0], v[:, 1], 1), np.polyfit(v[:, 0], v[:, 2], 1)
            print(f"{p}\t{len(v)}\t{sw:.2f}\t{iw:.3f}\t{ss:.2f}\t{is_:.3f}")


def t7():
    print("\n\nT7: U250 atlas with one added column. 'off' = methylated 5% of the time at every CpG; 'on' = 95% (control arm).")
    ids = [m["id"] for m in META] + BREAST
    base = lambda i: "results/breast" if i in BREAST else "results/step2"
    arms = {"original atlas": (base, 31), "plus off column": (lambda i: f"{O}/t7_off", 32), "plus on column": (lambda i: f"{O}/t7_on", 32)}
    for arm, (folder, n) in arms.items():
        pad = lambda m: np.concatenate([m, np.zeros(n - 31, bool)])
        print(f"\n--- {arm}")
        print("median share per group, unclipped\nscore\t" + "\t".join(f"{g} whole\t{g} single" for g in GROUPS))
        for name in ("epithelial16", "gut3", "nongut13", "placebo7"):
            print(f"{name}\t" + medians(scores(folder, pad(s1.MASK[name]), ids), GROUPS))
        if n == 32:
            print("added column\t" + medians(scores(folder, np.arange(32) == 31, ids), GROUPS))
        for clip in ("unclipped", "clipped"):
            print(f"epithelial score (16 types), {clip}, undiluted samples\n" + HEAD)
            print("\n".join(detection(scores(folder, pad(s1.MASK["epithelial16"]), ids, clip), np.random.default_rng(0))))
    print("\nwhere the added column's share came from: change of the median share against the original atlas (points of 100), unclipped")
    print("arm\tgroup\tversion\tadded column\tepithelial16\tplacebo7\tall other types")
    other = ~(s1.MASK["epithelial16"] | s1.MASK["placebo7"])
    for arm in ("off", "on"):
        for g, sel in GROUPS.items():
            for k, version in enumerate(("whole", "single")):
                o = np.array([pair(base(i), i)[k] for i in sel])
                a = np.array([pair(f"{O}/t7_{arm}", i)[k] for i in sel])
                d = [np.median(a[:, :31][:, m].sum(1)) - np.median(o[:, m].sum(1)) for m in (s1.MASK["epithelial16"], s1.MASK["placebo7"], other)]
                print(f"{arm}\t{g}\t{version}\t{100 * np.median(a[:, 31]):.1f}\t" + "\t".join(f"{100 * x:+.1f}" for x in d))


def t8():
    print("\n\nT8: the 22 controls on the pan-epithelial atlas, with and without its megakaryocyte column (unclipped)")
    ids = [m["id"] for m in CTRL]
    full = [h[:-5] for h in open(score_pair.ATLAS["max"]).readline().rstrip("\n").split("\t") if h.endswith("_METH")]
    cut = [h[:-5] for h in open(f"{O}/atlas/max_minus_megak.txt").readline().rstrip("\n").split("\t") if h.endswith("_METH")]
    groups = {k: GROUPS[k] for k in ("healthy HU012", "healthy HU083")}
    print("type\tatlas\t" + "\t".join(f"{g} whole\t{g} single" for g in groups))
    for n in full:
        print(f"{n}\twith megakaryocytes\t" + medians(scores(lambda i: "results/step7", np.array([x == n for x in full]), ids), groups))
        if n in cut:
            print(f"{n}\twithout\t" + medians(scores(lambda i: f"{O}/t8_no_megak", np.array([x == n for x in cut]), ids), groups))
    for tag, folder, names in (("with megakaryocytes", "results/step7", full), ("without", f"{O}/t8_no_megak", cut)):
        sc = scores(lambda i: folder, np.array([x == "Epithelial" for x in names]), ids)
        h = {b: np.array([sc[i] for i in groups[f"healthy {b}"]]) for b in ("HU012", "HU083")}
        print(f"healthy HU012 against HU083 on the epithelial score, {tag}: whole {auroc(h['HU012'][:, 0], h['HU083'][:, 0]):.3f}, "
              f"single {auroc(h['HU012'][:, 1], h['HU083'][:, 1]):.3f}")


def t9():
    print("\n\nT9: the four breast samples diluted into pooled healthy strings (U250, epithelial score, unclipped)")
    print("A demonstration, not a discriminating test: these samples carry 0-3% tumor by copy number, so a falling series is\n"
          "what a dilution series looks like when the thing being diluted is the run, not tumor.")
    e = s1.MASK["epithelial16"]
    print("sample\tversion\tundiluted\t" + "\t".join(f"x={x:g}" for x in XS) + "\tslope over the 5 levels\tcorrelation with dilution (6 points)")
    series = [(i, "results/breast", f"{O}/t9_breast_dil") for i in BREAST]
    series += [(m["id"], "results/step2", "results/step2") for m in CANC if m["tf"] < 0.02]
    for i, und, dil in series:
        for k, version in enumerate(("whole", "single")):
            y = [pair(und, i)[k][e].sum()] + [pair(dil, f"{i}__d{p}")[k][e].sum() for p in LEVELS]
            slope = np.polyfit(XS, y[1:], 1)[0]
            r = np.corrcoef(np.concatenate([[1], XS]), y)[0, 1]
            tag = "" if i in BREAST else " (colorectal, tumor share under 2%, from step 2)"
            print(f"{i}{tag}\t{version}\t" + "\t".join(f"{v:.3f}" for v in y) + f"\t{slope:+.3f}\t{r:+.2f}")


def t10():
    print("\n\nT10: simulated healthy controls (true epithelial share 0) with one switch on. Epithelial score = false signal.")
    e = s1.MASK["epithelial16"]
    ids = [m["id"] for m in CTRL]
    rng = np.random.default_rng(0)
    conds = [("none (D1 as stored)", "results/d1")] + [(c, f"{O}/t10_{c}") for c in
             ("scattered_1", "clumped_1", "scattered_3", "clumped_3", "symmetric_5", "symmetric_10")]
    for clip in ("unclipped", "clipped"):
        print(f"\n{clip}\ncondition\tmedian whole\tmedian single\tmean whole minus single\tci_lo\tci_hi\tmoved\t"
              "HU012-shaped files: whole\tsingle\tHU083-shaped files: whole\tsingle")
        for name, folder in conds:
            sc = scores(lambda i: folder, e, ids, clip)
            a = np.array([sc[i] for i in ids])
            d = a[:, 0] - a[:, 1]
            bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(s1.N_BOOT)]
            lo, hi = np.percentile(bs, [2.5, 97.5])
            b = [np.median([sc[m["id"]][k] for m in CTRL if m["batch"] == bt]) for bt in ("HU012", "HU083") for k in (0, 1)]
            print(f"{name}\t{np.median(a[:, 0]):.3f}\t{np.median(a[:, 1]):.3f}\t{d.mean():+.3f}\t{lo:+.3f}\t{hi:+.3f}\t"
                  f"{'yes' if lo > 0 or hi < 0 else 'no'}\t" + "\t".join(f"{x:.3f}" for x in b))
    real = scores(lambda i: "results/step2", e, ids)
    a = np.array([real[i] for i in ids])
    print(f"\nfor comparison, the real controls (unclipped): median whole {np.median(a[:, 0]):.3f}, single {np.median(a[:, 1]):.3f}, "
          f"mean whole minus single {np.mean(a[:, 0] - a[:, 1]):+.3f}")


if __name__ == "__main__":
    print("EXPLORATORY factor tests: pan-epithelial readout and stage 2. Plan: step_scripts/factors/PLAN.md\n")
    for name in sys.argv[1:] or ["t6", "t7", "t8", "t9", "t10"]:
        globals()[name]()
