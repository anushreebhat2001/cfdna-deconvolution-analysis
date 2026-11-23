import collections
import gzip
import os
import sys

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import score_pair

N_BOOT = 2000
LOW = 0.20
NAMES = [h[:-5] for h in open(score_pair.ATLAS["u250"]).readline().rstrip("\n").split("\t") if h.endswith("_METH")]
EPI = score_pair.EPI["u250"]
GUT = {"Colon-Ep", "Small-Int-Ep", "Gastric-Ep"}
PLACEBO = {"Neuron", "Oligodend", "Heart-Cardio", "Adipocytes", "Pancreas-Alpha", "Pancreas-Beta", "Pancreas-Delta"}
SETS = {"epithelial16": EPI, "gut3": GUT, "nongut13": EPI - GUT, "placebo7": PLACEBO}
assert all(s <= set(NAMES) for s in SETS.values()) and len(EPI) == 16 and not (EPI & PLACEBO)
MASK = {k: np.array([n in s for n in NAMES]) for k, s in SETS.items()}
BREAST = ["HU017.04.PL4596", "HU017and19.05.PL4434", "HU017and19.06.PL4518", "HU025.13.PL4704"]
PAIRS = {
    "A1_all_strings": ("results/step2/{id}/whole", "results/step2/{id}/cut"),
    "A2_2plus_strings": ("results/step2_2plus/{id}/whole", "results/step2_2plus/{id}/cut"),
    "A3_4plus_whole_vs_single": ("results/step4/{id}/whole", "results/step4/{id}/size1"),
    "A4_bucket1_whole_only": ("results/step3/{id}/b1_whole", "results/step3/{id}/b1_whole"),
    "A4_bucket2": ("results/step3/{id}/b2_whole", "results/step3/{id}/b2_cut"),
    "A4_bucket3": ("results/step3/{id}/b3_whole", "results/step3/{id}/b3_cut"),
    "D1_simulated": ("results/d1/{id}/whole", "results/d1/{id}/cut"),
}
META = []
for _line in list(open("meta/samples.tsv"))[1:]:
    _f = _line.rstrip("\n").split("\t")
    META.append(dict(id=_f[0], cancer=_f[1] == "crc", patient=_f[2], batch=_f[3], tf=float(_f[5])))
AGE = {f[0]: float(f[1]) for f in (line.rstrip("\n").split("\t") for line in list(open("meta/sample_age.tsv"))[1:])}
auroc = score_pair.auroc


def vec(path):
    v = np.loadtxt(path).ravel()
    return v / v.sum()


def load(a, b, clip, ids):
    return {i: (vec(a.format(id=i) + f".ish_{clip}.txt"), vec(b.format(id=i) + f".ish_{clip}.txt")) for i in ids}


def boot_gain(pos_by_cluster, neg, rng):
    pos = np.concatenate(pos_by_cluster)
    a, b = auroc(pos[:, 0], neg[:, 0]), auroc(pos[:, 1], neg[:, 1])
    bs = []
    for _ in range(N_BOOT):
        p = np.concatenate([pos_by_cluster[i] for i in rng.integers(0, len(pos_by_cluster), len(pos_by_cluster))])
        n = neg[rng.integers(0, len(neg), len(neg))]
        bs.append(auroc(p[:, 0], n[:, 0]) - auroc(p[:, 1], n[:, 1]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return a, b, a - b, lo, hi


def fmt(r):
    a, b, g, lo, hi = r
    return f"{a:.3f}\t{b:.3f}\t{g:+.3f}\t{lo:+.3f}\t{hi:+.3f}\t{'yes' if lo > 0 or hi < 0 else 'no'}"


def t1():
    print("T1: the same fits scored with four sets of cell types")
    print("positives: all20 = the 20 undiluted cancer samples; low6 = those with tumor share under 20%;\n"
          "breast4 = the four breast samples from the cancer flow cells (all strings only);\n"
          "healthy_HU012 = the HU012 controls counted as positives against the HU083 controls.")
    print("comparison\tclip\tscore\tpositives\tnegatives\tn_pos\tn_neg\tauroc_whole\tauroc_single\tgain\tci_lo\tci_hi\tmoved")
    ids = [m["id"] for m in META]
    for label, (a, b) in PAIRS.items():
        for clip in ("unclipped", "clipped"):
            fits = load(a, b, clip, ids)
            breast = load("results/breast/{id}/whole", "results/breast/{id}/cut", clip, BREAST) if label == "A1_all_strings" else {}
            for name, mask in MASK.items():
                rng = np.random.default_rng(0)
                sc = {i: np.array([v[0][mask].sum(), v[1][mask].sum()]) for i, v in {**fits, **breast}.items()}
                ctrl = {k: np.array([sc[m["id"]] for m in META if not m["cancer"] and k in ("all", m["batch"])])
                        for k in ("all", "HU012", "HU083")}
                groups = {"all20": lambda m: True, "low6": lambda m: m["tf"] < LOW}
                for gname, keep in groups.items():
                    by_pat = collections.defaultdict(list)
                    for m in META:
                        if m["cancer"] and keep(m):
                            by_pat[m["patient"]].append(sc[m["id"]])
                    clusters = [np.array(v) for v in by_pat.values()]
                    for k, neg in ctrl.items():
                        print(f"{label}\t{clip}\t{name}\t{gname}\t{k}\t{sum(len(c) for c in clusters)}\t{len(neg)}\t"
                              + fmt(boot_gain(clusters, neg, rng)))
                if breast:
                    clusters = [sc[i][None, :] for i in BREAST]
                    for k, neg in ctrl.items():
                        print(f"{label}\t{clip}\t{name}\tbreast4\t{k}\t4\t{len(neg)}\t" + fmt(boot_gain(clusters, neg, rng)))
                clusters = [x[None, :] for x in ctrl["HU012"]]
                print(f"{label}\t{clip}\t{name}\thealthy_HU012\tHU083\t11\t11\t" + fmt(boot_gain(clusters, ctrl["HU083"], rng)))

    print("\nT1b: median score per group (all strings, unclipped), and rank correlation with tumor share over the 20 cancers")
    fits = load(*PAIRS["A1_all_strings"], "unclipped", ids)
    fits.update(load("results/breast/{id}/whole", "results/breast/{id}/cut", "unclipped", BREAST))
    grp = {"healthy HU012": [m["id"] for m in META if m["batch"] == "HU012"],
           "healthy HU083": [m["id"] for m in META if m["batch"] == "HU083"],
           "breast4": BREAST,
           "cancer low6": [m["id"] for m in META if m["cancer"] and m["tf"] < LOW],
           "cancer 20%+": [m["id"] for m in META if m["cancer"] and m["tf"] >= LOW]}
    print("score\tversion\t" + "\t".join(grp) + "\trank_corr_with_tumor_share")
    canc = [m for m in META if m["cancer"]]
    for name, mask in MASK.items():
        for k, version in enumerate(("whole", "single")):
            med = [np.median([fits[i][k][mask].sum() for i in g]) for g in grp.values()]
            rho = spearmanr([m["tf"] for m in canc], [fits[m["id"]][k][mask].sum() for m in canc])[0]
            print(f"{name}\t{version}\t" + "\t".join(f"{x:.3f}" for x in med) + f"\t{rho:+.2f}")
    print("\nT1c: where the share sits, per single type (all strings, unclipped, group medians in points of 100)")
    print("type\tset\t" + "\t".join(f"{g} whole\t{g} single" for g in grp))
    for j, n in enumerate(NAMES):
        s = next((k for k in ("gut3", "nongut13", "placebo7") if n in SETS[k]), "other")
        print(f"{n}\t{s}\t" + "\t".join(f"{100 * np.median([fits[i][0][j] for i in g]):.1f}\t"
                                        f"{100 * np.median([fits[i][1][j] for i in g]):.1f}" for g in grp.values()))


def t2():
    print("\n\nT2: offset null (all strings, unclipped, epithelial score)")
    ids = [m["id"] for m in META]
    e = MASK["epithelial16"]
    fits = load(*PAIRS["A1_all_strings"], "unclipped", ids)
    sc = {i: (v[0][e].sum(), v[1][e].sum()) for i, v in fits.items()}
    br = load("results/breast/{id}/whole", "results/breast/{id}/cut", "unclipped", BREAST)
    gaps = np.array([v[0][e].sum() - v[1][e].sum() for v in br.values()])
    print("whole-minus-single gap of the breast samples: " + ", ".join(f"{g:+.3f}" for g in gaps) + f"; mean {gaps.mean():+.3f}")
    canc = [m for m in META if m["cancer"]]
    obs_gap = np.array([sc[m["id"]][0] - sc[m["id"]][1] for m in canc])
    low = np.array([m["tf"] < LOW for m in canc])
    print(f"whole-minus-single gap of the cancer samples: low6 mean {obs_gap[low].mean():+.3f}, 20%+ mean {obs_gap[~low].mean():+.3f}")
    offsets = [("none (single beads as they are)", 0.0), ("breast mean", gaps.mean())]
    offsets += [(f"breast mean without sample {k + 1}", np.delete(gaps, k).mean()) for k in range(4)]
    print("positives\tcontrols\tobserved_whole\t" + "\t".join(f"single + {n} ({o:+.3f})" for n, o in offsets))
    for gname, sel in (("all20", np.ones(len(canc), bool)), ("low6", low)):
        for k in ("all", "HU012", "HU083"):
            hw = np.array([sc[m["id"]][0] for m in META if not m["cancer"] and k in ("all", m["batch"])])
            w = np.array([sc[m["id"]][0] for m in canc])[sel]
            c = np.array([sc[m["id"]][1] for m in canc])[sel]
            print(f"{gname}\t{k}\t{auroc(w, hw):.3f}\t" + "\t".join(f"{auroc(c + o, hw):.3f}" for _, o in offsets))


def t3():
    print("\n\nT3: dilution series in HU012 and in HU083 background (B3 fits, epithelial score, unclipped)")
    print("A dilution file is a mixture: a share x of the cancer sample's strings and 1-x of healthy background.\n"
          "If estimates mix in proportion, a series is a straight line from the background's floor (x=0) to the\n"
          "undiluted sample's own estimate (x=1), whatever that estimate is made of. Then the slope is\n"
          "'undiluted estimate minus floor' and must differ between two backgrounds by exactly the difference of\n"
          "their floors. So a different slope in the two backgrounds does NOT separate tumor from run; this table\n"
          "only checks that the series behave as straight-line mixtures (line fitted to the 5 dilution levels).")
    e = MASK["epithelial16"]
    und = load(*PAIRS["A1_all_strings"], "unclipped", [m["id"] for m in META if m["cancer"]])
    xs = np.array([0.5, 0.25, 0.125, 0.0625, 0.03125])
    res = collections.defaultdict(dict)
    for m in META:
        if not m["cancer"]:
            continue
        for bg in ("HU012", "HU083"):
            ids = [f"{m['id']}__d{p}" for p in ("50", "25", "12.5", "6.25", "3.125")]
            f = load(f"results/b3/{bg}/{{id}}/whole", f"results/b3/{bg}/{{id}}/cut", "unclipped", ids)
            for k, version in enumerate(("whole", "single")):
                y = np.array([f[i][k][e].sum() for i in ids])
                slope, inter = np.polyfit(xs, y, 1)
                res[version][(m["id"], bg)] = (slope, inter, slope + inter, und[m["id"]][k][e].sum())
    print("version\tmean_floor_HU012_bg\tmean_floor_HU083_bg\tfloor_difference\tmean_slope_HU012_bg\tmean_slope_HU083_bg\t"
          "slope_difference\tmean |line at x=1 minus undiluted estimate| HU012_bg\tsame HU083_bg")
    for version, r in res.items():
        s = {bg: np.array([v for (i, b), v in r.items() if b == bg]) for bg in ("HU012", "HU083")}
        print(f"{version}\t{s['HU012'][:, 1].mean():.3f}\t{s['HU083'][:, 1].mean():.3f}\t"
              f"{s['HU012'][:, 1].mean() - s['HU083'][:, 1].mean():+.3f}\t{s['HU012'][:, 0].mean():.3f}\t"
              f"{s['HU083'][:, 0].mean():.3f}\t{s['HU083'][:, 0].mean() - s['HU012'][:, 0].mean():+.3f}\t"
              f"{np.abs(s['HU012'][:, 2] - s['HU012'][:, 3]).mean():.3f}\t{np.abs(s['HU083'][:, 2] - s['HU083'][:, 3]).mean():.3f}")
    print("strings per dilution file (the two backgrounds do not hold the same amount of data):")
    for bg in ("HU012", "HU083"):
        n = [int(open(f"results/b3/{bg}/{m['id']}__d{p}/whole.nlines").read()) for m in META if m["cancer"]
             for p in ("50", "25", "12.5", "6.25", "3.125")]
        print(f"  {bg} background: median {int(np.median(n))}, range {min(n)}-{max(n)}")


def marker_groups():
    spans = []
    for line in list(open("data/ref/Markers.U250.hg38.tsv"))[1:]:
        f = line.split("\t")
        g = "gut" if f[5] in GUT else "placebo" if f[5] in PLACEBO else "other_epi" if f[5] in EPI else None
        if g:
            spans.append((int(f[3]), int(f[4]), g))
    spans.sort()
    return np.array([s[0] for s in spans]), np.array([s[1] for s in spans]), [s[2] for s in spans]


def t4():
    print("\n\nT4: fully 'off' strings at gut markers and at placebo markers (strings with 2+ beads, as fitted in step 2)")
    starts, ends, groups = marker_groups()
    ids = [(m["id"], f"results/step2/{m['id']}/whole.pat.gz") for m in META] + [(i, f"results/breast/{i}/whole.pat.gz") for i in BREAST]
    counts = {}
    for i, path in ids:
        c = collections.defaultdict(lambda: [0, 0])
        with gzip.open(path, "rt") as fh:
            for line in fh:
                _, first, pattern, n = line.rstrip("\n").split("\t")
                beads = pattern.replace(".", "")
                if len(beads) < 2:
                    continue
                k = np.searchsorted(starts, int(first), side="right") - 1
                if k < 0 or int(first) >= ends[k]:
                    continue
                key = (groups[k], min(len(beads), 4))
                c[key][0] += int(n)
                c[key][1] += int(n) * (set(beads) == {"T"})
        counts[i] = c
    tot = collections.Counter()
    for c in counts.values():
        for (g, s), (n, _) in c.items():
            tot[s] += n
    w = {s: tot[s] / sum(tot.values()) for s in (2, 3, 4)}
    print("common string-length mix used for every file and marker group: " + ", ".join(f"{s}{'+' if s == 4 else ''} beads {w[s]:.2f}" for s in w))

    def rate(c, g):
        if any(c[(g, s)][0] == 0 for s in w):
            return None
        return sum(w[s] * c[(g, s)][1] / c[(g, s)][0] for s in w)

    def pooled(sel, g):
        c = collections.defaultdict(lambda: [0, 0])
        for i in sel:
            for key, (n, k) in counts[i].items():
                c[key][0] += n
                c[key][1] += k
        return rate(c, g), sum(c[(g, s)][0] for s in w)

    grp = {"healthy HU012": [m["id"] for m in META if m["batch"] == "HU012"],
           "healthy HU083": [m["id"] for m in META if m["batch"] == "HU083"],
           "breast4": BREAST,
           "cancer low6": [m["id"] for m in META if m["cancer"] and m["tf"] < LOW],
           "cancer 20%+": [m["id"] for m in META if m["cancer"] and m["tf"] >= LOW]}
    print("group (strings pooled over its files)\tgut markers: all-off share (strings)\tother epithelial markers\tplacebo markers\tgut minus placebo")
    for gname, sel in grp.items():
        r = {g: pooled(sel, g) for g in ("gut", "other_epi", "placebo")}
        print(f"{gname}\t" + "\t".join(f"{100 * r[g][0]:.2f}% ({r[g][1]})" for g in r) + f"\t{100 * (r['gut'][0] - r['placebo'][0]):+.2f} points")
    print("per file:\nsample\tgroup\ttumor_share\tgut strings (2+ beads)\tgut all-off share\tplacebo all-off share\tgut minus placebo (points)")
    tf = {m["id"]: m["tf"] for m in META}
    xs, ys = [], []
    for gname, sel in grp.items():
        for i in sel:
            c = counts[i]
            g, p = rate(c, "gut"), rate(c, "placebo")
            n = sum(c[("gut", s)][0] for s in w)
            d = None if g is None or p is None else 100 * (g - p)
            print(f"{i}\t{gname}\t{tf.get(i, float('nan')):.3f}\t{n}\t{'' if g is None else f'{100 * g:.2f}%'}\t"
                  f"{'' if p is None else f'{100 * p:.2f}%'}\t{'' if d is None else f'{d:+.2f}'}")
            if gname.startswith("cancer") and d is not None:
                xs.append(tf[i])
                ys.append(d)
    print(f"rank correlation of 'gut minus placebo' with tumor share over {len(xs)} cancer samples: {spearmanr(xs, ys)[0]:+.2f}")


def t5():
    print("\n\nT5: scores after removing the age line fitted in the healthy controls (all strings, unclipped, epithelial score)")
    ids = [m["id"] for m in META]
    e = MASK["epithelial16"]
    fits = load(*PAIRS["A1_all_strings"], "unclipped", ids)
    fits.update(load("results/breast/{id}/whole", "results/breast/{id}/cut", "unclipped", BREAST))
    sc = {i: np.array([v[0][e].sum(), v[1][e].sum()]) for i, v in fits.items()}
    ctrl = [m for m in META if not m["cancer"]]
    x = np.array([[AGE[m["id"]], m["batch"] == "HU012", m["batch"] == "HU083"] for m in ctrl], dtype=float)
    slope = np.array([np.linalg.lstsq(x, np.array([sc[m["id"]][k] for m in ctrl]), rcond=None)[0][0] for k in (0, 1)])
    print(f"age line in controls (points of score per 10 years): whole {1000 * slope[0]:+.1f}, single {1000 * slope[1]:+.1f}")
    print("ages: controls 20-79, cancer patients 22-86 (patient P1, 12 samples, is 22-23: at the young edge of the controls)")
    adj = {i: v - slope * (AGE[i] - 45) for i, v in sc.items()}
    print("scores\tpositives\tnegatives\tauroc_whole\tauroc_single\tgain\tci_lo\tci_hi\tmoved")
    for tag, s in (("as fitted", sc), ("age line removed", adj)):
        rng = np.random.default_rng(0)
        neg = {k: np.array([s[m["id"]] for m in ctrl if k in ("all", m["batch"])]) for k in ("all", "HU012", "HU083")}
        for gname, keep in (("all20", lambda m: True), ("low6", lambda m: m["tf"] < LOW)):
            by_pat = collections.defaultdict(list)
            for m in META:
                if m["cancer"] and keep(m):
                    by_pat[m["patient"]].append(s[m["id"]])
            for k in neg:
                print(f"{tag}\t{gname}\t{k}\t" + fmt(boot_gain([np.array(v) for v in by_pat.values()], neg[k], rng)))
        for k in neg:
            print(f"{tag}\tbreast4\t{k}\t" + fmt(boot_gain([s[i][None, :] for i in BREAST], neg[k], rng)))
        print(f"{tag}\thealthy_HU012\tHU083\t" + fmt(boot_gain([v[None, :] for v in neg["HU012"]], neg["HU083"], rng)))


if __name__ == "__main__":
    print("EXPLORATORY factor tests, stage 1 (no fits). Plan and outcomes fixed in advance: step_scripts/factors/PLAN.md\n")
    t1()
    t2()
    t3()
    t4()
    t5()
