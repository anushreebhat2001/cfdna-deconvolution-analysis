import argparse
import collections
import os
import re

import numpy as np

EPI = {"u250": {"Bladder-Ep", "Breast-Basal-Ep", "Breast-Luminal-Ep", "Colon-Ep", "Fallopian-Ep",
                "Gastric-Ep", "Head-Neck-Ep", "Kidney-Ep", "Lung-Ep-Alveo", "Lung-Ep-Bron",
                "Ovary-Ep", "Prostate-Ep", "Small-Int-Ep", "Thyroid-Ep", "Pancreas-Acinar",
                "Pancreas-Duct"},
       "max": {"Epithelial"}}
ATLAS = {"u250": "data/ref/31CT_U250_Netanel_atlas_over_regions.txt",
         "max": "data/ref/blood_epithelial_atlas_over_regions_max_markers.txt"}
DIL = re.compile(r"__d([\d.]+)$")
N_BOOT = 2000
TUMOR_BINS = [("tumor_below_1pct", 0, 0.01), ("tumor_1_to_5pct", 0.01, 0.05),
              ("tumor_5_to_20pct", 0.05, 0.20), ("tumor_20pct_up", 0.20, 2)]


def epi_share(path, mask):
    v = np.loadtxt(path).ravel()
    return float(v[mask].sum() / v.sum())


def auroc(pos, neg):
    d = pos[:, None] - neg[None, :]
    return ((d > 0).sum() + 0.5 * (d == 0).sum()) / d.size


def gains(ca, cb, exp, ha, hb):
    err_a = np.abs(ca - exp[:, None]).mean(0)
    err_b = np.abs(cb - exp[:, None]).mean(0)
    draws = range(ca.shape[1])
    auc_a = np.array([auroc(ca[:, d], ha[:, d]) for d in draws])
    auc_b = np.array([auroc(cb[:, d], hb[:, d]) for d in draws])
    mean = (err_a.mean(), err_b.mean(), (err_b - err_a).mean(),
            auc_a.mean(), auc_b.mean(), (auc_a - auc_b).mean())
    return mean, err_b - err_a, auc_a - auc_b


def summarise(rows, rng):
    by_pat = collections.defaultdict(list)
    ctrl = []
    for r in rows:
        (by_pat[r["patient"]] if r["is_cancer"] else ctrl).append(r)
    pats = [(np.array([r["a"] for r in v]), np.array([r["b"] for r in v]), np.array([r["exp"] for r in v]))
            for v in by_pat.values()]
    ha, hb = np.array([r["a"] for r in ctrl]), np.array([r["b"] for r in ctrl])

    def score(pick_p, pick_h):
        ca, cb, exp = (np.concatenate([pats[i][k] for i in pick_p]) for k in range(3))
        return gains(ca, cb, exp, ha[pick_h], hb[pick_h])

    point, err_draws, auc_draws = score(range(len(pats)), slice(None))
    boot = np.array([score(rng.integers(0, len(pats), len(pats)), rng.integers(0, len(ha), len(ha)))[0]
                     for _ in range(N_BOOT)])
    ci = {k: np.percentile(boot[:, k], [2.5, 97.5]) for k in (2, 5)}
    return point, ci, sum(len(p[2]) for p in pats), len(ha), len(pats), err_draws.std(), auc_draws.std()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--markers", choices=["u250", "max"], required=True)
    ap.add_argument("--a", required=True)
    ap.add_argument("--b", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--draws", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    header = open(ATLAS[args.markers]).readline().rstrip("\n").split("\t")
    names = [h[:-5] for h in header if h.endswith("_METH")]
    mask = np.array([n in EPI[args.markers] for n in names])

    meta = {}
    for line in list(open(args.meta))[1:]:
        f = line.rstrip("\n").split("\t")
        meta[f[0]] = dict(group=f[1], patient=f[2], batch=f[3], ichor=float(f[5]))
    ids = list(meta)
    for s in meta:
        if meta[s]["group"] == "crc":
            ids += [f"{s}__d{p}" for p in ("50", "25", "12.5", "6.25", "3.125")]

    out = []
    for clip in ("unclipped", "clipped"):
        rows, missing = [], 0
        for i in ids:
            pa = [args.a.format(id=i, draw=d) + f".ish_{clip}.txt" for d in range(1, args.draws + 1)]
            pb = [args.b.format(id=i, draw=d) + f".ish_{clip}.txt" for d in range(1, args.draws + 1)]
            if not all(os.path.exists(p) for p in pa + pb):
                missing += 1
                continue
            m = DIL.search(i)
            s, dil = (i[:m.start()], float(m.group(1)) / 100) if m else (i, 1.0)
            info = meta[s]
            cancer = info["group"] == "crc"
            rows.append(dict(id=i, patient=info["patient"], batch=info["batch"], is_cancer=cancer,
                             dil=dil if cancer else 0.0, exp=info["ichor"] * dil if cancer else 0.0,
                             a=np.array([epi_share(p, mask) for p in pa]),
                             b=np.array([epi_share(p, mask) for p in pb])))
        rng = np.random.default_rng(0)
        out.append(f"\n=== {args.label} | CelFiE-ISH {clip} | inputs scored {len(rows)}, missing {missing}, "
                   f"draws {args.draws}")
        out.append("subset\tcontrols_used\tcancer_files\tpatients\tcontrols\terror_A\terror_B\terror_gain\t"
                   "err_ci_lo\terr_ci_hi\terror_moved\tauroc_A\tauroc_B\tauroc_gain\tauc_ci_lo\tauc_ci_hi\t"
                   "auroc_moved\terror_gain_sd_over_draws\tauroc_gain_sd_over_draws")
        subsets = [("all", lambda r: True)]
        subsets += [(f"dilution_{d:g}", (lambda r, d=d: r["dil"] == d))
                    for d in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)]
        subsets += [(name, (lambda r, lo=lo, hi=hi: lo <= r["exp"] < hi)) for name, lo, hi in TUMOR_BINS]
        for name, batch, keep in [(n, b, k) for n, k in subsets for b in ("all", "HU012", "HU083")]:
            sub = [r for r in rows if (keep(r) if r["is_cancer"] else batch in ("all", r["batch"]))]
            if not any(r["is_cancer"] for r in sub) or all(r["is_cancer"] for r in sub):
                continue
            p, ci, n_c, n_h, n_p, sd_e, sd_a = summarise(sub, rng)
            moved_e = "yes" if ci[2][0] > 0 or ci[2][1] < 0 else "no"
            moved_a = "yes" if ci[5][0] > 0 or ci[5][1] < 0 else "no"
            out.append(f"{name}\t{batch}\t{n_c}\t{n_p}\t{n_h}\t{p[0]:.4f}\t{p[1]:.4f}\t{p[2]:+.4f}\t"
                       f"{ci[2][0]:+.4f}\t{ci[2][1]:+.4f}\t{moved_e}\t{p[3]:.3f}\t{p[4]:.3f}\t{p[5]:+.3f}\t"
                       f"{ci[5][0]:+.3f}\t{ci[5][1]:+.3f}\t{moved_a}\t{sd_e:.4f}\t{sd_a:.3f}")
        out.append("healthy controls: estimated epithelial share (should be about 0)")
        out.append("group\tcontrols\tmedian_A\tmedian_B\tmean_A_minus_B\tci_lo\tci_hi")
        for name, keep in (("all", lambda r: True), ("HU012", lambda r: r["batch"] == "HU012"),
                           ("HU083", lambda r: r["batch"] == "HU083")):
            h = [r for r in rows if not r["is_cancer"] and keep(r)]
            if not h:
                continue
            ha, hb = np.array([r["a"] for r in h]), np.array([r["b"] for r in h])
            d = (ha - hb).mean(1)
            bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(N_BOOT)]
            lo, hi = np.percentile(bs, [2.5, 97.5])
            out.append(f"{name}\t{len(h)}\t{np.median(ha, axis=0).mean():.4f}\t{np.median(hb, axis=0).mean():.4f}\t"
                       f"{d.mean():+.4f}\t{lo:+.4f}\t{hi:+.4f}")
        with open(f"{args.out}/{args.label}_{clip}_estimates.tsv", "w") as fh:
            fh.write("id\tpatient\tbatch\tis_cancer\tdilution\texpected\test_A\test_B\n")
            for r in rows:
                fh.write(f"{r['id']}\t{r['patient']}\t{r['batch']}\t{int(r['is_cancer'])}\t{r['dil']:g}\t"
                         f"{r['exp']:.5f}\t{r['a'].mean():.5f}\t{r['b'].mean():.5f}\n")
    text = "\n".join(out) + "\n"
    with open(f"{args.out}/{args.label}_summary.tsv", "w") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
