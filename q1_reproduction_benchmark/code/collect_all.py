import glob, os, zlib, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
P, RUN, R = "/scratch/mxf4642/code/project", "/scratch/mxf4642/data/katsman_run", "/scratch/mxf4642/results"
MIX = "/scratch/mxf4642/data/mix_a"
THRESH = 0.95   # detection-floor AUROC threshold; the sustained floor is the headline number
meta = pd.read_csv(f"{P}/meta/samples.tsv", sep="\t")
atl = {"u250": f"{RUN}/31CT_U250_Netanel_atlas_over_regions.txt",
       "u250cal": f"{RUN}/31CT_U250_Netanel_atlas_over_regions.txt",
       "max": f"{RUN}/blood_epithelial_atlas_over_regions_max_markers.txt"}
import glob as _g
for _k, _n in (("c250", "250"), ("c1000", "1000")):
    _f = _g.glob(f"{RUN}/blood_epithelial_atlas_over_regions*_{_n}.txt")
    if _f: atl[_k] = _f[0]
names = {k: [h[:-5] for h in open(v).readline().rstrip("\n").split("\t") if h.endswith("_METH")] for k, v in atl.items()}
EPI = {"Bladder-Ep","Breast-Basal-Ep","Breast-Luminal-Ep","Colon-Ep","Fallopian-Ep","Gastric-Ep","Head-Neck-Ep",
       "Kidney-Ep","Lung-Ep-Alveo","Lung-Ep-Bron","Ovary-Ep","Prostate-Ep","Small-Int-Ep","Thyroid-Ep",
       "Pancreas-Acinar","Pancreas-Duct","Epithelial"}
def epi_of(f, ms):
    v = np.loadtxt(f).ravel(); v = v / v.sum()
    return sum(x for n, x in zip(names[ms], v) if n in EPI)
rows = []
for d in ["m5", "m6", "m7", "m8", "m9"]:
    for f in glob.glob(f"{R}/{d}/est/*.txt"):
        method, ms, id_ = os.path.basename(f)[:-4].split("@")
        if not id_.startswith("CHECK"):
            rows.append(dict(method=method, markers=ms, id=id_, epi=epi_of(f, ms)))
for _, r in pd.read_csv(f"{R}/m4/m4_estimates.tsv", sep="\t").iterrows():
    rows.append(dict(method=f"celfieish-{r['clip']}", markers=r["markers"], id=r["sample"], epi=r["epi"]))
est = pd.DataFrame(rows)
mi = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"{MIX}/*.mixinfo.tsv")])
ids = pd.concat([pd.DataFrame(dict(id=meta["sample"], sample=meta["sample"], dilution=1.0)),
                 mi[["mix_id", "sample", "dilution"]].rename(columns={"mix_id": "id"})])
est = est.merge(ids, on="id").merge(meta[["sample", "group", "patient", "batch", "ichorCNA_TF", "filtered_frags"]], on="sample")
est["est_TF"] = np.where(est.group == "crc", est.ichorCNA_TF * est.dilution, 0.0)
est.to_csv(f"{R}/summary/all_estimates.tsv", sep="\t", index=False)

print("NOTE: uxm/celfeer 'max'     = 7 groups from real Loyfer reads + Megakaryocytes from SYNTHETIC independent-site reads")
print("      uxm/celfeer 'u250cal' = U250 with the Blood-Granul reference rebuilt from synthetic independent-site reads (calibration)\n")

def rng_for(*key):
    return np.random.default_rng(zlib.crc32("|".join(key).encode()))
def auc(pos, neg):
    return roc_auc_score(np.r_[np.ones(len(pos)), np.zeros(len(neg))], np.r_[pos, neg])
def boot_ci(posdf, neg, rng, B=1000):
    pats = posdf.patient.unique()
    byp = {p: posdf.epi[posdf.patient == p].values for p in pats}
    out = []
    for _ in range(B):
        pos = np.concatenate([byp[p] for p in rng.choice(pats, len(pats))])
        out.append(auc(pos, neg[rng.integers(0, len(neg), len(neg))]))
    return np.percentile(out, [2.5, 97.5])
groups = est.groupby(["method", "markers"])

print("=== 1. 42 ORIGINAL SAMPLES: CRC vs controls (95% CI = bootstrap over patients and controls)")
for (m, ms), g in groups:
    o = g[g.dilution == 1.0]; crc, ctl = o[o.group == "crc"], o[o.group == "control"]
    if len(crc) != 20 or len(ctl) != 22: print(f"  {m} {ms}: incomplete ({len(crc)} crc, {len(ctl)} ctl)"); continue
    lo, hi = boot_ci(crc, ctl.epi.values, rng_for(m, ms, "s1"))
    rmse = np.sqrt(np.mean((crc.epi - crc.ichorCNA_TF) ** 2)); r2 = np.corrcoef(crc.epi, crc.ichorCNA_TF)[0, 1] ** 2
    print(f"  {m:17s} {ms:7s} AUROC={auc(crc.epi, ctl.epi):.3f} [{lo:.3f}-{hi:.3f}]  RMSE={rmse:.3f}  R2={r2:.3f}  control median={ctl.epi.median():.4f}")

print(f"\n=== 2. LIMIT OF DETECTION per bin (AUROC [95% CI], #patients in bin)  paper: 1.7-3.0% 0.85-0.95; 3.1-4.3% 0.97-0.99")
lod, floors = [], []
for (m, ms), g in groups:
    ctl = g[(g.group == "control") & (g.dilution == 1.0)]; crc = g[g.group == "crc"].sort_values("est_TF")
    if len(ctl) != 22 or len(crc) != 120: continue
    rng = rng_for(m, ms, "s2")
    cells, bins = [], []
    for b in np.array_split(crc, 10):
        a = auc(b.epi, ctl.epi); lo, hi = boot_ci(b, ctl.epi.values, rng)
        lab = f"{100*b.est_TF.min():.1f}-{100*b.est_TF.max():.1f}%"
        cells.append(f"{lab}: {a:.2f} [{lo:.2f}-{hi:.2f}] ({b.patient.nunique()}p)")
        bins.append((lab, a))
        lod.append(dict(method=m, markers=ms, bin=lab, tf_lo=b.est_TF.min(), tf_hi=b.est_TF.max(),
                        auroc=a, ci_lo=lo, ci_hi=hi, n_patients=b.patient.nunique()))
    first = next((l for l, a in bins if a >= THRESH), "none")
    sust = next((bins[i][0] for i in range(len(bins)) if all(a >= THRESH for _, a in bins[i:])), "none")
    floors.append((m, ms, sust, first))
    print(f"  {m} | {ms}   floor (AUROC>={THRESH}): sustained={sust}  first crossing={first}"); [print("     " + c) for c in cells]
pd.DataFrame(lod).to_csv(f"{R}/summary/lod_bins.tsv", sep="\t", index=False)
print(f"\n=== 2b. DETECTION FLOOR SUMMARY (AUROC >= {THRESH}; sustained = headline, first crossing = secondary)")
for m, ms, s, f in floors:
    print(f"  {m:17s} {ms:7s} sustained={s:13s} first crossing={f}{'   <- differ' if s != f else ''}")

print("\n=== 3. DILUTION DIAGONAL RMSE (median over 20 CRC samples)  A=normalized, C=raw scale [paper best 0.034]")
for (m, ms), g in groups:
    a, c = [], []
    for s, gs in g[g.group == "crc"].groupby("sample"):
        base = gs.loc[gs.dilution == 1.0, "epi"]
        if len(base) == 0 or base.iloc[0] <= 0: continue
        b, d = base.iloc[0], gs[gs.dilution < 1]
        a.append(np.sqrt(np.mean((d.epi / b - d.dilution) ** 2))); c.append(np.sqrt(np.mean((d.epi - d.dilution * b) ** 2)))
    if a: print(f"  {m:17s} {ms:7s} A={np.median(a):.3f}  C={np.median(c):.3f}  (n={len(a)})")

print("\n=== 4. COVERAGE: controls' epithelial estimate vs sequencing depth")
for (m, ms), g in groups:
    ctl = g[(g.group == "control") & (g.dilution == 1.0)]
    if len(ctl) < 22: continue
    rho, p = spearmanr(ctl.filtered_frags, ctl.epi)
    med = ctl.groupby("batch").epi.median()
    print(f"  {m:17s} {ms:7s} spearman(depth, epi)={rho:+.2f} (p={p:.3f})  median HU012={med.get('HU012', np.nan):.4f} HU083={med.get('HU083', np.nan):.4f}")

print("\n=== 5. RUNTIME per run")
t = []
for d in ["m4", "m5", "m6", "m7", "m8", "m9"]:
    for f in glob.glob(f"{R}/{d}/logs/*.log"):
        for line in open(f):
            if line.startswith("TIMING"):
                _, name, sec, mem, rc = line.rstrip("\n").split("\t")
                step = "@".join(name.split("@")[:2]) if "@" in name else name.split("_")[0]
                if d == "m4" and not name.startswith("convert"): step = "m4:" + "_".join(name.split("_")[:2])
                t.append(dict(step=step, sec=float(sec[:-1]), gb=float(mem[:-2]), rc=rc))
t = pd.DataFrame(t); t.to_csv(f"{R}/summary/all_timings.tsv", sep="\t", index=False)
print(t.groupby("step").agg(n=("sec", "size"), median_min=("sec", lambda s: s.median() / 60), max_min=("sec", lambda s: s.max() / 60),
                            max_gb=("gb", "max"), failures=("rc", lambda r: (r != "rc=0").sum())).round(2).to_string())

print("\n=== 6. MARKER CpGs PER READ (U250 regions, reads with >=1 call in a region), share of reads")
for grp in ["control", "crc"]:
    h = sum(np.load(f"{R}/m6/feat/{s}.npz")["ncpg"] for s in meta[meta.group == grp]["sample"])
    tot = h.sum()
    print(f"  {grp:8s} 1:{h[1]/tot:.1%}  2:{h[2]/tot:.1%}  3:{h[3]/tot:.1%}  >=4:{h[4:].sum()/tot:.1%}  (reads={int(tot)})")
