import glob, os, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
P, RUN = "/scratch/mxf4642/code/project", "/scratch/mxf4642/data/katsman_run"
OUT, M4, MIX = "/scratch/mxf4642/results/m5", "/scratch/mxf4642/results/m4", "/scratch/mxf4642/data/mix_a"
meta = pd.read_csv(f"{P}/meta/samples.tsv", sep="\t")
atl = {"u250": f"{RUN}/31CT_U250_Netanel_atlas_over_regions.txt",
       "max": f"{RUN}/blood_epithelial_atlas_over_regions_max_markers.txt"}
names = {k: [h[:-5] for h in open(v).readline().rstrip("\n").split("\t") if h.endswith("_METH")] for k, v in atl.items()}
EPI = {"Bladder-Ep","Breast-Basal-Ep","Breast-Luminal-Ep","Colon-Ep","Fallopian-Ep","Gastric-Ep","Head-Neck-Ep",
       "Kidney-Ep","Lung-Ep-Alveo","Lung-Ep-Bron","Ovary-Ep","Prostate-Ep","Small-Int-Ep","Thyroid-Ep",
       "Pancreas-Acinar","Pancreas-Duct","Epithelial"}
def epi(f, ms):
    v = np.loadtxt(f).ravel(); v = v / v.sum()
    if len(v) != len(names[ms]):
        print("WARNING length mismatch", f, len(v), len(names[ms]))
    return sum(x for n, x in zip(names[ms], v) if n in EPI)

rows = []
for f in glob.glob(f"{OUT}/est/*.txt"):
    method, ms, id_ = os.path.basename(f)[:-4].split("@")
    if id_.startswith("CHECK"): continue
    rows.append(dict(method=method, markers=ms, id=id_, epi=epi(f, ms)))
m4 = pd.read_csv(f"{M4}/m4_estimates.tsv", sep="\t")      # CelFiE-ISH on the 42 originals
for _, r in m4.iterrows():
    rows.append(dict(method=f"celfieish-{r['clip']}", markers=r["markers"], id=r["sample"], epi=r["epi"]))
est = pd.DataFrame(rows)

mi = pd.concat([pd.read_csv(f, sep="\t") for f in glob.glob(f"{MIX}/*.mixinfo.tsv")])
ids = pd.concat([pd.DataFrame(dict(id=meta["sample"], sample=meta["sample"], dilution=1.0)),
                 mi[["mix_id", "sample", "dilution"]].rename(columns={"mix_id": "id"})])
est = est.merge(ids, on="id").merge(meta[["sample", "group", "patient", "ichorCNA_TF"]], on="sample")
est["est_TF"] = np.where(est.group == "crc", est.ichorCNA_TF * est.dilution, 0.0)
est.to_csv(f"{OUT}/m5_estimates.tsv", sep="\t", index=False)

print("\n=== DILUTION DIAGONAL: median over CRC samples of RMSE(epi_mix/epi_undiluted vs dilution)  [paper: CelFiE-ISH clip, max markers = 0.034 best]")
for (m, ms), g in est[est.group == "crc"].groupby(["method", "markers"]):
    r = []
    for s, gs in g.groupby("sample"):
        base = gs.loc[gs.dilution == 1.0, "epi"]
        if len(base) == 0 or base.iloc[0] <= 0: continue
        d = gs[gs.dilution < 1]
        r.append(np.sqrt(np.mean((d.epi / base.iloc[0] - d.dilution) ** 2)))
    print(f"  {m:17s} {ms:5s} n={len(r):2d}  median RMSE={np.median(r):.3f}")

print("\n=== LIMIT OF DETECTION: 120 virtual CRC samples in 10 equal bins vs 22 controls (AUROC)")
print("    paper (CelFiE-ISH/CelFEER): 1.7-3.0% -> 0.85-0.95 ; 3.1-4.3% -> 0.97-0.99 ; >4.3% -> >0.99")
tab = {}
for (m, ms), g in est.groupby(["method", "markers"]):
    ctl, crc = g[g.group == "control"], g[g.group == "crc"].sort_values("est_TF")
    if len(ctl) != 22 or len(crc) != 120:
        print(f"  skipping {m} {ms}: controls={len(ctl)} crc={len(crc)}"); continue
    col = []
    for b in np.array_split(crc, 10):
        y = np.r_[np.ones(len(b)), np.zeros(len(ctl))]
        col.append(roc_auc_score(y, np.r_[b.epi, ctl.epi]))
    tab[f"{m}|{ms}"] = col
    bins = [f"{100*b.est_TF.min():.1f}-{100*b.est_TF.max():.1f}%" for b in np.array_split(crc, 10)]
if tab:
    print(pd.DataFrame(tab, index=bins).round(3).to_string())

t = []
for f in glob.glob(f"{OUT}/logs/*.log"):
    for line in open(f):
        if line.startswith("TIMING"):
            _, name, sec, mem, rc = line.rstrip("\n").split("\t")
            t.append(dict(step="@".join(name.split("@")[:2]) if "@" in name else "mix", sec=float(sec[:-1]), gb=float(mem[:-2]), rc=rc))
t = pd.DataFrame(t); t.to_csv(f"{OUT}/m5_timings.tsv", sep="\t", index=False)
print("\n=== RUNTIME per run")
print(t.groupby("step").agg(n=("sec", "size"), median_s=("sec", "median"), max_s=("sec", "max"),
                             max_gb=("gb", "max"), failures=("rc", lambda r: (r != "rc=0").sum())).round(1))
