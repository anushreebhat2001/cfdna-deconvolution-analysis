import glob, os, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
P, OUT, RUN = "/scratch/mxf4642/code/project", "/scratch/mxf4642/results/m4", "/scratch/mxf4642/data/katsman_run"
meta = pd.read_csv(f"{P}/meta/samples.tsv", sep="\t")
atl = {"u250": f"{RUN}/31CT_U250_Netanel_atlas_over_regions.txt",
       "max": f"{RUN}/blood_epithelial_atlas_over_regions_max_markers.txt"}
names = {k: [h[:-5] for h in open(v).readline().rstrip("\n").split("\t") if h.endswith("_METH")] for k, v in atl.items()}
EPI = {"Bladder-Ep","Breast-Basal-Ep","Breast-Luminal-Ep","Colon-Ep","Fallopian-Ep","Gastric-Ep","Head-Neck-Ep",
       "Kidney-Ep","Lung-Ep-Alveo","Lung-Ep-Bron","Ovary-Ep","Prostate-Ep","Small-Int-Ep","Thyroid-Ep",
       "Pancreas-Acinar","Pancreas-Duct","Epithelial"}
rows = []
for f in glob.glob(f"{OUT}/est/*.txt"):
    ms, mode, s = os.path.basename(f)[:-4].split("_", 2)
    v = np.loadtxt(f).ravel(); v = v / v.sum()
    rows.append(dict(sample=s, markers=ms, clip=mode, epi=sum(x for n, x in zip(names[ms], v) if n in EPI)))
est = pd.DataFrame(rows).merge(meta, on="sample")
est.to_csv(f"{OUT}/m4_estimates.tsv", sep="\t", index=False)

print("\n=== ACCURACY (compare with paper)")
for (ms, mode), g in est.groupby(["markers", "clip"]):
    crc, ctl = g[g.group == "crc"], g[g.group == "control"]
    auc = roc_auc_score((g.group == "crc").astype(int), g.epi)
    rmse = np.sqrt(np.mean((crc.epi - crc.ichorCNA_TF) ** 2))
    r2 = np.corrcoef(crc.epi, crc.ichorCNA_TF)[0, 1] ** 2
    print(f"{ms:5s} {mode:7s} n={len(g):2d}  AUROC={auc:.3f}  RMSE={rmse:.3f}  R2={r2:.3f}  "
          f"control epi median={ctl.epi.median():.4f} max={ctl.epi.max():.4f}  crc epi median={crc.epi.median():.3f}")

print("\n=== CONTROL epithelial by batch (median)")
print(est[est.group == "control"].groupby(["markers", "clip", "batch"]).epi.median().unstack().round(4))

t = []
for f in glob.glob(f"{OUT}/logs/*.log"):
    for line in open(f):
        if line.startswith("TIMING"):
            _, name, sec, mem, rc = line.rstrip("\n").split("\t")
            t.append(dict(name=name, sec=float(sec[:-1]), gb=float(mem[:-2]), rc=rc))
t = pd.DataFrame(t)
t["step"] = t.name.str.split("_").str[:2].str.join("_").where(~t.name.str.startswith("convert"), "convert")
t.to_csv(f"{OUT}/m4_timings.tsv", sep="\t", index=False)
print("\n=== RUNTIME per run")
print(t.groupby("step").agg(n=("sec", "size"), median_s=("sec", "median"), max_s=("sec", "max"),
                             max_gb=("gb", "max"), failures=("rc", lambda r: (r != "rc=0").sum())).round(1))

print("\n=== PER SAMPLE epithelial fraction")
print(est.pivot_table(index=["group", "patient", "sample", "ichorCNA_TF"], columns=["markers", "clip"], values="epi").round(3).to_string())
