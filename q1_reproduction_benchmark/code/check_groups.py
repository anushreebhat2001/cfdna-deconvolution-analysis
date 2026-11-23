"""Which Loyfer 31-groups make up each column of an atlas?
Sums each group's (meth, cov) from the Loyfer .hg38.beta files at the atlas CpGs, then fits
atlas COV = sum_k w_k * group_k COV with non-negative least squares. w ~ 1 means 'member'.
usage: check_groups.py ATLAS GROUP_CSV LOYFER_DIR OUT_TSV"""
import sys, glob, numpy as np, pandas as pd
from scipy.optimize import nnls
atlas, grp_csv, loy, out = sys.argv[1:5]
A = pd.read_csv(atlas, sep="\t")
idx = A["START"].values - 1                      # CpG number -> 0-based row of .beta
names = [c[:-5] for c in A.columns if c.endswith("_METH")]
g = pd.read_csv(grp_csv)
groups = sorted(g.group.unique())
M = {k: np.zeros(len(A)) for k in groups}
C = {k: np.zeros(len(A)) for k in groups}
n_ok = 0
for s, grp in zip(g.name, g.group):
    f = glob.glob(f"{loy}/*_{s}.hg38.beta")
    if not f:
        print("missing beta:", s); continue
    b = np.memmap(f[0], dtype=np.uint8, mode="r").reshape(-1, 2)
    if n_ok == 0:
        print("beta rows:", b.shape[0], "| max atlas CpG:", idx.max() + 1)
    v = np.asarray(b[idx], dtype=float)
    M[grp] += v[:, 0]; C[grp] += v[:, 1]; n_ok += 1
print(f"beta files read: {n_ok} / {len(g)}")
X = np.column_stack([C[k] for k in groups])
rows = []
for a in names:
    y = A[f"{a}_COV"].values.astype(float)
    w, res = nnls(X, y)
    r2 = 1 - res**2 / max(((y - y.mean())**2).sum(), 1e-9)
    mem = [groups[i] for i in np.where(w > 0.5)[0]]
    pred = X[:, w > 0.5].sum(axis=1) if mem else np.zeros(len(y))
    exact = np.mean(pred == y)
    with np.errstate(all="ignore"):
        bm = np.where(pred > 0, sum(M[k] for k in mem) / pred, np.nan) if mem else np.full(len(y), np.nan)
        ba = np.where(y > 0, A[f"{a}_METH"].values / y, np.nan)
    ok = ~np.isnan(bm) & ~np.isnan(ba)
    corr = np.corrcoef(bm[ok], ba[ok])[0, 1] if ok.sum() > 10 else np.nan
    print(f"\n== {a}:  R2={r2:.3f}  members(w>0.5)={mem}  exact COV match={exact:.3f}  beta corr={corr:.3f}")
    for i in np.argsort(-w):
        if w[i] > 0.05:
            print(f"     {groups[i]:20s} w={w[i]:.3f}")
    rows.append({"atlas_group": a, "members": ";".join(mem), "R2": round(r2, 4),
                 "exact_cov_match": round(exact, 4), "beta_corr": round(corr, 4) if corr == corr else "NA"})
pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
print("\nwrote", out)
