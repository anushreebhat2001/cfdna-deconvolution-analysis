"""celfie.py output -> one line of proportions in atlas order (same format as the other methods).
usage: celfie_to_line.py OUTDIR ATLAS OUT.txt"""
import sys, glob, numpy as np, pandas as pd
outdir, atlas, out = sys.argv[1:4]
names = [h[:-5] for h in open(atlas).readline().rstrip("\n").split("\t") if h.endswith("_METH")]
f = sorted(glob.glob(f"{outdir}/*tissue_proportions*.txt"))
assert f, "no tissue_proportions file in " + outdir
df = pd.read_csv(f[0], sep="\t", index_col=0)
strip = lambda k: str(k)[:-5] if str(k).endswith("_METH") else str(k)
if all(n in [strip(c) for c in df.columns] for n in names):
    lab, how = {strip(k): v for k, v in df.iloc[0].items()}, "columns"
elif all(n in [strip(i) for i in df.index] for n in names):
    lab, how = {strip(k): v for k, v in df.iloc[:, 0].items()}, "index"
else:
    vals = df.iloc[0].values if df.shape[1] >= len(names) else df.iloc[:, 0].values
    assert len(vals) >= len(names), f"cannot match output {df.shape} to {len(names)} cell types"
    lab, how = dict(zip(names, vals[:len(names)])), "position"
v = np.array([lab[n] for n in names], dtype=float)
assert np.isfinite(v).all() and v.sum() > 0, f"bad proportions: {v}"
np.savetxt(out, (100 * v / v.sum()).reshape(1, -1), fmt="%.6f", delimiter="\t")
print(f"matched by {how}; file {f[0].split('/')[-1]}; output columns: {list(df.columns)[:4]}...")
