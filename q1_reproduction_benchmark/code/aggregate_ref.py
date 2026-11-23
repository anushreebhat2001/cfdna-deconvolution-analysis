import sys, numpy as np, pandas as pd
grp_csv, npzdir, atlas, regions, outdir = sys.argv[1:6]
tag = sys.argv[6] if len(sys.argv) > 6 else "u250"
names = [h[:-5] for h in open(atlas).readline().rstrip("\n").split("\t") if h.endswith("_METH")]
g = pd.read_csv(grp_csv)
g["group"] = g["group"].replace({"Ovary+Endom-Ep": "Ovary-Ep"})
print("atlas groups missing from csv:", set(names) - set(g.group), "| csv groups not in atlas:", set(g.group) - set(names))
regs = [l.split()[:3] for l in open(regions) if l.strip()]
R, G = len(regs), len(names)
U, CF = np.full((R, G), np.nan), np.zeros((R, G, 5))
for j, n in enumerate(names):
    fr = []
    for s in g[g.group == n]["name"]:
        d = np.load(f"{npzdir}/{s}.npz")
        fr.append(np.where(d["ux_tot"] > 0, d["ux_u"] / np.maximum(d["ux_tot"], 1), np.nan))
        CF[:, j] += d["cf"]
    with np.errstate(all="ignore"):
        U[:, j] = np.nanmean(np.vstack(fr), axis=0)
    print(f"  {n:20s} samples={len(fr):2d}  median UXM-eligible reads/region={np.median(sum(np.load(f'{npzdir}/{s}.npz')['ux_tot'] for s in g[g.group == n]['name'])):.0f}")
keep_u = ~np.isnan(U).any(axis=1)
keep_cf = (CF.sum(axis=2) > 0).all(axis=1)
print(f"regions: {R} | usable for UXM: {keep_u.sum()} | usable for CelFEER: {keep_cf.sum()}")
with open(f"{outdir}/uxm_{tag}_atlas.tsv", "w") as fo, open(f"{outdir}/uxm_{tag}_regions.bed", "w") as fb:
    fo.write("\t".join(["chr", "start", "end", "x1", "x2", "x3"] + names) + "\n")
    for k in np.where(keep_u)[0]:
        c, s, e = regs[k]
        fo.write("\t".join([c, s, e, ".", ".", "."] + [f"{x:.6f}" for x in U[k]]) + "\n")
        fb.write(f"{c}\t{s}\t{e}\n")
np.savez(f"{outdir}/celfeer_{tag}_ref.npz", cf=CF, keep=keep_cf, names=np.array(names), regs=np.array(regs))
