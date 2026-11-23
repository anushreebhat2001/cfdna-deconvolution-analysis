"""Per-CpG methylated / total read counts of one pat file at the atlas CpGs (site-level input for celfie.py).
Every read is counted once; reads starting up to SLOP CpGs before an atlas CpG are fetched.
usage: pat_to_sites.py SAMPLE.pat.gz ATLAS OUT.tsv"""
import sys, os, numpy as np, pandas as pd, pysam
pat, atlas, out = sys.argv[1:4]
SLOP = 150
A = pd.read_csv(atlas, sep="\t", usecols=["CHROM", "START", "END"])
M = np.zeros(len(A), dtype=np.int64)
D = np.zeros(len(A), dtype=np.int64)
row = {(c, int(s)): i for i, (c, s) in enumerate(zip(A.CHROM.values, A.START.values))}
idx = pat + ".tbi" if os.path.exists(pat + ".tbi") else pat + ".csi"
tb = pysam.TabixFile(pat, index=idx)
for c, g in A.groupby("CHROM", sort=False):
    wins = []
    for s in np.sort(g.START.values.astype(int)):
        a, b = max(0, s - 1 - SLOP), s
        if wins and a < wins[-1][1]:
            wins[-1][1] = max(wins[-1][1], b)
        else:
            wins.append([a, b])
    for a, b in wins:
        try:
            it = tb.fetch(c, a, b)
        except ValueError:
            continue
        for line in it:
            f = line.split("\t")
            rs, p, cnt = int(f[1]), f[2], int(f[3])
            for j, ch in enumerate(p):
                if ch == "C" or ch == "T":
                    k = row.get((c, rs + j))
                    if k is not None:
                        D[k] += cnt
                        if ch == "C":
                            M[k] += cnt
A["METH"], A["DEPTH"] = M, D
A.to_csv(out, sep="\t", index=False)
print(f"{os.path.basename(pat)}\tatlas_CpGs={len(A)}\tcovered={int(np.sum(D > 0))}\tcalls={int(D.sum())}\tfrac_meth={M.sum() / max(D.sum(), 1):.3f}")
