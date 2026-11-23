"""Per marker region: UXM counts (reads with >=4 CpGs in region; U = meth fraction <= 0.25),
CelFEER counts (reads with >=3 CpGs; 5 bins 0/.25/.5/.75/1), and a histogram of CpGs-in-region per read.
usage: count_regions.py file.pat.gz regions.bed out.npz"""
import sys, numpy as np, pysam
pat, regions, out = sys.argv[1:4]
regs = [l.split()[:3] for l in open(regions) if l.strip()]
R = len(regs)
ux_tot, ux_u, cf, hist = np.zeros(R), np.zeros(R), np.zeros((R, 5)), np.zeros(30)
import os
idx = pat + ".tbi" if os.path.exists(pat + ".tbi") else pat + ".csi"
tb = pysam.TabixFile(pat, index=idx)
SLOP = 80   # CpGs: reads starting up to this far before the region are fetched
for i, (c, s, e) in enumerate(regs):
    s, e = int(s), int(e)
    try:
        it = tb.fetch(c, max(0, s - 1 - SLOP), e)
    except ValueError:
        continue
    for line in it:
        f = line.split("\t")
        rs, p, cnt = int(f[1]), f[2], int(f[3])
        if rs + len(p) <= s or rs >= e:
            continue
        sub = p[max(0, s - rs): e - rs]
        m, u = sub.count("C"), sub.count("T")
        n = m + u
        if n == 0:
            continue
        hist[min(n, 29)] += cnt
        if n >= 4:
            ux_tot[i] += cnt
            if m / n <= 0.25:
                ux_u[i] += cnt
        if n >= 3:
            cf[i, int(np.floor(4 * m / n + 0.5))] += cnt
np.savez(out, ux_tot=ux_tot, ux_u=ux_u, cf=cf, ncpg=hist)
