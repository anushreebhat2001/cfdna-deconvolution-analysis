"""Synthetic reads under an independent-sites assumption (professor's 'simple Option 1').
Template reads (start CpG, length, which CpGs carry a call) come from a real pat file. Every call is
redrawn independently: methylated with probability METH/COV of one atlas column at that CpG.
CpGs not in the atlas, or with COV=0, become '.' (no call). Each output read has count 1.
usage: synth_reads.py TEMPLATE.pat.gz ATLAS COLUMN REGIONS OUT.pat SEED"""
import sys, os, numpy as np, pandas as pd, pysam
tpl, atlas, col, regions, out, seed = sys.argv[1:7]
rng = np.random.default_rng(int(seed))
A = pd.read_csv(atlas, sep="\t", usecols=["START", f"{col}_METH", f"{col}_COV"])
A = A[A[f"{col}_COV"] > 0]
pmeth = dict(zip(A["START"].tolist(), (A[f"{col}_METH"] / A[f"{col}_COV"]).tolist()))
print(f"atlas CpGs with {col} coverage: {len(pmeth)}  mean beta={np.mean(list(pmeth.values())):.3f}", flush=True)
SLOP = 80   # same as count_regions.py
win = {}
for l in open(regions):
    if not l.strip():
        continue
    c, s, e = l.split()[:3]
    win.setdefault(c, []).append((max(0, int(s) - 1 - SLOP), int(e)))
idx = tpl + ".tbi" if os.path.exists(tpl + ".tbi") else tpl + ".csi"
tb = pysam.TabixFile(tpl, index=idx)
n_in = n_out = calls = meth = 0
with open(out, "w") as fo:
    for c in win:
        ws = sorted(win[c])
        merged = [list(ws[0])]
        for a, b in ws[1:]:                      # merge overlapping windows -> no read written twice
            if a < merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        for a, b in merged:
            try:
                it = tb.fetch(c, a, b)
            except ValueError:
                continue
            for line in it:
                f = line.rstrip("\n").split("\t")
                rs, pat, cnt = int(f[1]), f[2], int(f[3])
                pos = [j for j, ch in enumerate(pat) if ch in "CT"]
                q = [pmeth.get(rs + j, -1.0) for j in pos]
                n_in += cnt
                for _ in range(cnt):
                    s_ = ["."] * len(pat)
                    draw = rng.random(len(pos))
                    for j, qq, d in zip(pos, q, draw):
                        if qq >= 0:
                            if d < qq:
                                s_[j] = "C"; meth += 1
                            else:
                                s_[j] = "T"
                            calls += 1
                    fo.write(f"{c}\t{rs}\t{''.join(s_)}\t1\n")
                    n_out += 1
print(f"synthetic reads written: {n_out} (template reads {n_in})  calls drawn: {calls}  "
      f"fraction methylated: {meth / max(calls, 1):.3f}", flush=True)
