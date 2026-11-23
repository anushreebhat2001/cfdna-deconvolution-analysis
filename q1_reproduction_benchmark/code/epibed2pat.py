"""Convert a Biscuit epiBED (Nanopore) to wgbstools-style pat (CpG-number coordinates).
Usage: epibed2pat.py in.epibed.gz hg38_cpg_index.bed.gz out.pat   (unsorted, uncompressed)"""
import sys, re, gzip
import numpy as np, pysam

epibed, cpg_index, out = sys.argv[1], sys.argv[2], sys.argv[3]
TOK = re.compile(r"([A-Za-z])(\d*)")
tb = pysam.TabixFile(cpg_index)
st = dict(reads=0, reads_with_calls=0, calls=0, off_cpg=0, dup=0)
cur, pos, idx = None, np.array([]), np.array([])

with gzip.open(epibed, "rt") as f, open(out, "w") as fo:
    for line in f:
        fld = line.split("\t", 8)
        c = fld[0]
        if c != cur:
            cur = c
            try:
                rows = [r.split("\t") for r in tb.fetch(c)]
            except ValueError:
                rows = []
            pos = np.array([int(r[1]) for r in rows], dtype=np.int64)
            idx = np.array([int(r[3]) for r in rows], dtype=np.int64)
        st["reads"] += 1
        if len(pos) == 0:
            continue
        p, calls = int(fld[1]), {}
        for ch, num in TOK.findall(fld[6]):
            n = int(num) if num else 1
            if ch == "i":                     # insertion: no reference positions
                continue
            if ch in "MU":
                for j in range(n):
                    q = p + j
                    k = int(np.searchsorted(pos, q))
                    if k < len(pos) and pos[k] == q:          # call on the C
                        i = int(idx[k])
                    elif k > 0 and pos[k - 1] == q - 1:       # call on the G
                        i = int(idx[k - 1])
                    else:
                        st["off_cpg"] += 1
                        continue
                    if i in calls:
                        st["dup"] += 1
                        continue
                    calls[i] = "C" if ch == "M" else "T"
                    st["calls"] += 1
            p += n
        if calls:
            st["reads_with_calls"] += 1
            lo, hi = min(calls), max(calls)
            fo.write(f"{c}\t{lo}\t{''.join(calls.get(i, '.') for i in range(lo, hi + 1))}\t1\n")

print("STATS", epibed.split("/")[-1], st, flush=True)
