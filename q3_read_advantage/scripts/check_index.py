import collections, gzip, re, sys
import numpy as np
idx = np.load(sys.argv[1]); markers = sys.argv[2]; epibed = sys.argv[3]
CH = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
pos = {c: idx[f"pos_{c}"].astype(np.int64) for c in CH}
off = {c: int(idx[f"off_{c}"]) for c in CH}

n = ok1 = ok0 = oklen = 0
for line in open(markers):
    if line.startswith("#"): continue
    f = line.split("\t"); c = f[0]
    if c not in pos: continue
    s, e, sc, ec = int(f[1]), int(f[2]), int(f[3]), int(f[4])
    n += 1
    cpos = int(pos[c][sc - off[c] - 1])
    ok1 += (cpos + 1 == s); ok0 += (cpos == s)
    lo, hi = np.searchsorted(pos[c], [s - 1, e - 1])
    oklen += (int(hi - lo) == ec - sc)
print(f"markers checked {n}: first CpG at 1-based start {ok1}, at 0-based start {ok0}, CpG count matches {oklen}")

TOKEN = re.compile(r"([A-Za-z])(\d*)")
calls = collections.defaultdict(list)
with gzip.open(epibed, "rt") as fh:
    for i, line in enumerate(fh):
        if i >= 400000: break
        f = line.split("\t", 7)
        if f[0] not in pos: continue
        p = int(f[1])
        for ch, num in TOKEN.findall(f[6]):
            k = int(num) if num else 1
            if ch in "MU":
                calls[f[0]].extend(range(p, p + k))
            if ch != "i": p += k
tot = sum(len(v) for v in calls.values()); res = {}
for d in (-1, 0, 1):
    h = 0
    for c, v in calls.items():
        q = np.asarray(v, dtype=np.int64) - d
        j = np.minimum(np.searchsorted(pos[c], q), len(pos[c]) - 1)
        h += int((pos[c][j] == q).sum())
    res[d] = h
print(f"calls {tot}: on C (offset 0) {res[0]}, on G (offset 1) {res[1]}, one before C {res[-1]}")
