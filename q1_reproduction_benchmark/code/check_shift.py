import sys, gzip, collections, numpy as np
atlas, pat = sys.argv[1], sys.argv[2]
hdr = open(atlas).readline().rstrip("\n").split("\t")
blood = [i for i, h in enumerate(hdr) if h.endswith("_METH") and (h.startswith("Blood") or h.startswith("Lymph"))]
A = {}
with open(atlas) as f:
    next(f)
    for line in f:
        x = line.rstrip("\n").split("\t")
        m = sum(float(x[i]) for i in blood); c = sum(float(x[i + 1]) for i in blood)
        if c > 0:
            A[(x[0], int(x[1]))] = m / c
pairs = {sh: ([], []) for sh in (-1, 0, 1)}
with gzip.open(pat, "rt") as f:
    for line in f:
        c, s, p = line.split("\t")[:3]
        s = int(s)
        for j, ch in enumerate(p):
            if ch == ".":
                continue
            for sh in (-1, 0, 1):
                b = A.get((c, s + j + sh))
                if b is not None:
                    pairs[sh][0].append(b); pairs[sh][1].append(1.0 if ch == "C" else 0.0)
for sh in (-1, 0, 1):
    a, y = np.array(pairs[sh][0]), np.array(pairs[sh][1])
    lo, hi = y[a < 0.3].mean(), y[a > 0.7].mean()
    print(f"{pat.split('/')[-1]}  shift {sh:+d}: n={len(a)}  r={np.corrcoef(a, y)[0,1]:.3f}  "
          f"meth% where blood<0.3: {100*lo:.1f}  where blood>0.7: {100*hi:.1f}")
