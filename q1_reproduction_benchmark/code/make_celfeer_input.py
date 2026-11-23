import sys, numpy as np
feat, ref, sid, out = sys.argv[1:5]
x, r = np.load(feat), np.load(ref)
cf_s, cf_r, keep, names, regs = x["cf"], r["cf"], r["keep"], list(r["names"]), r["regs"]
with open(out, "w") as fo:
    fo.write("\t".join(["chrom", "start", "end", sid, "chrom", "start", "end"] + names) + "\n")
    for k in np.where(keep)[0]:
        c, s, e = regs[k]
        row = [c, s, e] + [str(int(v)) for v in cf_s[k]] + [c, s, e] + [str(int(v)) for v in cf_r[k].ravel()]
        fo.write("\t".join(row) + "\n")
