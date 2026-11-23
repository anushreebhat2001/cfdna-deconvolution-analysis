"""In-silico dilution as in Katsman in_silico_mixing.sh, done on pat lines.
Each mixture has N_s reads in total (N_s = reads in the cancer sample): f*N_s from cancer, (1-f)*N_s from the pool.
Reads without CpG calls are not in pat files, so the number of informative reads drawn is hypergeometric."""
import sys, gzip, numpy as np, pandas as pd
meta_f, totals_f, pool_f, patdir, outdir, sample, seed = sys.argv[1:8]
DIL = [0.5, 0.25, 0.125, 0.0625, 0.03125]
tot = pd.read_csv(totals_f, sep="\t", index_col=0)
meta = pd.read_csv(meta_f, sep="\t")
pool = meta[meta.in_healthy_pool == "yes"]["sample"].tolist()
N_pool, M_pool = int(tot.loc[pool, "n_reads"].sum()), int(tot.loc[pool, "n_pat"].sum())
N_s, M_s = int(tot.loc[sample, "n_reads"]), int(tot.loc[sample, "n_pat"])
rng = np.random.default_rng(int(seed))

canc = [l.rstrip("\n") for l in gzip.open(f"{patdir}/{sample}.pat.gz", "rt")]
assert len(canc) == M_s, (len(canc), M_s)
mask = np.zeros(M_pool, dtype=np.uint8)
outs, info = [], []
for b, f in enumerate(DIL):
    k_c, k_h = int(f * N_s), int((1 - f) * N_s)
    c_inf = int(rng.hypergeometric(M_s, N_s - M_s, k_c))
    h_inf = int(rng.hypergeometric(M_pool, N_pool - M_pool, k_h))
    ci = rng.choice(M_s, c_inf, replace=False)
    hi = rng.choice(M_pool, h_inf, replace=False)
    mask[hi] |= np.uint8(1 << b)
    mix_id = f"{sample}__d{f*100:g}"
    fo = open(f"{outdir}/{mix_id}.pat", "w")
    for i in ci:
        fo.write(canc[i] + "\tT\n")
    outs.append(fo)
    info.append(dict(mix_id=mix_id, sample=sample, dilution=f, reads_total=k_c + k_h,
                     reads_cancer=k_c, pat_cancer=c_inf, pat_healthy=h_inf))
del canc
mb = mask.tobytes()
n = 0
with open(pool_f) as fp:
    for j, line in enumerate(fp):
        m = mb[j]
        if m:
            s = line.rstrip("\n") + "\tH\n"
            for b in range(len(DIL)):
                if (m >> b) & 1:
                    outs[b].write(s)
        n += 1
assert n == M_pool, (n, M_pool)
for fo in outs:
    fo.close()
pd.DataFrame(info).to_csv(f"{outdir}/{sample}.mixinfo.tsv", sep="\t", index=False)
print("MIXED", sample, info)
