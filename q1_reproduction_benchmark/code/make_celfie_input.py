"""Input file for the original celfie.py, one sample: sample per-CpG counts + atlas columns.
Rows with no sample coverage are dropped (they carry no information about the sample's proportions).
usage: make_celfie_input.py SITES.tsv ATLAS SAMPLE_ID OUT.txt"""
import sys, pandas as pd
sites, atlas, sid, out = sys.argv[1:5]
S = pd.read_csv(sites, sep="\t")
A = pd.read_csv(atlas, sep="\t")
assert len(S) == len(A) and (S.START.values == A.START.values).all(), "sites and atlas rows differ"
keep = S.DEPTH.values > 0
cols = [c for c in A.columns if c.endswith("_METH") or c.endswith("_COV")]
left = S.loc[keep, ["CHROM", "START", "END", "METH", "DEPTH"]].reset_index(drop=True)
left.columns = ["CHROM", "START", "END", f"{sid}_METH", f"{sid}_DEPTH"]
right = A.loc[keep, ["CHROM", "START", "END"] + cols].reset_index(drop=True)
pd.concat([left, right], axis=1).to_csv(out, sep="\t", index=False)
print(f"{sid}: rows kept {int(keep.sum())} / {len(S)}")
