import sys, numpy as np, pandas as pd
tp, ref, out = sys.argv[1:4]
names = list(np.load(ref)["names"])
df = pd.read_csv(tp, sep="\t", index_col=0)
v = df.iloc[0][names].values.astype(float)
np.savetxt(out, (100 * v / v.sum()).reshape(1, -1), fmt="%.4f", delimiter="\t")
