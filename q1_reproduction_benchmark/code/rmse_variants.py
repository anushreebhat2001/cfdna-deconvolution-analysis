import numpy as np, pandas as pd
est = pd.read_csv("/scratch/mxf4642/results/m5/m5_estimates.tsv", sep="\t")
crc = est[est.group == "crc"]
rows = []
for (m, ms), g in crc.groupby(["method", "markers"]):
    r = {k: [] for k in "ABCDE"}
    for s, gs in g.groupby("sample"):
        base = gs.loc[gs.dilution == 1.0, "epi"]
        if len(base) == 0 or base.iloc[0] <= 0:
            continue
        b, d, a = base.iloc[0], gs[gs.dilution < 1], gs
        norm = d.epi / b
        r["A"].append(np.sqrt(np.mean((norm - d.dilution) ** 2)))
        r["B"].append(np.sqrt(np.sum((norm - d.dilution) ** 2) / (len(d) + 1)))
        r["C"].append(np.sqrt(np.mean((d.epi - d.dilution * b) ** 2)))
        r["D"].append(np.sqrt(np.mean((a.epi - a.ichorCNA_TF * a.dilution) ** 2)))
        ok = norm > 0
        r["E"].append(np.sqrt(np.mean((np.log2(norm[ok]) - np.log2(d.dilution[ok])) ** 2)))
    rows.append(dict(method=m, markers=ms, **{k: np.median(v) for k, v in r.items()}))
out = pd.DataFrame(rows)
print(out.round(3).to_string(index=False))
print("\nmedian over methods, per marker set:")
print(out.groupby("markers")[list("ABCDE")].median().round(3))
print("\npaper: u250 median 0.312 | collapsed median 0.078 | best (celfieish-clip, max) 0.034")
