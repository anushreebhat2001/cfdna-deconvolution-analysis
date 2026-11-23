import sys, numpy as np
atlas, res = sys.argv[1], sys.argv[2]
names = [h[:-5] for h in open(atlas).readline().rstrip("\n").split("\t") if h.endswith("_METH")]
v = np.loadtxt(res).ravel()
v = v / v.sum()
EPI = {"Bladder-Ep","Breast-Basal-Ep","Breast-Luminal-Ep","Colon-Ep","Fallopian-Ep","Gastric-Ep",
       "Head-Neck-Ep","Kidney-Ep","Lung-Ep-Alveo","Lung-Ep-Bron","Ovary-Ep","Prostate-Ep",
       "Small-Int-Ep","Thyroid-Ep","Pancreas-Acinar","Pancreas-Duct","Epithelial"}
print(f"RESULT\t{res.split('/')[-1]}\tepithelial={sum(x for n, x in zip(names, v) if n in EPI):.4f}")
for n, x in sorted(zip(names, v), key=lambda t: -t[1])[:6]:
    print(f"   {n:20s} {x:.4f}")
