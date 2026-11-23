import sys

import numpy as np

EPI = {"Bladder-Ep", "Breast-Basal-Ep", "Breast-Luminal-Ep", "Colon-Ep", "Fallopian-Ep",
       "Gastric-Ep", "Head-Neck-Ep", "Kidney-Ep", "Lung-Ep-Alveo", "Lung-Ep-Bron", "Ovary-Ep",
       "Prostate-Ep", "Small-Int-Ep", "Thyroid-Ep", "Pancreas-Acinar", "Pancreas-Duct"}


def main():
    table, regions, atlas, out = sys.argv[1:5]

    n = 0
    with open(table) as fh, open(f"{out}/epithelial_u250.bed", "w") as oh:
        for line in fh:
            f = line.split("\t")
            if line.startswith("#") or f[0] in ("chrX", "chrY") or f[5] not in EPI:
                continue
            oh.write(f"{f[0]}\t{f[3]}\t{f[4]}\n")
            n += 1
    print(f"U250: {n} epithelial markers")

    header = open(atlas).readline().rstrip("\n").split("\t")
    names = [h[:-5] for h in header if h.endswith("_METH")]
    d = np.loadtxt(atlas, skiprows=1, usecols=range(1, len(header)))
    d = d[np.argsort(d[:, 0])]
    start, meth, cov = d[:, 0].astype(int), d[:, 2::2], d[:, 3::2]
    epi = names.index("Epithelial")
    n = total = 0
    with open(regions) as fh, open(f"{out}/epithelial_panepi.bed", "w") as oh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            lo, hi = np.searchsorted(start, [int(f[1]), int(f[2])])
            beta = meth[lo:hi].sum(0) / np.maximum(cov[lo:hi].sum(0), 1)
            total += 1
            if int(np.argmin(beta)) == epi:
                oh.write(f"{f[0]}\t{f[1]}\t{f[2]}\n")
                n += 1
    print(f"pan-epithelial: {n} of {total} regions have Epithelial as the least methylated group")


if __name__ == "__main__":
    main()
