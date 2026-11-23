import collections
import filecmp
import os
import sys

import numpy as np

EPI = {"Bladder-Ep", "Breast-Basal-Ep", "Breast-Luminal-Ep", "Colon-Ep", "Fallopian-Ep",
       "Gastric-Ep", "Head-Neck-Ep", "Kidney-Ep", "Lung-Ep-Alveo", "Lung-Ep-Bron", "Ovary-Ep",
       "Prostate-Ep", "Small-Int-Ep", "Thyroid-Ep", "Pancreas-Acinar", "Pancreas-Duct"}


def calls(path):
    c, lines, beads = collections.Counter(), 0, 0
    for line in open(path):
        f = line.split("\t")
        lines += 1
        for k, ch in enumerate(f[2]):
            if ch != ".":
                c[(int(f[1]) + k, ch)] += 1
                beads += 1
    return c, lines, beads


def load(path):
    v = np.loadtxt(path).ravel()
    return v / v.sum() * 100


def main():
    root, atlas = sys.argv[1], sys.argv[2]
    header = open(atlas).readline().rstrip("\n").split("\t")
    names = [h[:-5] for h in header if h.endswith("_METH")]
    epi = np.array([n in EPI for n in names])

    for s in sorted(os.listdir(root)):
        d = f"{root}/{s}"
        if not os.path.isdir(d):
            continue
        print(f"\n== {s}")
        cw, lw, bw = calls(f"{d}/whole.pat")
        cc, lc, bc = calls(f"{d}/cut.pat")
        print(f"A same beads     : {'PASS' if cw == cc and lc == bc == bw else 'FAIL'}"
              f"   ({lw} strings, {bw} beads; cut file has {lc} lines)")
        same = filecmp.cmp(f"{d}/whole_1bead.pat", f"{d}/cut_1bead.pat", shallow=False)
        n1 = sum(1 for _ in open(f"{d}/whole_1bead.pat"))
        print(f"B 1-bead bucket  : {'PASS' if same else 'FAIL'}   ({n1} one-bead strings, files identical: {same})")

        for tag, a, b in (("C strings (ISH unclipped)", "orig.ish_unclipped", "whole.ish_unclipped"),
                          ("C strings (ISH clipped)  ", "orig.ish_clipped", "whole.ish_clipped"),
                          ("D CelFiE whole vs cut    ", "whole.celfie", "cut.celfie")):
            diff = float(np.abs(load(f"{d}/{a}.txt") - load(f"{d}/{b}.txt")).max())
            print(f"{tag}: {'PASS' if diff < 0.01 else 'FAIL'}   (largest difference {diff:.4f} percentage points)")

        print("epithelial share (%):")
        for lab, f in (("  CelFiE-ISH unclipped, whole", "whole.ish_unclipped"),
                       ("  CelFiE-ISH unclipped, cut  ", "cut.ish_unclipped"),
                       ("  CelFiE-ISH clipped, whole  ", "whole.ish_clipped"),
                       ("  CelFiE-ISH clipped, cut    ", "cut.ish_clipped"),
                       ("  CelFiE, whole              ", "whole.celfie"),
                       ("  CelFiE, cut                ", "cut.celfie")):
            print(f"{lab}: {load(f'{d}/{f}.txt')[epi].sum():6.2f}")


if __name__ == "__main__":
    main()
