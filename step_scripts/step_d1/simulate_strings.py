import gzip
import os
import sys
import zlib

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from score_pair import ATLAS, EPI

OUT = "results/d1/sim"
TUMOR = "Colon-Ep"


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = [line.rstrip("\n").split("\t") for line in open(ATLAS["u250"])]
    names = [h[:-5] for h in rows[0] if h.endswith("_METH")]
    table = np.array([r[3:] for r in rows[1:]], dtype=float)
    with np.errstate(invalid="ignore"):
        rate = table[:, 0::2] / table[:, 1::2]
    known = ~np.isnan(rate)
    fill = np.where(known, rate, 0).sum(1) / np.maximum(known.sum(1), 1)
    rate = np.where(known, rate, fill[:, None])
    row_of = {int(r[1]): k for k, r in enumerate(rows[1:])}
    epi = np.array([n in EPI["u250"] for n in names])
    tumor = names.index(TUMOR)

    meta = [line.rstrip("\n").split("\t") for line in list(open("meta/samples.tsv"))[1:]]
    healthy = np.mean([np.loadtxt(f"results/step2/{f[0]}/whole.ish_unclipped.txt").ravel()
                       for f in meta if f[1] == "control"], axis=0)
    healthy[epi] = 0
    healthy /= healthy.sum()

    inputs = []
    for f in meta:
        share = float(f[5]) if f[1] == "crc" else 0.0
        inputs.append((f[0], share))
        if f[1] == "crc":
            inputs += [(f"{f[0]}__d{p}", share * float(p) / 100) for p in ("50", "25", "12.5", "6.25", "3.125")]

    for name, share in inputs:
        rng = np.random.default_rng(zlib.crc32(name.encode()))
        mix = healthy * (1 - share)
        mix[tumor] += share
        with gzip.open(f"results/step2/{name}/whole.pat.gz", "rt") as fh, open(f"{OUT}/{name}.pat", "w") as oh:
            for line in fh:
                chrom, first, pattern, count = line.rstrip("\n").split("\t")
                cell = rng.choice(len(mix), p=mix)
                new = [c if c == "." else ("C" if rng.random() < rate[row_of[int(first) + i], cell] else "T")
                       for i, c in enumerate(pattern)]
                oh.write(f"{chrom}\t{first}\t{''.join(new)}\t{count}\n")
    print(f"wrote {len(inputs)} simulated files to {OUT}")


if __name__ == "__main__":
    main()
