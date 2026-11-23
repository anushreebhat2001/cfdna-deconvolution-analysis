import gzip
import os
import subprocess
import sys
import zlib

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
from score_pair import ATLAS, EPI

OUT = "results/factors"
BREAST = {"HU017.04.PL4596": ("HU017", 974142, 0.01909), "HU017and19.05.PL4434": ("HU017", 3642454, 0.0331),
          "HU017and19.06.PL4518": ("HU017", 3493079, 0.02067), "HU025.13.PL4704": ("HU025", 3214388, 0.0)}
CONDITIONS = {"none": ("none", 0), "clumped_1": ("clumped", 0.01), "clumped_3": ("clumped", 0.03),
              "scattered_1": ("scattered", 0.01), "scattered_3": ("scattered", 0.03),
              "symmetric_5": ("symmetric", 0.05), "symmetric_10": ("symmetric", 0.10)}


def atlases():
    os.makedirs(f"{OUT}/atlas", exist_ok=True)
    rows = [line.rstrip("\n") for line in open(ATLAS["u250"])]
    for tag, name, meth in (("off", "Off-everywhere", 5), ("on", "On-everywhere", 95)):
        with open(f"{OUT}/atlas/u250_plus_{tag}.txt", "w") as fh:
            fh.write(f"{rows[0]}\t{name}_METH\t{name}_COV\n")
            fh.writelines(f"{r}\t{meth}\t100\n" for r in rows[1:])
    rows = [line.rstrip("\n").split("\t") for line in open(ATLAS["max"])]
    assert rows[0][-2:] == ["Megakaryocytes_METH", "Megakaryocytes_COV"]
    with open(f"{OUT}/atlas/max_minus_megak.txt", "w") as fh:
        fh.writelines("\t".join(r[:-2]) + "\n" for r in rows)


def breast_dilutions():
    lines = list(open("meta/samples.tsv"))
    meta = [lines[0]] + [x for x in lines[1:] if x.split("\t")[1] == "control"]
    meta += [f"{s}\tcrc\t{s}\t{b}\t{n}\t{tf}\tno\n" for s, (b, n, tf) in BREAST.items()]
    with open(f"{OUT}/meta_breast.tsv", "w") as fh:
        fh.writelines(meta)
    subprocess.run([sys.executable, "-I", "scripts/mix_dilutions.py", "--meta", f"{OUT}/meta_breast.tsv", "--totals",
                    "logs/convert_summary.tsv", "--marker-dir", "data/marker_pat", "--out", "data/mix_breast"], check=True)


def simulations():
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
    meta = [line.rstrip("\n").split("\t") for line in list(open("meta/samples.tsv"))[1:]]
    controls = [f[0] for f in meta if f[1] == "control"]
    healthy = np.mean([np.loadtxt(f"results/step2/{c}/whole.ish_unclipped.txt").ravel() for c in controls], axis=0)
    healthy[epi] = 0
    healthy /= healthy.sum()

    for cond in CONDITIONS:
        os.makedirs(f"{OUT}/sim/{cond}", exist_ok=True)
    for name in controls:
        rng = np.random.default_rng(zlib.crc32(name.encode()))
        knob = np.random.default_rng(zlib.crc32((name + "|switch").encode()))
        out = {cond: open(f"{OUT}/sim/{cond}/{name}.pat", "w") for cond in CONDITIONS}
        with gzip.open(f"results/step2/{name}/whole.pat.gz", "rt") as fh:
            for line in fh:
                chrom, first, pattern, count = line.rstrip("\n").split("\t")
                cell = rng.choice(len(healthy), p=healthy)
                base = [c if c == "." else ("C" if rng.random() < rate[row_of[int(first) + i], cell] else "T")
                        for i, c in enumerate(pattern)]
                u_string = knob.random()
                u_bead = [None if c == "." else knob.random() for c in base]
                for cond, (kind, p) in CONDITIONS.items():
                    if kind == "clumped" and u_string < p:
                        new = ["." if c == "." else "T" for c in base]
                    elif kind == "scattered":
                        new = ["T" if c == "C" and u < p else c for c, u in zip(base, u_bead)]
                    elif kind == "symmetric":
                        new = [c if c == "." or u >= p else ("T" if c == "C" else "C") for c, u in zip(base, u_bead)]
                    else:
                        new = base
                    out[cond].write(f"{chrom}\t{first}\t{''.join(new)}\t{count}\n")
        for fh in out.values():
            fh.close()
    print(f"wrote {len(controls)} simulated controls for each of {len(CONDITIONS)} conditions")


if __name__ == "__main__":
    atlases()
    breast_dilutions()
    simulations()
