import argparse
import gzip
import re
import sys

import numpy as np

TOKEN = re.compile(r"([A-Za-z])(\d*)")
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def calls(start, rle):
    pos = start
    for ch, num in TOKEN.findall(rle):
        n = int(num) if num else 1
        if ch == "M" or ch == "U":
            state = "C" if ch == "M" else "T"
            for k in range(n):
                yield pos + k, state
        if ch != "i":
            pos += n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("epibed")
    ap.add_argument("index")
    ap.add_argument("out")
    args = ap.parse_args()

    idx = np.load(args.index)
    cpos = {c: idx[f"pos_{c}"].astype(np.int64) for c in CHROMS}
    coff = {c: int(idx[f"off_{c}"]) for c in CHROMS}

    total = kept = n_calls = off_target = other_chrom = 0
    rows = []
    with gzip.open(args.epibed, "rt") as fh:
        for line in fh:
            total += 1
            f = line.split("\t", 7)
            chrom = f[0]
            if chrom not in cpos:
                other_chrom += 1
                continue
            pc = list(calls(int(f[1]), f[6]))
            if not pc:
                continue
            arr = cpos[chrom]
            p = np.fromiter((x[0] for x in pc), dtype=np.int64, count=len(pc))
            j = np.searchsorted(arr, p, side="right") - 1
            d = p - arr[np.maximum(j, 0)]
            ok = (j >= 0) & ((d == 0) | (d == 1))
            n_calls += len(pc)
            off_target += int((~ok).sum())
            if not ok.any():
                continue
            jj = j[ok]
            states = [pc[i][1] for i in np.nonzero(ok)[0]]
            first, last = int(jj[0]), int(jj[-1])
            pattern = ["."] * (last - first + 1)
            for k, s in zip(jj, states):
                pattern[int(k) - first] = s
            rows.append((coff[chrom] + first + 1, chrom, "".join(pattern)))
            kept += 1

    rows.sort(key=lambda r: r[0])
    with open(args.out, "w") as out:
        for gi, chrom, pattern in rows:
            out.write(f"{chrom}\t{gi}\t{pattern}\t1\n")
    sys.stderr.write(
        f"[epibed2pat] {args.epibed}\ttotal_reads={total}\treads_with_cpg={kept}\t"
        f"calls={n_calls}\tcalls_off_cpg={off_target}\tother_chrom_reads={other_chrom}\n")


if __name__ == "__main__":
    main()
