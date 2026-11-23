import argparse
import os
import zlib

import numpy as np


def load_regions(path):
    ivals = sorted(tuple(map(int, line.split("\t")[1:3])) for line in open(path))
    m = np.array(ivals, dtype=np.int64)
    return m[:, 0], m[:, 1]


def in_bucket(k, b):
    return k >= int(b[:-1]) if b.endswith("+") else k == int(b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pat")
    ap.add_argument("regions")
    ap.add_argument("out")
    ap.add_argument("--size", type=int, default=0)
    ap.add_argument("--bucket", choices=["1", "2", "3", "4+", "2+"])
    ap.add_argument("--match", type=int)
    args = ap.parse_args()
    if args.match is not None and args.bucket not in ("2", "3", "4+"):
        ap.error("--match needs --bucket 2, 3 or 4+")

    starts, ends = load_regions(args.regions)
    strings = []
    with open(args.pat) as fh:
        for line in fh:
            f = line.split("\t")
            chrom, first, pattern = f[0], int(f[1]), f[2]
            last = first + len(pattern)
            j = int(np.searchsorted(starts, last - 1, side="right")) - 1
            while j >= 0 and ends[j] > first:
                lo, hi = max(int(starts[j]), first), min(int(ends[j]), last)
                idx = [i for i in range(lo, hi) if pattern[i - first] != "."]
                j -= 1
                if idx:
                    strings.append((chrom, first, pattern, idx))
    if args.match is not None:
        target = min(sum(len(s[3]) for s in strings if in_bucket(len(s[3]), b)) for b in ("2", "3", "4+"))
    if args.bucket:
        strings = [s for s in strings if in_bucket(len(s[3]), args.bucket)]
    if args.match is not None:
        rng = np.random.default_rng([zlib.crc32(os.path.basename(args.pat).encode()), args.match])
        order = rng.permutation(len(strings))
        beads = np.cumsum([len(strings[i][3]) for i in order])
        strings = [strings[i] for i in order[beads <= target]]
    rows = []
    for chrom, first, pattern, idx in strings:
        step = args.size or len(idx)
        for a in range(0, len(idx), step):
            piece = idx[a:a + step]
            rows.append((piece[0], chrom, pattern[piece[0] - first:piece[-1] - first + 1]))
    rows.sort(key=lambda r: r[0])
    with open(args.out, "w") as oh:
        for start, chrom, pat in rows:
            oh.write(f"{chrom}\t{start}\t{pat}\t1\n")


if __name__ == "__main__":
    main()
