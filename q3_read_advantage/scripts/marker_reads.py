import gzip
import sys

import numpy as np


def load_regions(paths):
    ivals = []
    for path in paths:
        with open(path) as fh:
            for line in fh:
                f = line.split("\t")
                if f[0].startswith("#") or not f[1].isdigit():
                    continue
                ivals.append((int(f[1]), int(f[2])))
    ivals.sort()
    merged = [list(ivals[0])]
    for s, e in ivals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    m = np.array(merged, dtype=np.int64)
    return m[:, 0], m[:, 1]


def n_marker_calls(first, pattern, starts, ends):
    j = int(np.searchsorted(starts, first + len(pattern) - 1, side="right")) - 1
    n = 0
    while j >= 0 and ends[j] > first:
        lo, hi = max(int(starts[j]), first), min(int(ends[j]), first + len(pattern))
        seg = pattern[lo - first:hi - first]
        n += len(seg) - seg.count(".")
        j -= 1
    return n


def main():
    pat, out, regions = sys.argv[1], sys.argv[2], sys.argv[3:]
    starts, ends = load_regions(regions)
    total = kept = 0
    with gzip.open(pat, "rt") as fh, open(out, "w") as oh:
        for line in fh:
            total += 1
            f = line.split("\t")
            if n_marker_calls(int(f[1]), f[2], starts, ends) > 0:
                oh.write(line)
                kept += 1
    sys.stderr.write(f"[marker_reads] {pat}\treads_in={total}\tkept={kept}\n")


if __name__ == "__main__":
    main()
