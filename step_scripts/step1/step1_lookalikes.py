import argparse
import collections
import os

import numpy as np

BUCKETS = ["1", "2", "3", "4+"]
HIGH_TUMOR = 0.4
N_BOOT = 2000


def load_places(path):
    ivals = sorted(tuple(map(int, line.split("\t")[1:3])) for line in open(path))
    merged = [list(ivals[0])]
    for s, e in ivals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    m = np.array(merged, dtype=np.int64)
    return m[:, 0], m[:, 1]


def strings(pat_path, starts, ends):
    with open(pat_path) as fh:
        for line in fh:
            f = line.split("\t")
            first, pattern = int(f[1]), f[2]
            last = first + len(pattern)
            j = int(np.searchsorted(starts, last - 1, side="right")) - 1
            while j >= 0 and ends[j] > first:
                lo, hi = max(int(starts[j]), first), min(int(ends[j]), last)
                idx = [i for i in range(lo, hi) if pattern[i - first] != "."]
                if idx:
                    yield j, idx, [pattern[i - first] == "T" for i in idx]
                j -= 1


def bucket(k):
    return "4+" if k >= 4 else str(k)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--marker-dir", required=True)
    ap.add_argument("--places", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    controls, cancers = [], []
    for line in list(open(args.meta))[1:]:
        f = line.rstrip("\n").split("\t")
        if f[1] == "control" and args.batch in (None, f[3]):
            controls.append(f[0])
        elif f[1] == "crc" and float(f[5]) >= HIGH_TUMOR:
            cancers.append(f[0])

    starts, ends = load_places(args.places)
    data = {s: list(strings(f"{args.marker_dir}/{s}.pat", starts, ends)) for s in controls + cancers}

    n_cpg, u_cpg = collections.Counter(), collections.Counter()
    n_place, u_place = collections.Counter(), collections.Counter()
    for s in controls:
        for j, idx, off in data[s]:
            for i, o in zip(idx, off):
                n_cpg[i] += 1
                u_cpg[i] += o
            n_place[j] += len(idx)
            u_place[j] += sum(off)

    per = collections.defaultdict(lambda: [0, 0, 0.0, 0.0])
    skipped = collections.Counter()
    for group, samples, own in (("healthy", controls, 1), ("cancer_high_tumor", cancers, 0)):
        for s in samples:
            for j, idx, off in data[s]:
                k = len(idx)
                np_, up_ = n_place[j] - own * k, u_place[j] - own * sum(off)
                if np_ <= 0:
                    skipped[group] += 1
                    continue
                rate_place = up_ / np_
                p = 1.0
                for i, o in zip(idx, off):
                    n, u = n_cpg[i] - own, u_cpg[i] - own * o
                    p *= u / n if n > 0 else rate_place
                row = per[(group, s, bucket(k))]
                row[0] += 1
                row[1] += all(off)
                row[2] += p
                row[3] += rate_place ** k

    rng = np.random.default_rng(0)
    pooled = ["group\tbeads\tstrings\tall_off_seen\tshare_seen\tall_off_expected\tratio\t"
              "ratio_ci_lo\tratio_ci_hi\tall_off_expected_marker_rate\tratio_marker_rate"]
    for group, samples in (("healthy", controls), ("cancer_high_tumor", cancers)):
        for b in BUCKETS:
            t = np.array([per[(group, s, b)] for s in samples], dtype=float)
            n, seen, exp, exp_m = t.sum(0)
            boot = []
            for _ in range(N_BOOT):
                r = t[rng.integers(0, len(samples), len(samples))].sum(0)
                if r[2] > 0:
                    boot.append(r[1] / r[2])
            lo, hi = np.percentile(boot, [2.5, 97.5])
            pooled.append(f"{group}\t{b}\t{int(n)}\t{int(seen)}\t{seen / n:.4f}\t{exp:.1f}\t"
                          f"{seen / exp:.2f}\t{lo:.2f}\t{hi:.2f}\t{exp_m:.1f}\t{seen / exp_m:.2f}")

    lines = ["group\tsample\tbeads\tstrings\tall_off_seen\tall_off_expected\tall_off_expected_marker_rate"]
    for group, samples in (("healthy", controls), ("cancer_high_tumor", cancers)):
        for s in samples:
            for b in BUCKETS:
                n, seen, exp, exp_m = per[(group, s, b)]
                lines.append(f"{group}\t{s}\t{b}\t{n}\t{seen}\t{exp:.2f}\t{exp_m:.2f}")

    with open(f"{args.out}/{args.label}_pooled.tsv", "w") as fh:
        fh.write("\n".join(pooled) + "\n")
    with open(f"{args.out}/{args.label}_per_sample.tsv", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"[{args.label}] {len(controls)} healthy controls, {len(cancers)} cancer samples with tumor share >= {HIGH_TUMOR}")
    print(f"[{args.label}] strings that could not be predicted (no other healthy call in the place): {dict(skipped)}")
    print("\n".join(pooled))


if __name__ == "__main__":
    main()
