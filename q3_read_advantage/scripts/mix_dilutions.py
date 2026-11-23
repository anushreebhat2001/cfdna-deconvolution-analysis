import argparse
import os
import re
import zlib

import numpy as np

DILUTIONS = [50, 25, 12.5, 6.25, 3.125]


def read_lines(path):
    with open(path) as fh:
        return fh.readlines()


def write_sorted(path, lines):
    lines.sort(key=lambda l: int(l.split("\t", 2)[1]))
    with open(path, "w") as fh:
        fh.writelines(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--totals", required=True)
    ap.add_argument("--marker-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    total = {}
    for line in open(args.totals):
        s = re.search(r"epibed/(.+)\.epibed\.gz", line).group(1)
        total[s] = int(re.search(r"total_reads=(\d+)", line).group(1))

    crc, pool = [], []
    for line in list(open(args.meta))[1:]:
        f = line.rstrip("\n").split("\t")
        if f[1] == "crc":
            crc.append(f[0])
        elif f[1] == "control" and f[6] == "yes" and args.batch in (None, f[3]):
            pool.append(f[0])

    pool_lines = []
    for s in pool:
        pool_lines += read_lines(f"{args.marker_dir}/{s}.pat")
    pool_total = sum(total[s] for s in pool)
    print(f"healthy pool: {len(pool)} controls, {pool_total} reads, {len(pool_lines)} marker reads")

    log = ["sample\tdilution_pct\treads_total\tcancer_reads\thealthy_reads\tcancer_marker\thealthy_marker"]
    for s in crc:
        rng = np.random.default_rng(zlib.crc32(s.encode()))
        own = read_lines(f"{args.marker_dir}/{s}.pat")
        n = total[s]
        k0 = int(rng.hypergeometric(len(pool_lines), pool_total - len(pool_lines), n))
        mix0 = [pool_lines[i] for i in rng.choice(len(pool_lines), size=k0, replace=False)]
        for pct in DILUTIONS:
            m_c = int(pct / 100 * n)
            m_h = int((1 - pct / 100) * n)
            kc = int(rng.hypergeometric(len(own), n - len(own), m_c))
            kh = int(rng.hypergeometric(len(mix0), n - len(mix0), m_h))
            part = [own[i] for i in rng.choice(len(own), size=kc, replace=False)]
            part += [mix0[i] for i in rng.choice(len(mix0), size=kh, replace=False)]
            write_sorted(f"{args.out}/{s}__d{pct:g}.pat", part)
            log.append(f"{s}\t{pct:g}\t{n}\t{m_c}\t{m_h}\t{kc}\t{kh}")
    with open(f"{args.out}/mix_log.tsv", "w") as fh:
        fh.write("\n".join(log) + "\n")
    print(f"wrote {len(crc) * len(DILUTIONS)} dilutions for {len(crc)} cancer samples")


if __name__ == "__main__":
    main()
