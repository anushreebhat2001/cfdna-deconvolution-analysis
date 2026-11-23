#!/usr/bin/env python3
"""
Generate control splits for setups (b) and (c), plus the free re-scoring
sets for setup (a).

Writes:
  splits.tsv              one row per (condition, role, sample)
  pools/<condition>.txt   the pool sample names for that condition, one per
                          line, sorted. Feed each to the patched
                          mix_dilutions.py as its 8th argument, and build the
                          matching healthy_pool.pat by cat-ing the same
                          samples in the same order.

Conditions
----------
  a_rescore_<k>   pool = the paper's panel of normals (in_healthy_pool==yes).
                  score = a random stratified half of all controls.
                  NO NEW MIXTURES: re-scores the estimates already in
                  results/m5. Isolates the effect of halving the scoring set.

  b_split_<k>     pool = stratified half A. The same mixtures are scored twice:
                      score_leaky    = half A  (full overlap)
                      score_heldout  = half B  (no overlap)
                  The paired difference is leakage with pool size, pool
                  composition, scoring-set size and batch balance held fixed.

  c_batch_fwd     pool = batch 1, score_heldout = batch 2
  c_batch_rev     pool = batch 2, score_heldout = batch 1

Usage
-----
  python3 make_splits.py \
      --samples /scratch/mxf4642/code/project/meta/samples.tsv \
      --out    /scratch/$USER/mrd/data/splits.tsv \
      --pool-dir /scratch/$USER/mrd/data/pools
"""

import argparse
import csv
import os
import random
import sys
from collections import Counter, defaultdict


def load(path, sample_col, group_col, batch_col, pool_col, control_value):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"{path}: no rows")

    cols = set(rows[0])
    for col in (sample_col, group_col, batch_col, pool_col):
        if col not in cols:
            sys.exit(f"{path}: no column {col!r}. Columns: {sorted(cols)}")

    groups = sorted({r[group_col].strip() for r in rows})
    controls = [(r[sample_col].strip(), r[batch_col].strip())
                for r in rows if r[group_col].strip() == control_value]
    if not controls:
        sys.exit(f"{path}: no rows with {group_col} == {control_value!r}. "
                 f"Values present: {groups}. Pass --control-value.")

    paper_pool = sorted(r[sample_col].strip() for r in rows
                        if r[pool_col].strip().lower() == "yes")
    return controls, paper_pool, groups


def stratified_halves(controls, rng, flip):
    """Split into two halves, balanced within each batch.

    Odd-sized batches cannot split evenly; `flip` alternates which half gets
    the extra member so it evens out across splits.
    """
    by_batch = defaultdict(list)
    for s, b in controls:
        by_batch[b].append(s)

    a, b_ = [], []
    for i, batch in enumerate(sorted(by_batch)):
        members = sorted(by_batch[batch])
        rng.shuffle(members)
        n = len(members)
        take = (n + 1) // 2 if (flip + i) % 2 == 0 else n // 2
        a.extend(members[:take])
        b_.extend(members[take:])
    return sorted(a), sorted(b_)


def stratified_subset(controls, k, rng):
    """Draw k controls, spread across batches as evenly as possible."""
    by_batch = defaultdict(list)
    for s, b in controls:
        by_batch[b].append(s)
    batches = sorted(by_batch)
    for b in batches:
        rng.shuffle(by_batch[b])

    picked, i = [], 0
    while len(picked) < k and any(by_batch[x] for x in batches):
        b = batches[i % len(batches)]
        if by_batch[b]:
            picked.append(by_batch[b].pop())
        i += 1
    return sorted(picked)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--samples", required=True)
    p.add_argument("--out", default="splits.tsv")
    p.add_argument("--pool-dir", default=None,
                   help="directory for per-condition pool lists "
                        "(default: 'pools' beside --out)")
    p.add_argument("--n-splits", type=int, default=20)
    p.add_argument("--n-rescore", type=int, default=20)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--sample-col", default="sample")
    p.add_argument("--group-col", default="group")
    p.add_argument("--batch-col", default="batch")
    p.add_argument("--pool-col", default="in_healthy_pool")
    p.add_argument("--control-value", default="control")
    args = p.parse_args()

    pool_dir = args.pool_dir or os.path.join(
        os.path.dirname(os.path.abspath(args.out)) or ".", "pools")
    os.makedirs(pool_dir, exist_ok=True)

    controls, paper_pool, groups = load(
        args.samples, args.sample_col, args.group_col, args.batch_col,
        args.pool_col, args.control_value)
    rng = random.Random(args.seed)

    batch_counts = Counter(b for _, b in controls)
    print(f"groups in sheet: {groups}")
    print(f"{len(controls)} controls: "
          + ", ".join(f"{b} n={n}" for b, n in sorted(batch_counts.items())))
    print(f"paper pool ({args.pool_col}==yes): {len(paper_pool)} samples")
    not_in_pool = sorted({s for s, _ in controls} - set(paper_pool))
    print(f"  controls excluded from the paper pool: {not_in_pool}")

    rows = []
    pools = {}

    def emit(cond, role, samples, note=""):
        for s in samples:
            rows.append({"condition": cond, "role": role,
                         "sample": s, "note": note})
        if role == "pool":
            pools[cond] = sorted(samples)

    half_n = len(controls) // 2

    # setup (a) re-scoring: reuses existing mixtures, costs no compute
    for k in range(args.n_rescore):
        c = f"a_rescore_{k:02d}"
        emit(c, "pool", paper_pool, "paper pool; mixtures already exist")
        emit(c, "score_heldout", stratified_subset(controls, half_n, rng),
             f"random stratified {half_n} of {len(controls)}")

    # setup (b): one mixture set per split, scored two ways
    for k in range(args.n_splits):
        c = f"b_split_{k:02d}"
        ha, hb = stratified_halves(controls, rng, flip=k)
        assert not set(ha) & set(hb)
        emit(c, "pool", ha, "build mixtures from these")
        emit(c, "score_leaky", ha, "same controls as pool")
        emit(c, "score_heldout", hb, "disjoint from pool")

    # setup (c): batch swap
    by_batch = defaultdict(list)
    for s, b in controls:
        by_batch[b].append(s)
    batches = sorted(by_batch)
    if len(batches) == 2:
        x, y = batches
        emit("c_batch_fwd", "pool", sorted(by_batch[x]), f"pool = {x}")
        emit("c_batch_fwd", "score_heldout", sorted(by_batch[y]), f"score = {y}")
        emit("c_batch_rev", "pool", sorted(by_batch[y]), f"pool = {y}")
        emit("c_batch_rev", "score_heldout", sorted(by_batch[x]), f"score = {x}")
    else:
        print(f"  warning: expected 2 batches, found {batches}; skipping (c)")

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["condition", "role", "sample", "note"])
        w.writeheader()
        w.writerows(rows)

    for cond, members in pools.items():
        with open(os.path.join(pool_dir, f"{cond}.txt"), "w") as fh:
            fh.write("\n".join(members) + "\n")

    # Conditions needing new mixtures = everything except the (a) re-scores,
    # deduplicated by pool membership.
    needs_mix = {cond: tuple(m) for cond, m in pools.items()
                 if not cond.startswith("a_rescore")}
    distinct = len(set(needs_mix.values()))

    for c in needs_mix:
        pool = set(pools[c])
        held = {r["sample"] for r in rows
                if r["condition"] == c and r["role"] == "score_heldout"}
        assert not pool & held, f"{c}: pool and held-out scoring overlap"

    print(f"\nWrote {args.out} ({len(rows)} rows) and "
          f"{len(pools)} pool lists in {pool_dir}")
    print(f"Conditions needing new mixtures: {len(needs_mix)} "
          f"({distinct} distinct pools) -> {len(needs_mix) * 100:,} mixtures")
    print(f"Setup (a) re-scoring: 0 new mixtures, 0 deconvolution runs")
    print("Checked: no pool/held-out overlap in any (b) or (c) condition")


if __name__ == "__main__":
    main()
