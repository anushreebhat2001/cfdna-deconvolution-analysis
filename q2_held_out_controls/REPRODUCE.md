# Reproducing Q2

Everything below ran on NYU HPC Cloud Bursting (`ood.burst.hpc.nyu.edu`) under account
`csci_ga_3033_135_2026fa`. Paths assume the layout in `/scratch/$USER/mrd/`. Total compute is about
7 hours of wall time on 24 cores: 5 h 25 min for setup (b), 80 min for setup (c), and a few minutes
for the scoring steps, which run on the login node.

## 0. Prerequisites

The environment is Asgar's Singularity overlay, mounted read-only. Nothing needs to be installed.

```bash
SIF=/share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif
OVERLAY=/scratch/mxf4642/envs/mrd.ext3:ro
```

These inputs come from Q1 and must exist before anything here runs:

| Path | What it is |
| --- | --- |
| `/scratch/mxf4642/code/project/meta/samples.tsv` | sample sheet, including the `in_healthy_pool` column |
| `/scratch/mxf4642/code/project/meta/read_totals.tsv` | informative read counts per sample |
| `/scratch/mxf4642/results/summary/all_estimates.tsv` | Q1's estimates, needed for the rescore control |
| `/scratch/$USER/mrd/data/pat/` | the 42 project samples in pat format |

The first three live in Asgar's scratch, which is where Q1 produced them. That is fine for re-running
now, but it is a dependency on one person's scratch space and will break when that space is purged —
worth copying somewhere shared before this is archived.

Layout used below:

```
/scratch/$USER/mrd/
  code/     the scripts in code/ and patches/ of this repo
  data/     pat/, pools/, splits.tsv
  results/  b/, c/
  logs/     Slurm output
```

## 1. Generate the splits and pool lists

Deterministic — re-running gives byte-identical output.

```bash
cd /scratch/$USER/mrd/code
python3 make_splits.py \
    --samples  /scratch/mxf4642/code/project/meta/samples.tsv \
    --out      /scratch/$USER/mrd/data/splits.tsv \
    --pool-dir /scratch/$USER/mrd/data/pools
```

Writes `data/splits.tsv` (1,324 rows) and `data/pools/<condition>.txt` (42 lists). Expected
console output: 22 controls (HU012 n=11, HU083 n=11), paper pool of 20, exclusions
`['HU083.63', 'HU083.67']`, 2,200 mixtures planned.

The conditions it defines:

| Prefix | Count | What it is |
| --- | --- | --- |
| `a_rescore_00..19` | 20 | random stratified 11-control subsets, for the scoring-set-size control — no new mixtures |
| `b_split_00..19` | 20 | setup (b); pool is half A, scored against half A and half B |
| `c_batch_fwd`, `c_batch_rev` | 2 | setup (c); pool is one whole batch, scored against the other |

Copies of both outputs are in `meta/` in this repo, so the split assignment can be inspected without
re-running this step.

## 2. Setup (b)

`run_b.sh` is driven by environment variables: `CONDS` (space-separated condition names), `MARKERS`,
`METHODS`, `NPROC` (default 14), `NPROC_MIX` (default 6) and `KEEP_MIX`. For each condition it builds
the healthy pool by concatenating its pool-list members, mixes all 20 CRC samples, sorts, bgzips and
indexes, deconvolves, copies the mixinfo, and deletes the mixtures unless `KEEP_MIX=1`.

```bash
cd /scratch/$USER/mrd/code
sbatch run_b.SBATCH
```

**Do not pass a long `CONDS` through `--export`.** The header comments in `run_b.SBATCH` show that
form, and it fails two ways: a space inside the value gets split by the submitting shell, and a long
environment makes Slurm hold the job with `user env retrieval failed, requeued held` — which
`scontrol release` will not clear. Both `run_b.SBATCH` and `run_c.SBATCH` therefore build `CONDS`
themselves on the compute node from the contents of `data/pools/`. If you need a subset, edit the
`CONDS=` line inside the script rather than exporting it.

Runtime is 5 h 25 min on 24 cores for all 20 splits. Watch it with:

```bash
squeue -u $USER
tail -f /scratch/$USER/mrd/logs/b_<jobid>.out
```

Mixtures are deleted as each condition finishes, so peak disk stays modest. The per-mixture
manifests in `results/b/mix/*.mixinfo.tsv` record the seed, pool membership and read counts, and are
kept. Seeds are derived as `cksum` of `"<condition> <sample>"`, not Python's `hash()`, which is
salted per process and would not reproduce.

Jobs are resumable: `.done` markers mean a resubmitted job skips completed conditions.

## 3. Setup (c)

```bash
cd /scratch/$USER/mrd/code
unset CONDS METHODS MARKERS KEEP_MIX
sbatch run_c.SBATCH
```

`run_c.SBATCH` sets `METHODS="celfieish celfie"`, `MARKERS="max"` and `KEEP_MIX=1`, and builds its
condition list from `c_batch_*.txt`. Two conditions, 200 mixtures, 600 estimates, about 80 minutes.
Mixtures are kept here (roughly 1 GB) because the batch arms are worth spot-checking.

## 4. Scoring

All three steps run on the login node in a few minutes.

**The scoring-set-size control.** Standalone — reads Q1's `all_estimates.tsv` and `splits.tsv`, and
re-scores against the 20 random control subsets. No mixtures involved.

```bash
python rescore_controls.py
```

Expect a maximum median shift of 0.027 AUROC and an identical floor in 20 of 20 subsets on max
markers. It prints both `floor_first` and `floor_sustained` and flags any disagreement.

**Setup (b).** All 20 splits at once:

```bash
python score_b.py --conditions 'b_split_*'
```

**Setup (c): score each direction separately.** Running both at once reports medians across the two
conditions, which averages the asymmetry away and hides the entire result.

```bash
python score_b.py --conditions 'c_batch_fwd'
python score_b.py --conditions 'c_batch_rev'
```

`--threshold` defaults to 0.95 and sets the AUROC level for both floor definitions. Every floor in
the README used that default with the sustained definition.

## 5. What to expect

| Output | Value |
| --- | --- |
| Rescore control, max median shift | ≤ 0.027, floor unmoved 20/20 |
| Setup (b), median paired change, lowest three bins | −0.034, −0.011, −0.008 |
| Setup (b), floor unchanged / worse / better | 14 / 5 / 1 of 20 splits |
| Setup (c) fwd, AUROC lowest bin and floor | 0.871, 1.2–1.7% |
| Setup (c) rev, AUROC lowest bin and floor | 0.583, 1.7–3.0% |

## Known rough edges

- `score_b.py` prints its closing "paired delta" note even for conditions that have no leaky arm,
  such as the two `c_batch_*` runs. The number is absent; only the label is spurious.
- Reconnecting to a Jupyter terminal through Open OnDemand starts a *login* shell, which reads
  `~/.bash_profile` and not `~/.bashrc`, so `$WORK` and similar come back empty. Add
  `[ -f ~/.bashrc ] && . ~/.bashrc` to `~/.bash_profile` once.
- Paste whole command blocks into a plain shell, not into an interactive `singularity ... /bin/bash`
  session; lines pasted while the container shell is still starting are lost. The scripts use the
  `/bin/bash -c '...'` form for this reason.

## The one change to shared code

`patches/mix_dilutions.py` is Katsman et al.'s mixing script with an optional eighth positional
argument: a file listing which controls form the healthy pool. Without it the script recomputes pool
membership from the `in_healthy_pool` column exactly as the original does, so a 7-argument call is
byte-identical in behavior — which is what makes the setup (a) reproduction check in §2 a valid test
of the pipeline. It also records `seed`, `n_pool_samples`, `pool_list` and `pool_samples` into each
`mixinfo.tsv`.
