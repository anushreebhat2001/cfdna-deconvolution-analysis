# Held-out controls and batch swap (Q2)

Katsman et al. (2026) build every dilution's healthy background from 20 controls and then score the
mixtures against all 22, so 20 of the "negatives" also contributed reads to the "positives". We
measured what that overlap is worth, and then asked a second question the overlap test raised: how
much does the detection floor depend on *which* controls you score against at all?

Two answers. Removing the overlap costs a median of **0.034 AUROC** in the lowest tumor-fraction bin
and less above it, and leaves the detection floor unchanged in 14 of 20 matched splits — Katsman's
~1.7% floor is not materially inflated by the overlap. But swapping which control batch the mixtures
are scored against moves the same measure by **0.288** and shifts the floor a full bin, roughly eight
times the leakage effect and ten times the scoring-set-size artifact. On this evidence the reported
floor depends more on control batch composition than on the pool/scoring overlap in Katsman's design.

Floors throughout are the **sustained** floor: the lowest tumor-fraction bin from which AUROC stays at
or above 0.95. First-crossing floors are reported as secondary by the scoring script.

## Setup (b): held-out controls

Comparing setup (a) against a held-out version would change four things at once — pool composition,
pool size (20 to 11), scoring-set size (22 to 11) and overlap. Only the last is leakage. Instead,
each of 20 splits of the 22 controls (stratified so each half holds roughly equal numbers from HU012
and HU083) builds one set of mixtures from half A and scores it twice: against half A (leaky) and
against half B (held out). Pool size, pool composition, scoring-set size and the mixtures themselves
are identical between the two arms, so the paired difference is leakage and nothing else. It also
costs no extra compute, because it is one set of mixtures scored two ways.

100 mixtures per split, 2,000 in total, run as one Slurm job in 5 h 25 min on 24 cores. CelFiE-ISH on
the max-marker set, clipped and unclipped.

**Two checks before the sweep.** Re-running one setup (a) mixture through the patched pipeline
produced a byte-identical estimate file (epithelial 0.1197), confirming the pool-list argument does
not change behavior when it is absent. And re-scoring setup (a)'s existing estimates against 20
random stratified subsets of 11 controls shifted the median AUROC by at most **0.027** across all
eight method/marker combinations, with the floor unmoved in every one — so an effect in setup (b)
has to exceed roughly 0.03 to mean anything.

| Tumor fraction bin | Leaky | Held out | Median paired change | Worse / better / tied | Wilcoxon p |
| --- | --- | --- | --- | --- | --- |
| 0.0–0.5% | 0.799 | 0.746 | −0.034 | 11 / 9 / 0 | 0.083 |
| 0.5–1.1% | 0.955 | 0.939 | −0.011 | 12 / 6 / 2 | 0.073 |
| 1.2–1.7% | 0.970 | 0.966 | −0.008 | 12 / 6 / 2 | 0.112 |
| 1.7–3.0% | 1.000 | 1.000 | 0.000 | 3 / 1 / 16 | 0.438 |
| above 3.1% | 1.000 | 1.000 | 0.000 | 0 / 0 / 20 | — |

CelFiE-ISH with clipping, max markers. Unclipped behaves the same (−0.038, −0.015 and 0.000 in the
three lowest bins). Removing the overlap makes AUROC worse more often than better in every bin below
1.7%, in both clipping modes, but the magnitude is small, comparable to the scoring-set-size
artifact, and no bin reaches p < 0.05. Above 1.7% AUROC saturates at 1.000, so the experiment has no
resolution there. Compared pairwise within split, the floor is unchanged in 14 of 20 splits, worse in
5 and better in 1.

**Two caveats worth stating.** The 20 splits re-divide the same 22 people, so they are not 20
independent experiments; the Wilcoxon test is anti-conservative and the p-values should be read as
descriptive. The effect size is the output. Separately, an early 2-split pilot gave 0.114 — more than
three times the final answer. Two splits of a measurement this variable will tell you almost
anything, which is why the full run is twenty.

## Setup (c): batch swap

The 22 controls split into two sequencing batches of 11, HU012 and HU083, which differ in basecaller
pipeline, sequencing depth, and age and sex composition. Setup (c) asks whether a floor measured
against one batch holds when the comparison controls come from the other. Two conditions, 100
mixtures each, 600 estimates, about 80 minutes of wall time. CelFiE-ISH and CelFiE on max markers.

| Direction | AUROC 0–0.5% | 0.5–1.1% | 1.2–1.7% | Sustained floor |
| --- | --- | --- | --- | --- |
| Mix HU012 → score vs HU083 | 0.871 | 0.947 | 0.992 | 1.2–1.7% |
| Mix HU083 → score vs HU012 | 0.583 | 0.947 | 0.848 | 1.7–3.0% |

CelFiE-ISH, clipped, max markers. A gap of 0.288 in the lowest bin and 0.144 at 1.2–1.7%, from the
same 100 mixtures per direction scored against different control sets, with the floor differing by
one bin. The direction of the asymmetry was recorded as a prediction before the run.

**Mechanism.** The false epithelial baseline in healthy controls is 0.0147 in HU012 (sd 0.0051)
against 0.0117 in HU083 (sd 0.0040) — a location shift of about 0.003 with comparable spread,
Cohen's d ≈ 0.65. At a 1% detection floor the tumor signal being hunted is about 0.01, so a 0.003
shift in the baseline consumes roughly a third of it. Score against the batch with the lower baseline
and the mixtures stand out; score against the higher one and they are partly swallowed. Below roughly
1% tumor fraction the floor is therefore set by the control batch's baseline rather than by the
method's sensitivity.

**Depth is not the primary driver, but is not fully separated.** HU083 carries 45,197,910 informative
reads against HU012's 24,524,327, a factor of 1.84. On max markers — the set setup (c) ran on — depth
does not predict the control estimates (Spearman −0.21, p = 0.347), while on U250 it does strongly
(−0.76, −0.70 and −0.50, all p < 0.02). The batch medians still differ by 26% on max markers. A
downsample-and-rerun of the deeper batch, about 45 CPU-minutes, would settle it and has not been run.

**Incidental: CelFiE returns exact zeros.** On max markers CelFiE's control estimates are 0.0 at both
the median and the 25th percentile in both batches. It is not estimating a small epithelial fraction;
it is estimating none. That degeneracy is why its AUROC is erratic across bins and why the two floor
definitions disagree for it in four separate bins. Its apparent batch-robustness is an artifact of a
collapsed output distribution, not a strength, and should be read that way in the method comparison.

## How the two setups relate

Setup (b)'s held-out arm (0.746 in the lowest bin) falls between setup (c)'s two directions, fwd
(0.871) and rev (0.583). Stratifying splits so that each half holds both batches averages the two
regimes rather than avoiding them. That is a defensible design, but it means a stratified estimate is
a blend of two quite different situations, not a neutral baseline.

## Limitations

Twenty cancer samples come from five patients, twelve of them from one, so every result leans on P1
and the bootstrap resamples patients rather than samples. Low tumor fractions are reached by
in-silico dilution, not by real low-burden samples. Each bin holds 12 virtual samples from only 3–5
patients, which is why AUROC is not monotonic across bins and why the floor needed an explicit
definition. In the controls, depth and batch identity are tangled; sequencing chemistry is a third
confound — SRA lists Kit 14 (R10.4.1) for every sample, but the free-text description for HU012 reads
"lsk109, nbd114", and LSK109 is a ligation kit used with the older R9.4.1 pores, so HU012 was most
likely R9.4.1. Setup (c) is therefore a cross-batch and probably cross-chemistry check, and the three
effects cannot be separated in this cohort.

## Contents

```
code/      make_splits.py, run_b.sh, score_b.py, rescore_controls.py
patches/   mix_dilutions.py — Katsman et al.'s script with one added argument
slurm/     run_b.SBATCH, run_c.SBATCH
meta/      splits.tsv, pools/ — the exact split assignment every run used
logs/      Slurm job output
results/   b/, c/, rescore/, repro/
```

`REPRODUCE.md` has the command sequence.
