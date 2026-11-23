# Read-level vs site-level deconvolution of ultra-low-pass Nanopore cfDNA

A held-out test of the detection floor and the role of CpGs per read.

Asgar Fataymamode · Ronojoy Mitra · Anushree Bhat

Katsman et al. (bioRxiv 2026) report a 1.7–3.1% tumor-fraction detection floor for read-level
deconvolution (CelFiE-ISH) on 0.25× Nanopore cfDNA. Two things about that result invited a closer
look. Their mixing script builds every dilution's healthy background from 20 controls while the
paper scores against all 22, so mixing and scoring controls overlap. And CelFiE-ISH (Unterman et al.
2024) showed the read-level gain over site-level growing up to about 4 CpGs per read, but Katsman
et al. report no per-read figures. We reproduced the paper, benchmarked four methods on its own
data, and stress-tested both claims.

## The three questions

| | Question | Answer | Folder |
|---|---|---|---|
| **Q1** | Do the published results reproduce, and how do the four methods compare? | **Yes.** 6 of 7 published accuracy values match within 0.006 RMSE. Marker set matters more than method; CelFiE-ISH has the lowest floor in both marker sets. | [`q1_reproduction_benchmark/`](q1_reproduction_benchmark/) |
| **Q2** | Does the floor hold when mixing and scoring controls are kept separate? | **Yes.** Removing the overlap costs 0.034 AUROC in the lowest bin and leaves the floor unchanged in 14 of 20 splits. But swapping the control *batch* moves it by 0.288. | [`q2_held_out_controls/`](q2_held_out_controls/) |
| **Q3** | Does the read-level advantage grow with marker CpGs per read? | **Not shown.** The gain is there and grows with CpGs per read, but it cannot be separated from the sequencing run in this cohort. | [`q3_read_advantage/`](q3_read_advantage/) |

Each folder has its own `README.md` with the full result and a `REPRODUCE.md` with the run order.

## Shared setup

All three questions run on the same data and scoring rules, so they are stated once here.

- **Samples.** 20 colorectal cancer plasma samples from 5 patients (12 from one patient) and 22
  healthy controls in two sequencing batches, HU012 and HU083 (11 each). Ultra-low-pass Nanopore,
  Katsman et al. 2026, Zenodo 18100426. Atlas from Loyfer et al. 2023 (GSE186458).
- **Dilutions.** 100 in-silico mixtures: each CRC sample mixed with reads from the paper's
  20-control pool at 50, 25, 12.5, 6.25 and 3.125%.
- **Marker sets.** U250 (31 cell types) and max markers (8 groups, pan-epithelial); collapsed 250
  and 1,000 markers per group for the marker-count check in Q1.
- **Tumor score.** The summed epithelial share, since colon cancer is epithelial and healthy blood
  carries almost none.
- **Detection floor.** The lowest tumor-fraction bin from which AUROC stays at or above 0.95 — the
  *sustained* definition, confirmed with the course advisor. First crossing is reported as
  secondary throughout. The distinction matters because AUROC is not monotonic across bins: each
  bin holds 12 virtual samples drawn from only 3–5 patients.
- **Intervals.** Bootstrap resampling patients, not samples, so that P1's twelve samples do not
  count as twelve independent observations.

## Reading the three together

Three threads run across the questions that none of the individual write-ups states on its own.

**At low tumor fraction the control cohort, not the method, is the binding constraint.** Each
question found this independently and by a different route. Q1 found that shallower controls get
more false epithelial signal on U250 (Spearman −0.76). Q2 found that swapping which control batch
the mixtures are scored against moves AUROC by 0.288 and the floor by a full bin — eight times the
leakage effect it set out to measure — traced to a 0.003 difference in the healthy epithelial
baseline between batches, which is about a third of the signal at a 1% floor. Q3 found that whole
reads separate the two healthy batches from *each other* at AUROC 0.926, with no tumor involved at
all. A published detection floor is therefore a property of the method **and** of the cohort it was
measured against, and should be reported with that cohort named.

**The conditional-independence assumption is doing real damage.** CelFiE-ISH's likelihood treats
CpGs on a read as independent within a cell type. Q3 shows what that costs: fully unmethylated
reads are 4–8× more common in healthy blood than independence predicts, the model does not expect
them there, and it scores them as tumor. Adding an all-unmethylated atlas column absorbs them and
the spurious healthy-batch separation disappears (0.926 → 0.504). Q2 records a different pathology
in the same family — site-level CelFiE returns exact zeros for most controls on max markers, so its
apparent batch-robustness is a collapsed output distribution rather than a strength. A
haplotype-aware EM is the natural next step for both and is the clearest piece of future work this
project points to.

**Marker choice outweighs method choice.** Q1's floors differ more across marker sets than across
methods: CelFiE-ISH reaches 1.2–1.7% on max markers against 6.9–11.9% on U250, and the other three
methods land at 3.1–4.3% and 12.5–28% respectively. The marker-count sweep places the paper's own
curve between the 250- and 1,000-marker runs, which accounts for most of the gap between our floor
and theirs. Anyone applying these methods should spend their effort on the marker set first.

## Limitations that apply throughout

Twenty cancer samples come from five patients, twelve of them from one, so every result leans on
P1. Low tumor fractions are reached by in-silico dilution rather than by real low-burden samples.
Each tumor-fraction bin holds 12 virtual samples from 3–5 patients, which is why AUROC is
non-monotonic and the floor needed an explicit definition.

In the controls, three things are confounded and cannot be separated in this cohort: HU012 is the
shallower batch, it differs from HU083 in age and sex composition, and the SRA description for
HU012 reads "lsk109, nbd114" — LSK109 being a ligation kit used with the older R9.4.1 pores, where
every sample is otherwise listed as Kit 14 (R10.4.1). Separately, no healthy control was sequenced
in the same run as any cancer sample, which is what blocks the Q3 answer.

On max markers, UXM and CelFEER rely on synthetic megakaryocyte reads, since the real ones are not
public. The approximation was calibrated against a real cell type and costs at most 0.03 AUROC
there, but megakaryocyte markers are intermediately methylated, where an independent-sites
approximation is weakest.

## Repository map

| Folder | Contents |
|---|---|
| `q1_reproduction_benchmark/` | Full pipeline, four-method benchmark, marker-count sweep, figures fig1–fig5, chemistry check, runtime table |
| `q2_held_out_controls/` | Paired held-out design (setup b), batch swap (setup c), scoring-set-size control, split and pool manifests |
| `q3_read_advantage/` | Whole-read vs single-CpG comparison, CpG-matched strata, simulation, exploratory confound checks |

Each folder holds `code/`, `results/`, `logs/`, its own `README.md`, `REPRODUCE.md` and
`.gitignore`. Patched upstream code lives in `patches/` with a note on what changed.

**Not stored here** (roughly 11 GB): read data, converted pat files, mixtures, reference read
counts, and the conda/Singularity environments. Every one is regenerable from the scripts and the
recorded seeds; the `REPRODUCE.md` files say how.

## Deliverables

As proposed: a method comparison table (Q1), detection-floor curves for each method under each
setup (Q1, Q2), a figure of the read-level advantage against marker CpGs per read with and without
CpG-total matching (Q3), and this repository, which reproduces every figure.

## Sources

- Katsman et al., bioRxiv 2026 — the paper reproduced here. Data: Zenodo record 18100426.
- Unterman et al., *Genome Biology* 2024 — CelFiE-ISH.
- Loyfer et al., *Nature* 2023 — the methylation atlas (GSE186458).
