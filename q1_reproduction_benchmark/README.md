# Reproduction and four-method benchmark (Q1)

We reran Katsman et al. (2026) end to end and compared four deconvolution methods on the same
data: site-level CelFiE and read-level CelFiE-ISH (with and without clipping), CelFEER and UXM.
The paper reproduces: 6 of 7 published accuracy values match within 0.006 RMSE. The marker set
matters more than the method, and CelFiE-ISH has the lowest detection floor in both marker sets:
1.2–1.7% tumor fraction with max markers and 6.9–11.9% with U250. Our max-marker floor is lower
than the paper's curve (3.1–4.3% under our floor rule; the paper's text states 1.7–3.1%).
Rerunning with 250 and 1,000 markers per group places the paper between the two, so marker count
most likely explains most of the gap.

## Setup

- **Data:** 20 colorectal cancer (CRC) plasma samples from 5 patients (12 from one patient) and
  22 healthy controls in two batches (HU012, HU083), ultra-low-pass Nanopore (Katsman et al. 2026,
  Zenodo 18100426). 100 in-silico dilutions: each CRC sample mixed with reads from the paper's
  20-control pool at 50, 25, 12.5, 6.25 and 3.125%.
- **Marker sets:** U250 (31 cell types) and max markers (8 groups); collapsed 250 and 1,000
  markers per group for the marker-count check.
- **Tumor score:** the summed epithelial share.
- **Measures:** AUROC of cancer files against the 22 controls per tumor-fraction bin (10 bins of
  12 virtual samples) with patient-level bootstrap intervals; detection floor = lowest bin from
  which AUROC stays ≥ 0.95 (sustained; first crossing also reported); RMSE of the epithelial share
  against ichorCNA on the 20 undiluted CRC samples; dilution tracking; runtime and memory.

## Results

### Accuracy against the paper (RMSE vs ichorCNA, 20 CRC samples; lower is better)

| Method | Markers | Ours | Paper |
|---|---|---|---|
| CelFiE-ISH, clipping | U250 | 0.098 | 0.099 |
| CelFiE-ISH, no clipping | U250 | 0.132 | 0.126 |
| CelFEER | U250 | 0.126 | 0.125 |
| UXM | U250 | 0.195 | 0.195 |
| CelFiE | U250 | 0.094 | 0.104 |
| CelFiE-ISH, clipping | Collapsed 250 | 0.075 | 0.077 |
| CelFiE-ISH, no clipping | Collapsed 250 | 0.080 | 0.082 |

### Detection floor (sustained, AUROC ≥ 0.95)

| Markers | CelFiE-ISH | CelFiE, CelFEER*, UXM* |
|---|---|---|
| Max markers | 1.2–1.7% | 3.1–4.3% |
| U250 | 6.9–11.9% | 12.5–28% |

\* On max markers, UXM and CelFEER use synthetic megakaryocyte reads (see Checks).

### Marker count (CelFiE-ISH, clipping, same 8 groups)

| Markers per group | Floor |
|---|---|
| 250 | 4.4–6.7% |
| 1,000 | 1.7–3.0% |
| Max | 1.2–1.7% |
| Paper (Fig. 3F, our floor rule) | 3.1–4.3% |

In the lowest bins the paper's curve stays below all of ours, which the published code does not explain.

### Checks

- **Original CelFiE:** `celfie.py` (`-u 0`, one sample per run) agrees with the package CelFiE
  used above (AUROC 0.960 vs 0.959 on max markers, 0.850 vs 0.848 on U250; identical RMSE and
  floors). The remaining CelFiE difference from the paper is therefore not the implementation.
- **Synthetic megakaryocyte reads:** UXM and CelFEER need real reads for every reference group,
  and megakaryocyte reads are not public. Real granulocyte reads serve as templates for read
  positions and lengths, and each CpG is redrawn independently from the megakaryocyte methylation
  level. Rebuilding the U250 granulocyte reference the same way leaves UXM unchanged (AUROC 0.952
  vs 0.952) and costs CelFEER 0.03 (0.950 vs 0.920), with the same floor. The max-marker grouping
  matches the atlas (methylation correlation ≥ 0.996).
- **Control depth:** on U250, shallower controls get more false epithelial signal (Spearman
  −0.76); the effect is not significant on max markers.
- **Chemistry:** the HU012 control files most likely mix an R9.4.1 run and an R10.4.1 re-run (SRA
  description "lsk109, nbd114"); HU012 is also the shallower batch, so depth, batch and chemistry
  are confounded.
- **Cost:** every run takes 3 minutes or less and under 1 GB on one core, except CelFiE on U250
  (about 15 minutes).

## Deviations from the paper's tools

- CelFEER: best-restart selection crashes on exactly tied log-likelihoods; we compare on the
  log-likelihood only (`patches/celfeer_tiefix.patch`, two lines). Results are unaffected.

## Folder map

| Path | Content |
|---|---|
| `code/` | All scripts, as run on the cluster (flat, because run scripts call helpers as `$P/<name>`) |
| `code/meta/` | `samples.tsv` (42 samples, group, patient, batch, ichorCNA), `crc_samples.tsv`, `read_totals.tsv` |
| `results/summary.txt` | All methods and measures, as printed by `collect_all.py` |
| `results/all_estimates.tsv` | Every estimate (method, markers, sample or mixture, epithelial share) |
| `results/lod_bins.tsv` | AUROC per tumor-fraction bin, with floors |
| `results/all_timings.tsv` | Runtime and memory per run |
| `results/m4_summary.txt` | First full CelFiE-ISH run on the 42 samples |
| `results/figures/` | fig1–fig4 (`plot_part_a.py`), fig5 (marker count vs paper), method comparison table |
| `results/chem/` | HU012 chemistry check: per-sample error profile (`chem_check.tsv`) and the SRA run records (`runinfo.csv`, `es.xml`, `sra_full.xml`) |
| `results/reference_groups/` | Grouping files for max markers and the U250 calibration |
| `logs/slurm/` | Slurm job outputs |
| `patches/` | CelFEER fix |
| `tools.tsv` | External tools and the commits used |
| `REPRODUCE.md` | Run order |

Not stored (about 11 GB): read data, converted pat files, mixtures, reference read counts,
per-run estimate files and logs, environments.
