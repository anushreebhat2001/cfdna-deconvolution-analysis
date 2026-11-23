# Reproducing Part A (Q1)

All scripts were written for the NYU HPC Cloud Bursting cluster. Paths at the top of each
`code/run_*.sh` point to `/scratch/mxf4642/...`: change `P` (this `code/` folder), `REF`, `RUN`,
`PAT`, `MIX` and `OUT` for another machine. Each `.SBATCH` file runs its `.sh` inside the course
Singularity container.

## Requirements

- Tools (commits in `tools.tsv`): `deconvolution_models` (CelFiE, CelFiE-ISH, UXM), CelFEER,
  original CelFiE (`celfie.py`), epiread-tools; Python with numpy, pandas, scikit-learn, pysam,
  matplotlib; htslib (bgzip, tabix).
- CelFEER fix: `cd CelFEER/scripts && cp celfeer.py celfeer_tiefix.py && patch celfeer_tiefix.py < <this folder>/patches/celfeer_tiefix.patch`
- Data: Katsman et al. 2026 epiBED files, atlases and marker files (Zenodo 18100426); hg38;
  Loyfer et al. 2023 reference pat files (on the course cluster mount).

## Run order (submit with `sbatch`)

| Step | Script | What it does |
|---|---|---|
| 1 | `download_zenodo.SBATCH` | Downloads the paper's data |
| 2 | `build_cpg_index.SBATCH`, `first_runs.SBATCH` | CpG-number index from plain hg38; converts and fits two test samples |
| 3 | `count_reads.SBATCH` | Distinct reads per sample (unique read names in each epiBED file) |
| 4 | `run_m4.SBATCH` | Converts all 42 epiBED files to pat (`epibed2pat.py`), checks CpG offsets (`check_shift.py`), CelFiE-ISH with and without clipping on both marker sets; `collect_m4.py` |
| 5 | `run_m5.SBATCH` | Read totals, 20-control healthy pool, 100 mixtures (`mix_dilutions.py`), CelFiE and CelFiE-ISH fits on the mixtures; `collect_m5.py` |
| 6 | `run_m6.SBATCH` | Reads per marker region for the Loyfer references and the 142 inputs (`count_regions.py`), UXM and CelFEER references (`aggregate_ref.py`), UXM and CelFEER on U250 |
| 7 | `check_groups.SBATCH` | Checks the max-marker grouping against the atlas |
| 8 | `run_m7.SBATCH` | Max-marker grouping, synthetic megakaryocyte reads (`synth_reads.py`), U250 granulocyte calibration, UXM and CelFEER on max markers |
| 9 | `run_m8_sites.SBATCH`, then `run_m8.SBATCH` | Original `celfie.py` on per-CpG counts |
| 10 | `run_m9.SBATCH` | CelFiE-ISH on collapsed 250 and 1,000 markers (`celfieish_both.py` saves the clipped and unclipped answer of one fit) |
| 11 | `chem_check.SBATCH` | HU012 sequencing chemistry check |
| 12 | `python code/collect_all.py`, then `python code/plot_part_a.py` | Summary, `lod_bins.tsv`, figures 1–4 and the method table |

Steps 6–10 end by running `collect_all.py`.

## Helper scripts

| Script | Use |
|---|---|
| `runtimed.py` | Wraps every run and logs runtime and memory |
| `summarise.py` | Prints the epithelial share of one estimate |
| `celfieish_both.py` | One CelFiE-ISH fit, both clipped and unclipped outputs |
| `make_celfeer_input.py`, `celfeer_to_line.py` | CelFEER input and output conversion |
| `pat_to_sites.py`, `make_celfie_input.py`, `celfie_to_line.py` | Input and output conversion for original `celfie.py` |
| `check_epibed*.py`, `check_shift.py` | Early checks of the epiBED format and CpG numbering |
| `rmse_variants.py` | Compares dilution-error definitions against the paper's numbers |
