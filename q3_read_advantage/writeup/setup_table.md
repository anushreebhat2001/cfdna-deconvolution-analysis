# How each group of samples was sequenced

Source: Katsman et al. 2026, Supplementary Table 7 (run information), checked against
`meta/samples.tsv` for all 42 study samples.

| Group | Samples used here | Run | Library kit | Flow cell (pore chemistry) | Instrument | Base-calling model in the table |
|---|---|---|---|---|---|---|
| Healthy controls "HU012" | 11 (women, 20-50) | HU012 | SQK-LSK109 + EXP-NBD104 | FLO-MIN106D (R9.4.1) | MinION Mk1B | Guppy, R9.4.1 5mC model |
| The same 11 people, run again | (inside the "HU012" files) | HU028 | SQK-NBD114.24 | FLO-MIN114 (R10.4.1) | MinION Mk1C | Guppy, R10.4.1 5mC model |
| Healthy controls "HU083" | 11 (7 women, 4 men, 24-79; 5 with Alzheimer's) | HU083 | SQK-NBD114.96 | FLO-PRO114M (R10.4.1) | PromethION P2S | sup v4.3.0 |
| Colorectal cancer, patients P2-P5 | 8 | HU017, HU023, HU025 | SQK-NBD114.24 | FLO-MIN114 (R10.4.1) | MinION Mk1C | Guppy, R10.4.1 5mC model |
| Colorectal cancer, patient P1 | 12 (one patient, aged 22-23) | HU030 | SQK-NBD114.24 | FLO-MIN114 (R10.4.1) | MinION Mk1C | sup v4.2.0 |
| Breast cancer, almost no tumor by copy number (exploratory) | 4 | HU017, HU017/19, HU025 | SQK-NBD114.24 | FLO-MIN114 (R10.4.1) | MinION Mk1C | Guppy, R10.4.1 5mC model |

What this means for the comparison

- No healthy sample was sequenced in the same run as any cancer sample.
- The public "HU012" files appear to hold both runs of those 11 people (old R9.4.1 pores and the
  R10.4.1 re-run). This is inferred: the table lists each person in both runs, the archive has one
  file per person, and the archived read file of one of them names both base-calling models in
  its header. The share of old-chemistry reads is unknown and the reads cannot be separated with
  the files available.
- The re-run (HU028) matches the setup of the P2-P5 cancer samples and of the four breast samples.
- For patient P1 the table and the paper disagree: the table gives a "sup v4.2.0" model, the
  methods section says all colorectal samples were base-called with Guppy.
- HU083 differs from every other group in instrument, kit and model version.
- The two healthy groups also differ as people (age, sex, Alzheimer's disease).
- Copy-number tumor share (ichorCNA) reads up to 4.1% in healthy people, so values under about
  3-4% are inside the healthy range. Two cancer samples (1.4% and 1.8%) fall there.
