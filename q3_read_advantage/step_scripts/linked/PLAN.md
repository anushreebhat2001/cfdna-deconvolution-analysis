# Linked model: does reading whole strings help once the model knows beads go together?

EXPLORATORY. Chosen after seeing the frozen-plan results and the factor tests, in answer to the
professor's first point (CelFiE-ISH treats the beads of a string as unrelated). Written on
2026-11-20 after the code had been tried on six files (the ones `check_model.py` uses; for
five of them only rho and the likelihood were looked at, for control HU012.01.PL4406 also the
cell-type shares) and before any group of files was fitted or compared. Every test below is
reported, including nulls.

## The model

CelFiE-ISH plus one number, rho ("togetherness"), fitted per file together with the cell-type
shares. With chance rho all beads of a string share one state (on with the cell type's average
reference rate over that string's beads); otherwise the beads are drawn one by one as in
CelFiE-ISH. Code: `linked_model.py`. Formula:

    independent(t) = product over beads of  rate[t, bead]      if the bead is on
                                            1 - rate[t, bead]  if the bead is off

    P(string | t) = (1 - rho) * independent(t)
                    + rho * ( mean_rate[t]       if every bead is on
                              1 - mean_rate[t]   if every bead is off
                              0                  if the string is mixed )

mean_rate[t] is the average reference rate of cell type t over the string's own beads. Only
strings with 2 or more beads tell anything about rho. Data handling (empty strings, CpGs without
reference data, the pseudocount) is inherited unchanged from the authors' CelfieISH class.
 rho = 0 is exactly CelFiE-ISH. A 1-bead string
is scored exactly as in CelFiE-ISH whatever rho is, so the single-bead arm does not change and
is not refitted.

## Fixed in advance

- Three arms on the same stored string files: **linked whole** (new), **CelFiE-ISH whole** and
  **CelFiE-ISH single** (both as stored). No file is rebuilt; the new fit reads the stored file.
- Scores, control groups and the "moved" rule are those of the frozen plan (`score_pair.py`,
  unchanged). The unclipped answer is read first; the clipped one is shown beside it.
- Files fitted: every stored file that holds strings with 2+ beads: A1, A2, A3 (whole, 3- and
  2-bead pieces), A4 (2- and 3-bead buckets), A5 (10 draws), B3, D1, T10, the 4 breast samples,
  and the pan-epithelial run (step 7).
- Same-rho check: because rho is fitted per file, cancer and control files could differ only
  through rho. A1 and B3 are therefore fitted a second time with one rho for every file: the
  median of the fitted rho over the 42 study samples (labels not used).

## Code checks that must pass first (`check_model.py` -> `results/linked/model_checks.txt`)

rho = 0 gives back the stored CelFiE-ISH answer; single-bead files give the stored answer for
any rho; the likelihood never goes down; the start does not matter; the fitted rho is at the
top of the likelihood curve; made-up data with a known rho and tumor share gives both back.
In the full run every fit repeats the first check on its own file and records the difference.

Amendment, 2026-11-20, before any study file was compared. The made-up-data check as first
written (single real-sized files, fitted rho within 3 standard errors of the truth) failed in
one of 12 rows: true rho 0.10, fitted 0.081. Diagnosis: not a code error. With the cell-type
shares held at the truth rho is recovered, and with ten files stacked into one the joint fit
gives 0.098 and 0.101 for a true 0.10. At real file size (about 3,600 strings, 31 free shares)
rho reads about 0.01 low. The check was rewritten as 6a (stacked files, within 0.01) and 6b
(single files, within 0.02). Consequence for reading the results: a fitted rho is a slightly
low estimate, and differences in rho of about 0.01 to 0.02 between files mean nothing.

## Outcomes fixed in advance

| Test | Supports "the linked model helps" if | Refutes it if |
|---|---|---|
| L1 simulation where beads really are unrelated (D1) | fitted rho stays near 0 (median under 0.02) and the scores stay within 0.01 of CelFiE-ISH whole | rho is clearly above 0 or the scores get worse: the model invents linkage |
| L2 simulation with whole strings switched off (T10 clumped; scattered and symmetric are the control arms) | on clumped files the healthy false signal of linked whole falls below CelFiE-ISH whole (0.097 at 1%, 0.168 at 3%) to no more than single (0.075, 0.089); on the control arms rho stays near 0 | the false signal stays at the CelFiE-ISH whole level |
| L3 the two healthy runs (A1) | healthy HU012 against healthy HU083 falls from 0.926 to under 0.70, and the healthy false signal of linked whole is no higher than single in both runs | it stays above 0.85 |
| L4 detection (A1 per control run; B3 with matching background and controls) | linked whole beats single with a moved gain of the same sign against both control runs, in B3's matched rows ("all", below 1%, or 1 to 5%) | no moved gain, or the sign changes with the control run |
| L5 how much (A1, A3) | error of linked whole is lower than single and the gain moved | not moved, or single is better |
| L6 CpGs per read (A5, equal bead totals) | the detection gain of linked whole over single rises from 2 to 3 to 4+ beads and is moved | flat or not moved |

The project's headline ("we cannot show that reading whole reads helps") changes only if L4 is
supported. L1 to L3 say whether the model does what it was built for; L5 and L6 are secondary.

## Limits known before running

- One rho for all cell types and markers, and the same rho for all-on and all-off strings.
- The reference still holds only per-CpG averages, so the locked state uses the average rate.
  A reference built from reads (how often each cell type is all-off at a marker) would be a
  different, larger project.
- T10's clumped files switch strings off whatever the cell type; the linked model ties the
  locked state to the cell type. It is tested on data it was not built from.
