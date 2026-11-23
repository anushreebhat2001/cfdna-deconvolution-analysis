# Factor tests: what else could cause the whole-read advantage?

EXPLORATORY. Everything here was chosen after seeing the frozen-plan results. Written on
2026-11-18, before any of these tests was run. Every test is reported, including nulls.

Fixed in advance
- Fits used: all strings, U250 markers, unclipped answer (clipped shown beside it, never instead).
- Epithelial score: the 16 planned types. Gut types: Colon-Ep, Small-Int-Ep, Gastric-Ep.
  Non-gut epithelial: the other 13.
- Placebo types (should not be in blood at all): Neuron, Oligodend, Heart-Cardio, Adipocytes,
  Pancreas-Alpha, Pancreas-Beta, Pancreas-Delta.
- "Low-tumor" cancer samples: copy-number tumor share under 20% (6 samples).
- Every detection number is given against all controls, HU012 only and HU083 only, with the
  healthy HU012-against-HU083 number beside it.
- "Moved" = 95% bootstrap interval excludes zero (patients and controls redrawn, 2000 rounds).

Stage 1, no fits (step_scripts/factors/stage1_checks.py -> results/factors/stage1.txt)
| Test | Factor | Supports the factor if | Refutes it if |
|---|---|---|---|
| T1 placebo score | leftover spread over many columns / all-off strings | whole-over-single gain on the placebo and non-gut scores is about as large as on the epithelial score, in cancer-vs-controls and in healthy-vs-healthy | the gain appears only on the gut score |
| T2 offset null | a constant run offset is enough | single-bead score + the breast samples' whole-minus-single gap reproduces the whole-string detection score (within 0.03) | it stays well below the observed whole score |
| T3 background swap | dilution dilutes the run signature too | a sample's dilution slope differs between HU083 and HU012 background by about the difference of the two healthy floors | slopes are the same in both backgrounds |
| T4 gut-minus-placebo all-off rate | real tumor signal vs tissue-blind all-off strings | the contrast rises with tumor share but is about zero in low-tumor and breast samples | low-tumor samples show a clear gut excess |
| T5 age-adjusted score | age | the whole-over-single gain shrinks after removing the age line fitted in controls | gain unchanged |
| T6 pan-epithelial readout | missing megakaryocyte column | on the atlas that has the column the healthy floor is several times lower and healthy batches separate less | floors and batch separation as on U250 |

Stage 2, small fits (scripts in step_scripts/factors/, results in results/factors/)
| Test | Factor | Supports the factor if | Refutes it if |
|---|---|---|---|
| T7 extra "off everywhere" column (control arm: "on everywhere" column) | all-off strings tied to no tissue | with the off column the healthy whole floor drops to about the single-bead floor and healthy batches stop separating, and the on column does not do the same | nothing changes, or the on column does it equally (then it is generic leftover, not all-off strings) |
| T8 pan-epithelial atlas without megakaryocytes | missing cell type | the healthy epithelial floor rises several-fold, more for whole | unchanged |
| T9 breast samples diluted into healthy background | dilution dilutes the run signature | the whole-string estimate falls along the dilution series although there is almost no tumor | flat series |
| T10 simulation: bead flips vs whole strings switched off | calling error vs string-level events | only string switch-off gives whole > single on healthy files | bead flips alone give whole > single |

Fairness rules (checked by step_scripts/factors/check_inputs.py before any fit starts)
- Whole and cut versions of a file hold exactly the same beads.
- A modified atlas equals the original in every untouched column.
- Refitting with the unmodified atlas through the new runner reproduces the stored answer.
- Every arm of a test is fitted on the same input files as the arm it is compared with.
- Simulation variants keep every string's place and bead count; only bead values change, and
  the realised flip / switch-off rates are reported next to the nominal ones.
- Where two arms cannot hold the same amount of data (dilutions into different backgrounds),
  strings and beads per file are reported for both arms beside the result.

Amendments made on 2026-11-18 while writing the scripts, before any test was run
- T1, reworded. The factor says: under whole strings the leftover lands in many columns, so a
  score built from types that cannot be in blood should separate the groups about as well as
  the epithelial score does. Supported if, under whole strings, the placebo score separates
  low-tumor and breast samples from controls, and the healthy batches from each other, within
  0.05 of the epithelial score. Refuted if the placebo score sits near 0.5 where the epithelial
  score is high. (The first wording compared whole-over-single gains; that is the wrong
  quantity, because single beads may carry the same leftover in other columns.)
- T3, withdrawn as a test. A dilution file is a mixture, so a series is a straight line from the
  background's floor to the undiluted sample's estimate whatever that estimate is made of. A
  different slope in two backgrounds follows from their different floors and cannot separate
  tumor from run. T3 is kept only as a check that the series are straight lines.
- T9 for the same reason is a demonstration, not a discriminating test: it shows what a dilution
  series of an almost tumor-free sample looks like when scored like the cancer series.
- T10, arms fixed so that they are comparable: "clumped" switches whole strings off with
  probability 1% or 3%; "scattered" switches the same expected number of beads off, one bead at
  a time (each on-bead with probability 1% or 3%); "symmetric" flips every bead with probability
  5% or 10% (calling error in both directions). 22 simulated controls x 6 conditions x whole and
  single = 264 fits. With every switch at zero the generator must reproduce the stored D1 files.
