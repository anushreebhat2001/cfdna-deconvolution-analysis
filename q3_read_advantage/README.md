# Read-level advantage vs CpGs per read (Q3)

We asked whether CelFiE-ISH detects tumor DNA better when it sees each read whole than when it
sees the same CpG calls one at a time, and whether the benefit grows with the number of marker
CpGs on a read. In our pre-specified tests, whole reads separate the 20 undiluted cancer samples
from the controls better than single CpGs do (AUROC 0.991 against 0.857). The gain grows with
CpGs per read (+0.048, +0.073 and +0.150 for reads with 2, 3 and ≥4 marker CpGs), including when
each group is matched on its number of CpGs. Across all 120 cancer files, including the
dilutions, the gain is not significant, and whole reads estimate tumor fraction less accurately
(mean absolute error 0.182 against 0.140). The gain also shows up where there is no tumor: whole
reads separate the two healthy batches from each other (AUROC 0.93). Because no healthy control
was sequenced in the same run as a cancer sample, the gain cannot be credited to tumor DNA in
this cohort.

## Design

One model, CelFiE-ISH, is fitted twice on identical CpG calls: "whole" keeps the part of each
read inside a marker region intact, and "single" feeds the same calls one CpG per line. Any
difference therefore comes only from keeping CpGs grouped by read. Reads are grouped by the
number of marker CpGs they carry (1, 2, 3, ≥4); 1-CpG reads give identical files by
construction. To check that longer reads do not win simply by carrying more CpGs, each group is
also thinned to equal CpG totals (10 random draws). The plan, with all tests and scores, was
fixed before any result was seen (10,766 fits, none failed); anything chosen after seeing
results is labelled exploratory.

## Setup

U250 atlas (31 cell types, with the 16 epithelial types summed as the tumor score); 20 CRC
samples from 5 patients, 22 controls and the 100 setup (a) mixtures. Detection is the AUROC of
cancer files against the 22 controls; error is the mean absolute error against ichorCNA ×
dilution. A gain counts only if its 95% bootstrap interval (patients and controls redrawn, 2,000
rounds) excludes zero, and every detection number is also computed against each healthy batch
alone. Values are unclipped; clipped runs give smaller gains, often with intervals reaching
zero. Checks before fitting: cutting keeps every call (for example, 2,023 reads become 5,883
single-CpG lines), and site-level CelFiE gives the same estimate on both files (largest
difference 0.006 percentage points). One check failed because of a single read in one control,
which does not affect the comparison.

## Planned tests

| Planned test | Whole | Single | Gain [95% CI] |
|---|---|---|---|
| Detection, all reads, 20 undiluted CRC | 0.991 | 0.857 | +0.134 [0.045, 0.327] |
| Detection, all reads, all 120 CRC files | 0.821 | 0.780 | +0.041 [−0.021, 0.170] |
| 2-CpG reads, undiluted | 0.868 | 0.820 | +0.048 [0.010, 0.118] |
| 3-CpG reads, undiluted | 0.925 | 0.852 | +0.073 [0.018, 0.159] |
| ≥4-CpG reads, undiluted | 0.986 | 0.836 | +0.150 [0.038, 0.277] |
| 2 CpGs, equal CpG totals | 0.870 | 0.816 | +0.053 [0.017, 0.126] |
| 3 CpGs, equal CpG totals | 0.923 | 0.839 | +0.084 [0.027, 0.166] |
| ≥4 CpGs, equal CpG totals | 0.953 | 0.864 | +0.089 [0.052, 0.151] |
| Error against ichorCNA, all 120 files | 0.182 | 0.140 | −0.042 [−0.061, −0.036] |
| Dilutions into HU012 reads only | 0.846 | 0.625 | +0.222 [0.121, 0.358] |
| Dilutions into HU083 reads only | 0.773 | 0.854 | −0.081 [−0.199, −0.016] |
| Simulation: error | 0.034 | 0.059 | +0.025 [0.020, 0.041] |
| Simulation: healthy false signal (median) | 1.8% | 5.1% | |

AUROC unless stated. A positive gain means whole reads are better; for error, a negative gain
means whole reads are worse. The simulation assumes CpGs are independent within each cell type.

## Planned results

Most of the gain comes from the first 2–3 CpGs on a read. When ≥4-CpG reads are cut into
pieces, pieces of 2 CpGs already reach AUROC 0.970 and pieces of 3 reach 0.984, against 0.986
for the whole read. In real data, whole reads give a larger error because their false baseline
in healthy controls is higher (median 18.6% against 15.3%). After subtracting each version's
own control baseline, whole reads are slightly better but not significantly (0.047 against
0.059; +0.012 [−0.002, 0.020]; exploratory). When every dilution's healthy background comes
from one batch, the sign of the gain depends on the batch.

Among healthy reads with ≥4 CpGs, fully unmethylated reads are 4.1× more common than
independent CpGs predict in HU012 (129 of 7,267 against 31.8 expected; 95% CI 2.9–6.5) and 8.0×
in HU083 (44 of 15,172 against 5.5; 4.7–16.7). In cancers with ≥40% tumor the ratio is about
500× against the pooled controls. Part of the healthy excess is expected from the mixture of
blood cell types itself, as admixture creates LD. In the simulation, whole reads have lower
error and less false signal, but the AUROC gain is not established (+0.039 [−0.035, 0.104]).

## Exploratory checks

Chosen after seeing the planned results.

Whole reads also separate groups with little or no tumor:

- CRC samples under 20% tumor (6): 0.970 against 0.530.
- Four breast cancer samples with almost no tumor by copy number, from the same flow cells as
  the CRC samples: 0.955 against 0.500.
- The two healthy batches against each other: 0.926 against 0.430 (+0.496 [0.273, 0.719]).

No healthy control was sequenced in the same run as any cancer sample, so run and disease are
confounded.

Fully unmethylated reads make up the following share of reads with ≥2 CpGs:

- HU012: 4.7%; HU083: 2.3%.
- Breast samples: 5.4%.
- CRC under 20% tumor: 6.4%; CRC with ≥20% tumor: 9.6%.

Adding an atlas column that is unmethylated at every CpG absorbs these reads:

- The healthy batches no longer separate (0.926 → 0.504), nor do the breast samples
  (0.955 → 0.557).
- CRC under 20% tumor still separates (whole 0.894 against single 0.530).

In simulated healthy controls:

- Errors that flip whole reads add false signal for whole reads (+0.077 [0.065, 0.089] against
  single CpGs).
- Errors scattered over single CpGs do not (−0.040 [−0.050, −0.031]).
- The real controls behave like the first case (+0.041).

On the 8-group max markers the effect shrinks:

- CRC: 1.00 against 0.95 (+0.05 [0.00, 0.18]).
- Healthy batches: 0.63 against 0.53.
- Healthy false signal: 2.7–2.9% against 0–0.1%.

## Answer to Q3

Not shown. The planned test behaves as read-level theory predicts, but in this cohort the
whole-read gain cannot be separated from the sequencing run. The likely mechanism is the excess
of fully unmethylated reads in healthy blood. CelFiE-ISH assumes CpGs are independent within a
cell type, so it does not expect these reads in blood and scores them as tumor. A clean test
needs controls sequenced in the same runs as the cancers, or a model that lets CpGs on a read be
correlated within a cell type.

Run order, environment and known issues: [REPRODUCE.md](REPRODUCE.md).
