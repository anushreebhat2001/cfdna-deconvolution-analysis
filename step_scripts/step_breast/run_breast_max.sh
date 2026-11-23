#!/bin/bash
NP="${1:-2}"
for s in HU017.04.PL4596 HU017and19.05.PL4434 HU017and19.06.PL4518 HU025.13.PL4704; do
  echo "data/marker_pat/$s.pat results/breast_max/$s/whole max"
  echo "data/marker_pat/$s.pat results/breast_max/$s/cut max --size 1"
done | xargs -P "$NP" -L 1 bash step_scripts/common/run_one.sh
echo "breast max finished: $(ls results/breast_max/*/*.ish_clipped.txt | wc -l) fits, $(find results/breast_max -name '*.failed' | wc -l) failed"
