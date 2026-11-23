#!/bin/bash
NP="${1:-8}"
for s in HU017.04.PL4596 HU017and19.05.PL4434 HU017and19.06.PL4518 HU025.13.PL4704; do
  echo "data/marker_pat/$s.pat results/breast/$s/whole u250"
  echo "data/marker_pat/$s.pat results/breast/$s/cut u250 --size 1"
done | xargs -P "$NP" -L 1 bash step_scripts/common/run_one.sh
echo "breast finished: $(ls results/breast/*/*.ish_clipped.txt | wc -l) fits, $(ls results/breast/*/*.failed 2>/dev/null | wc -l) failed"
