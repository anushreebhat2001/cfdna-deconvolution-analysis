#!/bin/bash
NP="${1:-12}"
bash step_scripts/common/inputs.sh | while read f; do
  b="results/step3/$(basename "$f" .pat)"
  echo "$f $b/b1_whole u250 --bucket 1"
  for k in 2 3; do
    echo "$f $b/b${k}_whole u250 --bucket $k"
    echo "$f $b/b${k}_cut u250 --bucket $k --size 1"
  done
done | xargs -P "$NP" -L 1 bash step_scripts/common/run_one.sh
echo "step 3 unmatched finished: $(ls results/step3/*/*.ish_clipped.txt | wc -l) fits, $(ls results/step3/*/*.failed 2>/dev/null | wc -l) failed, $(ls results/step3/*/*.empty 2>/dev/null | wc -l) empty"
