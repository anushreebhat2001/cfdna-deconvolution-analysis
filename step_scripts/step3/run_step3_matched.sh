#!/bin/bash
NP="${1:-12}"; DRAWS="${2:-10}"
bash step_scripts/common/inputs.sh | while read f; do
  b="results/step3_matched/$(basename "$f" .pat)"
  for d in $(seq 1 "$DRAWS"); do
    for k in 2 3 4+; do
      echo "$f $b/d${d}_b${k}_whole u250 --bucket $k --match $d"
      echo "$f $b/d${d}_b${k}_cut u250 --bucket $k --match $d --size 1"
    done
  done
done | xargs -P "$NP" -L 1 bash step_scripts/common/run_one.sh
echo "step 3 matched finished: $(find results/step3_matched -name '*.ish_clipped.txt' | wc -l) fits, $(find results/step3_matched -name '*.failed' | wc -l) failed, $(find results/step3_matched -name '*.empty' | wc -l) empty"
