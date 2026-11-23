#!/bin/bash
NP="${1:-12}"
bash step_scripts/common/inputs.sh | while read f; do
  b="results/step7/$(basename "$f" .pat)"
  echo "$f $b/whole max"
  echo "$f $b/cut max --size 1"
done | xargs -P "$NP" -L 1 bash step_scripts/common/run_one.sh
echo "step 7 finished: $(ls results/step7/*/*.ish_clipped.txt | wc -l) fits, $(ls results/step7/*/*.failed 2>/dev/null | wc -l) failed"
