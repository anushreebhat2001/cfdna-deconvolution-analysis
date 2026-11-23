#!/bin/bash
NP="${1:-12}"
bash step_scripts/common/inputs.sh | while read f; do
  b="results/step2/$(basename "$f" .pat)"
  echo "$f $b/whole u250"
  echo "$f $b/cut u250 --size 1"
done | xargs -P "$NP" -L 1 bash step_scripts/common/run_one.sh
echo "step 2 finished: $(ls results/step2/*/*.ish_clipped.txt | wc -l) fits, $(ls results/step2/*/*.failed 2>/dev/null | wc -l) failed"
