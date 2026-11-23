#!/bin/bash
NP="${1:-12}"
bash step_scripts/common/inputs.sh | while read f; do
  b="results/step4/$(basename "$f" .pat)"
  echo "$f $b/whole u250 --bucket 4+"
  for n in 3 2 1; do
    echo "$f $b/size$n u250 --bucket 4+ --size $n"
  done
done | xargs -P "$NP" -L 1 bash step_scripts/common/run_one.sh
echo "step 4 finished: $(ls results/step4/*/*.ish_clipped.txt | wc -l) fits, $(ls results/step4/*/*.failed 2>/dev/null | wc -l) failed, $(ls results/step4/*/*.empty 2>/dev/null | wc -l) empty"
