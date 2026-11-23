#!/bin/bash
NP="${1:-12}"
bash step_scripts/common/inputs.sh | while read f; do
  b="results/step2_2plus/$(basename "$f" .pat)"
  echo "$f $b/whole u250 --bucket 2+"
  echo "$f $b/cut u250 --bucket 2+ --size 1"
done | xargs -P "$NP" -L 1 bash step_scripts/common/run_one.sh
echo "A2 finished: $(ls results/step2_2plus/*/*.ish_clipped.txt | wc -l) fits, $(find results/step2_2plus -name '*.failed' | wc -l) failed, $(find results/step2_2plus -name '*.empty' | wc -l) empty"
