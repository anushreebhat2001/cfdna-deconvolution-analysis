#!/bin/bash
NP="${1:-12}"
./env/bin/python -I step_scripts/step_d1/simulate_strings.py || exit 1
for f in results/d1/sim/*.pat; do
  o="results/d1/$(basename "$f" .pat)"
  echo "$f $o/whole u250"
  echo "$f $o/cut u250 --size 1"
done | xargs -P "$NP" -L 1 bash step_scripts/common/run_one.sh
echo "D1 finished: $(ls results/d1/*/*.ish_clipped.txt | wc -l) fits, $(find results/d1 -name '*.failed' | wc -l) failed"
