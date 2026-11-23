#!/bin/bash
NP="${1:-12}"
for b in HU012 HU083; do
  mkdir -p "results/b3/$b"
  awk -F'\t' '$2=="control"{print $1}' meta/samples.tsv | while read c; do
    ln -sfn "../../step2/$c" "results/b3/$b/$c"
  done
  for f in data/mix_$b/*__d*.pat; do
    o="results/b3/$b/$(basename "$f" .pat)"
    echo "$f $o/whole u250"
    echo "$f $o/cut u250 --size 1"
  done
done | xargs -P "$NP" -L 1 bash step_scripts/common/run_one.sh
echo "B3 finished: $(ls results/b3/*/*__d*/*.ish_clipped.txt | wc -l) fits, $(find results/b3 -name '*.failed' | wc -l) failed"
