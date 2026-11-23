#!/bin/bash
NP="${1:-12}"
ONE=step_scripts/linked/run_one_linked.sh
run() { xargs -P "$NP" -L 1 bash "$ONE"; }

for f in results/step2/*/whole.pat.gz; do echo "${f%.pat.gz} u250"; done | run
echo "$(date '+%m-%d %H:%M') A1 done"

if [ ! -f results/linked/same_rho.txt ]; then
  awk -F'\t' 'NR>1{print "results/linked/step2/"$1"/whole_linked.info.tsv"}' meta/samples.tsv \
    | xargs grep -h '^rho	' | cut -f2 | sort -n | awk '{v[NR]=$1} END{printf "%.4f\n", (v[NR/2]+v[NR/2+1])/2}' \
    > results/linked/same_rho.txt
fi
RHO=$(cat results/linked/same_rho.txt)
echo "same-rho check uses rho = $RHO"

{
  for f in results/step2_2plus/*/whole.pat.gz results/step4/*/whole.pat.gz results/step4/*/size3.pat.gz \
           results/step4/*/size2.pat.gz results/step3/*/b2_whole.pat.gz results/step3/*/b3_whole.pat.gz \
           results/b3/HU012/*__d*/whole.pat.gz results/b3/HU083/*__d*/whole.pat.gz results/d1/*/whole.pat.gz \
           results/factors/t10_*/*/whole.pat.gz results/breast/*/whole.pat.gz; do
    echo "${f%.pat.gz} u250"
  done
  for f in results/step7/*/whole.pat.gz; do echo "${f%.pat.gz} max"; done
  for f in results/step2/*/whole.pat.gz results/b3/HU012/*__d*/whole.pat.gz results/b3/HU083/*__d*/whole.pat.gz; do
    echo "${f%.pat.gz} u250 $RHO"
  done
  for f in results/step3_matched/*/d*_b*_whole.pat.gz; do echo "${f%.pat.gz} u250"; done
} | run

for b in HU012 HU083; do
  awk -F'\t' '$2=="control"{print $1}' meta/samples.tsv | while read c; do
    ln -sfn "../../step2/$c" "results/linked/b3/$b/$c"
  done
done
echo "$(date '+%m-%d %H:%M') all done: $(find results/linked -name '*_linked.info.tsv' | wc -l) linked fits," \
     "$(find results/linked -name '*_samerho.info.tsv' | wc -l) same-rho fits, $(find results/linked -name '*.failed' | wc -l) failed"
