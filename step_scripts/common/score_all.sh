#!/bin/bash
score() {
  local label="$1" dir="$2" a="$3" b="$4"; shift 4
  [ -d "$dir" ] || return 0
  ./env/bin/python -I step_scripts/common/score_pair.py --meta meta/samples.tsv --markers u250 \
    --out results/scores --label "$label" --a "$a" --b "$b" "$@" > /dev/null || echo "scoring failed: $label"
}
score A1_all_strings results/step2 "results/step2/{id}/whole" "results/step2/{id}/cut"
score A2_2plus_strings results/step2_2plus "results/step2_2plus/{id}/whole" "results/step2_2plus/{id}/cut"
score A3_ladder_whole_vs_single results/step4 "results/step4/{id}/whole" "results/step4/{id}/size1"
score A3_ladder_3bead_vs_single results/step4 "results/step4/{id}/size3" "results/step4/{id}/size1"
score A3_ladder_2bead_vs_single results/step4 "results/step4/{id}/size2" "results/step4/{id}/size1"
score A4_bucket1_whole_only results/step3 "results/step3/{id}/b1_whole" "results/step3/{id}/b1_whole"
score A4_bucket2 results/step3 "results/step3/{id}/b2_whole" "results/step3/{id}/b2_cut"
score A4_bucket3 results/step3 "results/step3/{id}/b3_whole" "results/step3/{id}/b3_cut"
for k in 2 3 4+; do
  score "A5_matched_bucket$k" results/step3_matched "results/step3_matched/{id}/d{draw}_b${k}_whole" \
    "results/step3_matched/{id}/d{draw}_b${k}_cut" --draws 10
done
for b in HU012 HU083; do
  score "B3_background_$b" "results/b3/$b" "results/b3/$b/{id}/whole" "results/b3/$b/{id}/cut"
done
score D1_simulated results/d1 "results/d1/{id}/whole" "results/d1/{id}/cut"
ls results/scores/*_summary.tsv
