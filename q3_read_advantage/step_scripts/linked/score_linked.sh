#!/bin/bash
L=results/linked
score() {
  local label="$1" a="$2" b="$3"; shift 3
  ./env/bin/python -I step_scripts/common/score_pair.py --meta meta/samples.tsv \
    --out $L/scores --label "$label" --a "$a" --b "$b" "$@" > /dev/null || echo "scoring failed: $label"
}
U="--markers u250"
score A1_linked_vs_single "$L/step2/{id}/whole_linked" "results/step2/{id}/cut" $U
score A1_linked_vs_ish    "$L/step2/{id}/whole_linked" "results/step2/{id}/whole" $U
score A1_samerho_vs_single "$L/step2/{id}/whole_samerho" "results/step2/{id}/cut" $U
score A2_linked_vs_single "$L/step2_2plus/{id}/whole_linked" "results/step2_2plus/{id}/cut" $U
score A2_linked_vs_ish    "$L/step2_2plus/{id}/whole_linked" "results/step2_2plus/{id}/whole" $U
score A3_linked_whole_vs_single "$L/step4/{id}/whole_linked" "results/step4/{id}/size1" $U
score A3_linked_whole_vs_ish    "$L/step4/{id}/whole_linked" "results/step4/{id}/whole" $U
score A3_linked_3bead_vs_single "$L/step4/{id}/size3_linked" "results/step4/{id}/size1" $U
score A3_linked_2bead_vs_single "$L/step4/{id}/size2_linked" "results/step4/{id}/size1" $U
score A4_bucket2_linked_vs_single "$L/step3/{id}/b2_whole_linked" "results/step3/{id}/b2_cut" $U
score A4_bucket3_linked_vs_single "$L/step3/{id}/b3_whole_linked" "results/step3/{id}/b3_cut" $U
for k in 2 3 4+; do
  score "A5_matched_bucket${k}_linked_vs_single" "$L/step3_matched/{id}/d{draw}_b${k}_whole_linked" \
    "results/step3_matched/{id}/d{draw}_b${k}_cut" $U --draws 10
done
for b in HU012 HU083; do
  score "B3_background_${b}_linked_vs_single" "$L/b3/$b/{id}/whole_linked" "results/b3/$b/{id}/cut" $U
  score "B3_background_${b}_linked_vs_ish" "$L/b3/$b/{id}/whole_linked" "results/b3/$b/{id}/whole" $U
  score "B3_background_${b}_samerho_vs_single" "$L/b3/$b/{id}/whole_samerho" "results/b3/$b/{id}/cut" $U
done
score D1_linked_vs_single "$L/d1/{id}/whole_linked" "results/d1/{id}/cut" $U
score D1_linked_vs_ish    "$L/d1/{id}/whole_linked" "results/d1/{id}/whole" $U
score PanEpi_linked_vs_single "$L/step7/{id}/whole_linked" "results/step7/{id}/cut" --markers max
score PanEpi_linked_vs_ish    "$L/step7/{id}/whole_linked" "results/step7/{id}/whole" --markers max
ls $L/scores/*_summary.tsv | wc -l
