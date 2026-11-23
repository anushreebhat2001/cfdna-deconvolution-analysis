#!/bin/bash
NP="${1:-12}"
for s in step4/run_step4.sh step3/run_step3.sh step2/run_step2_2plus.sh step3/run_step3_matched.sh \
         step_b3/run_b3.sh step_d1/run_d1.sh; do
  echo "$(date '+%m-%d %H:%M') start $s"
  bash "step_scripts/$s" "$NP"
  bash step_scripts/common/score_all.sh > logs/score_all.log 2>&1
  echo "$(date '+%m-%d %H:%M') scored"
done
echo "$(date '+%m-%d %H:%M') all done"
