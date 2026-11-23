#!/bin/bash
NP="${1:-12}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mkdir -p results/linked
echo "$(date '+%m-%d %H:%M') code checks"
./env/bin/python step_scripts/linked/check_model.py > results/linked/model_checks.txt 2> results/linked/model_checks.err
if ! grep -q "ALL CHECKS PASS" results/linked/model_checks.txt; then
  echo "code checks did not pass; no fits started. See results/linked/model_checks.txt"; exit 1
fi
echo "$(date '+%m-%d %H:%M') fits"
bash step_scripts/linked/run_linked.sh "$NP"
echo "$(date '+%m-%d %H:%M') scoring"
bash step_scripts/linked/score_linked.sh
./env/bin/python -I step_scripts/linked/readout.py > results/linked/readout.txt 2> results/linked/readout.err
echo "$(date '+%m-%d %H:%M') finished; read results/linked/readout.txt"
