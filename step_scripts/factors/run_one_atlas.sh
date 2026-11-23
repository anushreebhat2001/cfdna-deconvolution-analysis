#!/bin/bash
set -uo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
src="$1"; out="$2"; ATLAS="$3"; REGIONS="$4"; shift 4
[ -f "$out.ish_clipped.txt" ] || [ -f "$out.empty" ] && exit 0
mkdir -p "$(dirname "$out")"
rm -f "$out.failed"
./env/bin/python -I step_scripts/step0/cut_strings.py "$src" "$REGIONS" "$out.pat" "$@"
wc -l < "$out.pat" | tr -d ' ' > "$out.nlines"
if [ ! -s "$out.pat" ]; then rm -f "$out.pat"; touch "$out.empty"; exit 0; fi
bgzip -f "$out.pat" && tabix -f -s 1 -b 2 -e 2 "$out.pat.gz"
./env/bin/python scripts/run_ish.py "$out.ish_unclipped.tmp" "$out.ish_clipped.tmp" \
  -a "$ATLAS" -m "$out.pat.gz" --cpg_coordinates data/ref/hg38_pat_cpg.bed.gz -i "$REGIONS" -b \
  --epiformat pat --num_iterations 10000 --stop_criterion 0.0000001 --random_restarts 1 > "$out.log" 2>&1
if [ $? -eq 0 ] && [ -s "$out.ish_clipped.tmp" ]; then
  mv "$out.ish_unclipped.tmp" "$out.ish_unclipped.txt"; mv "$out.ish_clipped.tmp" "$out.ish_clipped.txt"
else
  rm -f "$out.ish_unclipped.tmp" "$out.ish_clipped.tmp"; touch "$out.failed"
fi
