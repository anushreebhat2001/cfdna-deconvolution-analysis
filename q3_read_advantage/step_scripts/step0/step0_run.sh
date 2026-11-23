#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

REGIONS=data/ref/31CT_U250_Netanel_merged_regions_file.bed
ATLAS=data/ref/31CT_U250_Netanel_atlas_over_regions.txt
CPG=data/ref/hg38_pat_cpg.bed.gz
OUT=results/step0
SAMPLES="data/marker_pat/HU012.01.PL4406 data/marker_pat/HU083.68 data/marker_pat/HU030.09.PL4937 data/mix/HU030.09.PL4937__d3.125"

prep() {
  local src="$1" d="$OUT/$(basename "$1")"
  mkdir -p "$d"
  cp "$src.pat" "$d/orig.pat"
  ./env/bin/python -I step_scripts/step0/cut_strings.py "$src.pat" "$REGIONS" "$d/whole.pat"
  ./env/bin/python -I step_scripts/step0/cut_strings.py "$src.pat" "$REGIONS" "$d/cut.pat" --size 1
  ./env/bin/python -I step_scripts/step0/cut_strings.py "$src.pat" "$REGIONS" "$d/whole_1bead.pat" --bucket 1
  ./env/bin/python -I step_scripts/step0/cut_strings.py "$src.pat" "$REGIONS" "$d/cut_1bead.pat" --bucket 1 --size 1
  for f in orig whole cut; do
    bgzip -k -f "$d/$f.pat"
    tabix -f -s 1 -b 2 -e 2 "$d/$f.pat.gz"
  done
}

run() {
  local d="$1" f="$2" model="$3"
  local common="-a $ATLAS -m $d/$f.pat.gz --cpg_coordinates $CPG -i $REGIONS -b --epiformat pat --num_iterations 10000 --stop_criterion 0.0000001 --random_restarts 1"
  if [ "$model" = ish ]; then
    ./env/bin/python scripts/run_ish.py "$d/$f.ish_unclipped.txt" "$d/$f.ish_clipped.txt" $common > "$d/$f.ish.log" 2>&1
  else
    ./env/bin/deconvolution --model celfie $common --outfile "$d/$f.celfie.txt" > "$d/$f.celfie.log" 2>&1
  fi
}
export -f run
export REGIONS ATLAS CPG OUT

for s in $SAMPLES; do prep "$s"; done

jobs=()
for s in $SAMPLES; do
  d="$OUT/$(basename "$s")"
  jobs+=("$d orig ish" "$d whole ish" "$d cut ish" "$d whole celfie" "$d cut celfie")
done
printf '%s\n' "${jobs[@]}" | xargs -P 12 -I{} bash -c 'run $0' {}
echo "step 0 runs finished"
