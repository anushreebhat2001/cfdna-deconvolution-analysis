#!/bin/bash
set -uo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
stored="$1"; markers="$2"; rho="${3:-}"
if [ "$markers" = max ]; then
  REGIONS=data/ref/blood_epithelial_merged_regions_file_max_markers.txt
  ATLAS=data/ref/blood_epithelial_atlas_over_regions_max_markers.txt
else
  REGIONS=data/ref/31CT_U250_Netanel_merged_regions_file.bed
  ATLAS=data/ref/31CT_U250_Netanel_atlas_over_regions.txt
fi
if [ -n "$rho" ]; then out="results/linked/${stored#results/}_samerho"; extra="--rho $rho"
else out="results/linked/${stored#results/}_linked"; extra=""; fi
[ -f "$out.info.tsv" ] && exit 0
mkdir -p "$(dirname "$out")"
rm -f "$out.failed"
./env/bin/python step_scripts/linked/fit_linked.py "$stored" "$out.tmp" "$ATLAS" "$REGIONS" $extra > "$out.log" 2>&1
if [ $? -eq 0 ] && [ -s "$out.tmp.info.tsv" ]; then
  for e in ish_unclipped.txt ish_clipped.txt info.tsv; do mv "$out.tmp.$e" "$out.$e"; done
else
  rm -f "$out.tmp."*; touch "$out.failed"
fi
