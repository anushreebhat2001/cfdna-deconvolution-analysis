#!/bin/bash
set -euo pipefail
NP="${1:-10}"
mkdir -p data/marker_pat logs/marker
one() {
  local f="$1" base
  base="$(basename "$f" .pat.gz)"
  ./env/bin/python -I scripts/marker_reads.py "$f" "data/marker_pat/$base.pat" \
    data/ref/31CT_U250_Netanel_merged_regions_file.bed \
    data/ref/blood_epithelial_merged_regions_file_max_markers.txt 2> "logs/marker/$base.log"
}
export -f one
ls data/pat/*.pat.gz | xargs -P "$NP" -I{} bash -c 'one "$@"' _ {}
cat logs/marker/*.log > logs/marker_summary.tsv
