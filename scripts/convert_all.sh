#!/bin/bash
set -euo pipefail
NP="${1:-10}"
mkdir -p data/pat logs/convert
one() {
  local f="$1" base
  base="$(basename "$f" .epibed.gz)"
  [ -f "data/pat/$base.pat.gz.tbi" ] && return 0
  ./env/bin/python -I scripts/epibed2pat.py "$f" data/ref_genome/cpg_index.npz "data/pat/$base.pat" 2> "logs/convert/$base.log"
  bgzip -f "data/pat/$base.pat"
  tabix -f -s 1 -b 2 -e 2 "data/pat/$base.pat.gz"
}
export -f one
ls data/epibed/*.epibed.gz | xargs -P "$NP" -I{} bash -c 'one "$@"' _ {}
cat logs/convert/*.log > logs/convert_summary.tsv
