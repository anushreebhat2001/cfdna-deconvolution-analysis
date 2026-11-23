set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export P=/scratch/mxf4642/code/project RUN=/scratch/mxf4642/data/katsman_run
export OUT=/scratch/mxf4642/results/m8 M6=/scratch/mxf4642/results/m6
mkdir -p $OUT/sites/u250 $OUT/sites/max $OUT/logs
echo "== 1. per-CpG counts at atlas CpGs (input for original celfie.py)"
sites() { ms=$1; id=$2; f=$3
  if [ $ms = max ]; then atl=$RUN/blood_epithelial_atlas_over_regions_max_markers.txt; else atl=$RUN/31CT_U250_Netanel_atlas_over_regions.txt; fi
  o=$OUT/sites/$ms/$id.tsv; [ -s $o ] && return
  python $P/runtimed.py sites@$ms@$id python $P/pat_to_sites.py $f $atl $o > $OUT/logs/sites@$ms@$id.log 2>&1; }
export -f sites
for ms in u250 max; do awk -v ms=$ms '{print ms, $1, $2}' $M6/inputs.txt; done | xargs -P 16 -L 1 bash -c 'sites $0 $1 $2'
for ms in u250 max; do echo "site files $ms: $(ls $OUT/sites/$ms/*.tsv 2>/dev/null | wc -l) / 142"; done
echo "examples:"; grep -h -E 'atlas_CpGs' $OUT/logs/sites@max@HU012.01.PL4406.log $OUT/logs/sites@max@HU030.01.PL3393.log $OUT/logs/sites@u250@HU012.01.PL4406.log
echo DONE
