set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export P=/scratch/mxf4642/code/project RUN=/scratch/mxf4642/data/katsman_run
export OUT=/scratch/mxf4642/results/m8 M6=/scratch/mxf4642/results/m6
export CELFIE=/scratch/mxf4642/code/celfie/scripts/celfie.py
export ATL_MAX=$RUN/blood_epithelial_atlas_over_regions_max_markers.txt ATL_U=$RUN/31CT_U250_Netanel_atlas_over_regions.txt
mkdir -p $OUT/est $OUT/celfie_tmp $OUT/logs
echo "== 2. original celfie.py (one sample per run, -u 0, other settings default)"
echo "-- relevant lines of celfie.py:"
grep -n -E 'tissue_proportions|to_csv|seed|def define_arrays|pseudo|header' $CELFIE | head -25
cf() { ms=$1; id=$2
  if [ $ms = max ]; then atl=$ATL_MAX; else atl=$ATL_U; fi
  o=$OUT/est/celfie-orig@$ms@$id.txt; [ -s $o ] && return
  d=$OUT/celfie_tmp/$ms/$id; rm -rf $d; mkdir -p $d
  lg=$OUT/logs/celfie-orig@$ms@$id.log
  python $P/make_celfie_input.py $OUT/sites/$ms/$id.tsv $atl $id $d/input.txt > $lg 2>&1 \
    && python $P/runtimed.py celfie-orig@$ms@$id python $CELFIE $d/input.txt $d 1 -u 0 >> $lg 2>&1 \
    && python $P/celfie_to_line.py $d $atl $o >> $lg 2>&1 && rm -f $d/input.txt; }
export -f cf
echo "-- test run on HU012.01.PL4406 (max markers)"
cf max HU012.01.PL4406
tail -8 $OUT/logs/celfie-orig@max@HU012.01.PL4406.log
ls $OUT/celfie_tmp/max/HU012.01.PL4406/
cat $OUT/est/celfie-orig@max@HU012.01.PL4406.txt
echo "-- full sweep"
for ms in max u250; do awk -v ms=$ms '{print ms, $1}' $M6/inputs.txt; done | xargs -P 16 -L 1 bash -c 'cf $0 $1'
for ms in max u250; do echo "celfie-orig $ms: $(ls $OUT/est/celfie-orig@$ms@* 2>/dev/null | wc -l)/142"; done
echo "== 3. summary of all methods"
python $P/collect_all.py | tee /scratch/mxf4642/results/summary/summary.txt
echo DONE
