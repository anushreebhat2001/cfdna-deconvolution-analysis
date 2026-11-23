set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export P=/scratch/mxf4642/code/project REF=/scratch/mxf4642/data/ref RUN=/scratch/mxf4642/data/katsman_run
export D=/scratch/mxf4642/data/katsman/epibed_biscuit_filtered_bams_links PAT=/scratch/mxf4642/data/pat
export MIX=/scratch/mxf4642/data/mix_a OUT=/scratch/mxf4642/results/m5 META=/scratch/mxf4642/code/project/meta/samples.tsv
mkdir -p $OUT/est $OUT/logs $MIX

echo "== 1. read totals (all reads vs reads with CpG calls)"
if [ ! -s $P/meta/read_totals.tsv ]; then
  tail -n +2 $META | cut -f1 | xargs -P 12 -I{} bash -c \
    'echo -e "{}\t$(zcat $D/{}.epibed.gz | wc -l)\t$(zcat $PAT/{}.pat.gz | wc -l)"' | sort > $OUT/totals.tmp
  (echo -e "sample\tn_reads\tn_pat"; cat $OUT/totals.tmp) > $P/meta/read_totals.tsv
fi
head -3 $P/meta/read_totals.tsv

echo "== 2. healthy pool (20 controls, as in the paper)"
POOL=$MIX/healthy_pool.pat
if [ ! -s $POOL.done ]; then
  for s in $(awk -F'\t' 'NR>1 && $7=="yes"{print $1}' $META); do zcat $PAT/$s.pat.gz; done > $POOL
  wc -l < $POOL > $POOL.done
fi
echo "pool lines: $(cat $POOL.done)   controls in pool: $(awk -F'\t' 'NR>1 && $7=="yes"' $META | wc -l)"

echo "== 3. mixtures"
mix() { s=$1; [ -f $MIX/$s.done ] && return
  seed=$(echo $s | cksum | cut -d' ' -f1)
  python $P/runtimed.py mix_$s python $P/mix_dilutions.py $META $P/meta/read_totals.tsv $MIX/healthy_pool.pat $PAT $MIX $s $seed > $OUT/logs/mix_$s.log 2>&1 || return
  for f in $MIX/${s}__d*.pat; do
    sort -k2,2n -S 500M -T $MIX $f | bgzip > $f.gz && rm $f && tabix -f -s 1 -b 2 -e 2 $f.gz
  done
  touch $MIX/$s.done
}
export -f mix
awk -F'\t' 'NR>1 && $2=="crc"{print $1}' $META | xargs -P 6 -I{} bash -c 'mix {}'
echo "mixtures: $(ls $MIX/*__d*.pat.gz | wc -l) / 100"

run() { method=$1; ms=$2; id=$3; f=$4
  if [ $ms = u250 ]; then A=$RUN/31CT_U250_Netanel_atlas_over_regions.txt; R=$RUN/31CT_U250_Netanel_merged_regions_file.bed
  else A=$RUN/blood_epithelial_atlas_over_regions_max_markers.txt; R=$RUN/blood_epithelial_merged_regions_file_max_markers.txt; fi
  COMMON="-a $A -m $f --cpg_coordinates $REF/hg38_pat_cpg.bed.gz -i $R -b --epiformat pat --num_iterations 10000 --stop_criterion 0.0000001 --random_restarts 1"
  if [ $method = celfieish ]; then
    o2=$OUT/est/celfieish-clip@$ms@$id.txt; [ -s $o2 ] && return
    python $P/runtimed.py celfieish@$ms@$id python $P/celfieish_both.py $OUT/est/celfieish-noclip@$ms@$id.txt $o2 0 \
      --model celfie-ish $COMMON --outfile /dev/null > $OUT/logs/celfieish@$ms@$id.log 2>&1
  else
    o=$OUT/est/celfie@$ms@$id.txt; [ -s $o ] && return
    python $P/runtimed.py celfie@$ms@$id deconvolution --model celfie $COMMON --outfile $o > $OUT/logs/celfie@$ms@$id.log 2>&1
  fi
}
export -f run

echo "== 4. wrapper check (compare with the earlier run of HU030.01, max markers)"
run celfieish max CHECK_HU030.01.PL3393 $PAT/HU030.01.PL3393.pat.gz
python $P/summarise.py $RUN/blood_epithelial_atlas_over_regions_max_markers.txt $OUT/est/celfieish-clip@max@CHECK_HU030.01.PL3393.txt | head -1
python $P/summarise.py $RUN/blood_epithelial_atlas_over_regions_max_markers.txt /scratch/mxf4642/results/m4/est/max_clip_HU030.01.PL3393.txt | head -1

echo "== 5. deconvolution runs"
{
  for m in $MIX/*__d*.pat.gz; do id=$(basename $m .pat.gz)
    for ms in max u250; do echo "celfieish $ms $id $m"; done; done
  for m in $MIX/*__d*.pat.gz; do id=$(basename $m .pat.gz)
    for ms in max u250; do echo "celfie $ms $id $m"; done; done
  for s in $(tail -n +2 $META | cut -f1); do
    for ms in max u250; do echo "celfie $ms $s $PAT/$s.pat.gz"; done; done
} | xargs -P 14 -L 1 bash -c 'run $0 $1 $2 $3'
echo "estimate files: $(ls $OUT/est/*.txt | wc -l)  (expected 400 celfie-ish + 284 celfie + 2 check)"
python $P/collect_m5.py
echo DONE
