set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export P=/scratch/mxf4642/code/project REF=/scratch/mxf4642/data/ref RUN=/scratch/mxf4642/data/katsman_run
export PAT=/scratch/mxf4642/data/pat MIX=/scratch/mxf4642/data/mix_a OUT=/scratch/mxf4642/results/m6
export C=/scratch/mxf4642/data/katsman_code/code LOY=/mnt/mrd/datasets/extracted CF=/scratch/mxf4642/code/CelFEER
export REG=$RUN/31CT_U250_Netanel_merged_regions_file.bed REFD=/scratch/mxf4642/data/readref
mkdir -p $OUT/est $OUT/logs $OUT/feat $OUT/celfeer_tmp $REFD/loyfer

echo "== 1. count reads in marker regions for each Loyfer reference sample"
cnt_loy() { n=$1; o=$REFD/loyfer/$n.npz; [ -s $o ] && return
  f=$(ls $LOY/*_${n}.hg38.pat.gz | head -1)
  python $P/runtimed.py refcount_$n python $P/count_regions.py $f $REG $o > $OUT/logs/refcount_$n.log 2>&1; }
export -f cnt_loy
tail -n +2 $C/group_31CT.csv | cut -d, -f1 | xargs -P 16 -I{} bash -c 'cnt_loy {}'
echo "reference samples counted: $(ls $REFD/loyfer/*.npz | wc -l) / $(tail -n +2 $C/group_31CT.csv | wc -l)"

echo "== 2. build UXM and CelFEER references"
python $P/aggregate_ref.py $C/group_31CT.csv $REFD/loyfer $RUN/31CT_U250_Netanel_atlas_over_regions.txt $REG $REFD

echo "== 3. count reads in marker regions for our 42 samples + 100 mixtures"
{ for s in $(tail -n +2 $P/meta/samples.tsv | cut -f1); do echo "$s $PAT/$s.pat.gz"; done
  for m in $MIX/*__d*.pat.gz; do echo "$(basename $m .pat.gz) $m"; done; } > $OUT/inputs.txt
feat() { o=$OUT/feat/$1.npz; [ -s $o ] || python $P/count_regions.py $2 $REG $o; }
export -f feat
xargs -P 16 -L 1 bash -c 'feat $0 $1' < $OUT/inputs.txt
echo "feature files: $(ls $OUT/feat/*.npz | wc -l) / 142"

echo "== 4. UXM and CelFEER runs"
run() { method=$1; id=$2; f=$3
  if [ $method = uxm ]; then
    o=$OUT/est/uxm@u250@$id.txt; [ -s $o ] && return
    python $P/runtimed.py uxm@u250@$id deconvolution --model uxm --percent_u $REFD/uxm_u250_atlas.tsv -m $f \
      --cpg_coordinates $REF/hg38_pat_cpg.bed.gz -i $REFD/uxm_u250_regions.bed -b --epiformat pat \
      --min_length 4 --u_threshold 0.25 --weights True --outfile $o > $OUT/logs/uxm@u250@$id.log 2>&1
  else
    o=$OUT/est/celfeer@u250@$id.txt; [ -s $o ] && return
    d=$OUT/celfeer_tmp/$id; mkdir -p $d
    python $P/make_celfeer_input.py $OUT/feat/$id.npz $REFD/celfeer_u250_ref.npz $id $d/input.txt
    python $P/runtimed.py celfeer@u250@$id python $CF/scripts/celfeer.py $d/input.txt $d 1 > $OUT/logs/celfeer@u250@$id.log 2>&1 \
      && python $P/celfeer_to_line.py $d/1_tissue_proportions.txt $REFD/celfeer_u250_ref.npz $o && rm -f $d/input.txt
  fi
}
export -f run
{ awk '{print "uxm", $1, $2}' $OUT/inputs.txt; awk '{print "celfeer", $1, $2}' $OUT/inputs.txt; } | xargs -P 16 -L 1 bash -c 'run $0 $1 $2'
echo "estimates: uxm $(ls $OUT/est/uxm@* 2>/dev/null | wc -l)/142  celfeer $(ls $OUT/est/celfeer@* 2>/dev/null | wc -l)/142"

echo "== 5. summary of all methods"
python $P/collect_all.py | tee /scratch/mxf4642/results/summary/summary.txt
echo DONE
