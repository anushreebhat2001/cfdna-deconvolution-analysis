set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export P=/scratch/mxf4642/code/project REF=/scratch/mxf4642/data/ref RUN=/scratch/mxf4642/data/katsman_run
export D=/scratch/mxf4642/data/katsman/epibed_biscuit_filtered_bams_links PAT=/scratch/mxf4642/data/pat
export OUT=/scratch/mxf4642/results/m4
mkdir -p $OUT/est $OUT/logs $PAT
SAMPLES=$(tail -n +2 $P/meta/samples.tsv | cut -f1)

convert() { s=$1
  [ -f $PAT/$s.pat.gz.tbi ] && return
  python $P/runtimed.py convert_$s python $P/epibed2pat.py $D/$s.epibed.gz $REF/hg38_cpg_index.bed.gz $PAT/$s.pat > $OUT/logs/convert_$s.log 2>&1
  sort -k2,2n -S 1G -T $PAT $PAT/$s.pat | bgzip > $PAT/$s.pat.gz && rm $PAT/$s.pat && tabix -f -s 1 -b 2 -e 2 $PAT/$s.pat.gz
}
deconv() { s=$1; ms=$2; mode=$3; name=${ms}_${mode}_$s
  [ -s $OUT/est/$name.txt ] && return
  if [ $ms = u250 ]; then A=$RUN/31CT_U250_Netanel_atlas_over_regions.txt; R=$RUN/31CT_U250_Netanel_merged_regions_file.bed
  else A=$RUN/blood_epithelial_atlas_over_regions_max_markers.txt; R=$RUN/blood_epithelial_merged_regions_file_max_markers.txt; fi
  clip=0.05; [ $mode = noclip ] && clip=0
  python $P/runtimed.py $name deconvolution --model celfie-ish -a $A -m $PAT/$s.pat.gz \
    --cpg_coordinates $REF/hg38_pat_cpg.bed.gz -i $R -b --epiformat pat \
    --num_iterations 10000 --stop_criterion 0.0000001 --random_restarts 1 \
    --read_max_clipping_cell_type_probabilty $clip --outfile $OUT/est/$name.txt > $OUT/logs/$name.log 2>&1
}
export -f convert deconv

echo "== converting"; echo "$SAMPLES" | xargs -P 8 -I{} bash -c 'convert {}'
echo "pat files: $(ls $PAT/*.pat.gz | wc -l)"

echo "== offset check"
for s in HU012.01.PL4406 HU083.64; do python $P/check_shift.py $RUN/31CT_U250_Netanel_atlas_over_regions.txt $PAT/$s.pat.gz; done

echo "== deconvolution"
for s in $SAMPLES; do for ms in u250 max; do for mode in clip noclip; do echo "$s $ms $mode"; done; done; done \
  | xargs -P 16 -L 1 bash -c 'deconv $0 $1 $2'
echo "estimates: $(ls $OUT/est/*.txt | wc -l) / 168"
python $P/collect_m4.py
echo DONE
