set -e
P=/scratch/mxf4642/code/project
REF=/scratch/mxf4642/data/ref
RUN=/scratch/mxf4642/data/katsman_run
D=/scratch/mxf4642/data/katsman/epibed_biscuit_filtered_bams_links
PAT=/scratch/mxf4642/data/pat
OUT=/scratch/mxf4642/results/first_runs

echo "== 1. CpG-number table from plain hg38"
python $P/build_cpg_index.py $REF/hg38.fa.gz $RUN/hg38_pat_cpg_from_netanel.bed.gz $REF/hg38_cpg_index.bed | tail -4
grep -q . $REF/hg38_cpg_index.bed
bgzip -f $REF/hg38_cpg_index.bed && tabix -f -p bed $REF/hg38_cpg_index.bed.gz

echo "== 2. index the authors' CpG-number file (needed by --cpg_coordinates)"
if [ ! -f $REF/hg38_pat_cpg.bed.gz.tbi ]; then
  zcat $RUN/hg38_pat_cpg_from_netanel.bed.gz | bgzip > $REF/hg38_pat_cpg.bed.gz
  tabix -f -s 1 -b 2 -e 2 $REF/hg38_pat_cpg.bed.gz
fi

echo "== 3. convert two samples"
for s in HU012.01.PL4406 HU030.01.PL3393; do
  python $P/runtimed.py convert_$s python $P/epibed2pat.py $D/$s.epibed.gz $REF/hg38_cpg_index.bed.gz $PAT/$s.pat
  sort -k2,2n -S 4G $PAT/$s.pat | bgzip > $PAT/$s.pat.gz && rm $PAT/$s.pat
  tabix -f -s 1 -b 2 -e 2 $PAT/$s.pat.gz
  echo "   $s: $(zcat $PAT/$s.pat.gz | wc -l) pat lines"; zcat $PAT/$s.pat.gz | head -3
done

echo "== 4. CelFiE-ISH on U250 (31 cell types), with and without clipping"
A=$RUN/31CT_U250_Netanel_atlas_over_regions.txt
R=$RUN/31CT_U250_Netanel_merged_regions_file.bed
for s in HU012.01.PL4406 HU030.01.PL3393; do
  for mode in noclip clip; do
    extra=""; [ $mode = clip ] && extra="--read_max_clipping_cell_type_probabilty 0.05"
    python $P/runtimed.py celfieish_${mode}_$s deconvolution --model celfie-ish -a $A -m $PAT/$s.pat.gz \
      --cpg_coordinates $REF/hg38_pat_cpg.bed.gz -i $R -b --epiformat pat \
      --num_iterations 10000 --stop_criterion 0.0000001 --random_restarts 1 $extra \
      --outfile $OUT/celfieish_${mode}_$s.txt
    python $P/summarise.py $A $OUT/celfieish_${mode}_$s.txt || true
  done
done
echo DONE
