set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export P=/scratch/mxf4642/code/project REF=/scratch/mxf4642/data/ref RUN=/scratch/mxf4642/data/katsman_run
export OUT=/scratch/mxf4642/results/m7 M6=/scratch/mxf4642/results/m6
export C=/scratch/mxf4642/data/katsman_code/code LOY=/mnt/mrd/datasets/extracted CF=/scratch/mxf4642/code/CelFEER
export REFD=/scratch/mxf4642/data/readref SYN=/scratch/mxf4642/data/readref/synth LMAX=/scratch/mxf4642/data/readref/loyfer_max
export ATL_MAX=$RUN/blood_epithelial_atlas_over_regions_max_markers.txt REG_MAX=$RUN/blood_epithelial_merged_regions_file_max_markers.txt
export ATL_U=$RUN/31CT_U250_Netanel_atlas_over_regions.txt REG_U=$RUN/31CT_U250_Netanel_merged_regions_file.bed
mkdir -p $OUT/est $OUT/logs $OUT/feat_max $OUT/celfeer_tmp $SYN $LMAX

echo "== 1. grouping files"
python - << 'PY'
import pandas as pd
C, REFD = "/scratch/mxf4642/data/katsman_code/code", "/scratch/mxf4642/data/readref"
g = pd.read_csv(f"{C}/group_31CT.csv")
EP = {"Bladder-Ep","Breast-Basal-Ep","Breast-Luminal-Ep","Colon-Ep","Fallopian-Ep","Gastric-Ep","Head-Neck-Ep",
      "Kidney-Ep","Lung-Ep-Alveo","Lung-Ep-Bron","Ovary+Endom-Ep","Prostate-Ep","Small-Int-Ep","Thyroid-Ep"}
SAME = {"Blood-Granul","Blood-Mono+Macro","Endothel","Eryth-prog","Liver-Hep"}
def to_max(x):
    if x in SAME: return x
    if x in EP: return "Epithelial"
    if x in ("Blood-B", "Blood-T"): return "Lymphocytes"
    return None
m = g.assign(group=g.group.map(to_max)).dropna()
granul = list(g[g.group == "Blood-Granul"].name)
mk = pd.DataFrame({"name": [f"MKsynth-{n}" for n in granul], "group": "Megakaryocytes"})
pd.concat([m, mk]).to_csv(f"{REFD}/group_max.csv", index=False)
cal = g.copy()
cal.loc[cal.group == "Blood-Granul", "name"] = [f"GranulSynth-{n}" for n in granul]
cal.to_csv(f"{REFD}/group_u250cal.csv", index=False)
print(pd.concat([m, mk]).group.value_counts().to_string())
PY

echo "== 2. count reads in max-marker regions for the Loyfer reference samples"
cnt_max() { n=$1; o=$LMAX/$n.npz; [ -s $o ] && return
  f=$(ls $LOY/*_${n}.hg38.pat.gz | head -1)
  python $P/runtimed.py refcountmax_$n python $P/count_regions.py $f $REG_MAX $o > $OUT/logs/refcountmax_$n.log 2>&1; }
export -f cnt_max
tail -n +2 $REFD/group_max.csv | cut -d, -f1 | grep -v '^MKsynth-' | xargs -P 16 -I{} bash -c 'cnt_max {}'
echo "max-marker reference samples counted: $(ls $LMAX/*.npz | grep -vc MKsynth) / $(tail -n +2 $REFD/group_max.csv | grep -vc '^MKsynth-')"

echo "== 3. synthetic independent-site reads: Megakaryocytes (max markers) and Blood-Granul (U250 calibration)"
synth() { kind=$1; n=$2; seed=$3
  if [ $kind = MK ]; then o=$SYN/MKsynth-$n; atl=$ATL_MAX; col=Megakaryocytes; reg=$REG_MAX; npz=$LMAX/MKsynth-$n.npz
  else o=$SYN/GranulSynth-$n; atl=$ATL_U; col=Blood-Granul; reg=$REG_U; npz=$REFD/loyfer/GranulSynth-$n.npz; fi
  [ -s $npz ] && return
  f=$(ls $LOY/*_${n}.hg38.pat.gz | head -1)
  python $P/runtimed.py synth_${kind}_$n python $P/synth_reads.py $f $atl $col $reg $o.pat $seed > $OUT/logs/synth_${kind}_$n.log 2>&1 \
    && sort -T $SYN -k1,1 -k2,2n $o.pat | bgzip > $o.pat.gz && tabix -f -s 1 -b 2 -e 2 $o.pat.gz && rm -f $o.pat \
    && python $P/count_regions.py $o.pat.gz $reg $npz; }
export -f synth
awk -F, '$2=="Blood-Granul"{print $1}' $C/group_31CT.csv | awk '{print "MK", $1, NR; print "GR", $1, NR+10}' | xargs -P 6 -L 1 bash -c 'synth $0 $1 $2'
grep -h -E '^(atlas|synthetic)' $OUT/logs/synth_*.log
echo "synthetic reference samples: MK $(ls $LMAX/MKsynth-*.npz 2>/dev/null | wc -l)/3  Granul $(ls $REFD/loyfer/GranulSynth-*.npz 2>/dev/null | wc -l)/3"

echo "== 4. build references"
python $P/aggregate_ref.py $REFD/group_max.csv $LMAX $ATL_MAX $REG_MAX $REFD max
python $P/aggregate_ref.py $REFD/group_u250cal.csv $REFD/loyfer $ATL_U $REG_U $REFD u250cal

echo "== 5. count reads in max-marker regions for our 42 samples + 100 mixtures"
feat() { o=$OUT/feat_max/$1.npz; [ -s $o ] || python $P/count_regions.py $2 $REG_MAX $o; }
export -f feat
xargs -P 16 -L 1 bash -c 'feat $0 $1' < $M6/inputs.txt
echo "feature files: $(ls $OUT/feat_max/*.npz | wc -l) / 142"

echo "== 6. UXM and CelFEER: max markers (synthetic MK) and U250 calibration"
run() { method=$1; ms=$2; id=$3; f=$4
  if [ $ms = max ]; then fe=$OUT/feat_max/$id.npz; else fe=$M6/feat/$id.npz; fi
  o=$OUT/est/$method@$ms@$id.txt; [ -s $o ] && return
  if [ $method = uxm ]; then
    python $P/runtimed.py uxm@$ms@$id deconvolution --model uxm --percent_u $REFD/uxm_${ms}_atlas.tsv -m $f \
      --cpg_coordinates $REF/hg38_pat_cpg.bed.gz -i $REFD/uxm_${ms}_regions.bed -b --epiformat pat \
      --min_length 4 --u_threshold 0.25 --weights True --outfile $o > $OUT/logs/uxm@$ms@$id.log 2>&1
  else
    d=$OUT/celfeer_tmp/$ms/$id; mkdir -p $d
    python $P/make_celfeer_input.py $fe $REFD/celfeer_${ms}_ref.npz $id $d/input.txt
    python $P/runtimed.py celfeer@$ms@$id python $CF/scripts/celfeer_tiefix.py $d/input.txt $d 1 > $OUT/logs/celfeer@$ms@$id.log 2>&1 \
      && python $P/celfeer_to_line.py $d/1_tissue_proportions.txt $REFD/celfeer_${ms}_ref.npz $o && rm -f $d/input.txt
  fi; }
export -f run
for ms in max u250cal; do for m in uxm celfeer; do awk -v m=$m -v ms=$ms '{print m, ms, $1, $2}' $M6/inputs.txt; done; done \
  | xargs -P 16 -L 1 bash -c 'run $0 $1 $2 $3'
for ms in max u250cal; do echo "estimates $ms: uxm $(ls $OUT/est/uxm@$ms@* 2>/dev/null | wc -l)/142  celfeer $(ls $OUT/est/celfeer@$ms@* 2>/dev/null | wc -l)/142"; done

echo "== 7. summary of all methods"
python $P/collect_all.py | tee /scratch/mxf4642/results/summary/summary.txt
echo DONE
