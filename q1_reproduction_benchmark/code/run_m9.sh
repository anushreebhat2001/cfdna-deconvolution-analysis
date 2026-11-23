set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export P=/scratch/mxf4642/code/project REF=/scratch/mxf4642/data/ref RUN=/scratch/mxf4642/data/katsman_run
export OUT=/scratch/mxf4642/results/m9 M6=/scratch/mxf4642/results/m6
mkdir -p $OUT/est $OUT/logs
echo "== 0. collapsed marker files"
ls $RUN | grep -E '_(250|1000)\.txt'
for n in 250 1000; do
  a=$(ls $RUN/blood_epithelial_atlas_over_regions*_$n.txt 2>/dev/null | head -1)
  r=$(ls $RUN/blood_epithelial_merged_regions_file*_$n.txt 2>/dev/null | head -1)
  echo "c$n atlas=$a regions=$r"
  if [ -z "$a" ] || [ -z "$r" ]; then echo "MISSING marker files for $n - stopping"; exit 1; fi
  export "ATL_$n=$a" "REG_$n=$r"
done
echo "== 1. CelFiE-ISH (clipped + unclipped from one fit) on collapsed 250 and 1000 markers"
ish() { n=$1; id=$2; f=$3; av=ATL_$n; rv=REG_$n
  oc=$OUT/est/celfieish-clip@c$n@$id.txt; on=$OUT/est/celfieish-noclip@c$n@$id.txt
  [ -s $oc ] && [ -s $on ] && return
  python $P/runtimed.py celfieish@c$n@$id python $P/celfieish_both.py $on $oc 0 --model celfie-ish \
    -a ${!av} -m $f --cpg_coordinates $REF/hg38_pat_cpg.bed.gz -i ${!rv} -b --epiformat pat \
    --num_iterations 10000 --stop_criterion 0.0000001 --random_restarts 1 --outfile /dev/null \
    > $OUT/logs/celfieish@c$n@$id.log 2>&1; }
export -f ish
echo "-- test run: HU030.01.PL3393, 250 markers"
ish 250 HU030.01.PL3393 /scratch/mxf4642/data/pat/HU030.01.PL3393.pat.gz
tail -3 $OUT/logs/celfieish@c250@HU030.01.PL3393.log
python $P/summarise.py $ATL_250 $OUT/est/celfieish-clip@c250@HU030.01.PL3393.txt
echo "-- full sweep"
for n in 250 1000; do awk -v n=$n '{print n, $1, $2}' $M6/inputs.txt; done | xargs -P 16 -L 1 bash -c 'ish $0 $1 $2'
for n in 250 1000; do echo "c$n: clip $(ls $OUT/est/celfieish-clip@c$n@* 2>/dev/null | wc -l)/142  noclip $(ls $OUT/est/celfieish-noclip@c$n@* 2>/dev/null | wc -l)/142"; done
echo "== 2. summary of all methods"
python $P/collect_all.py | tee /scratch/mxf4642/results/summary/summary.txt
echo DONE
