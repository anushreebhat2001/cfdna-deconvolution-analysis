set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Setup (b)/(c): build mixtures from a restricted control pool, deconvolve,
# then throw the mixtures away. Mirrors run_m5.sh; the only real differences
# are that the pool comes from a per-condition list and that each condition
# gets its own mixture and estimate directories.
#
# Driven by environment variables so the same script runs the pilot and the
# full sweep:
#   CONDS    space-separated condition names (default: two-split pilot)
#   MARKERS  "max", "u250", or "max u250"            (default: max)
#   METHODS  "celfieish", "celfie", or both          (default: celfieish)
#   NPROC    parallel deconvolution jobs             (default: 14)
#   KEEP_MIX 1 to keep mixtures on disk              (default: 0, delete)
#
# Reruns are safe: finished estimates are skipped and finished conditions
# are detected by their .done marker.

MATE=/scratch/mxf4642
WORK=/scratch/$USER/mrd
MP=$MATE/code/project                      # teammate's scripts, used as-is
REF=$MATE/data/ref
RUN=$MATE/data/katsman_run
PAT=$MATE/data/pat
META=$MP/meta/samples.tsv
TOTALS=$MP/meta/read_totals.tsv
MIXER=$WORK/code/mix_dilutions.py          # our patched copy
POOLS=$WORK/data/pools
MIXROOT=$WORK/data/mix_b
OUT=$WORK/results/b

CONDS=${CONDS:-"b_split_00 b_split_01"}
MARKERS=${MARKERS:-max}
METHODS=${METHODS:-celfieish}
NPROC=${NPROC:-14}
NPROC_MIX=${NPROC_MIX:-6}
KEEP_MIX=${KEEP_MIX:-0}

mkdir -p $MIXROOT $OUT/est $OUT/logs $WORK/logs
echo "[$(date -Is)] conditions: $CONDS"
echo "[$(date -Is)] markers: $MARKERS   methods: $METHODS   nproc: $NPROC"

# ---- sanity checks before anything expensive -------------------------------
for f in $META $TOTALS $MIXER; do
  [ -s "$f" ] || { echo "MISSING: $f" >&2; exit 1; }
done
for c in $CONDS; do
  [ -s "$POOLS/$c.txt" ] || { echo "MISSING pool list: $POOLS/$c.txt" >&2; exit 1; }
done

CRC=$(awk -F'\t' 'NR>1 && $2=="crc"{print $1}' $META)
NCRC=$(echo "$CRC" | wc -w)
echo "[$(date -Is)] $NCRC CRC samples"

# ---- one mixture set per cancer sample -------------------------------------
# $1 = sample. Reads COND, MIX, POOL from the environment.
mix() {
  s=$1
  [ -f $MIX/$s.done ] && return 0
  # Deterministic per (condition, sample). cksum rather than Python's hash(),
  # which is salted per process and would not reproduce.
  seed=$(echo "$COND $s" | cksum | cut -d' ' -f1)
  python $MP/runtimed.py mix_${COND}_$s \
    python $MIXER $META $TOTALS $POOL $PAT $MIX $s $seed $POOLS/$COND.txt \
    > $OUT/logs/mix_${COND}_$s.log 2>&1 || {
      echo "  MIX FAILED $COND $s (see $OUT/logs/mix_${COND}_$s.log)" >&2
      return 1; }
  for f in $MIX/${s}__d*.pat; do
    sort -k2,2n -S 500M -T $MIX "$f" | bgzip > "$f.gz" && rm "$f" \
      && tabix -f -s 1 -b 2 -e 2 "$f.gz"
  done
  touch $MIX/$s.done
}
export -f mix

# ---- one deconvolution ------------------------------------------------------
# $1 = method, $2 = marker set, $3 = estimate id, $4 = mixture file
run() {
  method=$1; ms=$2; id=$3; f=$4
  if [ $ms = u250 ]; then
    A=$RUN/31CT_U250_Netanel_atlas_over_regions.txt
    R=$RUN/31CT_U250_Netanel_merged_regions_file.bed
  else
    A=$RUN/blood_epithelial_atlas_over_regions_max_markers.txt
    R=$RUN/blood_epithelial_merged_regions_file_max_markers.txt
  fi
  COMMON="-a $A -m $f --cpg_coordinates $REF/hg38_pat_cpg.bed.gz -i $R -b \
          --epiformat pat --num_iterations 10000 --stop_criterion 0.0000001 \
          --random_restarts 1"
  if [ $method = celfieish ]; then
    o=$OUT/est/celfieish-clip@$ms@$id.txt
    [ -s $o ] && return 0
    python $MP/runtimed.py celfieish@$ms@$id \
      python $MP/celfieish_both.py $OUT/est/celfieish-noclip@$ms@$id.txt $o 0 \
      --model celfie-ish $COMMON --outfile /dev/null \
      > $OUT/logs/celfieish@$ms@$id.log 2>&1
  else
    o=$OUT/est/celfie@$ms@$id.txt
    [ -s $o ] && return 0
    python $MP/runtimed.py celfie@$ms@$id \
      deconvolution --model celfie $COMMON --outfile $o \
      > $OUT/logs/celfie@$ms@$id.log 2>&1
  fi
}
export -f run
export MP REF RUN PAT META TOTALS MIXER POOLS OUT

# ---- main loop --------------------------------------------------------------
for COND in $CONDS; do
  export COND
  export MIX=$MIXROOT/$COND
  export POOL=$MIX/healthy_pool.pat
  mkdir -p $MIX

  if [ -f $OUT/$COND.done ]; then
    echo "[$(date -Is)] $COND already complete, skipping"
    continue
  fi

  NPOOL=$(wc -l < $POOLS/$COND.txt)
  echo "[$(date -Is)] === $COND  (pool of $NPOOL controls)"

  # Build the pool by concatenating in the order the list gives, which
  # make_splits.py sorts. The mixer indexes into this file by line number,
  # so a different order with the same seed would select different reads.
  if [ ! -s $POOL.done ]; then
    while read -r s; do
      [ -s $PAT/$s.pat.gz ] || { echo "  MISSING pat: $s" >&2; exit 1; }
      zcat $PAT/$s.pat.gz
    done < $POOLS/$COND.txt > $POOL
    wc -l < $POOL > $POOL.done
  fi
  echo "[$(date -Is)]   pool lines: $(cat $POOL.done)"

  # Mixtures
  echo "$CRC" | tr ' ' '\n' | xargs -P $NPROC_MIX -I{} bash -c 'mix {}'
  NMIX=$(ls $MIX/*__d*.pat.gz 2>/dev/null | wc -l)
  echo "[$(date -Is)]   mixtures: $NMIX / $((NCRC * 5))"
  if [ "$NMIX" -ne "$((NCRC * 5))" ]; then
    echo "  INCOMPLETE mixtures for $COND, not deleting; stopping." >&2
    exit 1
  fi

  # Deconvolution. The id carries the condition so estimates never collide.
  {
    for m in $MIX/*__d*.pat.gz; do
      id=${COND}__$(basename $m .pat.gz)
      for ms in $MARKERS; do
        for meth in $METHODS; do echo "$meth $ms $id $m"; done
      done
    done
  } | xargs -P $NPROC -L 1 bash -c 'run $0 $1 $2 $3'

  NEST=$(ls $OUT/est/*@${COND}__*.txt 2>/dev/null | wc -l)
  echo "[$(date -Is)]   estimate files for $COND: $NEST"

  # Keep the manifests (they record seed and pool membership, so any mixture
  # can be regenerated), drop the mixtures themselves.
  # Copy straight to the final name. An earlier version copied first and
  # renamed in a second pass, whose glob also caught previous conditions'
  # files and gave them doubled prefixes.
  mkdir -p $OUT/mixinfo
  for t in $MIX/*.mixinfo.tsv; do
    [ -e "$t" ] || continue
    cp -f "$t" "$OUT/mixinfo/${COND}__$(basename "$t")"
  done

  if [ "$KEEP_MIX" != "1" ]; then
    rm -rf $MIX
    echo "[$(date -Is)]   deleted mixtures for $COND"
  fi
  touch $OUT/$COND.done
done

echo "[$(date -Is)] DONE. estimates: $(ls $OUT/est/*.txt 2>/dev/null | wc -l)"
echo "Disk in use under $WORK/data/mix_b: $(du -sh $MIXROOT 2>/dev/null | cut -f1)"
