#!/bin/bash
#SBATCH --job-name=mrd-splitb
#SBATCH --account=csci_ga_3033_135_2026fa
#SBATCH --partition=n2c48m24
#SBATCH --array=1-40%16
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/%u/mrd/logs/%x_%A_%a.out
#SBATCH --error=/scratch/%u/mrd/logs/%x_%A_%a.err

# One array task = one (condition, cancer sample) unit:
#   build the 5 dilutions from that condition's pool
#   run each method on each mixture
#   write result rows
#   DELETE the mixtures
#
# Mixtures are regenerable from (pool, sample, seed), so keeping 2,000 of
# them wastes hundreds of GB of scratch for no benefit. Only the estimates
# and the mixinfo manifests are kept.
#
# Submit with:  sbatch run_split_b.SBATCH
# Adjust --array to match the number of lines in tasks.tsv.
# The %16 caps concurrent tasks; raise it once the pilot is clean.

set -euo pipefail

MATE=/scratch/mxf4642
WORK=/scratch/$USER/mrd
TASKS=$WORK/data/tasks.tsv          # built by make_tasks.sh, see bottom
TMP=$WORK/scratch_mix/$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID

# ---- pick this task's work unit --------------------------------------
# tasks.tsv columns: condition<TAB>cancer_sample
read -r CONDITION SAMPLE < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKS")
if [[ -z "${CONDITION:-}" ]]; then
  echo "No task at line $SLURM_ARRAY_TASK_ID of $TASKS" >&2
  exit 1
fi
echo "[$(date -Is)] task=$SLURM_ARRAY_TASK_ID condition=$CONDITION sample=$SAMPLE"

OUT=$WORK/results/b/$CONDITION/$SAMPLE
mkdir -p "$OUT" "$TMP"
trap 'rm -rf "$TMP"' EXIT        # mixtures go away even if the job dies

# ---- the pool for this condition -------------------------------------
# splits.tsv columns: condition role sample note
awk -F'\t' -v c="$CONDITION" '$1==c && $2=="pool" {print $3}' \
    "$WORK/data/splits.tsv" > "$TMP/pool.txt"
echo "pool size: $(wc -l < "$TMP/pool.txt")"

# ---- run everything inside the container ------------------------------
# REPLACE the base image path with the one from your teammate's
# run_m5.SBATCH, and REPLACE the activation line if it differs.
singularity exec \
  --overlay "$MATE/envs/mrd.ext3:ro" \
  /scratch/work/public/singularity/<BASE-IMAGE>.sif \
  /bin/bash -eu -o pipefail -c "
    source /ext3/env.sh
    conda activate mrd

    cd $WORK

    # 1. Build the five dilutions for this sample from this pool.
    #    Your teammate's mix_dilutions.py needs a --pool argument added;
    #    it currently hardcodes the paper's 20 controls.
    python $WORK/code/mix_dilutions.py \
        --sample $SAMPLE \
        --pat-dir $MATE/data/pat \
        --pool-list $TMP/pool.txt \
        --out-dir $TMP \
        --seed \$(python -c \"print(abs(hash(('$CONDITION','$SAMPLE'))) % 10**8)\")

    # 2. Deconvolve each mixture. Max markers only for the pilot.
    for MIX in $TMP/*.pat.gz; do
      BASE=\$(basename \"\$MIX\" .pat.gz)
      python $WORK/code/runtimed.py --tag \"\$BASE\" -- \
        <CELFIE-ISH INVOCATION FROM run_m5.SBATCH, section 5.3 settings> \
          --epiformat pat \
          --cpg_coordinates $MATE/data/ref/hg38_pat_cpg.bed.gz \
          -i $MATE/data/katsman_run/blood_epithelial_merged_regions_file_max_markers.txt -b \
          --num_iterations 10000 --stop_criterion 1e-7 \
          --input \"\$MIX\" \
          --out $OUT/\$BASE.celfieish.txt
    done
  "

# ---- keep the manifests, drop the mixtures ----------------------------
cp "$TMP"/*.mixinfo.tsv "$OUT"/ 2>/dev/null || true
echo "[$(date -Is)] done; wrote $(ls "$OUT" | wc -l) files"


# =======================================================================
# Build tasks.tsv first (run this once, on the login/Jupyter terminal):
#
#   cd /scratch/$USER/mrd
#   awk -F'\t' '$2=="pool" && $1 ~ /^b_split_(00|01)$/ {print $1}' \
#       data/splits.tsv | sort -u > /tmp/conds.txt
#   awk -F'\t' 'NR>1 && $2=="crc" {print $1}' \
#       /scratch/mxf4642/code/project/meta/samples.tsv > /tmp/crc.txt
#   while read -r c; do while read -r s; do printf '%s\t%s\n' "$c" "$s";
#       done < /tmp/crc.txt; done < /tmp/conds.txt > data/tasks.tsv
#   wc -l data/tasks.tsv     # <- set --array=1-N to this
#
# (adjust the samples.tsv column test to whatever marks a CRC row)
# =======================================================================
