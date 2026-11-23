#!/bin/bash
#SBATCH --account=csci_ga_3033_135_2026fa
#SBATCH --partition=n2c48m24
#SBATCH --cpus-per-task=24
#SBATCH --mem=22GB
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --job-name=mrd-c
#SBATCH --output=/scratch/%u/mrd/logs/c_%j.out

# Setup (c): batch swap. Mix with one control batch, score against the other,
# then reverse.
#
# Submit with no environment at all:
#   cd /scratch/$USER/mrd/code && unset CONDS METHODS MARKERS KEEP_MIX
#   sbatch run_c.SBATCH
#
# Everything the run needs is set below rather than exported at submit time:
# a long CONDS in the submitting environment gets the job held on this
# cluster with "user env retrieval failed".
#
# Methods: CelFiE-ISH matches setup (b) so the two are comparable; CelFiE is
# added because max markers is where the advisor's answer puts the
# site-level vs read-level comparison (UXM and CelFEER lack a megakaryocyte
# read-level reference on this marker set).
#
# Cost: 2 conditions x 100 mixtures x (2.1 + 3.1) CPU-min ~= 17 CPU-hours,
# roughly 45 min of wall time on 24 cores.

export CONDS="$(cd /scratch/$USER/mrd/data/pools && \
  ls c_batch_*.txt 2>/dev/null | sed 's/\.txt$//' | tr '\n' ' ')"
export METHODS="celfieish celfie"
export MARKERS="max"
export KEEP_MIX=1          # keep mixtures: only 2 conditions, ~1 GB, and the
                           # batch arms are worth spot-checking before deletion

if [ -z "${CONDS// /}" ]; then
  echo "No c_batch_*.txt pool lists found — run make_splits.py first" >&2
  exit 1
fi
echo "[$(date -Is)] submitting conditions: $CONDS"

singularity exec --overlay /scratch/mxf4642/envs/mrd.ext3:ro \
  /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif \
  /bin/bash -c "source /ext3/env.sh; bash /scratch/$USER/mrd/code/run_b.sh"
