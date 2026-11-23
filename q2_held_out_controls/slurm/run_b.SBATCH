#!/bin/bash
#SBATCH --account=csci_ga_3033_135_2026fa
#SBATCH --partition=n2c48m24
#SBATCH --cpus-per-task=24
#SBATCH --mem=22GB
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --job-name=mrd-b
#SBATCH --output=/scratch/%u/mrd/logs/b_%j.out

# Submit with, for the pilot:
#   sbatch --export=ALL,CONDS="b_split_00 b_split_01" run_b.SBATCH
#
# Full setup (b), CelFiE-ISH on max markers:
#   sbatch --export=ALL,CONDS="$(cd /scratch/$USER/mrd/data/pools && \
#       ls b_split_*.txt | sed 's/\.txt$//' | tr '\n' ' ')" run_b.SBATCH
#
# Setup (c):
#   sbatch --export=ALL,CONDS="c_batch_fwd c_batch_rev" run_b.SBATCH
#
# Add methods or marker sets by also exporting METHODS and MARKERS, e.g.
#   --export=ALL,CONDS="...",METHODS="celfieish celfie",MARKERS="max u250"

export CONDS="$(cd /scratch/$USER/mrd/data/pools && ls b_split_*.txt | sed 's/\.txt$//' | tr '\n' ' ')"
singularity exec --overlay /scratch/mxf4642/envs/mrd.ext3:ro \
  /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif \
  /bin/bash -c "source /ext3/env.sh; bash /scratch/$USER/mrd/code/run_b.sh"
