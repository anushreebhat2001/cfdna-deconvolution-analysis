#!/bin/bash
#SBATCH --account=csci_ga_3033_135_2026fa
#SBATCH --partition=n2c48m24
#SBATCH --cpus-per-task=24
#SBATCH --mem=22GB
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --job-name=chem
#SBATCH --output=/scratch/mxf4642/results/chem_%j.out
mkdir -p /scratch/mxf4642/results/chem
singularity exec --overlay /scratch/mxf4642/envs/mrd.ext3:ro \
  /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif \
  /bin/bash -c "source /ext3/env.sh; python /scratch/mxf4642/code/project/chem_check.py \
    /scratch/mxf4642/data/katsman/epibed_biscuit_filtered_bams_links \
    /scratch/mxf4642/code/project/meta/samples.tsv \
    /scratch/mxf4642/results/chem/chem_check.tsv 200000 16"
