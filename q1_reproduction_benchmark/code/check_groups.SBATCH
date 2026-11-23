#!/bin/bash
#SBATCH --account=csci_ga_3033_135_2026fa
#SBATCH --partition=n2c48m24
#SBATCH --cpus-per-task=24
#SBATCH --mem=22GB
#SBATCH --time=1:00:00
#SBATCH --requeue
#SBATCH --job-name=grp
#SBATCH --output=/scratch/mxf4642/results/grp_%j.out
mkdir -p /scratch/mxf4642/results/m7
singularity exec --overlay /scratch/mxf4642/envs/mrd.ext3:ro \
  /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif \
  /bin/bash -c "source /ext3/env.sh; python /scratch/mxf4642/code/project/check_groups.py \
    /scratch/mxf4642/data/katsman_run/blood_epithelial_atlas_over_regions_max_markers.txt \
    /scratch/mxf4642/data/katsman_code/code/group_31CT.csv \
    /mnt/mrd/datasets/extracted \
    /scratch/mxf4642/results/m7/max_groups_check.tsv"
