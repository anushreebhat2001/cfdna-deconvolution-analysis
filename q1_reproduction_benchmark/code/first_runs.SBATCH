#!/bin/bash
#SBATCH --account=csci_ga_3033_135_2026fa
#SBATCH --partition=n2c48m24
#SBATCH --cpus-per-task=2
#SBATCH --mem=20GB
#SBATCH --time=04:00:00
#SBATCH --job-name=first_runs
#SBATCH --output=/scratch/mxf4642/results/first_runs_%j.out
cd /scratch/mxf4642/data/ref
wget -c -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
rm -f /scratch/mxf4642/data/ref/hg38_cpg_index.bed.gz /scratch/mxf4642/data/ref/hg38_cpg_index.bed.gz.tbi
singularity exec --overlay /scratch/mxf4642/envs/mrd.ext3:ro \
  /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif \
  /bin/bash -c "source /ext3/env.sh; bash /scratch/mxf4642/code/project/first_runs.sh"
