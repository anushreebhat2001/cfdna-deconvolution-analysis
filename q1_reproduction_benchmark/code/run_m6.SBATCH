#!/bin/bash
#SBATCH --account=csci_ga_3033_135_2026fa
#SBATCH --partition=n2c48m24
#SBATCH --cpus-per-task=24
#SBATCH --mem=22GB
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --job-name=m6
#SBATCH --output=/scratch/mxf4642/results/m6_%j.out
singularity exec --overlay /scratch/mxf4642/envs/mrd.ext3:ro \
  /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif \
  /bin/bash -c "source /ext3/env.sh; bash /scratch/mxf4642/code/project/run_m6.sh"
