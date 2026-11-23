#!/bin/bash
#SBATCH --account=csci_ga_3033_135_2026fa
#SBATCH --partition=n2c48m24
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --job-name=cpg_index
#SBATCH --output=/scratch/mxf4642/results/cpg_index_%j.out
cd /scratch/mxf4642/data/ref
wget -c -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz
ls -lh hg38.analysisSet.fa.gz
singularity exec --overlay /scratch/mxf4642/envs/mrd.ext3:ro \
  /share/apps/images/cuda12.1.1-cudnn8.9.0-devel-ubuntu22.04.2.sif /bin/bash -c "source /ext3/env.sh; \
  python /scratch/mxf4642/code/project/build_cpg_index.py \
    /scratch/mxf4642/data/ref/hg38.analysisSet.fa.gz \
    /scratch/mxf4642/data/katsman_run/hg38_pat_cpg_from_netanel.bed.gz \
    /scratch/mxf4642/data/ref/hg38_cpg_index.bed && \
  bgzip -f /scratch/mxf4642/data/ref/hg38_cpg_index.bed && \
  tabix -f -p bed /scratch/mxf4642/data/ref/hg38_cpg_index.bed.gz"
ls -lh /scratch/mxf4642/data/ref/
echo DONE
