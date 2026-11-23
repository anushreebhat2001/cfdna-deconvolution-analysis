#!/bin/bash
#SBATCH --account=csci_ga_3033_135_2026fa
#SBATCH --partition=n2c48m24
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --job-name=count_reads
#SBATCH --output=/scratch/mxf4642/results/count_reads_%j.out
D=/scratch/mxf4642/data/katsman/epibed_biscuit_filtered_bams_links
tail -n +2 /scratch/mxf4642/code/project/meta/samples.tsv | \
  xargs -P 8 -L 1 bash -c 'n=$(zcat '"$D"'/$0.epibed.gz | cut -f4 | sort -u -S 1G | wc -l); echo -e "$0\t$1\t$4\t$n"' | sort
