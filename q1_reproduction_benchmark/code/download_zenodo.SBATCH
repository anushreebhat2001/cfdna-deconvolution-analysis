#!/bin/bash
#SBATCH --account=csci_ga_3033_135_2026fa
#SBATCH --partition=n2c48m24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --job-name=zenodo_dl
#SBATCH --output=/scratch/mxf4642/data/zenodo_dl_%j.out
#SBATCH --error=/scratch/mxf4642/data/zenodo_dl_%j.err

cd /scratch/mxf4642/data
wget -c --tries=20 --timeout=60 -O epibed_files.zip "https://zenodo.org/records/18100426/files/epibed_files.zip?download=1"
ls -lh epibed_files.zip
unzip -o -q epibed_files.zip -d /scratch/mxf4642/data/katsman
du -sh /scratch/mxf4642/data/katsman
find /scratch/mxf4642/data/katsman -maxdepth 2 | head -60
echo DONE
