#!/bin/bash
#SBATCH --output=log_files/PCA-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=0:30:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files plink

module load plink/1.9

# PCA on pruned set
plink \
  --dog \
  --bfile plink/pruned \
  --pca 1 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --out plink/pca

# PCA on full set (all sites)
plink \
  --dog \
  --bfile plink/plink \
  --pca 1 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --out plink/all_site_pca