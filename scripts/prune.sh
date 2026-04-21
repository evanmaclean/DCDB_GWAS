#!/bin/bash
#SBATCH --output=log_files/prune_snps.out
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=0:15:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files plink

module load plink

plink \
  --dog \
  --bfile plink/cleaned \
  --indep-pairwise 1500 5 0.95 \
  --out plink/prune

plink \
  --dog \
  --bfile plink/cleaned \
  --extract plink/prune.prune.in \
  --make-bed \
  --out plink/pruned