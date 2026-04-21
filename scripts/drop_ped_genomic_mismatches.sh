#!/bin/bash
#SBATCH --output=log_files/drop_ped_genomic_mismatches.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files plink

module load plink/1.9

plink \
  --dog \
  --bfile plink/plink \
  --remove-fam plink/drop_dogs \
  --make-bed \
  --out plink/cleaned