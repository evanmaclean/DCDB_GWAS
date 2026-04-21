#!/bin/bash
#SBATCH --output=log_files/make_grm-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=3:00:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files gcta_files/grm

./gcta_files/gcta \
  --bfile plink/cleaned \
  --dog \
  --autosome-num 38 \
  --make-grm-bin \
  --thread-num "${SLURM_CPUS_PER_TASK}" \
  --out gcta_files/grm/grm