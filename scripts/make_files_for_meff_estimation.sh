#!/bin/bash
#SBATCH --output=log_files/thin_for_meff_estimation.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=5:00:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files plink meff/simpleM_mats

module load plink/1.9

# QC + remove missingness (and fill remaining)
plink \
  --dog \
  --bfile plink/cleaned \
  --maf 0.05 \
  --geno 0.05 \
  --fill-missing-a2 \
  --make-bed \
  --out plink/qc_nomiss_for_meff

module load plink/2.0

# LD pruning
plink2 \
  --dog \
  --bfile plink/qc_nomiss_for_meff \
  --indep-pairwise 50kb 1 0.5 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --out plink/pruned_50kb_r05

# Extract pruned set
plink2 \
  --dog \
  --bfile plink/qc_nomiss_for_meff \
  --extract plink/pruned_50kb_r05.prune.in \
  --make-bed \
  --out plink/qc_nomiss_pruned

module load plink/1.9

# Build per-chromosome matrices for simpleM
for CHR in $(seq 1 38); do
  plink \
    --dog \
    --bfile plink/qc_nomiss_pruned \
    --chr "${CHR}" \
    --recode A-transpose \
    --out meff/simpleM_mats/chr${CHR}

  tail -n +2 "meff/simpleM_mats/chr${CHR}.traw" \
    | cut -f7- \
    > "meff/simpleM_mats/chr${CHR}.matrix.txt"
done