#!/bin/bash
#SBATCH --output=log_files/cv_blup_ancestry_adjusted-%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=0:25:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS
#SBATCH --array=1-2

set -euo pipefail

mkdir -p log_files gcta_files/cv_blup_ancestry_adjusted

phen=$(sed -n "${SLURM_ARRAY_TASK_ID}p" gcta_files/pheno_list_ancestry_check)

[[ -f "gcta_files/qcovar_with_ancestry.qcovar" ]] || { echo "ERROR: missing gcta_files/qcovar_with_ancestry.qcovar" >&2; exit 1; }

# REML with cvBLUPs -- identical to cv_blup.sbatch, ancestry-augmented qcovar swapped in
./gcta_files/gcta \
  --reml \
  --grm gcta_files/grm/grm \
  --pheno "gcta_files/phenos/${phen}.phen" \
  --covar gcta_files/covar.covar \
  --qcovar gcta_files/qcovar_with_ancestry.qcovar \
  --autosome-num 38 \
  --cvblup \
  --thread-num "${SLURM_CPUS_PER_TASK}" \
  --out "gcta_files/cv_blup_ancestry_adjusted/${phen}"

# SNP BLUPs
./gcta_files/gcta \
  --bfile plink/cleaned \
  --autosome-num 38 \
  --thread-num "${SLURM_CPUS_PER_TASK}" \
  --blup-snp "gcta_files/cv_blup_ancestry_adjusted/${phen}.indi.cvblp" \
  --out "gcta_files/cv_blup_ancestry_adjusted/${phen}"

# Score individuals with SNP BLUPs
module load plink/1.9

plink \
  --bfile plink/cleaned \
  --autosome-num 38 \
  --score "gcta_files/cv_blup_ancestry_adjusted/${phen}.snp.blp" 1 2 3 \
  --out "gcta_files/cv_blup_ancestry_adjusted/${phen}"