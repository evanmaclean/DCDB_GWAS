#!/bin/bash
#SBATCH --output=log_files/n_per_snp.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0:10:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

# Ensure output directories exist
mkdir -p \
  log_files \
  magma_files/snp_sample_sizes \
  magma_files/snp_sample_sizes/missing_phenos

module load plink/1.9

: "${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set. Submit with sbatch --array=...}"

phen=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pheno_list)

# List missing individuals for this phenotype (FID IID where phenotype == NA)
awk '$3 == "NA" { print $1, $2 }' \
  "gcta_files/phenos/${phen}.phen" \
  > "magma_files/snp_sample_sizes/missing_phenos/${phen}"

# Compute per-SNP missingness after removing those individuals
plink \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --dog \
  --bfile plink/cleaned \
  --remove "magma_files/snp_sample_sizes/missing_phenos/${phen}" \
  --missing \
  --out "magma_files/snp_sample_sizes/${phen}"

# Drop per-individual missingness file (we only need .lmiss)
rm -f "magma_files/snp_sample_sizes/${phen}.imiss"