#!/bin/bash
#SBATCH --output=log_files/mendel_check.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files

module load plink/1.9

echo "[$(date)] checking inputs..."
[[ -f "plink/cleaned.bed" ]] || { echo "ERROR: missing plink/cleaned.bed" >&2; exit 1; }
[[ -f "plink/cleaned.bim" ]] || { echo "ERROR: missing plink/cleaned.bim" >&2; exit 1; }
[[ -f "plink/mendel_check.fam" ]] || { echo "ERROR: missing plink/mendel_check.fam" >&2; exit 1; }

# link the family-clustered .fam alongside the existing genotype fileset,
# without touching cleaned.fam (still needed elsewhere for GCTA)
[[ -e "plink/mendel_check.bed" ]] || ln -s cleaned.bed plink/mendel_check.bed
[[ -e "plink/mendel_check.bim" ]] || ln -s cleaned.bim plink/mendel_check.bim

echo "[$(date)] running PLINK --mendel..."
plink \
  --bfile plink/mendel_check \
  --dog \
  --mendel \
  --mendel-duos \
  --out plink/mendel_check_out

echo "[$(date)] PLINK Mendel check finished"