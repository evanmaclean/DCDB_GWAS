#!/bin/bash
#SBATCH --output=log_files/froh.out
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
[[ -f "plink/pruned.bed" ]] || { echo "ERROR: missing plink/pruned.bed" >&2; exit 1; }
[[ -f "plink/pruned.bim" ]] || { echo "ERROR: missing plink/pruned.bim" >&2; exit 1; }
[[ -f "plink/pruned.fam" ]] || { echo "ERROR: missing plink/pruned.fam" >&2; exit 1; }

echo "[$(date)] running PLINK --homozyg..."
# marker density ~2.35 kb/SNP (1,022,301 SNPs / ~2.4Gb autosomal genome);
# window-snp/homozyg-snp sized to span ~1.2Mb, matching the 1Mb minimum
# ROH 
plink \
  --bfile plink/pruned \
  --dog \
  --homozyg \
  --homozyg-window-snp 500 \
  --homozyg-snp 500 \
  --homozyg-kb 1000 \
  --homozyg-gap 1000 \
  --homozyg-density 50 \
  --homozyg-window-missing 5 \
  --homozyg-window-het 1 \
  --homozyg-window-threshold 0.05 \
  --out plink/froh

echo "[$(date)] PLINK ROH scan finished"