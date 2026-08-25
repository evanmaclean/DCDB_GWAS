#!/bin/bash
#SBATCH --output=log_files/lead_snp_matched_baseline-%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=00:30:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --array=1-5
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files gcta_files/lead_snp_matched_baseline

phen=$(sed -n "${SLURM_ARRAY_TASK_ID}p" gcta_files/pheno_list_lead_snp_check | tr -d '\r')
echo "[$(date)] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} phen='${phen}'"

[[ -n "${phen}" ]] || { echo "ERROR: phen is empty" >&2; exit 1; }
[[ -f "gcta_files/phenos/${phen}.phen" ]] || { echo "ERROR: missing gcta_files/phenos/${phen}.phen" >&2; exit 1; }
[[ -f "gcta_files/qcovar.qcovar" ]] || { echo "ERROR: missing gcta_files/qcovar.qcovar" >&2; exit 1; }
[[ -f "gcta_files/lead_snps_extract.txt" ]] || { echo "ERROR: run Rscripts/lead_snp_stats.R first" >&2; exit 1; }
[[ -x "./gcta_files/gcta" ]] || { echo "ERROR: ./gcta_files/gcta not executable" >&2; exit 1; }

echo "[$(date)] launching GCTA (matched baseline, no ancestry, no LOCO)..."
./gcta_files/gcta \
  --mlma \
  --covar "gcta_files/covar.covar" \
  --qcovar "gcta_files/qcovar.qcovar" \
  --grm "gcta_files/grm/grm" \
  --maf 0.05 \
  --bfile "plink/cleaned" \
  --extract "gcta_files/lead_snps_extract.txt" \
  --pheno "gcta_files/phenos/${phen}.phen" \
  --autosome-num 38 \
  --thread-num "${SLURM_CPUS_PER_TASK}" \
  --out "gcta_files/lead_snp_matched_baseline/${phen}"

echo "[$(date)] GCTA finished"
