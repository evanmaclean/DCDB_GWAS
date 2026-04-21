#!/bin/bash
#SBATCH --output=log_files/loco-%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=10:00:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --array=1-8
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files gcta_files/mlma_loco_results

phen=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pheno_list | tr -d '\r')
echo "[$(date)] SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} phen='${phen}'"

[[ -n "${phen}" ]] || { echo "ERROR: phen is empty" >&2; exit 1; }
[[ -f "gcta_files/phenos/${phen}.phen" ]] || { echo "ERROR: missing gcta_files/phenos/${phen}.phen" >&2; exit 1; }
[[ -x "./gcta_files/gcta" ]] || { echo "ERROR: ./gcta_files/gcta not executable" >&2; ls -l ./gcta_files/gcta; exit 1; }

echo "[$(date)] launching GCTA..."
./gcta_files/gcta \
  --mlma-loco \
  --covar "gcta_files/covar.covar" \
  --qcovar "gcta_files/qcovar.qcovar" \
  --grm "gcta_files/grm/grm" \
  --maf 0.05 \
  --bfile "plink/cleaned" \
  --pheno "gcta_files/phenos/${phen}.phen" \
  --autosome-num 38 \
  --thread-num "${SLURM_CPUS_PER_TASK}" \
  --out "gcta_files/mlma_loco_results/${phen}"

   


echo "[$(date)] GCTA finished"