#!/bin/bash
#SBATCH --job-name=ld_loop
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=log_files/ld_%j.out
#SBATCH --error=log_files/ld_%j.err
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files

module load plink/2.0

MANIFEST="locus_data/peak_manifest.tsv"
PLINK_PREFIX="plink/cleaned"
OUTDIR="locus_data/ld_results"

mkdir -p "${OUTDIR}"

while IFS=$'\t' read -r PEAK_ID CHR LEAD BP START END _; do
  [[ -z "${PEAK_ID:-}" ]] && continue
  [[ "${PEAK_ID}" == "PEAK_ID" ]] && continue

  # minimal sanity checks
  [[ -z "${CHR:-}" || -z "${LEAD:-}" || -z "${START:-}" || -z "${END:-}" ]] && {
    echo "Skipping PEAK_ID=${PEAK_ID}: missing required fields" >&2
    continue
  }

  echo "Processing peak ${PEAK_ID} (chr${CHR}:${START}-${END}, lead ${LEAD})"

  plink2 \
    --dog \
    --bfile "${PLINK_PREFIX}" \
    --chr "${CHR}" \
    --from-bp "${START}" \
    --to-bp "${END}" \
    --ld-snp "${LEAD}" \
    --ld-window-kb 5000 \
    --ld-window 999999 \
    --r2-unphased \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --out "${OUTDIR}/peak_${PEAK_ID}"
done < <(tail -n +2 "${MANIFEST}")