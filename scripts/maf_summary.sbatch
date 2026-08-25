#!/usr/bin/env bash
#SBATCH --job-name=maf_summary
#SBATCH --output=log_files/maf_summary.out
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail
module load plink/1.9

BFILE="plink/cleaned"
AUTOSOME_NUM=38
OUTDIR="gcta_files/qc"
OUT="${OUTDIR}/maf_summary"

mkdir -p log_files "${OUTDIR}"

# ---------------------------------------------------------
# SANITY CHECKS
# ---------------------------------------------------------

if [[ ! -f "${BFILE}.bed" ]]; then
  echo "ERROR: Missing genotype file: ${BFILE}.bed" >&2
  exit 1
fi

# ---------------------------------------------------------
# ALLELE FREQUENCY / MAF DISTRIBUTION
# ---------------------------------------------------------
# --freq writes a .frq file: CHR, SNP, A1, A2, MAF, NCHROBS (one row per
# variant). This is genome-wide across all tested SNPs, not restricted to
# the lead/significant loci - the point is to characterize how much
# standing variation exists on average in this population, addressing R2's
# "low genetic variation at causal loci" alternative explanation (we can't
# know the true causal loci, but we can show this population isn't
# depleted of variation in aggregate).

plink \
  --dog \
  --bfile "${BFILE}" \
  --freq \
  --out "${OUT}"

echo "Finished MAF summary."
echo "Allele-frequency table: ${OUT}.frq"
echo "Columns: CHR SNP A1 A2 MAF NCHROBS"
