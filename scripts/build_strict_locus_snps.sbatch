#!/usr/bin/env bash
#SBATCH --job-name=strict_locus_snps
#SBATCH --output=log_files/strict_locus_snps_%j.out
#SBATCH --error=log_files/strict_locus_snps_%j.err
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files

PROJECT_ROOT="${1:-.}"
ORIG_ROOT="${PROJECT_ROOT}/gcta_files/greml_estimates/locus_greml"
STRICT_ROOT="${PROJECT_ROOT}/gcta_files/greml_estimates/locus_greml_strict_r2"
LEADONLY_ROOT="${PROJECT_ROOT}/gcta_files/greml_estimates/locus_greml_leadonly"
STRICT_R2=0.9

if [[ ! -d "${ORIG_ROOT}" ]]; then
  echo "ERROR: ${ORIG_ROOT} not found. Point argument 1 at the directory containing" >&2
  echo "the original locus_greml.sbatch output." >&2
  exit 1
fi

mkdir -p "${STRICT_ROOT}" "${LEADONLY_ROOT}"

N_LOCI=0
while IFS= read -r LD_FILE; do
  LOCUS_DIR=$(dirname "${LD_FILE}")
  LEAD_SNP=$(basename "${LOCUS_DIR}")
  PHENO=$(basename "$(dirname "${LOCUS_DIR}")")
  LEAD_FILE="${LOCUS_DIR}/lead.txt"

  if [[ ! -f "${LEAD_FILE}" ]]; then
    echo "WARNING: missing lead.txt for ${PHENO}/${LEAD_SNP}, skipping" >&2
    continue
  fi
  LEAD_SNP_ID=$(cat "${LEAD_FILE}")

  # --- Stricter r2 >= 0.9 version: re-filter the EXISTING LD output ---
  STRICT_DIR="${STRICT_ROOT}/${PHENO}/${LEAD_SNP}"
  mkdir -p "${STRICT_DIR}"

  zcat "${LD_FILE}" | awk -v r2min="${STRICT_R2}" '
    NR==1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "SNP_A") a = i;
        if ($i == "SNP_B") b = i;
        if ($i == "R2")    r = i;
      }
      if (a == 0 || b == 0 || r == 0) {
        print "ERROR: could not find SNP_A/SNP_B/R2 columns in " FILENAME > "/dev/stderr";
        exit 2;
      }
      next
    }
    ($r + 0) >= r2min { print $a; print $b }
  ' | sort -u > "${STRICT_DIR}/locus_snps.txt"

  if ! grep -qxF "${LEAD_SNP_ID}" "${STRICT_DIR}/locus_snps.txt" 2>/dev/null; then
    echo "${LEAD_SNP_ID}" >> "${STRICT_DIR}/locus_snps.txt"
    sort -u "${STRICT_DIR}/locus_snps.txt" -o "${STRICT_DIR}/locus_snps.txt"
  fi

  N_STRICT=$(wc -l < "${STRICT_DIR}/locus_snps.txt" | tr -d ' ')

  # --- Lead-SNP-only version: trivial, just the lead SNP itself ---
  LEADONLY_DIR="${LEADONLY_ROOT}/${PHENO}/${LEAD_SNP}"
  mkdir -p "${LEADONLY_DIR}"
  echo "${LEAD_SNP_ID}" > "${LEADONLY_DIR}/locus_snps.txt"

  # Original count, for the comparison table
  N_ORIG=$(zcat "${LD_FILE}" | tail -n +2 | wc -l | tr -d ' ')

  echo "${PHENO}/${LEAD_SNP}: r2>=0.1 pairs=${N_ORIG} -> r2>=${STRICT_R2} SNPs=${N_STRICT} -> lead-only SNPs=1"
  N_LOCI=$((N_LOCI + 1))
done < <(find "${ORIG_ROOT}" -name "ld.ld.gz")

echo ""
echo "Done. Processed ${N_LOCI} trait x locus combinations."
echo "New locus_snps.txt trees:"
echo "  ${STRICT_ROOT}"
echo "  ${LEADONLY_ROOT}"
