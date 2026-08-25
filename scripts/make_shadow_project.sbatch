#!/usr/bin/env bash
#SBATCH --job-name=make_shadow_project
#SBATCH --output=log_files/make_shadow_project_%j.out
#SBATCH --error=log_files/make_shadow_project_%j.err
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files

PROJECT_ROOT="$(cd "${1:-.}" && pwd)"

for VARIANT in strict_r2 leadonly; do
  SHADOW="${PROJECT_ROOT}/shadow_${VARIANT}"
  mkdir -p "${SHADOW}/gcta_files/greml_estimates"

  if [[ ! -L "${SHADOW}/DoGA_data" ]]; then
    ln -s "${PROJECT_ROOT}/DoGA_data" "${SHADOW}/DoGA_data"
  fi

  TARGET="${PROJECT_ROOT}/gcta_files/greml_estimates/locus_greml_${VARIANT}"
  if [[ ! -d "${TARGET}" ]]; then
    echo "ERROR: ${TARGET} not found -- run build_strict_locus_snps.sbatch first." >&2
    exit 1
  fi

  LINK="${SHADOW}/gcta_files/greml_estimates/locus_greml"
  if [[ ! -L "${LINK}" ]]; then
    ln -s "${TARGET}" "${LINK}"
  fi

  echo "Built ${SHADOW} -> DoGA_data (real), locus_greml -> locus_greml_${VARIANT}"
done

echo ""
echo "Done. Now submit run_doga_intersect_variants.sbatch."
