#!/usr/bin/env bash
#SBATCH --job-name=overall_h2
#SBATCH --output=log_files/overall_h2_%a.out
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=28
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --array=1-8
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

# ---------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------

PHENO_LIST="pheno_list"

# Existing genome-wide GRM prefix
GW_GRM_PREFIX="gcta_files/grm/grm"

# Phenotypes / covariates
COVAR_FILE="gcta_files/covar.covar"
QCOVAR_FILE="gcta_files/qcovar.qcovar"

# Local GCTA binary
GCTA_BIN="gcta_files/gcta"

# Dog autosomes
AUTOSOME_NUM=38

# Output root
OUTROOT="gcta_files/greml_estimates/overall_h2"

mkdir -p log_files
mkdir -p "${OUTROOT}"

# ---------------------------------------------------------
# SANITY CHECKS
# ---------------------------------------------------------

if [[ ! -f "${PHENO_LIST}" ]]; then
  echo "ERROR: Phenotype list not found: ${PHENO_LIST}" >&2
  exit 1
fi

if [[ ! -x "${GCTA_BIN}" ]]; then
  echo "ERROR: GCTA binary not found or not executable: ${GCTA_BIN}" >&2
  exit 1
fi

if [[ ! -f "${GW_GRM_PREFIX}.grm.bin" || ! -f "${GW_GRM_PREFIX}.grm.N.bin" || ! -f "${GW_GRM_PREFIX}.grm.id" ]]; then
  echo "ERROR: Missing genome-wide GRM files for prefix ${GW_GRM_PREFIX}" >&2
  exit 1
fi

# ---------------------------------------------------------
# GET CURRENT PHENOTYPE
# ---------------------------------------------------------

IDX="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID not set}"
PHENO=$(awk -v i="${IDX}" 'NR == i {print; exit}' "${PHENO_LIST}")

if [[ -z "${PHENO}" ]]; then
  echo "ERROR: No phenotype found for array index ${IDX}" >&2
  exit 1
fi

PHENO_FILE="gcta_files/phenos/${PHENO}.phen"
if [[ ! -f "${PHENO_FILE}" ]]; then
  echo "ERROR: Missing phenotype file: ${PHENO_FILE}" >&2
  exit 1
fi

echo "Running overall SNP h2 for phenotype:"
echo "  PHENO = ${PHENO}"
echo "  PWD   = $(pwd)"

# ---------------------------------------------------------
# RUN SINGLE-GRM GREML
# ---------------------------------------------------------

OUTDIR="${OUTROOT}/${PHENO}"
mkdir -p "${OUTDIR}"

REML_OUT="${OUTDIR}/reml"

GCTA_ARGS=(
  --reml
  --grm "${GW_GRM_PREFIX}"
  --autosome-num "${AUTOSOME_NUM}"
  --thread-num "${SLURM_CPUS_PER_TASK}"
  --pheno "${PHENO_FILE}"
  --out "${REML_OUT}"
)

if [[ -f "${COVAR_FILE}" ]]; then
  GCTA_ARGS+=(--covar "${COVAR_FILE}")
fi

if [[ -f "${QCOVAR_FILE}" ]]; then
  GCTA_ARGS+=(--qcovar "${QCOVAR_FILE}")
fi

"${GCTA_BIN}" "${GCTA_ARGS[@]}"

echo "Finished ${PHENO}"
echo "Overall SNP heritability = V(G)/Vp in:"
echo "  ${REML_OUT}.hsq"