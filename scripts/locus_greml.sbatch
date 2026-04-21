#!/usr/bin/env bash
#SBATCH --job-name=locus_greml
#SBATCH --output=log_files/locus_greml_%A_%a.out
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=28
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --array=1-9
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

# ---------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------

# CSV with header, columns: SNP, pheno
MANIFEST="gcta_files/sig_snp_manifest.csv"

# Project paths
BFILE_PREFIX="plink/cleaned"
COVAR_FILE="gcta_files/covar.covar"
QCOVAR_FILE="gcta_files/qcovar.qcovar"

module load plink/1.9
PLINK_BIN="plink"

# Local GCTA binary
GCTA_BIN="gcta_files/gcta"

# Dog autosomes
AUTOSOME_NUM=38

# LD parameters
LD_R2=0.1
WINDOW_KB=1000

# Output root
OUTROOT="gcta_files/greml_estimates/locus_greml"

mkdir -p log_files
mkdir -p "${OUTROOT}"

# ---------------------------------------------------------
# SANITY CHECKS
# ---------------------------------------------------------
if [[ ! -f "${MANIFEST}" ]]; then
  echo "ERROR: Manifest not found: ${MANIFEST}" >&2
  echo "PWD: $(pwd)" >&2
  exit 1
fi

if [[ ! -x "${GCTA_BIN}" ]]; then
  echo "ERROR: GCTA binary not found or not executable: ${GCTA_BIN}" >&2
  exit 1
fi

if [[ ! -f "${BFILE_PREFIX}.bed" || ! -f "${BFILE_PREFIX}.bim" || ! -f "${BFILE_PREFIX}.fam" ]]; then
  echo "ERROR: Missing PLINK files for prefix ${BFILE_PREFIX} (.bed/.bim/.fam)" >&2
  exit 1
fi

# ---------------------------------------------------------
# GET CURRENT ROW FROM MANIFEST
# ---------------------------------------------------------

IDX="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID not set}"

# Read row IDX+1 (skip header). CSV assumed with no embedded commas in fields.
LINE=$(awk -v i=$((IDX+1)) 'NR==i{print; exit}' "${MANIFEST}")

if [[ -z "${LINE}" ]]; then
  echo "ERROR: No row found for array index ${IDX} in ${MANIFEST}" >&2
  exit 1
fi

# Parse CSV columns (SNP,pheno)
LEAD_SNP=$(echo "${LINE}" | awk -F',' '{gsub(/^"|"$/, "", $1); print $1}')
PHENO=$(echo "${LINE}"   | awk -F',' '{gsub(/^"|"$/, "", $2); print $2}')

if [[ -z "${LEAD_SNP}" || -z "${PHENO}" ]]; then
  echo "ERROR: Failed to parse SNP/pheno from line: ${LINE}" >&2
  exit 1
fi

PHENO_FILE="gcta_files/phenos/${PHENO}.phen"
if [[ ! -f "${PHENO_FILE}" ]]; then
  echo "ERROR: Missing phenotype file: ${PHENO_FILE}" >&2
  exit 1
fi

echo "Running locus GREML (exact complement) for:"
echo "  SNP   = ${LEAD_SNP}"
echo "  PHENO = ${PHENO}"
echo "  PWD   = $(pwd)"

# ---------------------------------------------------------
# PREP OUTPUT DIRECTORY
# ---------------------------------------------------------

OUTDIR="${OUTROOT}/${PHENO}/${LEAD_SNP}"
mkdir -p "${OUTDIR}"

# ---------------------------------------------------------
# 1) LD EXPANSION -> locus_snps.txt
# ---------------------------------------------------------

echo "${LEAD_SNP}" > "${OUTDIR}/lead.txt"

"${PLINK_BIN}" \
  --dog \
  --bfile "${BFILE_PREFIX}" \
  --ld-snp-list "${OUTDIR}/lead.txt" \
  --r2 gz \
  --ld-window-kb "${WINDOW_KB}" \
  --ld-window 999999 \
  --ld-window-r2 "${LD_R2}" \
  --out "${OUTDIR}/ld"

LD_FILE="${OUTDIR}/ld.ld.gz"
if [[ ! -f "${LD_FILE}" ]]; then
  echo "ERROR: LD output not found: ${LD_FILE}" >&2
  exit 1
fi

# Extract SNPs from LD file (expects SNP_A and SNP_B columns)
zcat "${LD_FILE}" \
  | awk '
      NR==1{
        for(i=1;i<=NF;i++){
          if($i=="SNP_A") a=i;
          if($i=="SNP_B") b=i;
        }
        if(a==0 || b==0){
          print "ERROR: Could not find SNP_A/SNP_B columns." > "/dev/stderr";
          exit 2;
        }
        next
      }
      { print $a; print $b }
    ' \
  | awk 'NF>0' \
  | sort -u > "${OUTDIR}/locus_snps.txt"

# Ensure lead SNP included
if ! grep -qxF "${LEAD_SNP}" "${OUTDIR}/locus_snps.txt"; then
  echo "${LEAD_SNP}" >> "${OUTDIR}/locus_snps.txt"
  sort -u "${OUTDIR}/locus_snps.txt" -o "${OUTDIR}/locus_snps.txt"
fi

N_LOCI_SNPS=$(wc -l < "${OUTDIR}/locus_snps.txt" | tr -d ' ')
echo "Locus SNP count: ${N_LOCI_SNPS}"

# ---------------------------------------------------------
# 2) Build EXACT COMPLEMENT SNP list (RoG) and GRMs
#    - Locus GRM: SNPs in locus_snps.txt
#    - RoG GRM: all SNPs EXCEPT locus_snps.txt
# ---------------------------------------------------------

# Get all SNP IDs from BIM
cut -f2 "${BFILE_PREFIX}.bim" | sort -u > "${OUTDIR}/all_snps.txt"

# RoG SNPs = all minus locus
# comm requires sorted inputs
sort -u "${OUTDIR}/locus_snps.txt" -o "${OUTDIR}/locus_snps.txt"
comm -23 "${OUTDIR}/all_snps.txt" "${OUTDIR}/locus_snps.txt" > "${OUTDIR}/rog_snps.txt"

N_ROG_SNPS=$(wc -l < "${OUTDIR}/rog_snps.txt" | tr -d ' ')
echo "RoG SNP count: ${N_ROG_SNPS}"

LOCUS_GRM_PREFIX="${OUTDIR}/grm_locus"
ROG_GRM_PREFIX="${OUTDIR}/grm_rog"

# Build locus GRM
"${GCTA_BIN}" \
  --bfile "${BFILE_PREFIX}" \
  --autosome-num "${AUTOSOME_NUM}" \
  --extract "${OUTDIR}/locus_snps.txt" \
  --make-grm \
  --thread-num "${SLURM_CPUS_PER_TASK}" \
  --out "${LOCUS_GRM_PREFIX}"

# Build RoG GRM (exact complement)
"${GCTA_BIN}" \
  --bfile "${BFILE_PREFIX}" \
  --autosome-num "${AUTOSOME_NUM}" \
  --extract "${OUTDIR}/rog_snps.txt" \
  --make-grm \
  --thread-num "${SLURM_CPUS_PER_TASK}" \
  --out "${ROG_GRM_PREFIX}"

# ---------------------------------------------------------
# 3) RUN TWO-GRM GREML (exact complement)
# ---------------------------------------------------------

MGRM_FILE="${OUTDIR}/mgrm.txt"
cat > "${MGRM_FILE}" <<EOF
${LOCUS_GRM_PREFIX}
${ROG_GRM_PREFIX}
EOF

REML_OUT="${OUTDIR}/reml"

GCTA_ARGS=(
  --reml
  --autosome-num "${AUTOSOME_NUM}"
  --thread-num "${SLURM_CPUS_PER_TASK}"
  --mgrm "${MGRM_FILE}"
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

echo "Finished ${PHENO} / ${LEAD_SNP}"
echo "Variance explained by locus = V(G1)/Vp in:"
echo "  ${REML_OUT}.hsq"