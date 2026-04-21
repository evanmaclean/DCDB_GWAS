#!/bin/bash
#SBATCH --output=log_files/gene-set-analysis-%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS
#SBATCH --array=1-8

set -euo pipefail

# Ensure key output dirs exist
mkdir -p log_files magma_files/magma_gene_set_analysis

phen=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pheno_list)
[[ -n "${phen}" ]] || { echo "ERROR: phen is empty" >&2; exit 1; }

MAGMA="magma_files/magma"
GENE_RESULTS="magma_files/magma_gene_analysis/${phen}.genes.raw"
GENESETS_DIR="magma_files/gene_sets"
OUTBASE="magma_files/magma_gene_set_analysis"

[[ -f "${GENE_RESULTS}" ]] || { echo "ERROR: missing ${GENE_RESULTS}" >&2; exit 1; }

run_set () {
  local set_name="$1"
  local set_path="${GENESETS_DIR}/${set_name}"
  local out_dir="${OUTBASE}/${set_name}"

  mkdir -p "${out_dir}"
  "${MAGMA}" \
    --gene-results "${GENE_RESULTS}" \
    --set-annot "${set_path}" col=2,1 \
    --out "${out_dir}/${phen}"
}

run_set epic_tissue_gene_sets
run_set DoGA_CNS_gene_sets
run_set zwir_temperament_2018
run_set big5_schwaba