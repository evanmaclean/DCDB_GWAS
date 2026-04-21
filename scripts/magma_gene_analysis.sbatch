#!/bin/bash
#SBATCH --output=log_files/magma-gene-analysis-%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=16:00:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS
#SBATCH --array=1-8

set -euo pipefail

mkdir -p log_files magma_files/magma_gene_analysis

phen=$(sed -n "${SLURM_ARRAY_TASK_ID}p" pheno_list)
[[ -n "${phen}" ]] || { echo "ERROR: phen is empty" >&2; exit 1; }

magma_files/magma \
  --bfile plink/cleaned \
  --gene-annot gene_info/anno.genes.annot \
  --pval "magma_files/pvals_for_magma/${phen}" ncol=3 \
  --gene-model snp-wise=mean \
  --out "magma_files/magma_gene_analysis/${phen}"