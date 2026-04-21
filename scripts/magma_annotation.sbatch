#!/bin/bash
#SBATCH --output=log_files/magma-annotate-%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files gene_info

MAGMA="magma_files/magma"

"${MAGMA}" \
  --annotate window=5,1.5 nonhuman \
  --snp-loc "plink/cleaned.bim" \
  --gene-loc "gene_info/magma_genes" \
  --out "gene_info/anno"