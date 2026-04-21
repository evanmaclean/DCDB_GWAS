#!/bin/bash
#SBATCH --output=log_files/pvals_for_magma.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6gb
#SBATCH --time=0:20:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS
# call using: sbatch scripts/call_pvals_for_magma.sh

mkdir -p log_files

module load R/4.4
module load gnu8

Rscript Rscripts/make_magma_pval_files.R