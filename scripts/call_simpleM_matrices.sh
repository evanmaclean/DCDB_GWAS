#!/bin/bash
#SBATCH --job-name=simpleM
#SBATCH --array=1-38
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=0-1
#SBATCH --output=log_files/simpleM_%A_%a.out
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files meff/simpleM_out

module load R/4.4

CHR="${SLURM_ARRAY_TASK_ID}"

Rscript Rscripts/simpleM_matrices.R "${CHR}" > "meff/simpleM_out/chr${CHR}.row.tsv"