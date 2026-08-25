#!/usr/bin/env bash
#SBATCH --job-name=doga_intersect_variants
#SBATCH --output=log_files/doga_intersect_variants_%j.out
#SBATCH --error=log_files/doga_intersect_variants_%j.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

set -euo pipefail

mkdir -p log_files

module load R
module load bedtools

Rscript Rscripts/doga_intersect.R \
  --project_root shadow_strict_r2 \
  --output_dir doga_intersections_strict_r2

Rscript Rscripts/doga_intersect.R \
  --project_root shadow_leadonly \
  --output_dir doga_intersections_leadonly

echo ""
echo "Done. Compare against the original doga_intersections/ output:"
echo "  doga_intersections/all_traits_combined_summary.tsv          (r2 >= 0.1, original)"
echo "  doga_intersections_strict_r2/all_traits_combined_summary.tsv (r2 >= 0.9)"
echo "  doga_intersections_leadonly/all_traits_combined_summary.tsv  (lead SNP only)"
