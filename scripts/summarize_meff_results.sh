#!/bin/bash
#SBATCH --output=log_files/summarize_meff
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=0-1   
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

echo -e "CHR\tFILE\tNIND\tNSNP_total\tNSNP_used\tMeff_chr" > simpleM_out/Meff_summary.tsv
cat simpleM_out/chr*.row.tsv >> simpleM_out/Meff_summary.tsv

awk -F'\t' 'NR>1{sum+=$6} END{
  printf("Meff_genomewide\t%.12f\n", sum);
  printf("GW_sig_0.05_over_Meff\t%.12e\n", 0.05/sum);
  printf("GW_suggestive_1_over_Meff\t%.12e\n", 1/sum);
}' simpleM_out/Meff_summary.tsv > simpleM_out/Meff_genomewide.txt