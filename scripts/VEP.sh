#!/bin/bash
#SBATCH --output=log_files/VEP
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:10:00
#SBATCH --partition=standard
#SBATCH --account=evanmaclean
#SBATCH --chdir=/groups/evanmaclean/GWAS

module load plink/2.0

plink2 --bfile plink/cleaned \
  --dog \
  --extract VEP/vep_prep/all_sig_snps_uniq.txt \
  --export vcf bgz \
  --out VEP/vep_prep/sig_snps
  
  
  # run VEP
  vep \
  --cache --offline \
  --species canis_lupus_familiarisgsd \
  --assembly UU_Cfam_GSD_1.0 \
  --dir_cache /groups/evanmaclean/GWAS/VEP \
  -i VEP/vep_prep/sig_snps.vcf.gz \
  -o VEP/vep_prep/sig_snps_vep.tsv \
  --tab \
  --symbol --canonical --nearest symbol \
  --pick \
  --force_overwrite
  
  # make a summary table
  grep -m1 "^#Uploaded_variation" VEP/vep_prep/sig_snps_vep.tsv
  cut -f1-7 VEP/vep_prep/sig_snps_vep.tsv | head -n 30
  
  
 # make table with results
  awk 'BEGIN{OFS="\t"}
/^#/ {next}
{
  snp=$1
  split($2,a,":")
  chr=a[1]
  pos=a[2]

  symbol=$18
  nearest=$22

  gene = (symbol != "-" ? symbol : nearest)

  consequence=$7
  impact=$14

  print snp, chr, pos, gene, consequence, impact
}' VEP/vep_prep/sig_snps_vep.tsv \
> VEP/vep_prep/sig_snps_summary.tsv

# add header row
(echo -e "SNP\tCHR\tPOS\tGENE\tCONSEQUENCE\tIMPACT" && \
 cat VEP/vep_prep/sig_snps_summary.tsv) \
> VEP/vep_prep/sig_snps_summary_with_header.tsv