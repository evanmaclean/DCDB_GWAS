library(tidyverse)
library(fs)

# ---- Settings ----
mlma_dir  <- "gcta_files/mlma_loco_results"
suffix    <- "loco.mlma"

pthresh   <- 1.486325802616e-06 # from meff
out_dir   <- "VEP/vep_prep"

dir_create(out_dir)

# ---- Read phenotype list ----
phenos <- read_lines("pheno_list") %>%
  discard(~ .x == "")

# ---- Collect significant SNPs ----
all_snps <- phenos %>%
  set_names() %>%
  map_dfr(function(pheno) {
    
    file <- path(mlma_dir, paste0(pheno, ".", suffix))
    
    if (!file_exists(file)) {
      warning("Missing file: ", file)
      return(tibble())
    }
    
    read_tsv(file, col_types = cols(
               Chr = col_double(),
             SNP = col_character(),
             bp = col_double(),
             A1 = col_character(),
             A2 = col_character(),
             Freq = col_double(),
             b = col_double(),
             se = col_double(),
             p = col_double()
    )) %>%
      filter(!is.na(p), p <= pthresh) %>%
      transmute(SNP, pheno)
    
  })

# ---- Unique SNP IDs ----
unique_snps <- all_snps %>%
  distinct(SNP)

# ---- Write output for PLINK --extract ----
write_lines(unique_snps$SNP,
            path(out_dir, "all_sig_snps_uniq.txt"))

cat("Total unique significant SNPs:",
    nrow(unique_snps), "\n")
