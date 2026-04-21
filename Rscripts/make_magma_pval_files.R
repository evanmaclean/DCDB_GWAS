#!/usr/bin/env Rscript
# run on cluster
library(tidyverse)
library(fs)
library(glue)

# for setting the working directory when not run interactively using an r project 
library(rprojroot)
root <- rprojroot::find_root(rprojroot::is_rstudio_project)
setwd(root)

res <- dir_ls('gcta_files/mlma_loco_results')
res <- res[str_ends(string = res, pattern = ".mlma")]

for (i in res) {
  dat <- read_delim(i, col_types = cols(
    Chr = col_double(),
    SNP = col_character(),
    bp = col_double(),
    A1 = col_character(),
    A2 = col_character(),
    Freq = col_double(),
    b = col_double(),
    se = col_double(),
    p = col_double()
  ))

my_var <- str_remove(i, "gcta_files/mlma_loco_results/") %>% str_remove(".loco.mlma")
dat <- select(dat, SNP, P = p)
dat <- filter(dat, !is.na(P))

# get the sample size per snp and add to the file
snps <- read_table(glue('magma_files/snp_sample_sizes/{my_var}.lmiss'), col_types = "ncnnn")
snps <- mutate(snps, N=N_GENO-N_MISS)
snps <- select(snps, SNP, N)

# join data and write
dat <- left_join(dat, snps)
write_delim(dat, file = glue('magma_files/pvals_for_magma/{my_var}'))
}
