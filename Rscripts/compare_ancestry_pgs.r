library(tidyverse)

ancestry <- read_csv('genomic_PCA/breed_composition.csv') |> select(iid = id, percent.lab)

compare_pheno <- function(phen) {
  orig <- read_table(glue::glue('gcta_files/cv_blup/{phen}.profile')) |> select(iid = IID, score_orig = SCORE)
  adj  <- read_table(glue::glue('gcta_files/cv_blup_ancestry_adjusted/{phen}.profile')) |> select(iid = IID, score_adj = SCORE)
  dat <- reduce(list(ancestry, orig, adj), full_join, by = "iid")
  
  tibble(
    phen = phen,
    r_ancestry_original = cor(dat$percent.lab, dat$score_orig, use = "complete.obs"),
    r_ancestry_adjusted  = cor(dat$percent.lab, dat$score_adj, use = "complete.obs")
  )
}

map_dfr(c("UT_looking_time", "UT_manipulate_time"), compare_pheno)