library(tidyverse)

qcovar <- read_delim('gcta_files/qcovar.qcovar', col_names = c('fid', 'iid', 'test_age_years'))
ancestry <- read_csv('genomic_PCA/breed_composition.csv') |> select(iid = id, percent.lab)

qcovar_with_ancestry <- qcovar |>
  left_join(ancestry, by = "iid")

stopifnot(sum(is.na(qcovar_with_ancestry$percent.lab)) == 0)  # every dog must have ancestry data

write_delim(qcovar_with_ancestry, "gcta_files/qcovar_with_ancestry.qcovar", col_names = FALSE)

# mini pheno list -- just the two flagged traits
# UT_looking_time = "gaze (unsolvable task)", UT_manipulate_time = "manipulation time (unsolvable task)"
# per gwas_dictionary.xlsx / the naming in your data-prep script -- please confirm these two codes
# are right before running
writeLines(c("UT_looking_time", "UT_manipulate_time"), "gcta_files/pheno_list_ancestry_check")