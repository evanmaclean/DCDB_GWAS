library(tidyverse)
library(fs)
library(gt)

lambda_tbl <- read_csv("gcta_files/lambda_gc_summary.csv", col_types = cols(
  trait = col_character(), n_snps = col_double(), lambda_gc = col_double()
))

lambda_gt <- lambda_tbl |>
  gt() |>
  cols_label(
    trait = "Trait",
    n_snps = "N SNPs tested",
    lambda_gc = "\u03BBGC"
  ) |>
  fmt_integer(columns = n_snps) |>
  fmt_number(columns = lambda_gc, decimals = 3)

lambda_gt

out_dir <- path("tables")
dir_create(out_dir)
gtsave(lambda_gt, path(out_dir, "lambda_gc_summary_table.docx"))
