library(tidyverse)
library(fs)
library(gt)

sig_threshold <- 1.486325802616e-06
dict <- readxl::read_excel("gwas_dictionary.xlsx")
peak_info <- read_csv("misc/peak_manifest_annotate.csv", col_types = cols(
  peak_id = col_double(), chr = col_double(), lead_snp = col_character(),
  bp = col_double(), start = col_double(), end = col_double(),
  trait = col_character(), lead_p = col_double(),
  minus = col_double(), plus = col_double()
))
sig_peaks <- peak_info |>
  filter(lead_p <= sig_threshold) |>
  mutate(trait = tolower(trait))
sig_traits <- sig_peaks |> distinct(trait) |> pull(trait)
trait_codes <- dict$old[match(sig_traits, tolower(dict$new))]
names(trait_codes) <- sig_traits
read_mlma_col <- cols(
  Chr = col_double(), SNP = col_character(), bp = col_double(),
  A1 = col_character(), A2 = col_character(), Freq = col_double(),
  b = col_double(), se = col_double(), p = col_double()
)

comparison <- map_dfr(sig_traits, function(tr) {
  code <- trait_codes[[tr]]
  this_snps <- filter(sig_peaks, trait == tr) |> pull(lead_snp)
  baseline <- read_delim(
    path("gcta_files", "lead_snp_matched_baseline", paste0(code, ".mlma")),
    col_types = read_mlma_col
  ) |>
    filter(SNP %in% this_snps) |>
    select(SNP, beta_baseline = b, se_baseline = se, p_baseline = p)
  adj <- read_delim(
    path("gcta_files", "lead_snp_ancestry_check", paste0(code, ".mlma")),
    col_types = read_mlma_col
  ) |>
    filter(SNP %in% this_snps) |>
    select(SNP, beta_ancestry_adj = b, se_ancestry_adj = se, p_ancestry_adj = p)
  baseline |>
    inner_join(adj, by = "SNP") |>
    mutate(
      trait = tr,
      beta_pct_change_ancestry = 100 * (beta_ancestry_adj - beta_baseline) / beta_baseline,
      sig_baseline = p_baseline <= sig_threshold,
      sig_ancestry_adj = p_ancestry_adj <= sig_threshold
    ) |>
    relocate(trait, .before = SNP)
})
print(comparison)

out_dir <- path("tables")
dir_create(out_dir)
write_csv(comparison, path(out_dir, "na6_ancestry_covariate_comparison.csv"))

message(
  "Ancestry covariate alone: median |beta % change| = ",
  round(median(abs(comparison$beta_pct_change_ancestry)), 2), "%. ",
  sum(comparison$sig_ancestry_adj), " of ", nrow(comparison), " stay significant."
)
message(
  "Correlation of beta, baseline vs. ancestry-adjusted: r = ",
  round(cor(comparison$beta_baseline, comparison$beta_ancestry_adj), 4)
)

## Build a formatted comparison table for the supplement ----
table_data <- comparison |>
  select(
    trait, SNP,
    beta_baseline, se_baseline, p_baseline,
    beta_ancestry_adj, se_ancestry_adj, p_ancestry_adj,
    beta_pct_change_ancestry
  )

na6_gt <- table_data |>
  gt() |>
  tab_header(
    title = "SNP effect size comparison: baseline vs. ancestry-adjusted models"
  ) |>
  tab_spanner(
    label = "Baseline (no ancestry covariate)",
    columns = c(beta_baseline, se_baseline, p_baseline)
  ) |>
  tab_spanner(
    label = "Ancestry-adjusted (+ percent Labrador ancestry)",
    columns = c(beta_ancestry_adj, se_ancestry_adj, p_ancestry_adj)
  ) |>
  cols_label(
    trait = "Trait",
    SNP = "SNP",
    beta_baseline = "Beta", se_baseline = "SE", p_baseline = "P",
    beta_ancestry_adj = "Beta", se_ancestry_adj = "SE", p_ancestry_adj = "P",
    beta_pct_change_ancestry = "% change"
  ) |>
  fmt_number(
    columns = c(starts_with("beta"), starts_with("se")),
    decimals = 3
  ) |>
  fmt_scientific(
    columns = starts_with("p_"),
    decimals = 2
  ) |>
  tab_source_note(
    source_note = "Primary models used a leave-one-chromosome-out relatedness matrix (--mlma-loco). For computational efficiency, the sensitivity analysis reported in this table instead uses --mlma with a single genome-wide genetic relatedness matrix common to both models."
  )
na6_gt

gtsave(na6_gt, path(out_dir, "na6_ancestry_covariate_comparison_table.docx"))