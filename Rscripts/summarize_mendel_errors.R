library(tidyverse)

mendel_dir <- "plink"
n_sites <- 6671897  # variants in mendel_check_out.lmendel; 
ref_panel_file <- "plink/ref_dogs.csv"  # columns: name, seq_group (Refpanel dogs only need to be listed)

# .mendel: one row per (child, SNP) Mendelian-incompatible call
mendel <- read_table(
  file.path(mendel_dir, "mendel_check_out.mendel"),
  col_types = cols(.default = "c")
)

# per-child error counts -- this is the correct, unambiguous quantity
child_errors <- mendel %>%
  count(KID, name = "n_errors") %>%
  mutate(
    n_sites = n_sites,
    error_rate = n_errors / n_sites
  ) %>%
  arrange(desc(n_errors))

stopifnot(child_errors$n_errors[child_errors$KID == "AGNES-II"] == 8818)
hollyn_kids_sum <- sum(child_errors$n_errors[child_errors$KID %in% c("NECTARINE", "REX-VII", "RITA-VII")])
stopifnot(hollyn_kids_sum == 29965)
message("Sanity checks passed: AGNES-II and HOLLYN-family sums match hand verification.")

# --- summary stats ---
n_offspring <- nrow(child_errors)
pooled_rate <- sum(child_errors$n_errors) / (n_offspring * n_sites)

summary_stats <- tibble(
  n_offspring = n_offspring,
  total_errors = sum(child_errors$n_errors),
  n_sites = n_sites,
  pooled_error_rate = pooled_rate,
  median_error_rate = median(child_errors$error_rate),
  mean_error_rate = mean(child_errors$error_rate),
  min_error_rate = min(child_errors$error_rate),
  max_error_rate = max(child_errors$error_rate)
)
print(summary_stats)

fam <- read_table(
  file.path(mendel_dir, "mendel_check.fam"),
  col_names = c("fid", "iid", "pat", "mat", "sex", "pheno"),
  col_types = "cccccc"
)

ref_panel <- read_csv(ref_panel_file, col_types = cols(.default = "c"))
ref_ids <- ref_panel$name

classified <- child_errors %>%
  left_join(fam, by = c("KID" = "iid")) %>%
  mutate(
    pat_known = pat != "0",
    mat_known = mat != "0",
    n_known_parents = as.integer(pat_known) + as.integer(mat_known),
    pat_refpanel = pat_known & pat %in% ref_ids,
    mat_refpanel = mat_known & mat %in% ref_ids,
    group = case_when(
      n_known_parents == 2 & pat_refpanel & mat_refpanel ~ "trio_both_refpanel",
      n_known_parents == 2 ~ "trio_one_or_more_lowcov",
      n_known_parents == 1 & ((pat_known & pat_refpanel) | (mat_known & mat_refpanel)) ~ "duo_refpanel_parent",
      n_known_parents == 1 ~ "duo_lowcov_parent",
      TRUE ~ "unclassified"
    )
  )

group_summary <- classified %>%
  group_by(group) %>%
  summarise(
    n = n(),
    median_error_rate = median(error_rate),
    mean_error_rate = mean(error_rate),
    min_error_rate = min(error_rate),
    max_error_rate = max(error_rate),
    .groups = "drop"
  ) %>%
  arrange(desc(median_error_rate))

print(group_summary)

outliers <- classified %>%
  group_by(group) %>%
  filter(error_rate > 2 * median(error_rate)) %>%
  ungroup() %>%
  select(KID, pat, mat, group, n_errors, error_rate) %>%
  arrange(desc(error_rate))

print(outliers)

# --- write outputs for documentation / re-use ---
write_csv(classified, file.path(mendel_dir, "mendel_check_per_dog_error_rates.csv"))
write_csv(summary_stats, file.path(mendel_dir, "mendel_check_summary_stats.csv"))
write_csv(group_summary, file.path(mendel_dir, "mendel_check_group_summary.csv"))
write_csv(outliers, file.path(mendel_dir, "mendel_check_outliers.csv"))