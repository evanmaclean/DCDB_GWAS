library(tidyverse)
library(glue)
library(readxl)
library(broom)

options(scipen = 999)

phens <- read_csv('pheno_list', col_names = "pheno") |> pull(pheno)
dict <- read_excel('gwas_dictionary.xlsx')

# --- shared covariate/ancestry/FROH tables (same across all traits) ---
covar <- read_delim('gcta_files/covar.covar',
                    col_names = c('iid', 'fid', 'sex', 'experimenter', 'region', 'battery_type', 'batch'))
qcovar <- read_delim('gcta_files/qcovar.qcovar',
                     col_names = c('iid', 'fid', 'test_age_years'))

ancestry <- read_csv('genomic_PCA/breed_composition.csv') |>
  select(iid = id, percent.lab) |>
  mutate(percent.lab_z = as.vector(scale(percent.lab)))

froh <- read_csv('plink/froh_per_dog.csv') |>
  select(iid = IID, FROH) |>
  mutate(FROH_z = as.vector(scale(FROH)))

stopifnot(nrow(ancestry) > 0, nrow(froh) > 0)
message(sprintf("ancestry table: %d dogs | FROH table: %d dogs", nrow(ancestry), nrow(froh)))

# --- per-trait correlations ---
results <- map_dfr(phens, function(pheno) {
  
  my_lab <- dict$new[match(pheno, dict$old)]
  
  phenos <- read_delim(glue('gcta_files/phenos/{pheno}.phen'), col_names = c("iid", "fid", "pheno"))
  
  score <- read_table(glue('gcta_files/cv_blup/{pheno}.profile')) |>
    select(iid = IID, score = SCORE) |>
    mutate(score_z = as.vector(scale(score)))
  
  mod_dat <- purrr::reduce(
    list(covar, qcovar, phenos, score, ancestry, froh),
    full_join, by = "iid"
  )
  
  base_mod <- lm(pheno ~ test_age_years + region + experimenter + battery_type + batch, data = mod_dat)
  mod_dat$resid <- augment(base_mod, newdata = mod_dat)$.resid
  
  test_pair <- function(x, y) {
    pear <- suppressWarnings(cor.test(x, y))
    spear <- suppressWarnings(cor.test(x, y, method = "spearman"))
    tibble(
      r = unname(pear$estimate), p = pear$p.value,
      rho = unname(spear$estimate), p_spearman = spear$p.value,
      n = sum(complete.cases(x, y))
    )
  }
  
  bind_rows(
    test_pair(mod_dat$percent.lab_z, mod_dat$score_z) |> mutate(predictor = "ancestry", outcome = "prs"),
    test_pair(mod_dat$percent.lab_z, mod_dat$pheno)   |> mutate(predictor = "ancestry", outcome = "raw_phenotype"),
    test_pair(mod_dat$percent.lab_z, mod_dat$resid)   |> mutate(predictor = "ancestry", outcome = "adjusted_phenotype"),
    test_pair(mod_dat$FROH_z, mod_dat$score_z)        |> mutate(predictor = "FROH", outcome = "prs"),
    test_pair(mod_dat$FROH_z, mod_dat$pheno)          |> mutate(predictor = "FROH", outcome = "raw_phenotype"),
    test_pair(mod_dat$FROH_z, mod_dat$resid)          |> mutate(predictor = "FROH", outcome = "adjusted_phenotype")
  ) |> mutate(task = my_lab, .before = 1)
})

print(results, n = Inf)

# --- joint models (standardized predictors -> directly comparable coefficients) ---
joint_results <- map_dfr(phens, function(pheno) {
  my_lab <- dict$new[match(pheno, dict$old)]
  phenos <- read_delim(glue('gcta_files/phenos/{pheno}.phen'), col_names = c("iid", "fid", "pheno"))
  score <- read_table(glue('gcta_files/cv_blup/{pheno}.profile')) |>
    select(iid = IID, score = SCORE) |> mutate(score_z = as.vector(scale(score)))
  mod_dat <- purrr::reduce(list(covar, qcovar, phenos, score, ancestry, froh), full_join, by = "iid")
  base_mod <- lm(pheno ~ test_age_years + region + experimenter + battery_type + batch, data = mod_dat)
  mod_dat$resid <- augment(base_mod, newdata = mod_dat)$.resid
  
  joint_prs <- lm(score_z ~ percent.lab_z + FROH_z, data = mod_dat) |> tidy() |> mutate(outcome = "prs")
  joint_raw <- lm(pheno ~ percent.lab_z + FROH_z, data = mod_dat) |> tidy() |> mutate(outcome = "raw_phenotype")
  joint_adj <- lm(resid ~ percent.lab_z + FROH_z, data = mod_dat) |> tidy() |> mutate(outcome = "adjusted_phenotype")
  bind_rows(joint_prs, joint_raw, joint_adj) |> mutate(task = my_lab, .before = 1)
})

print(joint_results, n = Inf)

write_csv(results, "gcta_files/froh_ancestry_prs_correlations.csv")
write_csv(joint_results, "gcta_files/froh_ancestry_prs_joint_models.csv")

# --- summary figures ---
plot_dat <- map_dfr(phens, function(pheno) {
  my_lab <- dict$new[match(pheno, dict$old)]
  phenos <- read_delim(glue('gcta_files/phenos/{pheno}.phen'), col_names = c("iid", "fid", "pheno"))
  score <- read_table(glue('gcta_files/cv_blup/{pheno}.profile')) |>
    select(iid = IID, score = SCORE) |> mutate(score_z = as.vector(scale(score)))
  purrr::reduce(list(phenos, score, ancestry, froh), full_join, by = "iid") |> mutate(task = my_lab)
})

dir_path <- "figures/froh_ancestry_prs"
if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)

ggplot(plot_dat, aes(x = percent.lab, y = score_z)) +
  geom_point(alpha = 1/8) + geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ task, scales = "free", nrow = 2) + theme_classic(8) +
  labs(x = "Percent Labrador ancestry", y = "Genomic prediction (z-scored)")
ggsave(file.path(dir_path, "ancestry_vs_prs.tiff"), width = 8, height = 4)

ggplot(plot_dat, aes(x = percent.lab, y = pheno)) +
  geom_point(alpha = 1/8) + geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ task, scales = "free", nrow = 2) + theme_classic(8) +
  labs(x = "Percent Labrador ancestry", y = "Raw (normalized, unadjusted) phenotype")
ggsave(file.path(dir_path, "ancestry_vs_raw_phenotype.tiff"), width = 8, height = 4)

ggplot(plot_dat, aes(x = FROH, y = score_z)) +
  geom_point(alpha = 1/8) + geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ task, scales = "free", nrow = 2) + theme_classic(8) +
  labs(x = "FROH", y = "Genomic prediction (z-scored)")
ggsave(file.path(dir_path, "froh_vs_prs.tiff"), width = 8, height = 4)

ggplot(plot_dat, aes(x = FROH, y = pheno)) +
  geom_point(alpha = 1/8) + geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ task, scales = "free", nrow = 2) + theme_classic(8) +
  labs(x = "FROH", y = "Raw (normalized, unadjusted) phenotype")
ggsave(file.path(dir_path, "froh_vs_raw_phenotype.tiff"), width = 8, height = 4)