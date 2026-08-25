library(tidyverse)
library(fs)
library(gt)
library(readxl)

# -----------------------------
# Locus GREML results
# -----------------------------
hsq_files <- dir_ls(
  "gcta_files/greml_estimates/locus_greml",
  recurse = TRUE,
  glob = "*.hsq"
)

parse_hsq_g1_g2 <- function(file) {
  lines <- readr::read_lines(file) |>
    stringr::str_squish()
  
  keep <- lines[stringr::str_detect(lines, "^(V\\(G1\\)/Vp|V\\(G2\\)/Vp)\\b")]
  
  tibble(line = keep) |>
    mutate(
      term = stringr::str_extract(line, "^(V\\(G1\\)/Vp|V\\(G2\\)/Vp)"),
      estimate = as.numeric(stringr::word(line, 2)),
      se = as.numeric(stringr::word(line, 3)),
      pheno = path_file(path_dir(path_dir(file))),
      snp = path_file(path_dir(file)),
      file = file,
      locus_snps_file = path(path_dir(file), "locus_snps.txt"),
      n_locus_snps = if_else(
        file_exists(locus_snps_file),
        length(readr::read_lines(locus_snps_file)),
        NA_integer_
      )
    ) |>
    select(pheno, snp, term, estimate, se, n_locus_snps, file)
}

results_long <- map_dfr(hsq_files, parse_hsq_g1_g2)

results <- results_long |>
  select(pheno, snp, n_locus_snps, term, estimate, se) |>
  pivot_wider(
    names_from = term,
    values_from = c(estimate, se)
  ) |>
  rename(
    locus_var = `estimate_V(G1)/Vp`,
    locus_se  = `se_V(G1)/Vp`,
    rog_var   = `estimate_V(G2)/Vp`,
    rog_se    = `se_V(G2)/Vp`
  ) |>
  relocate(pheno, snp, n_locus_snps)

# -----------------------------
# Overall SNP heritability
# -----------------------------
overall_h2_files <- dir_ls(
  "gcta_files/greml_estimates/overall_h2",
  recurse = TRUE,
  glob = "*.hsq"
)

parse_overall_h2 <- function(file) {
  lines <- readr::read_lines(file) |>
    stringr::str_squish()
  
  vg_line <- lines[stringr::str_detect(lines, "^V\\(G\\)/Vp\\b")]
  
  tibble(
    pheno = path_file(path_dir(file)),
    total_h2 = as.numeric(stringr::word(vg_line, 2)),
    total_h2_se = as.numeric(stringr::word(vg_line, 3))
  )
}

overall_h2 <- map_dfr(overall_h2_files, parse_overall_h2)

# -----------------------------
# Join overall h2 to locus table
# Keep traits with no significant loci
# -----------------------------
results_full <- overall_h2 |>
  left_join(results, by = "pheno") |>
  arrange(pheno, desc(locus_var)) |>
  mutate(
    locus_lower = locus_var - 1.96 * locus_se,
    locus_upper = locus_var + 1.96 * locus_se,
    rog_lower   = rog_var - 1.96 * rog_se,
    rog_upper   = rog_var + 1.96 * rog_se
  )

# -----------------------------
# Figures
# -----------------------------
plot_dat <- results |>
  mutate(
    locus_lower = pmax(0, locus_var - 1.96 * locus_se),
    locus_upper = locus_var + 1.96 * locus_se
  )

p_locus <- ggplot(
  plot_dat,
  aes(x = reorder(snp, locus_var), y = locus_var)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = locus_lower, ymax = locus_upper),
    width = .15
  ) +
  coord_flip() +
  labs(
    x = "Associated locus",
    y = "Variance explained (V(G1)/Vp)",
    title = "Variance explained by GWAS loci"
  ) +
  theme_bw()

p_compare <- plot_dat |>
  pivot_longer(
    cols = c(locus_var, rog_var),
    names_to = "component",
    values_to = "variance"
  ) |>
  mutate(
    component = recode(
      component,
      locus_var = "GWAS locus",
      rog_var = "Rest of genome"
    )
  ) |>
  ggplot(aes(component, variance)) +
  geom_boxplot(width = .5) +
  geom_jitter(width = .1, alpha = .5) +
  labs(
    x = "",
    y = "Variance explained",
    title = "Regional vs genome-wide SNP heritability"
  ) +
  theme_bw()

p_locus
p_compare

# -----------------------------
# Table
# -----------------------------
table_dat <- results_full |>
  mutate(
    total_h2_fmt = sprintf("%.3f (%.3f)", total_h2, total_h2_se),
    locus = if_else(
      is.na(locus_var),
      "",
      sprintf("%.3f (%.3f)", locus_var, locus_se)
    ),
    rog = if_else(
      is.na(rog_var),
      "",
      sprintf("%.3f (%.3f)", rog_var, rog_se)
    ),
    Locus = if_else(is.na(snp), "", snp),
    `SNPs in locus` = if_else(
      is.na(n_locus_snps),
      "",
      as.character(n_locus_snps)
    )
  ) |>
  group_by(pheno) |>
  mutate(
    `Total SNP h²` = if_else(row_number() == 1, total_h2_fmt, "")
  ) |>
  ungroup() |>
  select(
    pheno,
    Locus,
    `Total SNP h²`,
    `SNPs in locus`,
    `locus SNP h²` = locus,
    `rest of genome SNP h²` = rog
  )

dict <- read_excel("gwas_dictionary.xlsx")

table_dat <- left_join(
  table_dat,
  select(dict, Trait = new, pheno = old),
  by = "pheno"
) |>
  mutate(
    Trait = coalesce(Trait, pheno)
  ) |>
  select(-pheno) |>
  separate(
    Locus,
    into = c("Chromosome", "Position"),
    sep = ":",
    fill = "right",
    remove = TRUE
  ) |>
  mutate(
    Chromosome = replace_na(Chromosome, ""),
    Position = replace_na(Position, "")
  ) |> arrange(Trait)

gt_tbl <- table_dat |>
  gt(
    groupname_col = "Trait",
    row_group_as_column = TRUE
  ) |>
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_stub()
  ) |>
  opt_css(
    css = "
      .gt_stub_row_group {
        text-align: center !important;
        vertical-align: middle !important;
      }
    "
  ) |>
  tab_spanner(
    label = "Lead SNP",
    columns = c(Chromosome, Position)
  ) |>
  cols_align(
    align = "center",
    columns = c(
      Chromosome,
      Position,
      `Total SNP h²`,
      `SNPs in locus`,
      `locus SNP h²`,
      `rest of genome SNP h²`
    )
  ) |>
  cols_width(
    Chromosome ~ px(110),
    Position ~ px(110),
    `Total SNP h²` ~ px(140),
    `SNPs in locus` ~ px(110),
    `locus SNP h²` ~ px(125),
    `rest of genome SNP h²` ~ px(190)
  )

gt_tbl

gtsave(gt_tbl, "tables/locus_variance_explained_table.docx")

