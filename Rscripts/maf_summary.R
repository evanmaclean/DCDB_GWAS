library(tidyverse)
library(fs)
library(scales)

maf <- read_table(
  "gcta_files/qc/maf_summary.frq",
  col_types = cols(
    CHR = col_double(),
    SNP = col_character(),
    A1 = col_character(),
    A2 = col_character(),
    MAF = col_double(),
    NCHROBS = col_double()
  )
)
# Standard population-genetics MAF bins (rare / low-frequency / common).
# MAF is already coded as the minor-allele frequency (0-0.5 by definition),
# not the reference-allele frequency, so these thresholds apply directly.

maf_summary_stats <- tibble(
  n_snps = nrow(maf),
  mean_maf = mean(maf$MAF, na.rm = TRUE),
  median_maf = median(maf$MAF, na.rm = TRUE),
  pct_rare_under_1pct = mean(maf$MAF < 0.01, na.rm = TRUE) * 100,
  pct_low_freq_1_5pct = mean(maf$MAF >= 0.01 & maf$MAF < 0.05, na.rm = TRUE) * 100,
  pct_common_5pct_plus = mean(maf$MAF >= 0.05, na.rm = TRUE) * 100
)

print(maf_summary_stats)

qc_dir <- path("gcta_files", "qc")
dir_create(qc_dir)
write_csv(maf_summary_stats, path(qc_dir, "maf_summary_stats.csv"))

# Histogram of the full genome-wide MAF distribution. 

supp_dir <- path("figures", "supplemental")
dir_create(supp_dir)

maf_plot <- ggplot(maf, aes(x = MAF)) +
  geom_histogram(binwidth = 0.01, boundary = 0, fill = "steelblue", color = "white") +
  labs(
    x = "Minor allele frequency",
    y = "Number of SNPs",
    title = "Genome-wide MAF distribution of SNPs tested in GWAS"
  ) +
  scale_y_continuous(labels = label_comma()) +
  theme_bw()

ggsave(path(supp_dir, "FigureS_MAF_distribution.png"), maf_plot,
       width = 6, height = 4, dpi = 300, bg = "white")
