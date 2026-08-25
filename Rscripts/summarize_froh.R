
library(tidyverse)

plink_dir <- "plink"

# --- compute autosomal span directly from the marker panel ---
hom_summary <- read_table(
  file.path(plink_dir, "froh.hom.summary"),
  col_types = cols(.default = "c")
) %>%
  mutate(BP = as.numeric(BP))

chrom_span <- hom_summary %>%
  group_by(CHR) %>%
  summarise(span_bp = max(BP) - min(BP), .groups = "drop")

stopifnot(nrow(chrom_span) == 38)  # dog autosome count -- sanity check
genome_kb <- sum(chrom_span$span_bp) / 1000

message(sprintf("Autosomal span from marker panel: %.0f kb (%.3f Gb) across %d chromosomes",
                genome_kb, genome_kb / 1e6, nrow(chrom_span)))

# --- per-dog FROH ---
froh <- read_table(
  file.path(plink_dir, "froh.hom.indiv"),
  col_types = cols(.default = "c")
) %>%
  mutate(
    NSEG = as.integer(NSEG),
    total_kb = as.numeric(KB),
    FROH = total_kb / genome_kb
  ) %>%
  select(FID, IID, NSEG, total_kb, FROH)

summary_stats <- froh %>%
  summarise(
    n_dogs = n(),
    genome_kb = genome_kb,
    median_froh = median(FROH),
    mean_froh = mean(FROH),
    sd_froh = sd(FROH),
    min_froh = min(FROH),
    max_froh = max(FROH),
    p95_froh = quantile(FROH, 0.95),
    p99_froh = quantile(FROH, 0.99)
  )
print(summary_stats)

above_first_cousin <- froh %>%
  filter(FROH >= 0.0625) %>%
  arrange(desc(FROH))
print(above_first_cousin)

write_csv(froh, file.path(plink_dir, "froh_per_dog.csv"))
write_csv(summary_stats, file.path(plink_dir, "froh_summary_stats.csv"))
