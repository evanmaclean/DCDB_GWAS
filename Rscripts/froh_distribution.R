library(tidyverse)
library(fs)

autosome_length_kb <- 2227792

froh <- read_table("plink/froh.hom.indiv", col_types = cols(
  FID = col_character(), IID = col_character(), PHE = col_double(),
  NSEG = col_double(), KB = col_double(), KBAVG = col_double()
)) |>
  mutate(
    froh = KB / autosome_length_kb,
    froh_pct = froh * 100
  )

message(
  "N = ", nrow(froh),
  "; mean FROH = ", round(mean(froh$froh_pct), 2), "%",
  "; median FROH = ", round(median(froh$froh_pct), 2), "%",
  "; range = ", round(min(froh$froh_pct), 2), "-", round(max(froh$froh_pct), 2), "%"
)

p <- ggplot(froh, aes(x = froh_pct)) +
  geom_histogram(binwidth = 0.25, fill = "#2F5496", color = "white", boundary = 0) +
  geom_vline(xintercept = median(froh$froh_pct), linetype = "dashed", color = "grey30") +
  labs(
    x = expression(paste(F[ROH], " (% of autosomal genome in runs of homozygosity)")),
    y = "Number of dogs"
  ) +
  theme_minimal(base_size = 13)

print(p)

out_dir <- path("figures")
dir_create(out_dir)
ggsave(path(out_dir, "FigureS_FROH_distribution.png"), p, width = 7, height = 5, dpi = 300)
