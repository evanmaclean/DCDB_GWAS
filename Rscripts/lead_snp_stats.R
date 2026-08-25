library(tidyverse)
library(fs)
library(readxl)

sig_threshold <- 1.486325802616e-06

dict <- read_excel("gwas_dictionary.xlsx")

peak_info <- read_csv(
  "misc/peak_manifest_annotate.csv",
  col_types = cols(
    peak_id = col_double(), chr = col_double(), lead_snp = col_character(),
    bp = col_double(), start = col_double(), end = col_double(),
    trait = col_character(), lead_p = col_double(),
    minus = col_double(), plus = col_double()
  )
)

sig_peaks <- peak_info |>
  filter(lead_p <= sig_threshold) |>
  mutate(trait = tolower(trait))

stopifnot(nrow(sig_peaks) == 9) # sanity check against the "9 genomic regions" already reported in Results

writeLines(unique(sig_peaks$lead_snp), "gcta_files/lead_snps_extract.txt")

# Map each significant trait (dict$new, e.g. "reversal learning (cyl)") back
# to its .mlma file code (dict$old, e.g. "CYL.DT_C1"), same lookup direction
# manhattan_plot.R uses in reverse.
sig_traits <- sig_peaks |> distinct(trait) |> pull(trait)
trait_codes <- dict$old[match(sig_traits, tolower(dict$new))]
names(trait_codes) <- sig_traits

# Pheno-code list for the 5 significant-hit traits only, in the same
# line-per-phenotype format pheno_list uses, for
# sbatch_scripts/lead_snp_ancestry_check.sbatch's --array=1-N to index
# into (mirrors loco.sbatch's own pattern, just restricted to the traits
# that actually have a lead SNP to check).
writeLines(trait_codes, "gcta_files/pheno_list_lead_snp_check")

lead_snp_stats <- map_dfr(sig_traits, function(tr) {
  code <- trait_codes[[tr]]
  mlma_path <- path("gcta_files", "mlma_loco_results", paste0(code, ".loco.mlma"))

  mlma <- read_delim(mlma_path, col_types = cols(
    Chr = col_double(), SNP = col_character(), bp = col_double(),
    A1 = col_character(), A2 = col_character(), Freq = col_double(),
    b = col_double(), se = col_double(), p = col_double()
  ))

  this_trait_peaks <- filter(sig_peaks, trait == tr)

  # Join by SNP id (both peak_manifest's lead_snp and .mlma's SNP column use
  # "chr:pos" naming, confirmed against the same convention seen in
  # maf_summary.frq) - not by chr+bp separately, to avoid any off-by-one or
  # type-mismatch risk in a manual coordinate join.
  inner_join(this_trait_peaks, mlma, by = c("lead_snp" = "SNP")) |>
    transmute(
      trait = tr,
      chr = Chr,
      bp = bp.y,
      snp = lead_snp,
      a1 = A1,
      a2 = A2,
      freq_a1 = Freq,
      beta = b,
      se = se,
      p = p
    )
})

print(lead_snp_stats)

out_dir <- path("tables")
dir_create(out_dir)
write_csv(lead_snp_stats, path(out_dir, "lead_snp_summary_stats.csv"))

# Sanity check: the p-values here should match lead_p in peak_manifest_annotate.csv
# (same underlying .mlma row) - if these diverge meaningfully, something is
# off in the join 

check <- lead_snp_stats |>
  left_join(sig_peaks, by = c("trait", "snp" = "lead_snp")) |>
  mutate(p_diff = abs(p - lead_p))
print(select(check, trait, snp, p, lead_p, p_diff))
