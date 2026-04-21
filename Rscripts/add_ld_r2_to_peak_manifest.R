library(tidyverse)
library(readxl)
library(fs)
library(glue)
library(IRanges)
library(viridis)
library(patchwork)
library(here)

# set a window size to filter to
window = 5000000

# -----------------------------
# 1. Read peak manifest
# -----------------------------
manifest <- read_tsv("locus_data/peak_manifest.tsv", col_types = "nncnnnc")

# -----------------------------
# 2. Read all LD files
# -----------------------------
# List all files
ld_files <- dir_ls(path = "locus_data/ld_results/")
ld_files <- ld_files[str_ends(ld_files, pattern = ".vcor")]

# Helper to read a single LD file
read_ld_file <- function(ld_file) {
  read_tsv(ld_file, col_types = "nncnncn") %>%
    transmute(
      lead_snp = ID_A,
      SNP      = ID_B,
      BP       = POS_B,
      r2       = UNPHASED_R2
    )
}

# Read all LD files and bind into one table
ld_all <- map_dfr(ld_files, read_ld_file)


ld_all <- mutate(ld_all,
                 lead_bp = str_extract(lead_snp, "(?<=:).*") |> as.numeric(),
                 distance = abs(lead_bp-BP))

ld_all <- filter(ld_all, distance < window)
        
# -----------------------------
# 3. Merge LD with manifest
# -----------------------------
ld_with_peak <- ld_all %>%
  left_join(manifest %>% select(peak_id, lead_snp, chr, trait, lead_p), by = "lead_snp")

# I need to add the lead snps themselves
add_leads <- select(manifest, lead_snp = lead_snp, SNP = lead_snp, BP = bp, lead_p, peak_id, chr, trait) |> mutate(r2 = 1)
ld_with_peak <- bind_rows(ld_with_peak, add_leads)

# Now we have: peak_id, lead_snp, chr, SNP, BP, r2

# Now we need to join with the GWAS results.  Let's do this iteratively and create one file per task
dict <- read_excel('gwas_dictionary.xlsx')
task_names <- unique(ld_with_peak$trait)
dict <- filter(dict, new %in% task_names)

# now the 'old' column has the names we need to loop over to bring in the gwas results

# i <- dict$old[1]

for (i in dict$old) {
 
  dat <- read_delim(glue("gcta_files/mlma_loco_results/{i}.loco.mlma"), col_types = cols(
    Chr = col_double(),
    SNP = col_character(),
    bp = col_double(),
    A1 = col_character(),
    A2 = col_character(),
    Freq = col_double(),
    b = col_double(),
    se = col_double(),
    p = col_double()
  )
  )

  my_var <- str_remove(i, "gcta_files/mlma_loco_results/") %>% str_remove(".loco.mlma")
  my_lab <- dict$new[dict$old == my_var]
  
ld_task <- filter(ld_with_peak, trait == my_lab)
ld_gwas <- left_join(ld_task, select(dat, SNP, P = p)) |> mutate(P = replace_na(P,1))
ld_gwas <- ld_gwas |> mutate(task = my_lab)

ld_gwas <- mutate(ld_gwas, chr = as.character(chr))

# Optional (recommended): locus boundaries if present
# If not present, set them from the min/max BP per peak.
if (!all(c("start_bp", "end_bp") %in% names(ld_gwas))) {
  ld_gwas <- ld_gwas %>%
    group_by(peak_id) %>%
    mutate(start_bp = min(BP, na.rm = TRUE),
           end_bp   = max(BP, na.rm = TRUE)) %>%
    ungroup()
}

# add flag for lead SNP
ld_gwas <- mutate(ld_gwas, is_lead = SNP == lead_snp)
###

# remove columns I'm not using
ld_gwas <- select(ld_gwas, -lead_bp, -distance)

# I'm going to try the locuszoomr r package but am doing this locally.  Write the files here. 

outdir <- "locus_data/data_for_locus_plots"
dir.create(outdir, showWarnings = FALSE)

write_rds(ld_gwas, file = glue("{outdir}/{i}.rds"))
}

