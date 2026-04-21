library(tidyverse)
library(fs)
library(readxl)
library(bigsnpr)
library(viridis)
library(patchwork)


# set significant p threshold
sig = 1.486325802616e-06

dict <- read_excel('gwas_dictionary.xlsx')

res <- dir_ls('gcta_files/mlma_loco_results')
res <- res[str_ends(string = res, pattern = ".mlma")]
#load('nearest_gene/nearest_gene.RData')
#nearest.gene <- select(nearest.gene, SNP, distance, gene)

my_peaks <- map_dfr(res, function(i) {
 
  dat <- read_delim(i, col_types = cols(
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
  
  dat <- filter(dat, !is.na(p)) # for now drop any sites with NA
  my_var <- str_remove(i, "gcta_files/mlma_loco_results/") %>% str_remove(".loco.mlma")
  my_lab <- dict$new[dict$old == my_var]
  
  # find significant AND suggestive snps
  sig <- dat %>%
    filter(p < 2.97e-5) %>%
    arrange(Chr, bp)
  
  window <- 5000000  # 5 MB
  
  sig <- sig %>%
    group_by(Chr) %>%
    arrange(bp) %>%
    mutate(
      new_peak = c(TRUE, diff(bp) > window),
      peak_id = cumsum(new_peak)
    ) %>%
    ungroup()
  
  lead_snps <- sig %>%
    group_by(Chr, peak_id) %>%
    slice_min(p, with_ties = FALSE) %>%
    ungroup()
  
  peaks <- sig %>%
    group_by(Chr, peak_id) %>%
    summarise(
      peak_start = min(bp),
      peak_end   = max(bp),
      lead_pos   = bp[which.min(p)],
      lead_p     = min(p),
      .groups = "drop"
    )
  
  peaks$trait <- my_lab
  
  # get snp ID back 
  peaks <- left_join(peaks, select(sig, Chr, lead_pos = bp, SNP))
  
  peaks
  
})

# reformat file for input to next step

my_peaks <- select(my_peaks, 
                   peak_id,
                   chr = Chr,
                   lead_snp = SNP,
                   bp = lead_pos,
                   trait,
                   lead_p
                   )

# quick aside, save two columns of this for my locus GREML script
loc_greml <- select(my_peaks, SNP=lead_snp, trait = trait, p = lead_p)
loc_greml <- left_join(loc_greml, select(dict, trait = new, pheno = old))
loc_greml <- filter(loc_greml, p <= sig) |> select(SNP, pheno)
write_csv(loc_greml, file = "gcta_files/sig_snp_manifest.csv")
# end quick aside

my_peaks <- my_peaks |> mutate(peak_id = row_number(), start = bp-5000000, end = bp+5000000)
my_peaks <- mutate(my_peaks, start = if_else(start < 0, 0, false = start))

#put in order needed for next script that gets LD of nearby snps
my_peaks <- my_peaks |> select(peak_id, chr, lead_snp, bp, start, end, trait, lead_p)
write_tsv(my_peaks, 'locus_data/peak_manifest.tsv')
