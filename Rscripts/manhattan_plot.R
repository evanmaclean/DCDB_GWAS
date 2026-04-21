library(tidyverse)
library(fs)
library(readxl)
library(CMplot)

threshold <- 1.49e-6

dict <- read_excel('gwas_dictionary.xlsx')

res <- dir_ls('gcta_files/mlma_loco_results')
res <- res[str_ends(string = res, pattern = ".mlma")]

dat <- map_dfr(res, function(i) {
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
  ))

  dat <- filter(dat, !is.na(p)) # for now drop any sites with NA
  my_var <- str_remove(i, "gcta_files/mlma_loco_results/") %>% str_remove(".loco.mlma")
  my_lab <- dict$new[dict$old == my_var]
  dat <- mutate(dat, trait = my_lab) 
})

# Pivot to wide format - one p-value column per trait
wide <- dat |> 
  select(SNP, Chr, bp, p, trait) %>%
  pivot_wider(names_from = trait, values_from = p)

thinning_threshold <- 1e-2  # relaxed - preserves full peak shape
plot_threshold <- 1.49e-6   # actual significance line

sig <- wide |>
  filter(if_any(-c(SNP, Chr, bp), ~ . < thinning_threshold))

# move prop back up to 0.15 or higher when ready for final plot
nonsig <- wide |>
  filter(if_all(-c(SNP, Chr, bp), ~ . >= thinning_threshold | is.na(.))) |>
  slice_sample(prop = 0.15)

downsampled <- bind_rows(sig, nonsig)

# multi track plot ----
# trait_colors <- c(
#   "#390099", "#390099",  # positions 1-2
#   "#9e0059", "#9e0059",  # positions 3-4
#   "#ff0054",             # position 5
#   "#ff5400",             # position 6
#   "#ffbd00", "#ffbd00"   # positions 7-8
# )
# 
# # Build matrix: 8 traits x 1 color each (CMplot will recycle across chromosomes)
# col_matrix <- matrix(trait_colors, nrow = 8, ncol = 2, byrow = TRUE)
# 
# CMplot(downsampled,
#        plot.type = "m",
#        type = "p",
#        multracks = TRUE,
#        threshold = 1.49e-6,
#        threshold.lty = 1,
#        threshold.lwd = 1,
#        threshold.col = "red",
#        amplify = FALSE,
#        bin.size = 1e6,
#        chr.den.col = NULL,
#        col = col_matrix,
#        file = "png",
#        dpi = 300,
#        verbose = TRUE,
#        file.output = TRUE,
#        file.name = "multi_trait_manhattan",
#        width = 10, height = 2, 
#        cex=c(0.5,0.1,1),
#        axis.cex = 0.5,
#        main.cex = 0.1,
#        lab.cex = 0.5,
#        mar.between = 1,
#        ylim = c(2,7)
#        ) # second number is for manhattan point size

## single track plot ----

col_matrix <- matrix(
  c("#023e8a", "#0077b6",
    "#34a0a4", "#76c893",
    "#e9c46a",
    "#f4a261",
    "#ff6d00", "#e5383b"),
  nrow = 8, ncol = 2, byrow = TRUE
)

CMplot(downsampled,
       plot.type = "m",
       type = "p",
       multracks = FALSE,
       multraits = TRUE,
       threshold = 1.49e-6,
       amplify = TRUE,
       bin.size = 1e6,
       chr.den.col = NULL,
       col = col_matrix,
       points.alpha = 150L,
       file = "tiff",
       dpi = 300,
       verbose = TRUE,
       file.output = TRUE,
       file.name = "manhattan",
       width = 11, height = 4,
       cex = 0.5,
       main.cex = 0.1,
       signal.cex = 1,
       lab.cex = 2,
       axis.cex = 1.5,
       mar.between = 1,
       ylim = c(2, 7),
       legend.ncol = 4,
       legend.pos = "middle",
       legend.cex = 0.8,
       threshold.col = "black",
       threshold.lty = 3) # dotted

file_move(path = 'Multi-traits_Manhtn.manhattan.tiff',
          new_path = 'figures/manhattan.tiff')

