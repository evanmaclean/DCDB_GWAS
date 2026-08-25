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
# NOTE: this is the combined, color-overlaid figure that is currently Figure 2.
# Reviewer 3 (minor point 2) found this hard to read because all 8 traits'
# points share the same track/coordinate space and visually overlap ("stacked
# on top of each other"). The commented-out multi-track version above avoids
# overlap but stacks 8 short tracks into one cramped image, which is the
# version you said you tried and didn't like. The new loop below (NA17
# section, further down) gives each trait a fully separate, standalone plot
# instead, which is what the reviewer explicitly asked for. Left this section
# as-is in case you still want the combined overlay figure for another
# purpose (e.g. a compact overview panel).

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
       lab.cex = 1,
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

supp_dir <- path("figures", "supplemental")
dir_create(supp_dir)

# Toggle for the two per-trait loops (Manhattan below, Q-Q further down).

run_individual_plots <- FALSE

traits <- setdiff(names(wide), c("SNP", "Chr", "bp"))

trait_colors <- c(
  "#023e8a", "#0077b6", "#34a0a4", "#76c893",
  "#e9c46a", "#f4a261", "#ff6d00", "#e5383b"
)

trait_colors <- rep_len(trait_colors, length(traits))
names(trait_colors) <- traits

if (run_individual_plots) {
  for (i in seq_along(traits)) {
    
    tr <- traits[i]
    tr_color <- trait_colors[[tr]]
    tr_slug <- make.names(tr) # safe-for-filenames version of the trait label
    
    trait_dat <- wide |>
      select(SNP, Chr, bp, all_of(tr)) |>
      filter(!is.na(.data[[tr]]))
    
    sig_tr <- filter(trait_dat, .data[[tr]] < thinning_threshold)
    nonsig_tr <- trait_dat |>
      filter(.data[[tr]] >= thinning_threshold) |>
      slice_sample(prop = 0.15)
    trait_downsampled <- bind_rows(sig_tr, nonsig_tr)
    
    CMplot(trait_downsampled,
           plot.type = "m",
           type = "p",
           multracks = FALSE,
           threshold = plot_threshold,
           threshold.col = "black",
           threshold.lty = 3,
           amplify = TRUE,
           bin.size = 1e6,
           chr.den.col = NULL,
           col = c(tr_color, "grey60"),
           points.alpha = 150L,
           file = "png", # was tiff — the full-page version came back at 200MB
           # as tiff, so switching all of these to png for the
           # same reason (lossless but properly compressed,
           # much better suited to scatter plots this size)
           dpi = 300,
           verbose = TRUE,
           file.output = TRUE,
           file.name = paste0("manhattan_", tr_slug),
           width = 10, height = 3,
           cex = 0.5,
           main.cex = 0.1,
           signal.cex = 1,
           lab.cex = 1.5,
           axis.cex = 1.2,
           mar.between = 1)
    
    produced_name <- paste0("Rect_Manhtn.manhattan_", tr_slug, ".png")
    if (file_exists(produced_name)) {
      file_move(produced_name, path(supp_dir, paste0("FigureS_Manhattan_", tr_slug, ".png")))
    } else {
      message("Could not find expected CMplot output for '", tr,
              "' — check the working directory for the actual filename and adjust produced_name.")
    }
  }
} # end run_individual_plots (Manhattan)

col_matrix_full <- matrix(
  as.vector(rbind(trait_colors, "grey60")),
  nrow = length(traits), ncol = 2, byrow = TRUE
)

page_height <- 18
per_track_height <- page_height / length(traits)

CMplot(downsampled,
       plot.type = "m",
       type = "p",
       multracks = TRUE,
       threshold = plot_threshold,
       threshold.lty = 3,
       threshold.lwd = 1,
       threshold.col = "black",
       amplify = FALSE,
       bin.size = 1e6,
       chr.den.col = NULL,
       col = col_matrix_full,
       file = "png",
       dpi = 300,
       verbose = TRUE,
       file.output = TRUE,
       file.name = "manhattan_full_page",
       width = 13, height = per_track_height,
       cex = 0.3,
       signal.cex = 0.5,
       axis.cex = 0.6,
       main.cex = 0.1,
       lab.cex = 0.38,
       mar.between = 0.3,
       ylim = c(2, 7))
)
produced_name_full <- "Multi-tracks_Manhtn.manhattan_full_page.png"
if (file_exists(produced_name_full)) {
  file_move(produced_name_full, path(supp_dir, "FigureS_Manhattan_all_traits.png"))
} else {
  message("Could not find expected CMplot output for the full-page multi-panel plot — ",
          "check the working directory for the actual filename and adjust produced_name_full.")
}

compute_lambda <- function(p) {
  p <- p[!is.na(p)]
  chisq_obs <- qchisq(1 - p, df = 1)
  median(chisq_obs, na.rm = TRUE) / qchisq(0.5, df = 1)
}

lambda_tbl <- map_dfr(traits, function(tr) {
  p <- wide[[tr]]
  p <- p[!is.na(p)]
  tibble(
    trait = tr,
    n_snps = length(p),
    lambda_gc = compute_lambda(p)
  )
})

print(lambda_tbl)
write_csv(lambda_tbl, "gcta_files/lambda_gc_summary.csv")

library(patchwork)
library(ggfastman)

if (run_individual_plots) {
  for (i in seq_along(traits)) {
    
    tr <- traits[i]
    tr_slug <- make.names(tr)
    
    p <- wide[[tr]]
    p <- p[!is.na(p)]
    
    qq_panel <- fast_qq(pvalue = p, speed = "fast", pointsize = 2) + ggtitle(tr)
    ggsave(path(supp_dir, paste0("FigureS_QQ_", tr_slug, ".png")), qq_panel,
           width = 6, height = 6, dpi = 300, bg = "white")
  }
}

qq_panels <- map(traits, function(tr) {
  p <- wide[[tr]]
  p <- p[!is.na(p)]
  # pointsize bumped up from the default 1.2, same as the individual-plot
  # loop above.
  fast_qq(pvalue = p, speed = "fast", pointsize = 4) + ggtitle(tr)
})

qq_grid <- wrap_plots(qq_panels, ncol = 4)

ggsave(path(supp_dir, "FigureS_QQ_all_traits.png"), qq_grid,
       width = 16, height = 8, dpi = 300, bg = "white")
