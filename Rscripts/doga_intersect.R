#!/usr/bin/env Rscript
# =============================================================================
# gwas_doga_intersect.R
#
# Intersect GWAS LD-expanded SNP lists with DoGA promoter and enhancer tracks.
#
# SLURM USAGE
# -------------------------------------------------------------------------
# Create a separate wrapper script (e.g. run_doga_intersect.sbatch):
#
#   #!/bin/bash
#   #SBATCH --job-name=doga_intersect
#   #SBATCH --output=logs/doga_intersect_%j.out
#   #SBATCH --error=logs/doga_intersect_%j.err
#   #SBATCH --ntasks=1
#   #SBATCH --cpus-per-task=1
#   #SBATCH --mem=8G
#   #SBATCH --time=01:00:00
#   #SBATCH --partition=your_partition
#
#   module load R/4.x.x
#   module load bedtools
#
#   Rscript /path/to/gwas_doga_intersect.R \
#     --project_root /path/to/project
#
# INTERACTIVE USAGE (for debugging)
# -------------------------------------------------------------------------
# Option 1 — source() from an R session after pre-defining opt:
#
#   opt <- list(project_root = ".",
#               snp_file     = "locus_snps.txt",
#               window       = 0L,
#               output_dir   = "doga_intersections")
#   source("gwas_doga_intersect.R")
#
# Option 2 — run Rscript with no arguments; parse_args() uses defaults
#   and the script runs from the current working directory.
#
# Option 3 — source() to load the functions, then call them manually.
#   This is the most flexible for step-by-step debugging.
#
# PROJECT STRUCTURE ASSUMED:
#   <project_root>/
#     DoGA_data/
#       enhancers/
#         enhancers_with_annotation_column_TE_bed/   <- tissue-specific BEDs
#         bed/_Enhancer_robust.bed                   <- all-tissue master BEDs
#         bed/_Enhancer_comp.bed
#       promoters/
#         Promoters_TE_bed/bed/                      <- tissue-specific BEDs
#         Promoters_with_annotation_column_TE_txt/   <- same regions + gene names
#     gcta_files/greml_estimates/locus_greml/
#       <trait>/
#         <locus>/
#           locus_snps.txt     <- LD-expanded SNPs, one per line as CHR:POS
#
# OUTPUT:  <project_root>/doga_intersections/
#   <trait>__<locus>__snp_hits.tsv         - per SNP x DoGA element hit detail
#   <trait>__<locus>__summary.tsv          - per track hit counts
#   all_traits_combined_summary.tsv        - long-format across all traits/loci
#   brain_endocrine_hits.tsv               - hits in brain/neuroendocrine tissues
#   brain_endocrine_gene_hits.tsv          - promoter hits collapsed to gene level
#
# DEPENDENCIES: bedtools on PATH; R packages: optparse, tidyverse, fs
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(fs)
})

# Set path to local bedtools
Sys.setenv(PATH = paste("/opt/ohpc/pub/apps/bedtools2/2.29.2/bin",
                        Sys.getenv("PATH"), sep = ":"))

# ── Argument parsing ───────────────────────────────────────────────────────────
# Only runs parse_args() when not interactive OR when `opt` hasn't been
# pre-defined (so interactive users can set opt manually and source() the file).
if (!interactive() || !exists("opt")) {
  option_list <- list(
    make_option("--project_root", type = "character", default = ".",
                help = "Path to project root [default: current directory]"),
    make_option("--snp_file", type = "character", default = "locus_snps.txt",
                help = "SNP filename inside each locus directory [default: locus_snps.txt]"),
    make_option("--window", type = "integer", default = 0L,
                help = "bp window around each SNP (0 = exact position) [default: 0]"),
    make_option("--output_dir", type = "character", default = "doga_intersections",
                help = "Output directory name at project root [default: doga_intersections]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
}

project_root <- path_real(opt$project_root)
output_dir   <- path(project_root, opt$output_dir)
dir_create(output_dir, recurse = TRUE)

# ── Check bedtools ─────────────────────────────────────────────────────────────
if (!nzchar(Sys.which("bedtools"))) {
  stop("bedtools not found on PATH. Load the module or install bedtools.")
}
message("bedtools: ", Sys.which("bedtools"))

# ── Tissue categories for priority filtering ───────────────────────────────────
BRAIN_TISSUES <- c(
  "frontal cortex_parietal cortex_temporal cortex_occipital cortex",
  "hippocampal formation",
  "caudate nucleus_amygdala",
  "substantia nigra_pons_medulla oblongata",
  "cerebellum hemisphere_cerebellum vermis",
  "corpus callosum",
  "olfactory bulb",
  "olfactory bulb_vomeronasal organ_nasal cavity olfactory mucosa",
  "piriform lobe",
  "spinal cord",
  "spinal ganglion",
  "diencephalon (ent nucleus reuniens)_ lateral geniculate nucleus",
  "neurohypophysis",
  "adenohypophysis"
)

ENDOCRINE_TISSUES <- c(
  "adrenal gland", "thyroid gland", "parathyroid gland",
  "neurohypophysis", "adenohypophysis",
  "ovary", "testis", "prostate gland", "uterus"
)

PRIORITY_TISSUES <- unique(c(BRAIN_TISSUES, ENDOCRINE_TISSUES))

# ══════════════════════════════════════════════════════════════════════════════
# FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

#' Parse a one-column SNP file (CHR:POS format) into a BED tibble.
#'
#' BED coordinates are 0-based half-open, so SNP at POS becomes [POS-1, POS].
#' The optional window expands this symmetrically.
parse_snp_file <- function(snp_file, window = 0L) {
  lines <- read_lines(snp_file) |>
    str_trim() |>
    discard(~ str_starts(.x, "#") | !nzchar(.x))

  if (length(lines) == 0L) return(NULL)

  parsed <- str_match(lines, "^(?:chr)?(\\w+):(\\d+)$")
  bad    <- which(is.na(parsed[, 1L]))

  if (length(bad) > 0L) {
    warning(sprintf("Skipping %d unrecognised lines in %s", length(bad), snp_file))
    parsed <- parsed[-bad, , drop = FALSE]
    lines  <- lines[-bad]
  }
  if (nrow(parsed) == 0L) return(NULL)

  tibble(
    chrom  = str_c("chr", parsed[, 2L]),
    start  = pmax(0L, as.integer(parsed[, 3L]) - 1L - window),
    end    = as.integer(parsed[, 3L]) + window,
    snp_id = lines
  )
}


#' Collect metadata for every DoGA BED track under doga_dir.
collect_doga_tracks <- function(doga_dir) {
  bed_files <- dir_ls(doga_dir, glob = "*.bed", recurse = TRUE) |>
    keep(~ file_size(.x) > 100)

  map_dfr(bed_files, function(p) {
    name <- path_file(p)

    element_type <- case_when(
      str_detect(str_to_lower(name), "enhancer") ~ "enhancer",
      str_detect(str_to_lower(name), "promoter") ~ "promoter",
      .default = "unknown"
    )

    specificity <- case_when(
      str_detect(str_to_lower(name), "robust") ~ "robust",
      str_detect(str_to_lower(name), "comp")   ~ "comp",
      .default = "unknown"
    )

    # Tissue from filename e.g. Enhancer_robust_hippocampal formation_enriched.bed
    tissue_m <- str_match(name, "(?:robust|comp)_(.+?)(?:_enriched)?\\.bed$")
    tissue <- if (!is.na(tissue_m[1L, 2L])) {
      tissue_m[1L, 2L]
    } else if (str_starts(name, "_")) {
      "ALL_TISSUES"
    } else {
      "unknown"
    }

    tibble(path = p, element_type, specificity, tissue, filename = name)
  })
}


#' Build a lookup table from the promoter TXT files that carries gene names.
#'
#' TXT files have an extra leading column (promoter_name e.g. "p6\@NETO1"),
#' then chr/start/end/coords/score/strand. The coords field matches the
#' track_name (name field) in the BED files, so we use it as the join key.

build_promoter_gene_lookup <- function(doga_dir) {
  txt_dir <- path(doga_dir, "promoters",
                  "Promoters_with_annotation_column_TE_txt")

  if (!dir_exists(txt_dir)) {
    message("  [NOTE] Promoter TXT directory not found; gene names unavailable.")
    return(NULL)
  }

  txt_files <- dir_ls(txt_dir, glob = "*.txt")
  if (length(txt_files) == 0L) return(NULL)

  message(sprintf("  Building gene-name lookup from %d promoter TXT files...",
                  length(txt_files)))

  map_dfr(txt_files, function(f) {
    tryCatch(
      read_tsv(f,
               col_names = c("promoter_name", "chrom", "start",
                             "end", "coords", "score", "strand"),
               col_types  = "cciicdc",
               skip       = 1L,
               progress   = FALSE) |>
        filter(!is.na(chrom)) |>
        transmute(
          coords_key    = coords,
          promoter_name = promoter_name,
          gene_name     = str_extract(promoter_name, "(?<=@).+")
        ),
      error = function(e) NULL
    )
  }) |>
    # Same promoter appears across multiple tissue files; keep one row
    distinct(coords_key, .keep_all = TRUE)
}


#' Run bedtools intersect; return a tibble of hits or NULL if none.
bedtools_intersect <- function(snp_bed, track_bed) {
  cmd <- sprintf(
    "bedtools intersect -a %s -b %s -wa -wb 2>/dev/null",
    shQuote(snp_bed), shQuote(track_bed)
  )
  raw <- system(cmd, intern = TRUE)
  if (length(raw) == 0L) return(NULL)

  # SNP cols  : chrom(1) start(2) end(3) snp_id(4)
  # Track cols: chrom(5) start(6) end(7) name(8)  score(9) strand(10)
  str_split_fixed(raw, "\t", n = 10L) |>
    as_tibble(.name_repair = ~ c("snp_chrom", "snp_start", "snp_end", "snp_id",
                                 "track_chrom", "track_start", "track_end",
                                 "track_name", "track_score", "track_strand")) |>
    mutate(
      snp_pos     = as.integer(snp_end),     # BED end = 1-based SNP position
      track_start = as.integer(track_start),
      track_end   = as.integer(track_end)
    ) |>
    select(-snp_start, -snp_end)
}


#' Discover all SNP files by walking the locus_greml directory tree.
discover_snp_files <- function(project_root, snp_filename = "locus_snps.txt") {
  locus_root <- path(project_root, "gcta_files",
                     "greml_estimates", "locus_greml")

  if (!dir_exists(locus_root)) {
    stop("locus_greml directory not found: ", locus_root)
  }

  snp_paths <- dir_ls(locus_root,
                      regexp   = str_c(snp_filename, "$"),
                      recurse  = TRUE)

  if (length(snp_paths) == 0L) {
    stop("No '", snp_filename, "' files found under ", locus_root)
  }

  # Path components relative to locus_greml/ are: <trait>/<locus>/<file>
  map_dfr(snp_paths, function(p) {
    parts <- path_split(path_rel(p, locus_root))[[1L]]
    tibble(trait = parts[1L], locus = parts[2L], snp_path = p)
  })
}


#' Core function: intersect one trait x locus against all DoGA tracks.
run_intersection <- function(trait, locus, snp_path,
                             tracks, promoter_gene_lookup,
                             window, output_dir) {

  tag <- str_c(trait, "__", locus)
  message(sprintf("\n── %s ──", tag))

  # Parse SNPs
  snp_df <- parse_snp_file(snp_path, window = window)
  if (is.null(snp_df) || nrow(snp_df) == 0L) {
    message("  No valid SNPs; skipping.")
    return(invisible(NULL))
  }
  message(sprintf("  %d SNPs loaded.", nrow(snp_df)))

  # Temporary BED file, deleted on function exit
  tmp_bed <- file_temp(ext = "bed")
  on.exit(file_delete(tmp_bed), add = TRUE)
  write_tsv(snp_df, tmp_bed, col_names = FALSE)

  # Intersect every track
  raw_results <- pmap(tracks, function(path, element_type, specificity,
                                       tissue, filename) {
    hits <- bedtools_intersect(tmp_bed, path)
    list(
      summary = tibble(
        trait, locus, tissue, element_type, specificity,
        n_hits        = if (!is.null(hits)) nrow(hits) else 0L,
        n_unique_snps = if (!is.null(hits)) n_distinct(hits$snp_id) else 0L,
        track_file    = filename
      ),
      hits = if (!is.null(hits) && nrow(hits) > 0L)
        mutate(hits, trait, locus, tissue, element_type, specificity,
               track_file = filename)
      else NULL
    )
  })

  summary_df <- map_dfr(raw_results, "summary") |> arrange(desc(n_hits))
  hits_df    <- map_dfr(raw_results, "hits")

  # Report hits to console
  summary_df |>
    filter(n_hits > 0L) |>
    pwalk(function(tissue, element_type, specificity, n_hits, n_unique_snps, ...) {
      message(sprintf("  HIT  %-10s | %-7s | %-50s | %d hits (%d SNPs)",
                      element_type, specificity,
                      str_trunc(tissue, 50L, ellipsis = ""),
                      n_hits, n_unique_snps))
    })

  # Attach gene names for promoter hits
  if (nrow(hits_df) > 0L) {
    if (!is.null(promoter_gene_lookup) && nrow(promoter_gene_lookup) > 0L) {
      hits_df <- left_join(
        hits_df,
        select(promoter_gene_lookup, coords_key, promoter_name, gene_name),
        by = c("track_name" = "coords_key")
      )
    } else {
      hits_df <- mutate(hits_df,
                        promoter_name = NA_character_,
                        gene_name     = NA_character_)
    }

    write_tsv(hits_df, path(output_dir, str_c(tag, "__snp_hits.tsv")))
    message(sprintf("  Hit detail written (%d rows).", nrow(hits_df)))
  } else {
    message("  No hits with any DoGA track.")
  }

  write_tsv(summary_df, path(output_dir, str_c(tag, "__summary.tsv")))

  list(hits = hits_df, summary = summary_df)
}

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
message("\n=== DoGA x GWAS Intersection Pipeline ===")
message("  Project root : ", project_root)
message("  Output dir   : ", output_dir)
message("  SNP file     : ", opt$snp_file)
message("  Window       : ", opt$window, " bp\n")

doga_dir <- path(project_root, "DoGA_data")
if (!dir_exists(doga_dir)) stop("DoGA_data not found under ", project_root)

message("Collecting DoGA tracks...")
tracks <- collect_doga_tracks(doga_dir)
message(sprintf("  %d track files found.", nrow(tracks)))

promoter_gene_lookup <- build_promoter_gene_lookup(doga_dir)
if (!is.null(promoter_gene_lookup)) {
  message(sprintf("  Gene lookup: %d unique promoters.\n",
                  nrow(promoter_gene_lookup)))
}

message("Discovering SNP files...")
snp_manifest <- discover_snp_files(project_root, snp_filename = opt$snp_file)
message(sprintf("  %d trait x locus combinations found:", nrow(snp_manifest)))
print(select(snp_manifest, trait, locus))

# Run all intersections
all_results <- pmap(snp_manifest, run_intersection,
                    tracks               = tracks,
                    promoter_gene_lookup = promoter_gene_lookup,
                    window               = opt$window,
                    output_dir           = output_dir)

# ── Combined outputs ───────────────────────────────────────────────────────────
message("\n── Writing combined outputs ──")

all_summaries <- map_dfr(all_results, "summary")
if (nrow(all_summaries) > 0L) {
  write_tsv(all_summaries,
            path(output_dir, "all_traits_combined_summary.tsv"))
  message("  all_traits_combined_summary.tsv written.")
}

all_hits <- map_dfr(all_results, "hits")
if (nrow(all_hits) > 0L) {
  priority_hits <- filter(all_hits, tissue %in% PRIORITY_TISSUES)

  if (nrow(priority_hits) > 0L) {
    write_tsv(priority_hits, path(output_dir, "brain_endocrine_hits.tsv"))
    message(sprintf("  brain_endocrine_hits.tsv: %d rows, %d unique SNPs.",
                    nrow(priority_hits), n_distinct(priority_hits$snp_id)))

    if ("gene_name" %in% names(priority_hits)) {
      gene_summary <- priority_hits |>
        filter(!is.na(gene_name) | !is.na(promoter_name)) |>   # keep LOC entries too
        mutate(
          gene_label = if_else(!is.na(gene_name), gene_name, promoter_name)
        ) |>
        summarise(
          genes      = str_c(sort(unique(gene_label)), collapse = "; "),
          n_elements = n(),
          .by = c(trait, locus, snp_id, tissue, element_type, specificity)
        ) |>
        arrange(trait, locus, tissue)

      write_tsv(gene_summary,
                path(output_dir, "brain_endocrine_gene_hits.tsv"))
      message("  brain_endocrine_gene_hits.tsv written.")
    }
  } else {
    message("  No hits in brain/endocrine tissues.")
  }
}

message("\nDone.")
