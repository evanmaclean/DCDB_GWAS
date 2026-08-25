library(tidyverse)
library(fs)
library(glue)

collections <- c(
  "zwir_temperament_2018",
  "big5_schwaba",
  "DoGA_CNS_gene_sets",
  "epic_tissue_gene_sets"
)

collection_labels <- c(
  "zwir_temperament_2018" = "Human temperament (Zwir et al., 2020)",
  "big5_schwaba" = "Big Five personality (Schwaba et al., 2025)",
  "DoGA_CNS_gene_sets" = "CNS regional expression, DoGA (Hörtenhuber et al., 2024)",
  "epic_tissue_gene_sets" = "Tissue-specific expression, EPIC (Son et al., 2023)"
)

reference_task <- "HI_gaze" # arbitrary -- NGENES doesn't vary by task 

# --- full membership, as originally defined ---
# distinct() because the raw definition files contain literal duplicate (set, gene)
# rows (e.g. zwir_temperament_2018's "sensitive" set lists THEMIS twice) 

full_membership <- map_dfr(collections, function(coll) {
  path_in <- path("magma_files", "gene_sets", coll)
  read_tsv(path_in, col_names = c("set", "gene"), col_types = "cc") |>
    mutate(
      collection = coll,
      collection_label = collection_labels[[coll]],
      .before = 1
    ) |>
    distinct()
})

# --- NGENES actually analyzed by MAGMA, from one representative task per collection ---
analyzed_counts <- map_dfr(collections, function(coll) {
  gsa_path <- path(
    "magma_files",
    "magma_gene_set_analysis",
    coll,
    glue("{reference_task}.gsa.out")
  )
  gsa_raw <- read_table(gsa_path, comment = "#", col_types = cols(.default = "c"))
  # Use FULL_NAME, not VARIABLE, when it's present -- MAGMA truncates VARIABLE
  # with "..." for long set names (e.g. "caudate-nucleus_Tissue-Enric..."),
  # which silently fails to join against the untruncated set names from
  # full_membership. MAGMA only emits the FULL_NAME column at all when at
  # least one VARIABLE in that file needed truncating, so short-name
  # collections (Zwir, Big Five, and apparently EPIC) don't have it -- fall
  # back to VARIABLE for those.
  name_col <- if ("FULL_NAME" %in% names(gsa_raw)) "FULL_NAME" else "VARIABLE"
  gsa_raw |>
    transmute(
      collection = coll,
      set = .data[[name_col]],
      n_genes_analyzed = as.integer(NGENES)
    )
})

# --- combine: defined-set size vs. analyzed size, per set ---
set_sizes <- full_membership |>
  count(collection, collection_label, set, name = "n_genes_defined") |>
  left_join(analyzed_counts, by = c("collection", "set")) |>
  arrange(collection_label, desc(n_genes_defined))

print(set_sizes, n = Inf)

out_dir <- path("tables")
dir_create(out_dir)
write_excel_csv(set_sizes, path(out_dir, "gene_set_sizes_summary.csv"))

## Build a formatted table for the supplement (Table S3) ----
library(gt)

gene_set_gt <- set_sizes |>
  # EPIC: drop "_overall" and "_Group-Enriched" sets -- not interpreted/used in the paper
  filter(
    !(collection == "epic_tissue_gene_sets" &
        str_detect(set, "_overall$|_Group-Enriched$"))
  ) |>
  # DoGA CNS: set names are already tissue-restricted by definition, so the
  # "_Tissue-Enriched" suffix is redundant in the table
  mutate(
    set = if_else(
      collection == "DoGA_CNS_gene_sets",
      str_remove(set, "_Tissue-Enriched$"),
      set
    )
  ) |>
  select(collection_label, set, n_genes_defined, n_genes_analyzed) |>
  gt(groupname_col = "collection_label") |>
  cols_label(
    set = "Set",
    n_genes_defined = "N genes (defined)",
    n_genes_analyzed = "N genes (analyzed by MAGMA)"
  ) |>
  fmt_integer(columns = c(n_genes_defined, n_genes_analyzed)) |>
  sub_missing(columns = n_genes_analyzed, missing_text = "—")

gtsave(gene_set_gt, path(out_dir, "gene_set_sizes_table.docx"))

# --- full gene lists, deposited as supplementary data:

supp_dir <- path("supplemental")
dir_create(supp_dir)
full_membership_export <- full_membership |>
  filter(
    !(collection == "epic_tissue_gene_sets" &
        str_detect(set, "_overall$|_Group-Enriched$"))
  ) |>
  select(-collection) 
  write_tsv(
    full_membership_export,
    path(supp_dir, "Supplementary Data 2.tsv")
  )

# Sanity check: analyzed count should never exceed the defined count -- if it does,
# something is off in the join 
  
bad <- filter(set_sizes, n_genes_analyzed > n_genes_defined)
if (nrow(bad) > 0) {
  warning(
    "n_genes_analyzed > n_genes_defined for: ",
    paste(bad$set, collapse = ", ")
  )
} else {
  message(
    "Sanity check passed: n_genes_analyzed <= n_genes_defined for all sets."
  )
}
