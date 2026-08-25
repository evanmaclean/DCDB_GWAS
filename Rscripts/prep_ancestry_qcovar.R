library(tidyverse)

qcovar <- read_delim(
  "gcta_files/qcovar.qcovar",
  delim = " ",
  col_names = c("FID", "IID", "age"),
  col_types = cols(FID = col_character(), IID = col_character(), age = col_double())
)

pca <- read_delim(
  "genomic_PCA/pca.eigenvec",
  delim = " ",
  col_names = c("FID", "IID", paste0("PC", 1:20)), # extra PC columns harmless if unused
  col_types = cols(.default = col_double(), FID = col_character(), IID = col_character())
)

qcovar_ancestry <- qcovar |>
  left_join(select(pca, FID, IID, PC1), by = c("FID", "IID"))

# How many dogs would be dropped/have missing PC1 if any - worth checking
# before trusting the run, since GCTA silently excludes rows with missing
# covariate values rather than erroring.
n_missing <- sum(is.na(qcovar_ancestry$PC1))
if (n_missing > 0) {
  message(n_missing, " of ", nrow(qcovar_ancestry),
          " dogs have no PC1 match - these will be dropped from the ancestry-adjusted run.")
}

write_delim(
  qcovar_ancestry,
  "gcta_files/qcovar_with_ancestry.qcovar",
  delim = " ",
  col_names = FALSE
)

message("Wrote gcta_files/qcovar_with_ancestry.qcovar (FID IID age PC1).")
