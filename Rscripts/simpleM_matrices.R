#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(readr))

# for setting the working directory when not run interactively using an r project 
library(rprojroot)
root <- rprojroot::find_root(rprojroot::is_rstudio_project)
setwd(root)

# ---------------------------
# Parameters
# ---------------------------
PCA_cutoff <- 0.995
mat_dir <- "meff/simpleM_mats"
out_dir <- file.path("meff", "simpleM_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# chromosome from arg1 or SLURM_ARRAY_TASK_ID
args <- commandArgs(trailingOnly = TRUE)
chr <- if (length(args) >= 1) args[1] else Sys.getenv("SLURM_ARRAY_TASK_ID")
if (is.na(chr) || chr == "") stop("No chromosome provided. Use arg1 or SLURM_ARRAY_TASK_ID.")
chr <- as.character(chr)

fn <- file.path(mat_dir, paste0("chr", chr, ".matrix.txt"))
if (!file.exists(fn)) stop("Matrix file not found: ", fn)

# ---------------------------
# Functions
# ---------------------------
Meff_PCA <- function(eigenValues, percentCut){
  totalEigenValues <- sum(eigenValues)
  myCut <- percentCut * totalEigenValues
  myEigenSum <- 0
  index_Eigen <- 0
  for (i in seq_along(eigenValues)) {
    if (myEigenSum <= myCut) {
      myEigenSum <- myEigenSum + eigenValues[i]
      index_Eigen <- i
    } else break
  }
  index_Eigen
}

read_snp_matrix <- function(fn){
  df <- read_table(
    file = fn,
    col_names = FALSE,
    col_types = cols(.default = col_integer()),
    progress = FALSE
  )
  X <- as.matrix(df)          # SNPs x individuals
  storage.mode(X) <- "integer"
  X
}

# ---------------------------
# Main
# ---------------------------
X <- read_snp_matrix(fn)  # rows=SNPs, cols=individuals
if (any(!is.finite(X))) stop("Matrix contains non-finite values: ", fn)

nSNP_total <- nrow(X)
nInd <- ncol(X)

# Drop zero-variance SNPs (prevents NA correlations)
var_snp <- apply(X, 1, var)
keep <- var_snp > 0
if (!all(keep)) X <- X[keep, , drop = FALSE]
nSNP_used <- nrow(X)

if (nSNP_used < 2) {
  Meff_chr <- 0L
} else {
  # cor() expects variables in columns; we want SNP-SNP correlation
  # so use individuals x SNPs
  G <- t(X)  # individuals x SNPs
  
  # Optional: scale columns for numerical stability (cor does this internally)
  # Compute SNP correlation matrix
  CLD <- cor(G)
  
  # Eigenvalues
  ev <- eigen(CLD, symmetric = TRUE, only.values = TRUE)$values
  ev <- abs(ev)
  
  # Meff via PCA cutoff
  Meff_chr <- Meff_PCA(ev, PCA_cutoff)
}

# Write key/value output
out_chr <- file.path(out_dir, paste0("chr", chr, ".out"))
writeLines(c(
  paste("FILE", fn, sep="\t"),
  paste("CHR", chr, sep="\t"),
  paste("PCA_cutoff", PCA_cutoff, sep="\t"),
  paste("NIND", nInd, sep="\t"),
  paste("NSNP_total", nSNP_total, sep="\t"),
  paste("NSNP_used", nSNP_used, sep="\t"),
  paste("Meff_chr", format(Meff_chr, digits=12), sep="\t")
), con = out_chr)

# Emit a single TSV row for easy concatenation
cat(paste(chr, fn, nInd, nSNP_total, nSNP_used, format(Meff_chr, digits=12), sep="\t"), "\n")