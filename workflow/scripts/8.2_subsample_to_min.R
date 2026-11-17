# For ONE cell type: equalize the number of cells across phenotypes.
#   N = number of usable cells of the smallest phenotype (cells with a donor).
#   The smallest phenotype is kept complete; the others are cut down to N.
#   The cut is proportional by donor: every donor of a phenotype loses the same
#   fraction (donor with 50 cells and donor with 30, remove 10% -> 5 and 3).
#   Rounding uses the largest-remainder method, so the total is exactly N.
# Harmony/PCA/UMAP embeddings of the kept cells are not recomputed.

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

in_rds     <- unlist(snakemake@input[["integrated"]])
out_rds    <- unlist(snakemake@output[["subsampled"]])
out_tsv    <- snakemake@output[["report"]]
phenotypes <- unlist(snakemake@params[["phenotypes"]])
donor_col  <- snakemake@params[["donor_col"]]
seed       <- as.integer(snakemake@params[["seed"]])
min_cells_warn <- snakemake@params[["min_cells_warn"]]

# Inputs and outputs come from expand() over the same phenotype list: check the pairing
stopifnot(length(in_rds) == length(phenotypes), length(out_rds) == length(phenotypes))
stopifnot(basename(dirname(dirname(in_rds)))  == phenotypes)
stopifnot(basename(dirname(dirname(out_rds))) == phenotypes)

set.seed(seed)

# 1. Read every phenotype and find its usable cells
objs <- list()
donors_by_pheno <- list()   # named vector: cell -> donor, only cells with a donor
n_dropped <- integer(0)

for (i in seq_along(phenotypes)) {
  p <- phenotypes[i]
  message("Reading: ", in_rds[i])
  obj <- readRDS(in_rds[i])
  if (!donor_col %in% colnames(obj@meta.data)) {
    stop("Column '", donor_col, "' not found in meta.data of ", in_rds[i])
  }
  donor <- as.character(obj@meta.data[[donor_col]])
  names(donor) <- rownames(obj@meta.data)

  # Cells without donor cannot be allocated: drop them BEFORE computing N
  bad <- is.na(donor)
  n_dropped[p] <- sum(bad)
  if (any(bad)) warning(p, ": ", sum(bad), " cells without '", donor_col, "' are dropped.")

  objs[[p]] <- obj
  donors_by_pheno[[p]] <- donor[!bad]
}

n_usable <- vapply(donors_by_pheno, length, integer(1))
message("Usable cells per phenotype before subsampling:")
print(n_usable)

n_target <- min(n_usable)
message("N target (smallest phenotype) = ", n_target)
if (n_target < min_cells_warn) {
  warning("N target is low (", n_target, " cells per phenotype). Check this cell type.")
}

#  2. Proportional allocation helper 
# Distribute 'target' cells among donors in proportion to their counts.
# Largest-remainder method; ties are broken at random (uses the seed above).
alloc_proportional <- function(counts, target) {
  q    <- counts * target / sum(counts)
  base <- floor(q)
  left <- target - sum(base)
  if (left > 0) {
    frac <- q - base
    ord  <- order(-frac, runif(length(frac)))
    add  <- ord[seq_len(left)]
    base[add] <- base[add] + 1
  }
  base
}

#  3. Subsample each phenotype and save
report <- list()

for (i in seq_along(phenotypes)) {
  p      <- phenotypes[i]
  donor  <- donors_by_pheno[[p]]
  counts <- table(donor)

  if (length(donor) <= n_target) {
    alloc <- as.integer(counts)                       # smallest phenotype: keep all
  } else {
    alloc <- alloc_proportional(as.integer(counts), n_target)
  }
  names(alloc) <- names(counts)

  keep <- character(0)
  for (d in names(counts)) {
    cells_d <- names(donor)[donor == d]
    n_keep  <- alloc[[d]]
    keep <- c(keep, cells_d[sample.int(length(cells_d), n_keep)])
  }

  # Checks
  stopifnot(length(keep) == n_target, !anyDuplicated(keep))
  if (any(alloc == 0)) {
    warning(p, ": donor(s) with 0 cells after subsampling: ",
            paste(names(alloc)[alloc == 0], collapse = ", "))
  }

  report[[p]] <- tibble(
    phenotype    = p,
    donor        = names(counts),
    cells_before = as.integer(counts),
    cells_after  = as.integer(alloc[names(counts)]),
    cells_dropped_no_donor = n_dropped[[p]]
  ) %>%
    mutate(fraction_before = cells_before / sum(cells_before),
           fraction_after  = cells_after  / sum(cells_after),
           fraction_kept   = cells_after  / cells_before)

  obj_sub <- subset(objs[[p]], cells = keep)
  dir.create(dirname(out_rds[i]), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj_sub, out_rds[i], compress = "gzip")
  message("Saved: ", out_rds[i], " (", ncol(obj_sub), " cells)")
  objs[[p]] <- NULL   # free memory
}

report <- bind_rows(report)
message("Max change in a donor's share within its phenotype: ",
        round(max(abs(report$fraction_after - report$fraction_before)), 4))

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
write_tsv(report, out_tsv)
message("Saved report: ", out_tsv)

sink(type = "message"); sink()
