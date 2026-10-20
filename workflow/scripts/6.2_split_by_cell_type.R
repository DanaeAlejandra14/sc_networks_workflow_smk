log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

phenotype_wc  <- snakemake@wildcards[["phenotype"]]
cell_type_wc  <- snakemake@wildcards[["cell_type"]]
phenotype_dir <- snakemake@input[["phenotype_dir"]]
out_dir       <- snakemake@output[[1]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# There is one run for each phenotype x cell_type combination
message("Split por cell_type: ", cell_type_wc, " | phenotype: ", phenotype_wc)

rds_files <- list.files(phenotype_dir, pattern = "\\.rds$", full.names = TRUE)
message("Found ", length(rds_files), " specimen(s) en este phenotype")

# Number of specimens that contain at least one cell of cell_type_wc
# (ideally, every specimen of this phenotype has at least one cell of this type)
n_saved <- 0 

for (rds_file in rds_files) { # iterate over each specimen
  specimen_id <- tools::file_path_sans_ext(basename(rds_file))
  obj <- readRDS(rds_file) 

  # Count the cells of this specimen annotated as cell_type_wc
  n_match <- sum(obj$cell_type == cell_type_wc, na.rm = TRUE) 

  if (n_match == 0) { # this specimen has no cells of the requested type
    message("  ", specimen_id, ": 0 celulas de ", cell_type_wc, " -> se omite")
    next
  }

  # Keep only the cells annotated as cell_type_wc
  obj_sub <- subset(obj, subset = cell_type == cell_type_wc)

  message("  ", specimen_id, ": ", ncol(obj_sub), " celulas de ", cell_type_wc)

  # Save one object per specimen, named after the specimen ID
  saveRDS(obj_sub, file.path(out_dir, paste0(specimen_id, ".rds")), compress = "gzip")
  n_saved <- n_saved + 1
}

message("Phenotype ", phenotype_wc, " | cell_type ", cell_type_wc,
        " | specimens saved: ", n_saved, " of ", length(rds_files))

sink(type = "message"); sink()