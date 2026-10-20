log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

phenotype_wc <- snakemake@wildcards[["phenotype"]]
cell_type_wc <- snakemake@wildcards[["cell_type"]]
input_dirs   <- unlist(snakemake@input)   # one folder per library_batch
out_file     <- snakemake@output[[1]]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

message("Merging phenotype: ", phenotype_wc, " | cell_type: ", cell_type_wc)

obj_list     <- list()   # list that stores the Seurat object of each specimen
specimen_ids <- c()      # specimen IDs, in the same order as obj_list

for (d in input_dirs) {
  rds_files <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  for (rds_file in rds_files) {
    specimen_id <- tools::file_path_sans_ext(basename(rds_file))
    message("  Loading: ", specimen_id, " (from ", d, ")")
    obj_list[[specimen_id]] <- readRDS(rds_file)
    specimen_ids <- c(specimen_ids, specimen_id)
  }
}

message("Total specimens found: ", length(obj_list))

# Stop if no Seurat object was found
if (length(obj_list) == 0) {
  stop("No specimen found for phenotype=", phenotype_wc,
       " cell_type=", cell_type_wc, ". Check that this combination exists in the data.")
}

# Specimen IDs must be unique, otherwise objects would be overwritten in obj_list
stopifnot("Duplicated specimen IDs across batches" = !anyDuplicated(specimen_ids))

# Objects are merged with Seurat's merge(). Cell names are prefixed with the
# specimen ID to avoid barcode collisions between different specimens.
# In Seurat v5 the merged object keeps one layer per specimen; layers are
# joined later (JoinLayers) when the expression matrix is needed.
if (length(obj_list) == 1) {
  # A single object cannot be merged: only rename its cells with the specimen ID prefix
  merged_obj <- RenameCells(obj_list[[1]], add.cell.id = specimen_ids[1])
} else {
  merged_obj <- merge(
    x = obj_list[[1]],          # first Seurat object
    y = obj_list[-1],           # all the remaining objects
    add.cell.ids = specimen_ids # prefix per specimen, so cell names are unique
  )
}

message("Final object: ", ncol(merged_obj), " cells, ",
        length(unique(merged_obj$individualID)), " distinct individualID(s)")

saveRDS(merged_obj, out_file, compress = "gzip")

message("Merge complete: ", out_file)
sink(type = "message"); sink()