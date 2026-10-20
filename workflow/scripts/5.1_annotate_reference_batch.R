log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

batch_name <- snakemake@wildcards[["library_batch"]]
norm_dir <- snakemake@input[["norm_dir"]]
annotated_dir <- snakemake@output[["annotated_dir"]]

reference_csv <- snakemake@config[["annotation"]][["reference_csv"]]

dir.create(annotated_dir, recursive = TRUE, showWarnings = FALSE)

message("Running reference annotation for batch: ", batch_name)

cell_annotation <- readr::read_csv(reference_csv, show_col_types = FALSE)


# We construct a lookup table (annotation_map) from the reference CSV file, which maps each cell ID to its corresponding cell type and state. This will be used to annotate the cells in the Seurat objects for the current batch.
#Converts the CSV into a lookup table: cell_id -> cell.type/state
annotation_map <- cell_annotation %>%
  dplyr::select(cell, cell.type, state) %>%
  dplyr::distinct() %>%
  tibble::column_to_rownames("cell") # los convierte en el nombre de la fila del dataframe, para poder acceder a ellos por el nombre de la celula

rds_files <- list.files(norm_dir, pattern = "\\.rds$", full.names = TRUE)
message("Found ", length(rds_files), " specimen(s) in this batch")

for (rds_file in rds_files) {

  specimen_id <- tools::file_path_sans_ext(basename(rds_file))
  message("Processing specimen: ", specimen_id)

  obj <- readRDS(rds_file)

  # Reconstruct the cell keys for the current Seurat object by concatenating the batch name and cell barcodes. This will be used to look up the corresponding cell type and state in the annotation_map.
  cell_key <- paste0(batch_name, "_", colnames(obj))

   
   # The match is done by using the cell_key to index into the annotation_map, 
   # which returns the corresponding cell type and state for each cell. If a cell_key does not exist in the annotation_map,
   # R will return NA automatically, resulting in "unassigned" cells.
  obj$cell_type  <- annotation_map[cell_key, "cell.type", drop = TRUE]
  obj$cell_state <- annotation_map[cell_key, "state", drop = TRUE]

  message("  Cells: ", ncol(obj),
          " | Annotated: ", sum(!is.na(obj$cell_type)),
          " | Unassigned: ", sum(is.na(obj$cell_type)))

  saveRDS(obj, file.path(annotated_dir, paste0(specimen_id, ".rds")), compress = "gzip")
}

message("Reference annotation complete for batch: ", batch_name)
sink(type = "message"); sink()