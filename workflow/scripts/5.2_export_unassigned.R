log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)
library(SeuratDisk) # To convert a H5Seurat object to h5ad format for Python compatibility

annotated_dirs <- unlist(snakemake@input) 
h5ad_path <- snakemake@output[["h5ad"]]

dir.create(dirname(h5ad_path), recursive = TRUE, showWarnings = FALSE)

message("Collecting unassigned cells from ", length(annotated_dirs), " batch folder(s)")

unassigned_objs <- list() # For loop to process each batch and within each batch, process each RDS file
## another for loop to process each RDS file in the batch folder, extract unassigned cells, rename them with specimen_id prefix, and store them in the unassigned_objs list

for (batch_dir in annotated_dirs) {
  rds_files <- list.files(batch_dir, pattern = "\\.rds$", full.names = TRUE)

# For loop to process each RDS file in the batch folder, extract unassigned cells, rename them with specimen_id prefix, 
# and store them in the unassigned_objs list
  for (rds_file in rds_files) {
    specimen_id <- tools::file_path_sans_ext(basename(rds_file))
    obj <- readRDS(rds_file)
    un_cells <- colnames(obj)[is.na(obj$cell_type)]

    if (length(un_cells) > 0) { #some specimens may not have unassigned cells, so we only process those that do
      message("  ", specimen_id, ": ", length(un_cells), " unassigned cell(s)")

      obj_un <- subset(obj, cells = un_cells) #subet the Seurat object to only include unassigned cells
      
      # Prefix the cell barcodes with specimen_id to avoid barcode collisions across different specimens (same 10x barcode, different real cell)
    
      obj_un <- RenameCells(obj_un, add.cell.id = specimen_id)
      unassigned_objs[[specimen_id]] <- obj_un
    }
  }
}

if (length(unassigned_objs) == 0) {
  stop("No unassigned cells found across any batch. Nothing to export.")
}

message("Merging unassigned cells into a single object...") 
merged_unassigned <- if (length(unassigned_objs) == 1) { # the merge function requires at least 2 objects, so if there is only one object, we just take it as is
  unassigned_objs[[1]]
} else {
  merge(unassigned_objs[[1]], unassigned_objs[-1])
}

# Join layers containing counts and other assays into a single layer to ensure compatibility with SeuratDisk and h5ad export. This is necessary because Seurat v5 uses multiple layers for different assays, while h5ad expects a single layer.
merged_unassigned[["RNA"]] <- JoinLayers(merged_unassigned[["RNA"]]) 
message("Total unassigned cells to export: ", ncol(merged_unassigned))

print(sapply(merged_unassigned@meta.data, class))
saveRDS(merged_unassigned, "results/unassigned_export/merged_debug.rds")

# Convert the RNA assay to the "Assay" class to ensure compatibility with SeuratDisk and h5ad export. This is necessary because Seurat v5 uses a different class for assays, while h5ad expects a specific structure.
merged_unassigned[["RNA"]] <- as(object = merged_unassigned[["RNA"]], Class = "Assay")

#convertir archivo a h5ad para poder abrirlo en Python
h5seurat_path <- sub("\\.h5ad$", ".h5Seurat", h5ad_path)
message("Writing h5Seurat: ", h5seurat_path)
SaveH5Seurat(merged_unassigned, filename = h5seurat_path, overwrite = TRUE)

message("Converting to h5ad: ", h5ad_path)
Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)

message("Export complete.")
sink(type = "message"); sink()