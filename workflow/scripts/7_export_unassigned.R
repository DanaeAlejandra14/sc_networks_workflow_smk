log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)
library(SeuratDisk) #Convertir a h5ad

annotated_dirs <- unlist(snakemake@input) 
h5ad_path <- snakemake@output[["h5ad"]]

dir.create(dirname(h5ad_path), recursive = TRUE, showWarnings = FALSE)

message("Collecting unassigned cells from ", length(annotated_dirs), " batch folder(s)")

unassigned_objs <- list() # es un lopp de batches , no sabemos cuantos batches hay, por eso es una lista, para ir guardando los objetos de cada batch
##loop anidado para recorrer cada batch y dentro de cada batch recorrer cada rds

for (batch_dir in annotated_dirs) {
  rds_files <- list.files(batch_dir, pattern = "\\.rds$", full.names = TRUE)

#carga cada rds y extrae las células unassigned, las renombra con el prefijo del specimen_id y las guarda en la lista unassigned_objs
  for (rds_file in rds_files) {
    specimen_id <- tools::file_path_sans_ext(basename(rds_file))
    obj <- readRDS(rds_file)
    un_cells <- colnames(obj)[is.na(obj$cell_type)]

    if (length(un_cells) > 0) { #algunos specimens pueden no tener células unassignedn
      message("  ", specimen_id, ": ", length(un_cells), " unassigned cell(s)")

      obj_un <- subset(obj, cells = un_cells) #subset de Seurat para extraer las células unassigned
      
      # Prefijo con specimen_id antes de fusionar, para evitar colisiones
      # de barcodes entre specimens distintos (mismo barcode 10x, distinta célula real)
      obj_un <- RenameCells(obj_un, add.cell.id = specimen_id)
      unassigned_objs[[specimen_id]] <- obj_un
    }
  }
}

if (length(unassigned_objs) == 0) {
  stop("No unassigned cells found across any batch. Nothing to export.")
}

message("Merging unassigned cells into a single object...") 
merged_unassigned <- if (length(unassigned_objs) == 1) { #se fusionan los objetos en una lista en uno solo 
  unassigned_objs[[1]]
} else {
  merge(unassigned_objs[[1]], unassigned_objs[-1])
}

# Join layers los combina en una sola capa de data 
# Join layers los combina en una sola capa de data 
merged_unassigned[["RNA"]] <- JoinLayers(merged_unassigned[["RNA"]]) 
message("Total unassigned cells to export: ", ncol(merged_unassigned))

## let see 
print(sapply(merged_unassigned@meta.data, class))
saveRDS(merged_unassigned, "results/unassigned_export/merged_debug.rds")

# Convertir Assay5 (Seurat v5) a Assay clásico, compatible con SeuratDisk
merged_unassigned[["RNA"]] <- as(object = merged_unassigned[["RNA"]], Class = "Assay")

#convertir archivo a h5ad para poder abrirlo en Python
h5seurat_path <- sub("\\.h5ad$", ".h5Seurat", h5ad_path)
message("Writing h5Seurat: ", h5seurat_path)
SaveH5Seurat(merged_unassigned, filename = h5seurat_path, overwrite = TRUE)

message("Converting to h5ad: ", h5ad_path)
Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)

message("Export complete.")
sink(type = "message"); sink()