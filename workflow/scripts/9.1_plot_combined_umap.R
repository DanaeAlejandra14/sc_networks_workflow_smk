log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

annotated_dirs <- snakemake@input[["annotated_dirs"]]

message("Cargando specimens de ", length(annotated_dirs), " batch(es)")

obj_list <- list()

for (dir in annotated_dirs) {
  library_batch <- basename(dir)
  rds_files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)

  for (rds_file in rds_files) {
    specimen_id <- tools::file_path_sans_ext(basename(rds_file))
    message("  Cargando ", library_batch, " / ", specimen_id)

    obj <- readRDS(rds_file)
    obj$library_batch <- library_batch
    obj$specimen_id <- specimen_id

    key <- paste0(library_batch, "_", specimen_id)
    obj_list[[key]] <- obj
  }
}

message("Total specimens cargados: ", length(obj_list))

# Merge sin corrección de batch: solo para inspección visual, no es integración formal
merged <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
merged <- JoinLayers(merged)

merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:30, verbose = FALSE)

p_celltype <- DimPlot(merged, group.by = "cell_type", label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("Todos los specimens - cell_type")

p_source <- DimPlot(merged, group.by = "annotation_source") +
  ggtitle("Todos los specimens - annotation_source")

p_batch <- DimPlot(merged, group.by = "library_batch") +
  ggtitle("Todos los specimens - library_batch (chequeo de batch effect)")

ggsave(snakemake@output[["umap_celltype"]], p_celltype, width = 8, height = 7, dpi = 150)
ggsave(snakemake@output[["umap_source"]], p_source, width = 8, height = 7, dpi = 150)
ggsave(snakemake@output[["umap_batch"]], p_batch, width = 8, height = 7, dpi = 150)

message("UMAPs combinados guardados en results/annotation_plots/")
sink(type = "message"); sink()