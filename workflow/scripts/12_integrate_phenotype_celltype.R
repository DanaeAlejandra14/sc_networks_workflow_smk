log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(harmony)
  library(patchwork)
})

phenotype_wc <- snakemake@wildcards[["phenotype"]]
cell_type_wc <- snakemake@wildcards[["cell_type"]]

merged_path <- snakemake@input[["merged"]]
output_path <- snakemake@output[["integrated"]]
plot_path   <- snakemake@output[["umap_plot"]]

obj <- readRDS(merged_path)

message("Phenotype: ", phenotype_wc, " | Cell type: ", cell_type_wc)
message("Celulas: ", ncol(obj), " | Individuos: ", length(unique(obj$individualID)))

# ---- ANTES: PCA normal (sin corregir nada) ----
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap_before")

if (length(unique(obj$individualID)) < 2) {
  message("ADVERTENCIA: solo hay 1 individuo. Se omite Harmony.")
  obj$harmony_skipped <- TRUE
  obj[["umap_after"]] <- obj[["umap_before"]]  # no hay "despues" real que mostrar
} else {
  message("Ejecutando integracion con Harmony, agrupando por individualID...")
  obj <- RunHarmony(obj, group.by.vars = "individualID")
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap_after")
  obj$harmony_skipped <- FALSE
}

# ---- Plot antes/despues ----
p_antes <- DimPlot(obj, reduction = "umap_before", group.by = "individualID") +
  ggtitle(paste0("ANTES (PCA)\n", phenotype_wc, " | ", cell_type_wc)) +
  theme(legend.position = "none")

p_despues <- DimPlot(obj, reduction = "umap_after", group.by = "individualID") +
  ggtitle(paste0("DESPUES (Harmony)\n", phenotype_wc, " | ", cell_type_wc)) +
  theme(legend.position = "none")

p_final <- p_antes | p_despues

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, p_final, width = 12, height = 6, dpi = 150, bg = "white")

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, output_path, compress = "gzip")

message("Objeto guardado en: ", output_path)
message("Plot guardado en: ", plot_path)

sink(type = "message"); sink()