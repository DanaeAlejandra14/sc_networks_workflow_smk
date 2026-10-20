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
message("Cells: ", ncol(obj), " | Individuals: ", length(unique(obj$individualID)))

# iNTEGRATION 
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap_before")

if (length(unique(obj$individualID)) < 2) {
  message("WARNING: There is only 1 individual. Harmony is omitted.")
  obj$harmony_skipped <- TRUE
  obj[["umap_after"]] <- obj[["umap_before"]]  
} else {
  message("Executing integration with Harmony, grouping by individualID...")
  obj <- RunHarmony(obj, group.by.vars = "individualID")
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap_after")
  obj$harmony_skipped <- FALSE
}

#Plots before and after
p_antes <- DimPlot(obj, reduction = "umap_before", group.by = "individualID") +
  ggtitle(paste0("Before (PCA)\n", phenotype_wc, " | ", cell_type_wc)) +
  theme(legend.position = "none")

p_despues <- DimPlot(obj, reduction = "umap_after", group.by = "individualID") +
  ggtitle(paste0("After (Harmony)\n", phenotype_wc, " | ", cell_type_wc)) +
  theme(legend.position = "none")

p_final <- p_antes | p_despues

dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)
ggsave(plot_path, p_final, width = 12, height = 6, dpi = 150, bg = "white")

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, output_path, compress = "gzip")

message(" Save in ", output_path)
message("Plot save in: ", plot_path)

sink(type = "message"); sink()