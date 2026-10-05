#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

option_list <- list(
  make_option(c("-p", "--phenotype"), type = "character",
              help = "Fenotipo a revisar (ej. 'AD-NC_SYM', 'MCI', 'control')"),
  make_option(c("-c", "--cell_type"), type = "character",
              help = "Tipo celular a revisar (ej. 'Excitatory Neurons')"),
  make_option(c("-o", "--out"), type = "character", default = NULL,
              help = "Ruta del PNG de salida (default: results/checks/<phenotype>_<cell_type>_before_after.png)")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$phenotype) || is.null(opt$cell_type)) {
  stop("Debes especificar --phenotype y --cell_type. Ejemplo:\n",
       "Rscript scripts/check_integration_before_after.R -p 'AD-NC_SYM' -c 'Excitatory Neurons'")
}

merged_path     <- file.path("results/phenotype", opt$phenotype, opt$cell_type, "merged.rds")
integrated_path <- file.path("results/phenotype", opt$phenotype, opt$cell_type, "integrated.rds")

if (!file.exists(merged_path))     stop("No existe: ", merged_path)
if (!file.exists(integrated_path)) stop("No existe: ", integrated_path)

out_png <- if (!is.null(opt$out)) {
  opt$out
} else {
  dir.create("results/checks", showWarnings = FALSE)
  file.path("results/checks",
             paste0(opt$phenotype, "_", gsub(" ", "_", opt$cell_type), "_before_after.png"))
}

message("Phenotype: ", opt$phenotype, " | Cell type: ", opt$cell_type)

# ---- ANTES: objeto fusionado, sin integrar ----
message("Cargando merged.rds y calculando PCA/UMAP normal...")
merged <- readRDS(merged_path)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)
merged <- RunUMAP(merged, reduction = "pca", dims = 1:30)

p_antes <- DimPlot(merged, reduction = "umap", group.by = "individualID") +
  ggtitle(paste0("ANTES (PCA normal)\n", opt$phenotype, " | ", opt$cell_type)) +
  theme(legend.position = "bottom")

# ---- DESPUES: objeto ya integrado con Harmony ----
message("Cargando integrated.rds...")
integrated <- readRDS(integrated_path)

if (isTRUE(integrated$harmony_skipped[1])) {
  message("ADVERTENCIA: esta combinacion tuvo Harmony omitido (solo 1 individuo). ",
          "El plot de 'despues' va a verse igual al PCA normal.")
}

p_despues <- DimPlot(integrated, reduction = "umap", group.by = "individualID") +
  ggtitle(paste0("DESPUES (Harmony)\n", opt$phenotype, " | ", opt$cell_type)) +
  theme(legend.position = "bottom")

p_final <- (p_antes | p_despues)

ggsave(out_png, p_final, width = 14, height = 7, dpi = 150, bg = "white")

message("Listo. Plot guardado en: ", out_png)