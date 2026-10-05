log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(tidyverse)
library(Seurat)
library(harmony)

merged_path <- snakemake@input[[1]]
output_path <- snakemake@output[[1]]

message("Cargando objeto fusionado: ", merged_path)
obj <- readRDS(merged_path)

message("Numero de celulas totales: ", ncol(obj))
message("Numero de individuos (individualID) distintos: ", length(unique(obj$individualID)))

if (length(unique(obj$individualID)) < 2) {
  message("ADVERTENCIA: solo hay 1 individuo en esta combinacion phenotype x cell_type. ",
          "Se omite Harmony (no hay variacion inter-individuo que corregir) y se usa PCA directamente.")
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj$harmony_skipped <- TRUE
} else {
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)

  message("Ejecutando integracion con Harmony, agrupando por individualID...")
  obj <- RunHarmony(obj, group.by.vars = "individualID")
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30)
  obj$harmony_skipped <- FALSE
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(obj, output_path, compress = "gzip")

message("Objeto integrado guardado en: ", output_path)

sink(type = "message")
sink()