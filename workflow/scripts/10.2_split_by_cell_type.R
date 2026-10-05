log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

phenotype_wc  <- snakemake@wildcards[["phenotype"]]
cell_type_wc  <- snakemake@wildcards[["cell_type"]]
phenotype_dir <- snakemake@input[["phenotype_dir"]]
out_dir       <- snakemake@output[[1]]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#Ahora hay una corrida para cada combinacion phenotype x cell_type
message("Split por cell_type: ", cell_type_wc, " | phenotype: ", phenotype_wc)

rds_files <- list.files(phenotype_dir, pattern = "\\.rds$", full.names = TRUE)
message("Found ", length(rds_files), " specimen(s) en este phenotype")

n_saved <- 0 # cuantos rds coinciden con el cell_type_wc , idealmente todos los specimens 
#de phenotype_wc deberian tener al menos 1 celula de cell_type_wc

for (rds_file in rds_files) { #recorre cada rds
  specimen_id <- tools::file_path_sans_ext(basename(rds_file))
  #carga cada rds y extrae las células del cell_type_wc, 
  #las renombra con el prefijo del specimen_id y las guarda en la carpeta de salida
  obj <- readRDS(rds_file) 

  #compara cada una de las celulas del rds con el cell_type_wc, y devuelve un vector logico de TRUE/FALSE
  n_match <- sum(obj$cell_type == cell_type_wc, na.rm = TRUE) 

  if (n_match == 0) { #si ese specimen no tiene ninguna celula del cell_type_wc
    message("  ", specimen_id, ": 0 celulas de ", cell_type_wc, " -> se omite")
    next
  }

  # subset de Seurat para extraer las células del cell_type_wc, se queda con True cells
  obj_sub <- subset(obj, subset = cell_type == cell_type_wc)

  message("  ", specimen_id, ": ", ncol(obj_sub), " celulas de ", cell_type_wc)

  saveRDS(obj_sub, file.path(out_dir, paste0(specimen_id, ".rds")), compress = "gzip") #lo guarda con el nombre de specimen Id
  n_saved <- n_saved + 1
}

message("Phenotype ", phenotype_wc, " | cell_type ", cell_type_wc,
        " | specimens guardados: ", n_saved, " de ", length(rds_files))

sink(type = "message"); sink()