log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

phenotype_wc   <- snakemake@wildcards[["phenotype"]]
annotated_dirs <- unlist(snakemake@input[["annotated_dirs"]])
out_dir        <- snakemake@output[[1]]

# creamos el directorio de salida si no existe
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Split por phenotype: ", phenotype_wc, " | buscando en ", length(annotated_dirs), " batch(es)")

# Iteramos sobre los directorios de annotated_dirs y buscamos los archivos .rds
n_matched <- 0 # cuantos rds coinciden con el phenotype_wc
n_total   <- 0 # del  total

#un loop dentro de otro loop
for (annotated_dir in annotated_dirs) {  
#recibe una lista de carpetas 
  batch_name <- basename(annotated_dir) 
  rds_files <- list.files(annotated_dir, pattern = "\\.rds$", full.names = TRUE)
  message("Batch ", batch_name, ": ", length(rds_files), " specimen(s)")
  n_total <- n_total + length(rds_files)

  for (rds_file in rds_files) {  #recorre cada archivo rds dentro de la carpeta
    specimen_id <- tools::file_path_sans_ext(basename(rds_file))
    obj <- readRDS(rds_file)

    # Obtenemos el phenotype_label del specimen 
    #phenotype_label es uniforme por specimen (un specimen = un individualID)
    specimen_label <- unique(obj$phenotype_label)  #phenotype_label nos dice el fenotipo de cada celula
    specimen_label <- specimen_label[!is.na(specimen_label)]

    if (length(specimen_label) > 1) { # si hay mas de un phenotype_label no-NA para un specimen, se lanza una advertencia
      warning("  Specimen ", specimen_id,
              " (batch ", batch_name, ") tiene mas de un phenotype_label no-NA: ",
              paste(specimen_label, collapse = ", "))
    }

    specimen_label <- if (length(specimen_label) >= 1) specimen_label[1] else NA_character_

    # Si el phenotype_label del specimen coincide con el wildcard phenotype, lo guardamos en out_dir
    if (!is.na(specimen_label) && specimen_label == phenotype_wc) {
      dest_file <- file.path(out_dir, paste0(specimen_id, ".rds"))

      if (file.exists(dest_file)) {
        warning("  COLISION: ", specimen_id,
                " ya existe en ", out_dir,
                " (viene de otro batch con el mismo specimenID) — se va a sobreescribir.")
      }

      message("  Match: ", specimen_id, " (batch ", batch_name, ") -> ", phenotype_wc)
      saveRDS(obj, dest_file, compress = "gzip")
      n_matched <- n_matched + 1
    }

    rm(obj) #borramos el objeto de memoria para liberar espacio, ya que estamos procesando muchos archivos
  }
}

message("Phenotype ", phenotype_wc,
        " | specimens asignados: ", n_matched, " de ", n_total, " totales revisados")

sink(type = "message"); sink()