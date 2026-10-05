writeLines(paste("EMPECE A CORRER:", Sys.time()), "/tmp/debug_attach_filter.txt")

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

batch_name <- snakemake@wildcards[["library_batch"]]
specimens_dir <- snakemake@input[["specimens_dir"]]
phenotype_path <- snakemake@input[["phenotype"]]
filtered_dir <- snakemake@output[["filtered_dir"]]
excluded_dir <- snakemake@output[["excluded_dir"]]

dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(excluded_dir, recursive = TRUE, showWarnings = FALSE)

message("Adjuntando y filtrando fenotipo para batch: ", batch_name)

phenotype <- read_csv(phenotype_path, show_col_types = FALSE)
pheno_map <- setNames(phenotype$phenotype_label, phenotype$individualID)

rds_files <- list.files(specimens_dir, pattern = "\\.rds$", full.names = TRUE)
message("Found ", length(rds_files), " specimen(s) in this batch")

n_kept <- 0
n_excluded <- 0

for (rds_file in rds_files) {
  specimen_id <- tools::file_path_sans_ext(basename(rds_file))
  message("Processing specimen: ", specimen_id)

  obj <- readRDS(rds_file)

  # individualID ya existe en cada celula desde el paso 1 (demultiplex_batch.R)
  obj$phenotype_label <- unname(pheno_map[as.character(obj$individualID)])

  # un specimen = un individualID -> el phenotype_label debe ser uniforme en el objeto
  unique_labels <- unique(obj$phenotype_label)
  non_na_labels <- unique_labels[!is.na(unique_labels)]

  if (length(non_na_labels) > 1) {
    warning("  Specimen ", specimen_id,
            " tiene mas de un phenotype_label no-NA (individualID mixto?): ",
            paste(non_na_labels, collapse = ", "))
  }

  specimen_label <- if (length(non_na_labels) >= 1) non_na_labels[1] else NA_character_

  message("  Cells: ", ncol(obj),
          " | Specimen phenotype_label: ", specimen_label,
          " | Cells with NA phenotype: ", sum(is.na(obj$phenotype_label)))

  if (is.na(specimen_label)) {
    message("  -> Excluido (sin fenotipo valido): ", specimen_id)
    saveRDS(obj, file.path(excluded_dir, paste0(specimen_id, ".rds")), compress = "gzip")
    n_excluded <- n_excluded + 1
  } else {
    message("  -> Conservado: ", specimen_id)
    saveRDS(obj, file.path(filtered_dir, paste0(specimen_id, ".rds")), compress = "gzip")
    n_kept <- n_kept + 1
  }
}

message("Batch ", batch_name, " completo. Conservados: ", n_kept,
        " | Excluidos: ", n_excluded)

sink(type = "message"); sink()