log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

batch_name <- snakemake@wildcards[["library_batch"]]
annotated_dir <- snakemake@input[["annotated_dir"]]
phenotype_path <- snakemake@input[["phenotype"]]
phenotyped_dir <- snakemake@output[["phenotyped_dir"]]

dir.create(phenotyped_dir, recursive = TRUE, showWarnings = FALSE)

message("Adjuntando fenotipo para batch: ", batch_name)

phenotype <- read_csv(phenotype_path, show_col_types = FALSE)
pheno_map <- setNames(phenotype$phenotype_label, phenotype$individualID)

rds_files <- list.files(annotated_dir, pattern = "\\.rds$", full.names = TRUE)
message("Found ", length(rds_files), " specimen(s) in this batch")

for (rds_file in rds_files) {
  specimen_id <- tools::file_path_sans_ext(basename(rds_file))
  message("Processing specimen: ", specimen_id)

  obj <- readRDS(rds_file)

  # individualID ya existe en cada celula desde el paso 1 (demultiplex_batch.R)
  obj$phenotype_label <- unname(pheno_map[as.character(obj$individualID)])

  message("  Cells: ", ncol(obj),
          " | Phenotype found: ", sum(!is.na(obj$phenotype_label)),
          " | Phenotype NA: ", sum(is.na(obj$phenotype_label)))

  saveRDS(obj, file.path(phenotyped_dir, paste0(specimen_id, ".rds")), compress = "gzip")
}

message("Fenotipo adjuntado para batch: ", batch_name)
sink(type = "message"); sink()