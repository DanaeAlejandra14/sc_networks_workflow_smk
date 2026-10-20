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

message("Attaching and filtering phenotype for batch: ", batch_name)

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

  # individualID alredy exists in the Seurat object metadata, we just need to map it to phenotype_label
  obj$phenotype_label <- unname(pheno_map[as.character(obj$individualID)])

  # 
  # an specimen= an individualID -> the phenotype_label should be uniform in the object
  unique_labels <- unique(obj$phenotype_label)
  non_na_labels <- unique_labels[!is.na(unique_labels)]

  if (length(non_na_labels) > 1) {
    warning("  Specimen ", specimen_id,
            " has more than one non-NA phenotype_label (mixed individualID?): ",
            paste(non_na_labels, collapse = ", "))
  }

  specimen_label <- if (length(non_na_labels) >= 1) non_na_labels[1] else NA_character_

  message("  Cells: ", ncol(obj),
          " | Specimen phenotype_label: ", specimen_label,
          " | Cells with NA phenotype: ", sum(is.na(obj$phenotype_label)))

  if (is.na(specimen_label)) {
    message("  -> Excluded: has no valid phenotype label: ", specimen_id)
    saveRDS(obj, file.path(excluded_dir, paste0(specimen_id, ".rds")), compress = "gzip")
    n_excluded <- n_excluded + 1
  } else {
    message("  -> Included: ", specimen_id)
    saveRDS(obj, file.path(filtered_dir, paste0(specimen_id, ".rds")), compress = "gzip")
    n_kept <- n_kept + 1
  }
}

message("Batch ", batch_name, " complete. Kept: ", n_kept,
        " |excluded: ", n_excluded)

sink(type = "message"); sink()