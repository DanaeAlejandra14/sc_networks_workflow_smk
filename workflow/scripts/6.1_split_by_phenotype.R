log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

phenotype_wc   <- snakemake@wildcards[["phenotype"]]
annotated_dirs <- unlist(snakemake@input[["annotated_dirs"]])
out_dir        <- snakemake@output[[1]]

# Craete the output directory if it doesn't exist
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Split por phenotype: ", phenotype_wc, " | search in ", length(annotated_dirs), " batch(es)")

# We will count how many specimens match the phenotype_wc and how many total specimens we process
n_matched <- 0 # rds matched with the phenotype_wc
n_total   <- 0 # out of the total specimens processed (across all batches)

# Loop through each annotated_dir (each batch) and process the RDS files
for (annotated_dir in annotated_dirs) {  
# annotated_dir is a directory containing RDS files for a specific batch of specimens. We will iterate through each directory to find and process the RDS files.
  batch_name <- basename(annotated_dir) 
  rds_files <- list.files(annotated_dir, pattern = "\\.rds$", full.names = TRUE)
  message("Batch ", batch_name, ": ", length(rds_files), " specimen(s)")
  n_total <- n_total + length(rds_files)

  for (rds_file in rds_files) {  # Iterate through each RDS file in the batch directory to process individual specimens
    specimen_id <- tools::file_path_sans_ext(basename(rds_file))
    obj <- readRDS(rds_file)

    # We have the phenotype_label del specimen 
    #phenotype_label its consistent (un specimen = un individualID)
    specimen_label <- unique(obj$phenotype_label)  #phenotype_label nos dice el fenotipo de cada celula
    specimen_label <- specimen_label[!is.na(specimen_label)]

    if (length(specimen_label) > 1) { # If there is more than one phenotype to one SpecimenId
      warning("  Specimen ", specimen_id,
              " (batch ", batch_name, ") has more than one Pheotype, no-NA: ",
              paste(specimen_label, collapse = ", "))
    }

    specimen_label <- if (length(specimen_label) >= 1) specimen_label[1] else NA_character_

    # If  phenotype_label of the specimen mayches with the wildcard phenotype, we save it in out_dir
  
    if (!is.na(specimen_label) && specimen_label == phenotype_wc) {
      dest_file <- file.path(out_dir, paste0(specimen_id, ".rds"))

      if (file.exists(dest_file)) {
      stop("COLLISION: ", specimen_id, " already exists in ", out_dir,
       " (it comes from another batch with the same specimenID). ",
       "Refusing to overwrite; check the library names.")
}

      message("  Match: ", specimen_id, " (batch ", batch_name, ") -> ", phenotype_wc)
      saveRDS(obj, dest_file, compress = "gzip")
      n_matched <- n_matched + 1
    }

    rm(obj) #Erase the previous rds of the memory
  }   # end of loop over specimens
}     # end of loop over batches

message("Phenotype ", phenotype_wc,
        " | assigned specimens ", n_matched, " out of ", n_total, " totales checked")

sink(type = "message"); sink()