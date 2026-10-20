log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

batch_name <- snakemake@wildcards[["library_batch"]]
doublet_dir <- snakemake@input[["doublet_dir"]]
norm_dir <- snakemake@output[["norm_dir"]]

#The parameter scale_factor is obtained from the configuration file and is used in the NormalizeData function to scale the normalized data.
scale_factor <- snakemake@config[["normalization"]][["scale_factor"]]

dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)

message("Running normalization for batch: ", batch_name)
message("Using scale.factor: ", scale_factor)

rds_files <- list.files(doublet_dir, pattern = "\\.rds$", full.names = TRUE)
message("Found ", length(rds_files), " specimen(s) in this batch")

for (rds_file in rds_files) {

  specimen_id <- tools::file_path_sans_ext(basename(rds_file))
  message("Processing specimen: ", specimen_id)

  obj <- readRDS(rds_file)

# git 

# method “LogNormalize” that normalizes the feature expression measurements 
#for each cell by the total expression, multiplies this by a scale factor 
#(10,000 by default)

  DefaultAssay(obj) <- "RNA" # The normalization is applied to the RNA assay, 
  #which contains the raw counts of gene expression for each cell. This is important because Seurat objects 
  #can contain multiple assays (e.g., RNA, protein, etc.), and we want to ensure that we are normalizing the correct data.
  

  obj <- NormalizeData(  #Take the raw counts and aplly normalization
  
    obj,
    normalization.method = "LogNormalize",
    scale.factor = scale_factor,
    verbose = FALSE
  )

  saveRDS(
    obj,
    file = file.path(norm_dir, paste0(specimen_id, ".rds")),
    compress = "gzip"
  )
}

message("Normalization complete for batch: ", batch_name)

sink(type = "message"); sink()