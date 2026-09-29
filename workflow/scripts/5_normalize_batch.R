log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

# Define inputs and outputs
batch_name <- snakemake@wildcards[["library_batch"]]
doublet_dir <- snakemake@input[["doublet_dir"]]
norm_dir <- snakemake@output[["norm_dir"]]

#El parámetro scale_factor se obtiene del archivo de configuración y se utiliza en la función NormalizeData para escalar los datos normalizados.
scale_factor <- snakemake@config[["normalization"]][["scale_factor"]]

#creamos la carpea de salida si no existe
dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)

message("Running normalization for batch: ", batch_name)
message("Using scale.factor: ", scale_factor)

rds_files <- list.files(doublet_dir, pattern = "\\.rds$", full.names = TRUE)
message("Found ", length(rds_files), " specimen(s) in this batch")

for (rds_file in rds_files) {

  specimen_id <- tools::file_path_sans_ext(basename(rds_file))
  message("Processing specimen: ", specimen_id)

  obj <- readRDS(rds_file)

# de git 

# "m"ethod “LogNormalize” that normalizes the feature expression measurements 
#for each cell by the total expression, multiplies this by a scale factor 
#(10,000 by default)

  DefaultAssay(obj) <- "RNA" # cuando quiera assay saber que es el RNA
  #donde se realizara la normalización

  obj <- NormalizeData(  #toma los count crudos y les aplica la formula de normalizacion
  
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