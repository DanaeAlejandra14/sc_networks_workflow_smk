#scDoble Finder Script 

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse) 
library(Seurat)
library(scDblFinder)

batch_name <- snakemake@wildcards[["library_batch"]]
qc_dir <- snakemake@input[["qc_dir"]]
doublet_dir <- snakemake@output[["doublet_dir"]]
summary_path <- snakemake@output[["summary"]]

seed <- snakemake@config[["doublet_detection"]][["seed"]]
set.seed(seed) #una semilla para la reproducibilidad 

dir.create(doublet_dir, recursive = TRUE, showWarnings = FALSE)

message("Running doublet detection for batch: ", batch_name)
message("Using seed: ", seed)

rds_files <- list.files(qc_dir, pattern = "\\.rds$", full.names = TRUE)
message("Found ", length(rds_files), " specimen(s) in this batch")

doublet_summary <- list() # list to save the summary of doublet detection for each specimen_id, that will be saved in a csv at the end
doublet_cells   <- list()   # list to save the per-cell doublet metrics for each specimen_id, that will be saved in a csv at the end

for (rds_file in rds_files) {

  specimen_id <- tools::file_path_sans_ext(basename(rds_file))
  message("Processing specimen: ", specimen_id)

  obj <- readRDS(rds_file)

# Trandsform the Seurat object to a SingleCellExperiment object, which is the input format required by scDblFinder
  sce <- as.SingleCellExperiment(obj, assay = "RNA") 
  sce <- scDblFinder(sce) # Take the SingleCellExperiment object and run scDblFinder, which simulates doublets and classifies cells as "singlet" or "doublet" based on their expression profile. The results are added to the metadata of the SingleCellExperiment object.
  #This tool estimates the probability of each cell being a doublet based on its expression profile and the expected doublet rate, which is inferred from the data. The results are added to the metadata of the SingleCellExperiment object.
  # This information is save in metadata columns "scDblFinder.score" (the probability of being a doublet) and "scDblFinder.class" (the classification as "singlet" or "doublet").

# Add the doublet detection results back to the Seurat object, 
  obj$scDblFinder.class <- sce$scDblFinder.class #singlet/doublet
  obj$scDblFinder.score <- sce$scDblFinder.score #probability of being a doublet

  n_cells <- ncol(obj) 
  n_doublet <- sum(obj$scDblFinder.class == "doublet")
  n_singlet <- sum(obj$scDblFinder.class == "singlet")

  cell_df <- tibble(
  library_batch      = batch_name,
  specimenID         = specimen_id,
  cell_barcode       = colnames(obj),
  nCount_RNA         = obj$nCount_RNA,
  nFeature_RNA       = obj$nFeature_RNA,
  scDblFinder.score  = obj$scDblFinder.score,
  scDblFinder.class  = obj$scDblFinder.class
)

 doublet_cells[[specimen_id]] <- cell_df   
#Eliminate doublets from the Seurat object, keeping only singlets, and save the filtered object to a new .rds file in the doublet_dir. 
  #The summary of doublet detection (number of cells, number of singlets, number of doublets, percentage of doublets) is saved in a tibble for each specimen_id, which will be combined into a single summary table at the end.

obj_filtered <- subset(obj, subset = scDblFinder.class == "singlet")

saveRDS(
  obj_filtered,
  file = file.path(doublet_dir, paste0(specimen_id, ".rds")),
  compress = "gzip"
)

doublet_summary[[specimen_id]] <- tibble(
  specimenID  = specimen_id,
  n_cells     = n_cells,
  n_singlet   = n_singlet,
  n_doublet   = n_doublet,
  pct_doublet = round(100 * n_doublet / n_cells, 2)
)
}   

doublet_cells_df <- bind_rows(doublet_cells)
write_csv(doublet_cells_df, snakemake@output[["cells"]])
message("Per-cell doublet metrics written to: ", snakemake@output[["cells"]])

# The summary of doublet detection (number of cells, number of singlets, number of doublets, percentage of doublets) 
#is saved in a tibble for each specimen_id, which will be combined into a single summary table at the end.
doublet_summary_df <- bind_rows(doublet_summary)
write_csv(doublet_summary_df, summary_path)

message("Doublet detection complete for batch: ", batch_name)
message("Summary written to: ", summary_path)

sink(type = "message"); sink()