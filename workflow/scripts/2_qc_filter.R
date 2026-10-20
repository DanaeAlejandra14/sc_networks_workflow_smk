log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

batch_name <- snakemake@wildcards[["library_batch"]]
specimens_dir <- snakemake@input[["specimens_dir"]] 
qc_dir <- snakemake@output[["qc_dir"]]
summary_path <- snakemake@output[["summary"]]

min_features <- snakemake@config[["qc"]][["min_features_per_cell"]]
min_counts   <- snakemake@config[["qc"]][["min_counts_per_cell"]]
max_pct_mito <- snakemake@config[["qc"]][["max_pct_mito"]]

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

message("Running QC for batch: ", batch_name)
message("Thresholds -> min_features: ", min_features,
        ", min_counts: ", min_counts,
        ", max_pct_mito: ", max_pct_mito)

rds_files <- list.files(specimens_dir, pattern = "\\.rds$", full.names = TRUE) # regex( cualquier cosa que termine con .rds )
message("Found ", length(rds_files), " specimen(s) in this batch")

qc_summary <- list() # The list to save the QC summary for each specimen_id, that will be saved in a csv at the end
qc_cells   <- list()   # <- NEW list

for (rds_file in rds_files) { # foor loop to process each specimen_id, that is, each Seurat object

  specimen_id <- tools::file_path_sans_ext(basename(rds_file)) # wich specimen_id is this rds file
  message("Processing specimen: ", specimen_id)

  obj <- readRDS(rds_file)  # read the Seurat object for this specimen_id
  n_cells_before <- ncol(obj) # number of cells in this Seurat object before QC filtering

  # The QC filter starts :)

# Calculate teh percentage pf mitocondrial genes, using the Seurat function PercentageFeatureSet, that takes the object 
#and a search pattern for mitochondrial genes (in this case, those that start with "MT-")
  
  # The result is the percentage of counts that are from mitochondrial genes for each cell, and it is stored in the metadata of the Seurat object
  # it is asigned to "percent.mt" in the metadata of the Seurat object
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-") 
  
# report
  n_fail_features <- sum(obj$nFeature_RNA <= min_features)  
  n_fail_counts   <- sum(obj$nCount_RNA   <= min_counts)
  n_fail_mito     <- sum(obj$percent.mt   >= max_pct_mito, na.rm = TRUE)

  
  # Using "subset" of seurat to filter the cells that meet the 3 conditions
  # In this case, we filter the cells that meet all 3 conditions 
  obj_filtered <- subset(
    obj,
    subset = nFeature_RNA > min_features &
      nCount_RNA > min_counts &
      percent.mt < max_pct_mito
  )

  n_cells_after <- ncol(obj_filtered)   #  After the filtering, we count the number of cells that remain in the Seurat object

  saveRDS(
    obj_filtered,
    file = file.path(qc_dir, paste0(specimen_id, ".rds")), 
    compress = "gzip"
  )

# A summary of QC for each specimen_id is made and saved in a tibble, that will be saved in a csv at the end
  qc_summary[[specimen_id]] <- tibble(
    specimenID = specimen_id,
    n_cells_before = n_cells_before,
    n_cells_after = n_cells_after,
    pct_cells_kept = round(100 * n_cells_after / n_cells_before, 2),
    n_fail_min_features = n_fail_features,
    n_fail_min_counts = n_fail_counts,
    n_fail_max_pct_mito = n_fail_mito
  )

  # Metrics FOR EACH CELL are saved in a tibble, that will be saved in a csv at the end
  cell_df <- tibble(
    library_batch = batch_name,
    specimenID    = specimen_id,
    cell_barcode  = colnames(obj),
    nFeature_RNA  = obj$nFeature_RNA,
    nCount_RNA    = obj$nCount_RNA,
    percent_mt    = obj$percent.mt
  ) %>%
    mutate(
      keep_cf  = if_else(nFeature_RNA > min_features & nCount_RNA > min_counts, "keep", "remove"),
      keep_mt  = if_else(percent_mt < max_pct_mito, "keep", "remove"),
      keep_all = if_else(keep_cf == "keep" & keep_mt == "keep", "keep", "remove")
    )
  qc_cells[[specimen_id]] <- cell_df
  

}

qc_cells_df <- bind_rows(qc_cells)
write_csv(qc_cells_df, snakemake@output[["cells"]])
message("Per-cell QC metrics written to: ", snakemake@output[["cells"]])

# A summary of QC for each specimen_id is made and saved in a tibble, that will be saved in a csv at the end 
qc_summary_df <- bind_rows(qc_summary) # takes the list of tibbles and binds them into a single tibble, that will be saved in a csv at the end
write_csv(qc_summary_df, summary_path) 

message("QC complete for batch: ", batch_name)
message("Summary written to: ", summary_path)

sink(type = "message"); sink()