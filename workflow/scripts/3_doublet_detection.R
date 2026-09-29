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

doublet_summary <- list() # empieza la lista para almacenar los resultados de cada specimen
doublet_cells   <- list()   # <- nuevo

for (rds_file in rds_files) {

  specimen_id <- tools::file_path_sans_ext(basename(rds_file))
  message("Processing specimen: ", specimen_id)

  obj <- readRDS(rds_file)

#Transforma el objeto Seurat a un objeto SingleCellExperiment para usar scDblFinder
  sce <- as.SingleCellExperiment(obj, assay = "RNA") 
  sce <- scDblFinder(sce) #Toma el objeto SingleCellExperiment y ejecuta la detección de dobletes, 
  #simulando dobletes y clasificando las células como "singlet" o "doublet" según su perfil de expresión.
  #se agrega la columna dentro de la misma Metadata interna 

# Agrega los resultados de scDblFinder al objeto Seurat original
  obj$scDblFinder.score <- sce$scDblFinder.score #probabilidad de que una célula sea un doblete
  obj$scDblFinder.class <- sce$scDblFinder.class #singlet/doublet

  n_cells <- ncol(obj) 
  n_doublet <- sum(obj$scDblFinder.class == "doublet")
  n_singlet <- sum(obj$scDblFinder.class == "singlet")

  # Ahora sí se eliminan los doublets del objeto que se guarda
  obj_filtered <- subset(obj, subset = scDblFinder.class == "singlet")

  saveRDS(
    obj_filtered,                                      # <- antes era "obj"
    file = file.path(doublet_dir, paste0(specimen_id, ".rds")),
    compress = "gzip"
  )

  doublet_summary[[specimen_id]] <- tibble( 
    specimenID = specimen_id,
    n_cells = n_cells,
    n_singlet = n_singlet,
    n_doublet = n_doublet,
    pct_doublet = round(100 * n_doublet / n_cells, 2)
  )

  
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

n_cells   <- ncol(obj)
n_doublet <- sum(obj$scDblFinder.class == "doublet")
n_singlet <- sum(obj$scDblFinder.class == "singlet")

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
}   # <- aquí SÍ cierra el for loop

doublet_cells_df <- bind_rows(doublet_cells)
write_csv(doublet_cells_df, snakemake@output[["cells"]])
message("Per-cell doublet metrics written to: ", snakemake@output[["cells"]])

# Los dobletes se eliminan en este mismo paso; el .rds guardado en doublet_dir
# solo contiene singlets. El resumen (n_doublet, pct_doublet) se calcula
# ANTES de esa eliminación, sobre todas las células.
doublet_summary_df <- bind_rows(doublet_summary)
write_csv(doublet_summary_df, summary_path)

message("Doublet detection complete for batch: ", batch_name)
message("Summary written to: ", summary_path)

sink(type = "message"); sink()