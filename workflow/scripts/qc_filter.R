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

qc_summary <- list() # crea la lista vacia para guardar el resumen de QC de cada specimen_id

for (rds_file in rds_files) { # ciclo "for" en cada vuelta de el loop 

  specimen_id <- tools::file_path_sans_ext(basename(rds_file))  #solo es para decirme que specimen se esta procesando, quitando la extension del archivo y el path
  message("Processing specimen: ", specimen_id)

  obj <- readRDS(rds_file)  #lee el archivo y lo guarda en obj 
  n_cells_before <- ncol(obj) # número de celulas de el specimen Id antes de los filtros de QC

  #Inician los filtros de QC

#Calculamos el porcentane de genes mitocondriales, usando la funcion PercentageFeatureSet de Seurat, que toma el objeto y un patron de busqueda para los genes mitocondriales (en este caso, los que empiezan con "MT-")
  #Busca el patrón en los rownames el patron de MT
  #resultado: porcentaje de cada celula de cada uno de los specimen_Id 
  #se asigna a "percent.mt" como metadata de mi objeto Seurat
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-") 
  
#solo es un reporte diagnositco , solo es para saber cuantas celulas no cumplen con  cada uno los filtros de QC
#no estoy segura de que sea necesario 
  n_fail_features <- sum(obj$nFeature_RNA <= min_features)  
  n_fail_counts   <- sum(obj$nCount_RNA   <= min_counts)
  n_fail_mito     <- sum(obj$percent.mt   >= max_pct_mito, na.rm = TRUE)

  
  #El filtro real , usando la funcion subset de Seurat, que toma el objeto y un conjunto de condiciones para filtrar las celulas
  #En este caso, se filtran las celulas que cumplen con las 3 condiciones 
  obj_filtered <- subset(
    obj,
    subset = nFeature_RNA > min_features &
      nCount_RNA > min_counts &
      percent.mt < max_pct_mito
  )

  n_cells_after <- ncol(obj_filtered)   # número de celulas de el specimen Id déspues de los filtros de QC

  saveRDS(
    obj_filtered,
    file = file.path(qc_dir, paste0(specimen_id, ".rds")), # se guarda en el path que define qc_dir y con el nombre de specimen_id.rds
    compress = "gzip"
  )

#Hacer un resumen de QC para cada specimen_id y guardarlo en un tibble, que luego se va a guardar en un csv
  qc_summary[[specimen_id]] <- tibble(
    specimenID = specimen_id,
    n_cells_before = n_cells_before,
    n_cells_after = n_cells_after,
    pct_cells_kept = round(100 * n_cells_after / n_cells_before, 2),
    n_fail_min_features = n_fail_features,
    n_fail_min_counts = n_fail_counts,
    n_fail_max_pct_mito = n_fail_mito
  )
}

# Déspues de el lopp por cada specimen_id 
qc_summary_df <- bind_rows(qc_summary) #toma cada una de las filas por specimenID y las junta en un solo tibble
write_csv(qc_summary_df, summary_path) # guarda esa tabla en un csb y en el path que se define en el Snakefile 

message("QC complete for batch: ", batch_name)
message("Summary written to: ", summary_path)

sink(type = "message"); sink()