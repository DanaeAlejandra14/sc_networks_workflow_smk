log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

phenotype_wc <- snakemake@wildcards[["phenotype"]]
cell_type_wc <- snakemake@wildcards[["cell_type"]]
input_dirs   <- snakemake@input  # una carpeta por library_batch
out_file     <- snakemake@output[[1]]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

message("Merge para phenotype: ", phenotype_wc, " | cell_type: ", cell_type_wc)

obj_list     <- list() # una lista para almacenar los objetos Seurat de cada specimen
specimen_ids <- c() # un vector para almacenar los specimen_ids correspondientes a cada objeto Seurat

for (d in input_dirs) {
  rds_files <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  for (rds_file in rds_files) {
    specimen_id <- tools::file_path_sans_ext(basename(rds_file))
    message("  Cargando: ", specimen_id, " (de ", d, ")")
    obj_list[[specimen_id]] <- readRDS(rds_file)
    specimen_ids <- c(specimen_ids, specimen_id)
  }
}

message("Total de specimens encontrados: ", length(obj_list))

if (length(obj_list) == 0) { #se detiene si la lista de objetos Seurat esta vacia
  stop("No se encontro ningun specimen para phenotype=", phenotype_wc,
       " cell_type=", cell_type_wc, ". Revisar si esta combinacion existe en los datos.")
}
# La fusion de objetos Seurat se realiza con la funcion merge, que combina los objetos en uno solo, 
# agregando un prefijo a los nombres de las celulas para evitar colisiones de barcodes 
#entre specimens distinto
if (length(obj_list) == 1) { #si solo hay un objeto Seurat, se renombra con el prefijo del specimen_id
  merged_obj <- RenameCells(obj_list[[1]], add.cell.id = specimen_ids[1])
} else {
  merged_obj <- merge( # si hay mas de un objeto Seurat, se fusionan todos en uno solo con Merge de seurat
    x = obj_list[[1]], # extrae el contenido de el primer objeto Seurat
    y = obj_list[-1], # todos los demas objetos Seurat
    add.cell.ids = specimen_ids # los IDs de los specimens  , importante para 
    # evitar colisiones de barcodes entre specimens distintos , 
    #para hacer la integracion por individal ID
  )
}

message("Objeto final: ", ncol(merged_obj), " celulas, ",
        length(unique(merged_obj$individualID)), " individualID(s) distintos")

saveRDS(merged_obj, out_file, compress = "gzip") 

message("Merge completo: ", out_file)
sink(type = "message"); sink()