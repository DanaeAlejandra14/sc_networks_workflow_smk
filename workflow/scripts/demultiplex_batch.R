log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

batch_name <- snakemake@wildcards[["library_batch"]] #input: batch name, o sea el nombre del directorio que contiene los archivos de la matriz, features y barcodes
out_dir <- snakemake@output[[1]] 
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading and demultiplexing batch: ", batch_name) #leer los archvos que yo tenga

#paso 1: construccion de la matriz de conteos y creacion del objeto seurat
counts <- ReadMtx( #funcion de seurat para leer los archivos de la matriz, features y barcodes
  mtx = snakemake@input[["matrix"]],
  features = snakemake@input[["features"]],
  cells = snakemake@input[["barcodes"]]
)


#se crea el objeto seurat con la matriz de conteos y se le asigna el nombre del batch
#"estas celulas vienen de este batch"
seurat_obj <- CreateSeuratObject(counts = counts, project = batch_name)
seurat_obj$sample <- batch_name

#paso 2: Hacer el subset del mapping
demultiplex_df <- read_csv(snakemake@input[["mapping"]], show_col_types = FALSE) #leer el mapping_demultiplexing, este viene desde la funcion de el Snakefile
demux_subset <- demultiplex_df %>%
  filter(libraryBatch == batch_name) #¿el valor de la columna libraryBatch de esta fila es exactamente igual al texto '190403-B4-A'?".
  #solo se queda con las filas que cumplen esa condicion, es decir, con el batch que estoy procesando en este momento

#paso 3- Demutiplexing

#definimos "meta", luego esta la uniremos a el seurat objet con addMetadta  
meta <- tibble(cell = colnames(seurat_obj)) %>% #crea una fila por celula
  #realiza un join con el subset del mapping
  #cell barcode es igual a cell -> Pega la columna individualID del mapping a cada celula del objeto seurat
  #Si NO lo encuentra, deja esas columnas vacías (NA).
  left_join(demux_subset, by = c("cell" = "cellBarcode")) %>%
  filter(!is.na(individualID)) %>% #esta celda esta vacia?
  column_to_rownames("cell") #se transforma a una columna de nombres de fila, para que se pueda agregar al objeto seurat

seurat_obj <- AddMetaData(seurat_obj, metadata = meta) # unimos la metadata de cada celula al objeto seurat, ahora cada celula tiene su individualID
seurat_ind <- SplitObject(seurat_obj, split.by = "individualID")
names(seurat_ind) <- paste0(batch_name, "_", names(seurat_ind))

message("Processed individuals: ", length(seurat_ind))

# guarda cada specimenID como su propio .rds, en vez de un solo archivo combinado
for (specimen_id in names(seurat_ind)) {
  saveRDS(
    seurat_ind[[specimen_id]],
    file = file.path(out_dir, paste0(specimen_id, ".rds")),
    compress = "gzip"
  )
}

sink(type = "message"); sink()