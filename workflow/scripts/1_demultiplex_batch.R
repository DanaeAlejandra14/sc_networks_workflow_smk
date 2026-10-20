log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

batch_name <- snakemake@wildcards[["library_batch"]] #input from the Snakefile, this is the name of the batch that we are processing in this script
out_dir <- snakemake@output[[1]] 
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading and demultiplexing batch: ", batch_name) #leer los archvos que yo tenga

#Step 1 : Construct the count matrix and create the Seurat object
counts <- ReadMtx( #Function from Seurat to read the matrix, features and barcodes files and create a count matrix
  mtx = snakemake@input[["matrix"]],
  features = snakemake@input[["features"]],
  cells = snakemake@input[["barcodes"]]
)

# Creation of the Seurat object with the count matrix and assign the batch name
#"This cells came from this batch"
seurat_obj <- CreateSeuratObject(counts = counts, project = batch_name)
seurat_obj$batch <- batch_name

#Step 2: Read the demultiplexing mapping file and filter for the current batch
demultiplex_df <- read_csv(snakemake@input[["mapping"]], show_col_types = FALSE) #leer el mapping_demultiplexing, este viene desde la funcion de el Snakefile
demux_subset <- demultiplex_df %>%
  filter(libraryBatch == batch_name) #The value of the libraryBatch column in this row is exactly equal to the text '190403-B4-A'?
  
  # We just keep the rows that meet that condition, that is, with the batch that I am processing at this moment

#Step 3: Add the demultiplexing information to the Seurat object and split by individualID
# we define "meta", to then join it to the seurat object with addMetadata
meta <- tibble(cell = colnames(seurat_obj)) %>% #crea una fila por celula

  #Cell barcode is the same as ceell -> Join the individualID column from the mapping to each cell of the Seurat object
  #If it does not find it, it leaves those columns empty (NA).
  left_join(demux_subset, by = c("cell" = "cellBarcode")) %>%
  filter(!is.na(individualID)) %>% # This cell its empty 
  column_to_rownames("cell") # Transform the cell column to row names

seurat_obj <- AddMetaData(seurat_obj, metadata = meta) # add de metadata to the Seurat object, now each cell has its individualID
seurat_ind <- SplitObject(seurat_obj, split.by = "individualID") # we divide by indiviadualID, so that each individual has its own Seurat object, and we store them in a list
names(seurat_ind) <- paste0(batch_name, "_", names(seurat_ind))# Now we have Specimen   IDs as names of the list, so that we can save them later with the correct name

message("Processed individuals: ", length(seurat_ind))

#Save each individual Seurat object as a separate .rds file, instead of a single combined file
for (specimen_id in names(seurat_ind)) {
  saveRDS(
    seurat_ind[[specimen_id]],
    file = file.path(out_dir, paste0(specimen_id, ".rds")),
    compress = "gzip"
  )
}

sink(type = "message"); sink()