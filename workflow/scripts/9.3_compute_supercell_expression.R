log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(Seurat)
library(SuperCell)
library(Matrix)

integrated_rds <- snakemake@input[["integrated"]] # integrated Seurat object (expression)
sc_rds         <- snakemake@input[["sc_rds"]]     # SuperCell object for this gamma
output_rds     <- snakemake@output[[1]]

message("Reading integrated object: ", integrated_rds)
obj <- readRDS(integrated_rds)

message("Reading SuperCell object: ", sc_rds)
sc <- readRDS(sc_rds)

# Raw expression matrix (counts) of the original cells
obj <- JoinLayers(obj)
GE <- GetAssayData(obj, layer = "counts")

# Make sure the cells match between the integrated object and the membership vector
common_cells <- intersect(colnames(GE), names(sc$membership))
message("Cells in common: ", length(common_cells), " of ", ncol(GE), " in the integrated object")

GE <- GE[, common_cells]
membership <- sc$membership[common_cells]


message("Computing expression per metacell...")
# supercell_GE() takes the gene-by-cell matrix and, using the membership vector
# (which says which cells belong to each metacell), computes the average
# expression of each gene within each group.

SC.GE <- supercell_GE(ge = GE, groups = membership, mode = "average")

message("Metacell matrix: ", nrow(SC.GE), " genes x ", ncol(SC.GE), " metacells")

saveRDS(SC.GE, output_rds)
message("Saved: ", output_rds)

sink(type = "message"); sink()