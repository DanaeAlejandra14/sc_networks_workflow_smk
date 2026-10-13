log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(Seurat)
library(SuperCell)
library(Matrix)

integrated_rds <- snakemake@input[["integrated"]]
sc_rds         <- snakemake@input[["sc_rds"]]
output_rds     <- snakemake@output[[1]]

message("Leyendo objeto integrado: ", integrated_rds)
obj <- readRDS(integrated_rds)

message("Leyendo objeto SuperCell: ", sc_rds)
sc <- readRDS(sc_rds)

# Matriz de expresion cruda (counts) de las celulas originales
GE <- GetAssayData(obj, slot = "counts")

# Asegurar que coincidan las celulas entre el objeto integrado y el membership
common_cells <- intersect(colnames(GE), names(sc$membership))
message("Celulas en comun: ", length(common_cells), " de ", ncol(GE), " en el objeto integrado")

GE <- GE[, common_cells]
membership <- sc$membership[common_cells]

message("Calculando expresion por metacelula...")
SC.GE <- supercell_GE(ge = GE, groups = membership, mode = "average")

message("Matriz de metacelulas: ", nrow(SC.GE), " genes x ", ncol(SC.GE), " metacelulas")

saveRDS(SC.GE, output_rds)
message("Guardado: ", output_rds)

sink(type = "message"); sink()