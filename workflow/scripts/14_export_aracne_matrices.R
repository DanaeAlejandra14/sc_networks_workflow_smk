log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(Matrix)

ge_rds    <- snakemake@input[["ge"]]
out_mat   <- snakemake@output[["matrix"]]
out_cent  <- snakemake@output[["centroids"]]
normalize <- snakemake@params[["normalize"]]
min_frac  <- snakemake@params[["min_frac_expressed"]]

message("Leyendo: ", ge_rds)
ge <- as.matrix(readRDS(ge_rds))
message("Dimensiones originales (genes x metacelulas): ", paste(dim(ge), collapse = " x "))

# Filtro por fraccion de metacelulas con expresion > 0 (sobre valores crudos)
keep <- rowMeans(ge > 0) >= min_frac
ge <- ge[keep, , drop = FALSE]
message("Genes tras filtro (>= ", min_frac, " de metacelulas): ", nrow(ge))

# Normalizacion - Si este paso se hace es despues de el filtro de genes
if (normalize == "log_cpm") {
  ge <- log2(t(t(ge) / colSums(ge)) * 1e6 + 1)
} else if (normalize != "none") {
  stop("normalize debe ser 'log_cpm' o 'none'")
}

# ARACNe-AP necesita nombres de muestra
colnames(ge) <- paste0("MC", seq_len(ncol(ge)))

out <- data.frame(gene = rownames(ge), ge, check.names = FALSE)
dir.create(dirname(out_mat), recursive = TRUE, showWarnings = FALSE)
write.table(out, out_mat, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(rownames(ge), out_cent)

message("Escrito: ", out_mat, " y ", out_cent)
sink(type = "message"); sink()