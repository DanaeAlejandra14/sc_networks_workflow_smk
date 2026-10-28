# Write all stdout and messages/warnings/errors to the Snakemake log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(Matrix)

# Inputs, outputs and parameters provided by the Snakemake rule
ge_rds    <- snakemake@input[["ge"]]
out_mat   <- snakemake@output[["matrix"]]
out_cent  <- snakemake@output[["centroids"]]
normalize <- snakemake@params[["normalize"]]
min_frac  <- snakemake@params[["min_frac_expressed"]]

message("Reading: ", ge_rds)
ge <- as.matrix(readRDS(ge_rds))
message("Original dimensions (genes x metacells): ", paste(dim(ge), collapse = " x "))

# Filter genes by the fraction of metacells with expression > 0 (computed on raw values)
#Es necesario?
keep <- rowMeans(ge > 0) >= min_frac
ge <- ge[keep, , drop = FALSE]
message("Genes after filtering (>= ", min_frac, " of metacells): ", nrow(ge))

# Normalization (if requested) is applied after the gene filter
if (normalize == "log_cpm") {
  ge <- log2(t(t(ge) / colSums(ge)) * 1e6 + 1)
} else if (normalize != "none") {
  stop("normalize must be 'log_cpm' or 'none'")
}

# ARACNe-AP requires sample names
colnames(ge) <- paste0("MC", seq_len(ncol(ge))) 

# Write the expression matrix (first column = gene)  - Matrices 
out <- data.frame(gene = rownames(ge), ge, check.names = FALSE)
dir.create(dirname(out_mat), recursive = TRUE, showWarnings = FALSE)
write.table(out, out_mat, sep = "\t", quote = FALSE, row.names = FALSE)

# list of retained genes - Centroids 
writeLines(rownames(ge), out_cent)

message("Escrito: ", out_mat, " y ", out_cent)
sink(type = "message"); sink()