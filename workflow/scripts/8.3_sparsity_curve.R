# 8.3_sparsity_curve.R  (Snakemake script)
#
# For ONE phenotype x cell type: how does the sparsity (fraction of zeros) of the matrix
# that ARACNe receives change with the SuperCell gamma?
#
# It reproduces 9.1 (SCimplify_from_embedding) and 9.3 (supercell_GE, layer "data",
# mode "average") IN MEMORY with the settings of config/config.yaml. It never saves
# SuperCell objects or matrices; it only writes a small table (one row per gamma).
#
# "ARACNe input" = metacells x genes detected in >= min_frac_expressed of the metacells
# (the same gene filter as 10_export_aracne_matrices.R).

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(Seurat)
  library(SuperCell)
  library(Matrix)
})

# Settings from config/config.yaml
phenotype <- snakemake@wildcards[["phenotype"]]
cell_type <- snakemake@wildcards[["cell_type"]]
rds       <- snakemake@input[["subsampled"]]
out_curve <- snakemake@output[["curve"]]

k_knn     <- snakemake@params[["k_knn"]]
n_pc      <- snakemake@params[["n_pc"]]
do_approx <- snakemake@params[["do_approx"]]
approx_n  <- snakemake@params[["approx_n"]]
seed      <- snakemake@params[["seed"]]
reduction <- snakemake@params[["reduction"]]

min_frac   <- as.numeric(snakemake@params[["min_frac_expressed"]])
gamma_min  <- as.numeric(snakemake@params[["gamma_min"]])
gamma_max  <- as.numeric(snakemake@params[["gamma_max"]])
gamma_step <- as.numeric(snakemake@params[["gamma_step"]])

gammas <- sort(unique(c(seq(gamma_min, gamma_max, by = gamma_step), gamma_max)))
message("Phenotype: ", phenotype, " | cell type: ", cell_type)
message("Gammas: ", paste(gammas, collapse = ", "), " | min_frac_expressed: ", min_frac)

# Input 
obj <- readRDS(rds)

# Same reduction logic as 9.1
reduction_use <- if (reduction %in% Reductions(obj)) reduction else "pca" # use Harmony as reduction if available, else PCA
emb_all   <- Embeddings(obj, reduction = reduction_use)
n_pc_use  <- min(n_pc, ncol(emb_all))
embedding <- emb_all[, 1:n_pc_use, drop = FALSE] 

# Same expression matrix as 9.3 (log-normalized "data")
obj <- JoinLayers(obj)
GE  <- GetAssayData(obj, layer = "data")
rm(obj); invisible(gc())

n_cells <- ncol(GE)
message("Cells: ", n_cells, " | genes: ", nrow(GE),
        " | reduction: ", reduction_use, " (", n_pc_use, " PCs)")

# Genes never detected are all-zero rows at every gamma and can never pass the filter,
# so they are dropped now to save memory.
GE_e <- GE
rm(GE); invisible(gc())

#  One row per gamma 
rows <- lapply(gammas, function(g) {
  message("gamma = ", g)
  set.seed(seed)  # as in 9.1
  SC <- SCimplify_from_embedding(
    X         = embedding,
    gamma     = g,
    k.knn     = k_knn,
    n.pc      = n_pc_use,
    do.approx = do_approx,
    approx.N  = approx_n,
    seed      = seed
  )
  membership <- SC$membership[colnames(GE_e)]
  stopifnot(!anyNA(membership))
  SC.GE <- as.matrix(supercell_GE(ge = GE_e, groups = membership, mode = "average"))

  keep <- rowMeans(SC.GE > 0) >= min_frac
  data.frame(
    cell_type       = cell_type,
    phenotype       = phenotype,
    gamma           = g,
    n_cells         = n_cells,
    n_metacells     = SC$N.SC,
    mean_size       = mean(SC$supercell_size),
    n_genes_aracne  = sum(keep),
    sparsity_aracne = if (any(keep)) mean(SC.GE[keep, , drop = FALSE] == 0) else NA_real_
  )
})
curve <- do.call(rbind, rows)
print(curve, digits = 4)

dir.create(dirname(out_curve), recursive = TRUE, showWarnings = FALSE)
write.table(curve, out_curve, sep = "\t", quote = FALSE, row.names = FALSE)
message("Done: ", out_curve)

sink(type = "message"); sink()