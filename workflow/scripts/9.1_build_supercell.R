log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# Install in the env
if (!requireNamespace("SuperCell", quietly = TRUE)) {
  message("SuperCell is not in the env, need be installed from github..")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("GfellerLab/SuperCell", upgrade = "never")
}
suppressPackageStartupMessages(library(SuperCell))

phenotype_wc <- snakemake@wildcards[["phenotype"]]
cell_type_wc <- snakemake@wildcards[["cell_type"]]
gamma_wc     <- as.numeric(snakemake@wildcards[["gamma"]])

integrated_path <- snakemake@input[["integrated"]]
out_path        <- snakemake@output[[1]]

# All af this parameteres are set in the config.yaml file
k_knn          <- snakemake@params[["k_knn"]]
n_pc           <- snakemake@params[["n_pc"]]
do_approx      <- snakemake@params[["do_approx"]]
approx_n       <- snakemake@params[["approx_n"]]
seed           <- snakemake@params[["seed"]]
reduction_pref <- snakemake@params[["reduction"]]

set.seed(seed) # important for the reproducibility of the supercell construction

message("Phenotype: ", phenotype_wc, " | Cell type: ", cell_type_wc, " | Gamma: ", gamma_wc)

obj <- readRDS(integrated_path)
message("Cells: ", ncol(obj))

# (combinations with a single individual do not have a harmony reduction)
reduction_use <- if (reduction_pref %in% Reductions(obj)) {
  reduction_pref
} else {
  message("WARNING: reduction '", reduction_pref, "' not found, using 'pca'")
  "pca"
}

n_pc_use  <- min(n_pc, ncol(Embeddings(obj, reduction = reduction_use))) 
# Number of PCs used to build the metacell graph; it cannot exceed the number
# of PCs available in the Seurat object
embedding <- Embeddings(obj, reduction = reduction_use)[, 1:n_pc_use]
# Create the SuperCell object from the dimensionality-reduction embedding and the gamma parameter

message("Using reduction: ", reduction_use, " | dimensions: ", n_pc_use)
message("Building metacells with gamma = ", gamma_wc, "...")

SC <- SCimplify_from_embedding(
  X         = embedding,
  gamma     = gamma_wc,
  k.knn     = k_knn,
  n.pc      = n_pc_use,
  do.approx = do_approx,
  approx.N  = approx_n,
  seed      = seed
)

message("Metacells built: ", SC$N.SC)
message("Mean size: ", round(mean(SC$supercell_size), 1), " cells/metacell")


dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(SC, out_path)

message("Saved: ", out_path)

# Summary of this run
summary_row <- data.frame(
  phenotype       = phenotype_wc,
  cell_type       = cell_type_wc,
  gamma           = gamma_wc,
  reduction_used  = reduction_use,
  n_cells_input   = ncol(obj),
  n_metacells     = SC$N.SC,
  mean_size       = round(mean(SC$supercell_size), 2),
  min_size        = min(SC$supercell_size),
  max_size        = max(SC$supercell_size)
)

write.csv(summary_row, snakemake@output[["summary"]], row.names = FALSE)
message("Run summary saved in: ", snakemake@output[["summary"]])

sink(type = "message"); sink()