log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

if (!requireNamespace("SuperCell", quietly = TRUE)) {
  message("SuperCell no esta instalado en este ambiente, instalando desde GitHub...")
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

set.seed(seed) # impotant for reproducibility of the supercell construction

message("Phenotype: ", phenotype_wc, " | Cell type: ", cell_type_wc, " | Gamma: ", gamma_wc)

obj <- readRDS(integrated_path)
message("Celulas: ", ncol(obj))

# Usar la reduccion preferida (harmony) si existe; si no, caer a pca
# (las combinaciones con 1 solo individuo no tienen reduccion harmony)
reduction_use <- if (reduction_pref %in% Reductions(obj)) {
  reduction_pref
} else {
  message("ADVERTENCIA: reduccion '", reduction_pref, "' no encontrada, usando 'pca'")
  "pca"
}

n_pc_use  <- min(n_pc, ncol(Embeddings(obj, reduction = reduction_use)))
embedding <- Embeddings(obj, reduction = reduction_use)[, 1:n_pc_use]

message("Usando reduccion: ", reduction_use, " | dimensiones: ", n_pc_use)
message("Construyendo metacelulas con gamma = ", gamma_wc, "...")

SC <- SCimplify_from_embedding(
  X         = embedding,
  gamma     = gamma_wc,
  k.knn     = k_knn,
  n.pc      = n_pc_use,
  do.approx = do_approx,
  approx.N  = approx_n,
  seed      = seed
)

message("Metacelulas construidas: ", SC$N.SC)
message("Tamano promedio: ", round(mean(SC$supercell_size), 1), " celulas/metacelula")

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(SC, out_path)

message("Guardado: ", out_path)

# resumen individual de esta corrida 
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
message("Resumen individual guardado en: ", snakemake@output[["summary"]])

sink(type = "message"); sink()