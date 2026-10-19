log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(Seurat)
  library(SuperCell)
})

integrated_rds <- snakemake@input[["integrated"]]
sc_rds         <- snakemake@input[["sc_rds"]]
outdir         <- snakemake@output[[1]]
individual_col <- snakemake@params[["individual_col"]]
batch_col      <- snakemake@params[["batch_col"]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Leyendo objetos...")
obj <- readRDS(integrated_rds)
sc  <- readRDS(sc_rds)
N   <- sc$N.SC
nw  <- sc$graph.supercells

cells <- names(sc$membership)
stopifnot(all(cells %in% colnames(obj)))
meta <- obj@meta.data[cells, , drop = FALSE]

# Funciones auxiliares 
has_col <- function(col) {
  ok <- !is.null(col) && col %in% colnames(meta)
  if (!ok) message("Columna ausente en la metadata: ", col)
  ok
}

save_png <- function(file, plot_fun, width = 1200, height = 1000) {
  path <- file.path(outdir, file)
  png(path, width = width, height = height, res = 150)
  ok <- tryCatch({ plot_fun(); TRUE },
                 error = function(e) {
                   message("ERROR en ", file, ": ", conditionMessage(e)); FALSE
                 })
  dev.off()
  if (!ok && file.exists(path)) unlink(path)
  invisible(ok)
}

grad_plot <- function(file, values, title) {
  save_png(file, function() {
    supercell_plot_GE(nw, ge = values, main = title, seed = 42)
  })
}

# 00 Red base (nodo = metacelula, tamano segun numero de celulas)
save_png("00_network.png", function() {
  supercell_plot(nw, seed = 42,
                 main = paste0("Metacell network, gamma = ", sc$gamma))
})

# 01 Tamano de metacelulas 
grad_plot("01_size.png", log10(sc$supercell_size), "log10 metacell size (cells)")

# 02-03 Pureza: que tan mezclados estan individuos y batches 
pur_ind <- NULL; pur_batch <- NULL

if (has_col(individual_col)) {
  pur_ind <- as.numeric(supercell_purity(
    clusters = as.character(meta[[individual_col]]),
    supercell_membership = sc$membership, method = "max_proportion"))
  grad_plot("02_purity_individual.png", pur_ind,
            "Purity by individual (max proportion)")
}

if (has_col(batch_col)) {
  pur_batch <- as.numeric(supercell_purity(
    clusters = as.character(meta[[batch_col]]),
    supercell_membership = sc$membership, method = "max_proportion"))
  grad_plot("03_purity_batch.png", pur_batch,
            "Purity by batch (max proportion)")
}

# 04 Histogramas resumen
save_png("04_histograms.png", function() {
  op <- par(mfrow = c(1, 3)); on.exit(par(op))
  hist(sc$supercell_size, breaks = 20, main = "Metacell size", xlab = "cells")
  if (!is.null(pur_ind))
    hist(pur_ind, breaks = seq(0, 1, 0.05),
         main = "Purity by individual", xlab = "max proportion")
  if (!is.null(pur_batch))
    hist(pur_batch, breaks = seq(0, 1, 0.05),
         main = "Purity by batch", xlab = "max proportion")
}, width = 1800, height = 600)

n_png <- length(list.files(outdir, pattern = "\\.png$"))
message("PNG generados: ", n_png)
if (n_png == 0) stop("No se genero ningun PNG; revisa los errores arriba")

sink(type = "message"); sink()