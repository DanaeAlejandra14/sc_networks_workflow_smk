#!/usr/bin/env Rscript
# choose_n_pc.R
#
# Read-only analysis: how many principal components (n_pc) are worth using to build
# the metacells? It looks at the PCA stored in every subsampled.rds and applies three
# common criteria. It does NOT change the pipeline or any object.
#
#   n_elbow : the "elbow" of the variance curve (point farthest from the straight line
#             that joins the first and the last PC, after scaling both axes to 0-1)
#   n_hbc   : rule of the Harvard Chan Bioinformatics Core tutorial:
#             min( first PC with cumulative variance > 90% and variance < 5%,
#                  last PC whose variance drop to the next one is > 0.1% )
#   n_var80 / n_var90 : first PC where the cumulative variance reaches 80% / 90%
#
# Note: the variance is relative to the PCs that were computed (e.g. 50), as in the
# Seurat tutorials, not to the total variance of all genes.
#
# Usage (from the project root, in the 'seurat' env):
#   Rscript choose_n_pc.R [phenotype_dir] [out_dir]
# Defaults: results/phenotype  results/qc/npc

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args      <- commandArgs(trailingOnly = TRUE)
base_dir  <- if (length(args) >= 1) args[1] else "results/phenotype"
out_dir   <- if (length(args) >= 2) args[2] else "results/qc/npc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


file_name <- if (length(args) >= 3) args[3] else "subsampled.rds"
files <- Sys.glob(file.path(base_dir, "*", "*", file_name))
if (length(files) == 0) stop("No subsampled.rds found under ", base_dir)
message("Objects found: ", length(files))

# ---- Elbow: point of the curve farthest from the line first-last ----------------------
elbow_pc <- function(y) {
  x  <- seq_along(y)
  xn <- (x - min(x)) / (max(x) - min(x))
  yn <- (y - min(y)) / (max(y) - min(y))
  # distance from each point to the line that joins the first and the last point
  d <- abs((yn[length(yn)] - yn[1]) * xn - (xn[length(xn)] - xn[1]) * yn +
           xn[length(xn)] * yn[1] - yn[length(yn)] * xn[1])
  which.max(d)
}

rows   <- list()
curves <- list()

for (f in files) {
  cell_type <- basename(dirname(f))
  phenotype <- basename(dirname(dirname(f)))
  message("== ", phenotype, " / ", cell_type)

  obj <- readRDS(f)
  if (!"pca" %in% Reductions(obj)) {
    message("   no 'pca' reduction, skipped")
    next
  }

  sdev <- Stdev(obj[["pca"]])
  pct  <- sdev^2 / sum(sdev^2) * 100
  cum  <- cumsum(pct)
  n    <- length(pct)

  harmony_dims <- if ("harmony" %in% Reductions(obj)) ncol(Embeddings(obj, "harmony")) else NA_integer_
  n_cells <- ncol(obj)
  rm(obj); invisible(gc())

  # HBC rule
  co1 <- which(cum > 90 & pct < 5)[1]
  co2 <- sort(which((pct[-n] - pct[-1]) > 0.1), decreasing = TRUE)[1] + 1
  vals <- c(co1, co2); vals <- vals[!is.na(vals)]
  n_hbc <- if (length(vals) > 0) min(vals) else NA_integer_

  rows[[length(rows) + 1]] <- data.frame(
    phenotype = phenotype, cell_type = cell_type, n_cells = n_cells,
    pcs_computed = n, harmony_dims = harmony_dims,
    n_elbow = elbow_pc(pct), n_hbc = n_hbc,
    n_var80 = which(cum >= 80)[1], n_var90 = which(cum >= 90)[1],
    cum_var_at_30 = round(cum[min(30, n)], 1)
  )
  curves[[length(curves) + 1]] <- data.frame(
    phenotype = phenotype, cell_type = cell_type, pc = seq_len(n), pct = pct, cum = cum
  )
}

res <- do.call(rbind, rows)
cv  <- do.call(rbind, curves)
write.csv(res, file.path(out_dir, "n_pc_summary.csv"), row.names = FALSE)
print(res)

# ---- Summary across all objects ----------------------------------------------------------
message("\nMedian / max over all objects:")
for (m in c("n_elbow", "n_hbc", "n_var80", "n_var90")) {
  message(sprintf("  %-8s median = %s | max = %s", m,
                  median(res[[m]], na.rm = TRUE), max(res[[m]], na.rm = TRUE)))
}
message("  smallest number of Harmony dimensions in any object: ",
        min(res$harmony_dims, na.rm = TRUE))

# ---- Figure: variance per PC, same layout as the sparsity figure ---------------------------
theme_file <- "workflow/scripts/R/theme_pub.R"
have_theme <- file.exists(theme_file)
if (have_theme) source(theme_file)

ph_order <- if (have_theme && exists("pal_phenotype")) names(pal_phenotype) else character(0)
ct_order <- if (have_theme && exists("pal_celltype"))  names(pal_celltype)  else character(0)
ph_lev <- c(intersect(ph_order, unique(cv$phenotype)), setdiff(unique(cv$phenotype), ph_order))
ct_lev <- c(intersect(ct_order, unique(cv$cell_type)),  setdiff(unique(cv$cell_type),  ct_order))
cv$phenotype  <- factor(cv$phenotype,  levels = ph_lev)
cv$cell_type  <- factor(cv$cell_type,  levels = ct_lev)
res$phenotype <- factor(res$phenotype, levels = ph_lev)
res$cell_type <- factor(res$cell_type, levels = ct_lev)

cols <- if (have_theme && exists("pal_phenotype")) pal_phenotype else character(0)
cols[setdiff(ph_lev, names(cols))] <- "grey40"

p <- ggplot(cv, aes(pc, pct, color = phenotype)) +
  geom_vline(data = res, aes(xintercept = n_elbow), color = "grey50", linetype = "dotted") +
  geom_vline(xintercept = 30, color = "grey80", linetype = "dashed") +   # current n_pc
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.9) +
  facet_grid(phenotype ~ cell_type, labeller = labeller(cell_type = label_wrap_gen(10))) +
  scale_color_manual(values = cols[ph_lev], guide = "none") +
  labs(x = "Principal component", y = "% of variance (of the computed PCs)",
       subtitle = "Dotted: elbow of each curve | light dashed: n_pc = 30") +
  (if (have_theme) theme_pub() else theme_bw())



# ---- Second figure: cumulative variance ----
# Cumulative % of the variance of the computed PCs. It always ends at 100%, so use it as
# support for the elbow, not as the only criterion.
p_cum <- ggplot(cv, aes(pc, cum, color = phenotype)) +
  geom_hline(yintercept = c(80, 90), color = "grey70", linetype = "dotted") +
  geom_vline(xintercept = 30, color = "grey80", linetype = "dashed") +   # current n_pc
  geom_line(linewidth = 0.5) +
  geom_point(size = 0.9) +
  facet_grid(phenotype ~ cell_type, labeller = labeller(cell_type = label_wrap_gen(10))) +
  scale_color_manual(values = cols[ph_lev], guide = "none") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Principal component", y = "Cumulative % of variance (of the computed PCs)",
       subtitle = "Dotted: 80% and 90% | light dashed: n_pc = 30") +
  (if (have_theme) theme_pub() else theme_bw())

ggsave(file.path(out_dir, "n_pc_cumulative.png"), p_cum, width = 14, height = 7, dpi = 300)

ggsave(file.path(out_dir, "n_pc_elbow.png"), p, width = 14, height = 7, dpi = 300)
message("Done. Outputs in ", out_dir)