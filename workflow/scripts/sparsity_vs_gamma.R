# sparsity_vs_gamma.R
#
# Exploratory, READ-ONLY analysis: how does the sparsity (fraction of zeros) of the
# metacell expression matrix change with the SuperCell gamma, and is it linear?
#
# Scope: AD-NC_ASYM x Astrocyte, gammas 3, 5, 8, 10, 15, 20, 30.
#
# It reproduces scripts 9.1 (SCimplify_from_embedding) and 9.3 (supercell_GE, layer
# "data", mode "average") in memory, with the parameters of config/config.yaml, and it
# never writes SuperCell objects or matrices. Only the outputs listed below are written,
# all in results/qc/sparsity_gamma/.
#
# Usage (from the project root, in an env that has Seurat + SuperCell, i.e. supercell.yaml):
#   Rscript workflow/scripts/sparsity_vs_gamma.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SuperCell)
  library(Matrix)
  library(ggplot2)
  library(yaml)
})

# ---- Settings ---------------------------------------------------------------------
phenotype <- "AD-NC_ASYM"
cell_type <- "Astrocyte"
gammas    <- c(3, 5, 8, 10, 15, 20, 30)

cfg        <- read_yaml("config/config.yaml")
sc_cfg     <- cfg$supercell
min_frac   <- cfg$aracne$matrix$min_frac_expressed  # same gene filter as 10_export_aracne_matrices.R
integrated <- file.path("results/phenotype", phenotype, cell_type, "integrated.rds")
out_dir    <- "results/qc/sparsity_gamma"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source("workflow/scripts/R/theme_pub.R")

stopifnot(file.exists(integrated))

# ---- Input ------------------------------------------------------------------------
message("Reading: ", integrated)
obj <- readRDS(integrated)
message("Cells: ", ncol(obj))

# Same reduction logic as 9.1
reduction_use <- if (sc_cfg$reduction %in% Reductions(obj)) sc_cfg$reduction else "pca"
n_pc_use <- min(sc_cfg$n_pc, ncol(Embeddings(obj, reduction = reduction_use)))
embedding <- Embeddings(obj, reduction = reduction_use)[, 1:n_pc_use]
message("Reduction: ", reduction_use, " | PCs: ", n_pc_use)

# Same expression matrix as 9.3 (log-normalized "data")
obj <- JoinLayers(obj)
GE  <- GetAssayData(obj, layer = "data")

# Reference: sparsity of the single-cell matrix itself
expressed_genes <- Matrix::rowSums(GE > 0) > 0
cell_sparsity <- c(
  all_genes      = 1 - Matrix::nnzero(GE) / prod(dim(GE)),
  expressed_only = 1 - Matrix::nnzero(GE[expressed_genes, ]) / prod(dim(GE[expressed_genes, ]))
)
message(sprintf("Single-cell sparsity: all genes = %.4f | expressed genes only = %.4f",
                cell_sparsity[["all_genes"]], cell_sparsity[["expressed_only"]]))

# ---- Sparsity per gamma ------------------------------------------------------------
# Three definitions, because the denominator matters:
#   all_genes      every gene of the object (includes genes never detected)
#   expressed_only genes detected in >= 1 cell (fixed gene set across gammas)
#   aracne_input   genes kept by min_frac_expressed, i.e. the matrix ARACNe-AP receives
zero_frac <- function(m) mean(m == 0)

rows <- lapply(gammas, function(g) {
  message("gamma = ", g)
  set.seed(sc_cfg$seed)  # as in 9.1
  SC <- SCimplify_from_embedding(
    X         = embedding,
    gamma     = g,
    k.knn     = sc_cfg$k_knn,
    n.pc      = n_pc_use,
    do.approx = sc_cfg$do_approx,
    approx.N  = sc_cfg$approx_n,
    seed      = sc_cfg$seed
  )
  membership <- SC$membership[colnames(GE)]
  SC.GE <- as.matrix(supercell_GE(ge = GE, groups = membership, mode = "average"))

  keep <- rowMeans(SC.GE > 0) >= min_frac
  data.frame(
    phenotype      = phenotype,
    cell_type      = cell_type,
    gamma          = g,
    n_metacells    = SC$N.SC,
    mean_size      = mean(SC$supercell_size),
    n_genes        = nrow(SC.GE),
    n_genes_aracne = sum(keep),
    all_genes      = zero_frac(SC.GE),
    expressed_only = zero_frac(SC.GE[expressed_genes[rownames(SC.GE)], , drop = FALSE]),
    aracne_input   = zero_frac(SC.GE[keep, , drop = FALSE])
  )
})
res <- do.call(rbind, rows)

# Sanity check against the existing gamma 10 run (read-only)
sm10 <- file.path("results/supercell", phenotype, cell_type, "summary_gamma10.csv")
if (file.exists(sm10)) {
  ref <- read.csv(sm10)
  message("Check gamma 10: metacells here = ", res$n_metacells[res$gamma == 10],
          " | existing run = ", ref$n_metacells)
}

write.table(res, file.path(out_dir, "sparsity_vs_gamma.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(res, digits = 4)

# ---- Is the relationship linear? ---------------------------------------------------
# For each definition: linear in gamma, linear in log(gamma), and a quadratic term test
# (anova linear vs quadratic; small p => curvature). With 7 points, treat as descriptive.
metrics <- c("all_genes", "expressed_only", "aracne_input")

fit_tab <- do.call(rbind, lapply(metrics, function(m) {
  d <- data.frame(y = res[[m]], g = res$gamma)
  lin  <- lm(y ~ g, d)
  lg   <- lm(y ~ log(g), d)
  quad <- lm(y ~ g + I(g^2), d)
  data.frame(
    metric            = m,
    slope_per_gamma   = coef(lin)[["g"]],
    r2_linear         = summary(lin)$r.squared,
    r2_log            = summary(lg)$r.squared,
    aic_linear        = AIC(lin),
    aic_log           = AIC(lg),
    p_quadratic_term  = anova(lin, quad)$`Pr(>F)`[2],
    spearman_rho      = cor(d$g, d$y, method = "spearman")
  )
}))
write.table(fit_tab, file.path(out_dir, "sparsity_vs_gamma_fits.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(fit_tab, digits = 4)

# ---- Plot --------------------------------------------------------------------------
metric_labels <- c(
  all_genes      = "All genes",
  expressed_only = "Genes detected in >= 1 cell",
  aracne_input   = paste0("ARACNe input (>= ", min_frac * 100, "% metacells)")
)

long <- do.call(rbind, lapply(metrics, function(m) {
  data.frame(gamma = res$gamma, sparsity = res[[m]], metric = unname(metric_labels[m]))
}))
long$metric <- factor(long$metric, levels = unname(metric_labels))

# Same data twice: linear x axis and log x axis. A straight line on the log panel
# means sparsity is linear in log(gamma) rather than in gamma.
long_both <- rbind(
  transform(long, xscale = "gamma"),
  transform(long, xscale = "log(gamma)")
)
long_both$x <- ifelse(long_both$xscale == "gamma", long_both$gamma, log(long_both$gamma))
long_both$xscale <- factor(long_both$xscale, levels = c("gamma", "log(gamma)"))

p <- ggplot(long_both, aes(x, sparsity)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              color = "grey55", linetype = "dashed", linewidth = 0.6) +
  geom_line(color = pal_phenotype[[phenotype]], linewidth = 0.7) +
  geom_point(color = pal_phenotype[[phenotype]], size = 2.4) +
  facet_grid(metric ~ xscale, scales = "free", switch = "x",
             labeller = labeller(metric = label_wrap_gen(16))) +
  labs(x = NULL, y = "Fraction of zeros") +
  theme_pub() +
  theme(strip.placement = "outside",
        strip.text.y = element_text(angle = 0))

save_pub(p,
         file.path(out_dir, "sparsity_vs_gamma.png"),
         file.path(out_dir, "sparsity_vs_gamma.pdf"),
         width = 7, height = 7)

message("Done. Outputs in ", out_dir)
