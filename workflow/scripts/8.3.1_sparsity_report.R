
#
# Joins the sparsity curves of every phenotype x cell type (made by 8.3) and writes:
#   - one figure: rows = phenotypes, columns = cell types, one curve per panel,
#     dotted horizontal line at the target sparsity OF THAT CELL TYPE, dotted vertical line
#     at the reference gamma OF THAT CELL TYPE (visual only), dashed linear trend
#   - sparsity_vs_gamma_all.csv        every gamma of every phenotype x cell type
#   - gamma_for_sparsity_levels.csv    gamma needed to reach 0.2, 0.3, 0.4 ... of sparsity
#   - suggested_gamma_per_cell_type.csv  gamma that satisfies ALL phenotypes of a cell type
#   - sparsity_vs_gamma_fits.csv       how linear is each curve (R2, curvature test)
# It only reports: it does not choose or change any gamma of the workflow.

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(ggplot2))

files    <- unlist(snakemake@input[["curves"]])
target_raw <- snakemake@params[["target_sparsity"]]   # list: default + optional cell types (or one number)
reference_raw <- snakemake@params[["reference_gamma"]]   # list: default + optional cell types (or one number); visual only

# Target sparsity of a cell type: its own value if the config has one, otherwise "default"
target_for <- function(ct) {
  v <- if (is.list(target_raw)) {
    if (!is.null(target_raw[[ct]])) target_raw[[ct]] else target_raw[["default"]]
  } else {
    target_raw
  }
  if (is.null(v) || is.na(as.numeric(v))) {
    stop("No target_sparsity for '", ct, "' and no 'default' in the config.")
  }
  as.numeric(v)
}
# Reference gamma of a cell type (same rule: own value, otherwise "default")
reference_for <- function(ct) {
  v <- if (is.list(reference_raw)) {
    if (!is.null(reference_raw[[ct]])) reference_raw[[ct]] else reference_raw[["default"]]
  } else {
    reference_raw
  }
  if (is.null(v) || is.na(as.numeric(v))) {
    stop("No reference_gamma for '", ct, "' and no 'default' in the config.")
  }
  as.numeric(v)
}
levels_cfg <- as.numeric(unlist(snakemake@params[["sparsity_levels"]]))
min_frac <- as.numeric(snakemake@params[["min_frac_expressed"]])

out_png     <- snakemake@output[["png"]]
out_pdf     <- snakemake@output[["pdf"]]
out_all     <- snakemake@output[["curves_csv"]]
out_levels  <- snakemake@output[["levels_csv"]]
out_suggest <- snakemake@output[["suggested_csv"]]
out_fits    <- snakemake@output[["fits_csv"]]
for (f in c(out_png, out_pdf, out_all, out_levels, out_suggest, out_fits)) {
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
}

curves <- do.call(rbind, lapply(files, function(f) read.delim(f, stringsAsFactors = FALSE)))
message("Curves read: ", length(files), " files | rows: ", nrow(curves))

# ---- Gamma at which a curve reaches a sparsity level ---------------------------------
# Sparsity goes down when gamma goes up. First gamma of the scanned range where the
# sparsity is <= level, with linear interpolation inside the segment, rounded UP.
cross_gamma <- function(g, s, level) {
  o <- order(g); g <- g[o]; s <- s[o]
  d <- !is.na(s); g <- g[d]; s <- s[d]
  ok <- which(s <= level)
  if (length(ok) == 0) {
    return(list(interp = NA_real_, up = NA_integer_, status = "not_reached_by_gamma_max"))
  }
  i <- ok[1]
  if (i == 1) {
    return(list(interp = g[1], up = as.integer(ceiling(g[1])),
                status = "already_reached_at_gamma_min"))
  }
  gs <- g[i - 1] + (s[i - 1] - level) / (s[i - 1] - s[i]) * (g[i] - g[i - 1])
  list(interp = gs, up = as.integer(ceiling(gs)), status = "interpolated")
}

combos <- unique(curves[, c("cell_type", "phenotype")])

# ---- Table: gamma needed for each sparsity level -------------------------------------
level_rows <- list()
for (i in seq_len(nrow(combos))) {
  d <- curves[curves$cell_type == combos$cell_type[i] & curves$phenotype == combos$phenotype[i], ]
  # the configured levels plus the target of this cell type
  for (L in sort(unique(c(levels_cfg, target_for(combos$cell_type[i]))))) {
    cg <- cross_gamma(d$gamma, d$sparsity_aracne, L)
    level_rows[[length(level_rows) + 1]] <- data.frame(
      cell_type = combos$cell_type[i], phenotype = combos$phenotype[i],
      sparsity_level = L, gamma_interpolated = round(cg$interp, 2),
      gamma_rounded_up = cg$up, status = cg$status
    )
  }
}
levels_tab <- do.call(rbind, level_rows)

# ---- Table: gamma that satisfies every phenotype of a cell type ----------------------
is_target <- mapply(function(ct, L) abs(L - target_for(ct)) < 1e-9,
                    levels_tab$cell_type, levels_tab$sparsity_level)
tg <- levels_tab[is_target, ]
suggest <- do.call(rbind, lapply(split(tg, tg$cell_type), function(x) {
  data.frame(
    cell_type = x$cell_type[1],
    target_sparsity = target_for(x$cell_type[1]),
    suggested_gamma = if (all(is.na(x$gamma_rounded_up))) NA_integer_ else max(x$gamma_rounded_up, na.rm = TRUE),
    limiting_phenotype = if (all(is.na(x$gamma_rounded_up))) NA_character_ else x$phenotype[which.max(x$gamma_rounded_up)],
    n_phenotypes_not_reached = sum(x$status == "not_reached_by_gamma_max"),
    gamma_by_phenotype = paste0(x$phenotype, "=", x$gamma_rounded_up, collapse = "; ")
  )
}))
rownames(suggest) <- NULL
print(suggest)

# ---- Table: how linear is each curve ---------------------------------------------------
fit_rows <- list()
for (i in seq_len(nrow(combos))) {
  d <- curves[curves$cell_type == combos$cell_type[i] & curves$phenotype == combos$phenotype[i], ]
  d <- d[!is.na(d$sparsity_aracne), ]
  if (nrow(d) < 4) next
  dd <- data.frame(y = d$sparsity_aracne, g = d$gamma)
  lin  <- lm(y ~ g, dd)
  lg   <- lm(y ~ log(g), dd)
  quad <- lm(y ~ g + I(g^2), dd)
  fit_rows[[length(fit_rows) + 1]] <- data.frame(
    cell_type = combos$cell_type[i], phenotype = combos$phenotype[i],
    slope_per_gamma = coef(lin)[["g"]],
    r2_linear = summary(lin)$r.squared,
    r2_log = summary(lg)$r.squared,
    p_quadratic_term = anova(lin, quad)$`Pr(>F)`[2],
    spearman_rho = cor(dd$g, dd$y, method = "spearman")
  )
}
fits <- if (length(fit_rows) > 0) do.call(rbind, fit_rows) else data.frame()

write.csv(curves, out_all, row.names = FALSE)
write.csv(levels_tab, out_levels, row.names = FALSE)
write.csv(suggest, out_suggest, row.names = FALSE)
write.csv(fits, out_fits, row.names = FALSE)

# ---- Figure: rows = phenotypes, columns = cell types ---------------------------------
theme_file <- "workflow/scripts/R/theme_pub.R"
have_theme <- file.exists(theme_file)
if (have_theme) source(theme_file)

# Order of rows / columns and colors come from theme_pub.R
ph_present <- unique(curves$phenotype)
ct_present <- unique(curves$cell_type)
ph_order <- if (have_theme && exists("pal_phenotype")) names(pal_phenotype) else character(0)
ct_order <- if (have_theme && exists("pal_celltype"))  names(pal_celltype)  else character(0)
ph_levels <- c(intersect(ph_order, ph_present), setdiff(ph_present, ph_order))
ct_levels <- c(intersect(ct_order, ct_present), setdiff(ct_present, ct_order))
curves$phenotype <- factor(curves$phenotype, levels = ph_levels)
curves$cell_type <- factor(curves$cell_type, levels = ct_levels)

cols <- if (have_theme && exists("pal_phenotype")) pal_phenotype else character(0)
cols[setdiff(ph_levels, names(cols))] <- "grey40"

# One target per column (cell type)
hline_df <- data.frame(cell_type = factor(ct_levels, levels = ct_levels),
                       yintercept = unname(sapply(ct_levels, target_for)))
vline_df <- data.frame(cell_type = factor(ct_levels, levels = ct_levels),
                       xintercept = unname(sapply(ct_levels, reference_for)))

p <- ggplot(curves, aes(gamma, sparsity_aracne, color = phenotype)) +
  geom_hline(data = hline_df, aes(yintercept = yintercept), inherit.aes = FALSE,
             color = "grey50", linetype = "dotted", linewidth = 0.5) +
  geom_vline(data = vline_df, aes(xintercept = xintercept), inherit.aes = FALSE,
             color = "grey50", linetype = "dotted", linewidth = 0.5) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.3) +
  facet_grid(phenotype ~ cell_type, labeller = labeller(cell_type = label_wrap_gen(10))) +
  scale_color_manual(values = cols[ph_levels], guide = "none") +
  labs(x = "Gamma",
       y = paste0("Fraction of zeros (genes in >= ", min_frac * 100, "% of metacells)"),
       subtitle = paste0("Horizontal dotted: target sparsity of each cell type | vertical dotted: reference gamma of each cell type (visual only)")) +
  (if (have_theme) theme_pub() else theme_bw())

if (have_theme && exists("save_pub")) {
  save_pub(p, out_png, out_pdf, width = 14, height = 7)
} else {
  ggsave(out_png, p, width = 14, height = 7, dpi = 300)
  ggsave(out_pdf, p, width = 14, height = 7)
}

message("Done.")
sink(type = "message"); sink()