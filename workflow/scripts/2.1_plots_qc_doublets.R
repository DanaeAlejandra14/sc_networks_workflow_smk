log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(cowplot)
})

## ---- 1. Agregar QC + doublets en una sola tabla ----------------------

qc_files      <- unlist(snakemake@input[["qc"]])
doublet_files <- unlist(snakemake@input[["doublets"]])
summary_path  <- snakemake@output[["summary"]]

message("Aggregating ", length(qc_files), " QC summary file(s) and ",
        length(doublet_files), " doublet summary file(s)")

# library_batch se obtiene del NOMBRE del archivo, no de specimenID,
# porque Snakemake garantiza el patrón "{library_batch}_qc_summary.csv"
# / "{library_batch}_doublet_summary.csv" vía el wildcard. Evita parsear
# specimenID con regex (frágil, p.ej. para "191122-B6-R1710143-alone").
read_with_batch <- function(path, suffix) {
  batch <- basename(path) %>% str_remove(fixed(suffix))
  read_csv(path, show_col_types = FALSE) %>%
    mutate(library_batch = batch, .before = 1)
}

qc_all      <- map_dfr(qc_files, read_with_batch, suffix = "_qc_summary.csv")
doublet_all <- map_dfr(doublet_files, read_with_batch, suffix = "_doublet_summary.csv")

message("QC rows: ", nrow(qc_all), " | Doublet rows: ", nrow(doublet_all))

df <- full_join(qc_all, doublet_all, by = c("library_batch", "specimenID"))

n_unmatched <- sum(is.na(df$n_cells) | is.na(df$n_cells_after))
if (n_unmatched > 0) {
  warning(n_unmatched, " specimen(s) no coincidieron entre QC y doublet reports — revisar specimenID.")
}

write_csv(df, summary_path)
message("Combined QC + doublet summary written to: ", summary_path)
message("Total specimens: ", nrow(df), " across ", n_distinct(df$library_batch), " library batch(es)")

## ---- 2. Plots globales (todas las muestras/batches juntas) -----------

base_theme <- theme_cowplot() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

save_plot_png <- function(plot, path, width = 8, height = 6, dpi = 300) {
  ggsave(filename = path, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
}

bar_height <- max(6, 0.28 * nrow(df))

# 01 - pct_cells_kept por specimen, coloreado por library_batch
p_pct_kept <- df %>%
  ggplot(aes(x = reorder(specimenID, pct_cells_kept), y = pct_cells_kept, fill = library_batch)) +
  geom_col(colour = "black") +
  coord_flip() +
  ggtitle("Porcentaje de células retenidas por specimen (QC)") +
  labs(x = "specimen", y = "pct_cells_kept (%)", fill = "library_batch") +
  base_theme
save_plot_png(p_pct_kept, snakemake@output[["pct_kept"]], height = bar_height)

# 02 - pct_doublet por specimen, coloreado por library_batch
p_pct_doublet <- df %>%
  ggplot(aes(x = reorder(specimenID, pct_doublet), y = pct_doublet, fill = library_batch)) +
  geom_col(colour = "black") +
  coord_flip() +
  ggtitle("Porcentaje de doublets por specimen") +
  labs(x = "specimen", y = "pct_doublet (%)", fill = "library_batch") +
  base_theme
save_plot_png(p_pct_doublet, snakemake@output[["pct_doublet"]], height = bar_height)

# 03 - n_cells_before vs n_cells_after por specimen
p_cells_before_after <- df %>%
  select(specimenID, n_cells_before, n_cells_after) %>%
  pivot_longer(c(n_cells_before, n_cells_after), names_to = "stage", values_to = "cells") %>%
  mutate(stage = recode(stage, n_cells_before = "before", n_cells_after = "after")) %>%
  ggplot(aes(x = reorder(specimenID, cells), y = cells, fill = stage)) +
  geom_col(position = "dodge", colour = "black") +
  coord_flip() +
  ggtitle("Células por specimen: antes vs después del filtro QC") +
  labs(x = "specimen", y = "n células", fill = "etapa") +
  base_theme
save_plot_png(p_cells_before_after, snakemake@output[["cells_before_after"]], height = bar_height + 1)

# 04 - desglose de razones de falla de QC por specimen (stacked bar)
p_fail_reasons <- df %>%
  select(specimenID, n_fail_min_features, n_fail_min_counts, n_fail_max_pct_mito) %>%
  pivot_longer(
    c(n_fail_min_features, n_fail_min_counts, n_fail_max_pct_mito),
    names_to = "criterio", values_to = "n_celulas_fallidas"
  ) %>%
  mutate(criterio = recode(criterio,
    n_fail_min_features = "min_features",
    n_fail_min_counts   = "min_counts",
    n_fail_max_pct_mito = "max_pct_mito"
  )) %>%
  ggplot(aes(x = reorder(specimenID, n_celulas_fallidas), y = n_celulas_fallidas, fill = criterio)) +
  geom_col(colour = "black") +
  coord_flip() +
  ggtitle("Células que fallan cada criterio de QC, por specimen") +
  labs(x = "specimen", y = "n células (puede fallar más de un criterio)", fill = "criterio") +
  base_theme
save_plot_png(p_fail_reasons, snakemake@output[["qc_failure_reasons"]], height = bar_height + 1)

# 05 - scatter: pct_cells_kept vs pct_doublet, coloreado por library_batch
p_scatter <- df %>%
  ggplot(aes(x = pct_cells_kept, y = pct_doublet, colour = library_batch)) +
  geom_point(size = 2.5, alpha = 0.85) +
  ggtitle("Relación entre retención de QC y fracción de doublets") +
  labs(x = "pct_cells_kept (%)", y = "pct_doublet (%)", colour = "library_batch") +
  base_theme
save_plot_png(p_scatter, snakemake@output[["scatter_kept_vs_doublet"]], width = 7, height = 6)

# 06 - pct_doublet por batch (comparación entre batches)
p_box_doublet <- df %>%
  ggplot(aes(x = library_batch, y = pct_doublet)) +
  geom_boxplot(fill = "grey85", colour = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, colour = "cornflowerblue") +
  ggtitle("Distribución de pct_doublet por library_batch") +
  labs(x = "library_batch", y = "pct_doublet (%)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot_png(p_box_doublet, snakemake@output[["doublet_by_batch"]], width = 8, height = 6)

# 07 - pct_cells_kept por batch (comparación entre batches)
p_box_kept <- df %>%
  ggplot(aes(x = library_batch, y = pct_cells_kept)) +
  geom_boxplot(fill = "grey85", colour = "black", outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, colour = "cornflowerblue") +
  ggtitle("Distribución de pct_cells_kept por library_batch") +
  labs(x = "library_batch", y = "pct_cells_kept (%)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot_png(p_box_kept, snakemake@output[["kept_by_batch"]], width = 8, height = 6)

## ---- 3. PDF combinado, 2 plots por página -----------------------------

all_plots <- list(p_pct_kept, p_pct_doublet, p_cells_before_after,
                   p_fail_reasons, p_scatter, p_box_doublet, p_box_kept)

pdf(snakemake@output[["pdf_all"]], width = 8.5, height = 11)
for (i in seq(1, length(all_plots), by = 2)) {
  p1 <- all_plots[[i]]
  p2 <- if (i + 1 <= length(all_plots)) all_plots[[i + 1]] else NULL
  page <- if (!is.null(p2)) {
    cowplot::plot_grid(p1, p2, nrow = 2, ncol = 1, align = "v")
  } else {
    cowplot::plot_grid(p1, nrow = 1, ncol = 1)
  }
  print(page)
}
dev.off()

message("Plots written to: ", dirname(snakemake@output[["pct_kept"]]))
message("Total specimens plotted: ", nrow(df))

sink(type = "message"); sink()


## ---- 2b. Cargar datos por CÉLULA (todos los specimens juntos) --------

qc_cell_files      <- unlist(snakemake@input[["qc_cells"]])
doublet_cell_files <- unlist(snakemake@input[["doublet_cells"]])

message("Reading per-cell QC data from ", length(qc_cell_files), " batch file(s)")
qc_cells_all <- map_dfr(qc_cell_files, read_csv, show_col_types = FALSE)

message("Reading per-cell doublet data from ", length(doublet_cell_files), " batch file(s)")
doublet_cells_all <- map_dfr(doublet_cell_files, read_csv, show_col_types = FALSE)

message("Total cells (QC, pre-filtro): ", nrow(qc_cells_all))
message("Total cells (doublet, post-QC): ", nrow(doublet_cells_all))

# Umbrales de QC (los mismos que usó qc_batch, vía config)
min_feat  <- snakemake@config[["qc"]][["min_features_per_cell"]]
min_count <- snakemake@config[["qc"]][["min_counts_per_cell"]]
mt_max    <- snakemake@config[["qc"]][["max_pct_mito"]]

qc_cells_all <- qc_cells_all %>%
  mutate(
    keep_cf  = factor(keep_cf,  levels = c("keep", "remove")),
    keep_mt  = factor(keep_mt,  levels = c("keep", "remove")),
    keep_all = factor(keep_all, levels = c("keep", "remove"))
  )

# Downsample por specimen SOLO para los scatters (velocidad de render)
downsample_per_specimen <- function(df, cap = 10000L) {
  parts <- split(df, df$specimenID, drop = TRUE)
  bind_rows(lapply(parts, function(dd) {
    m <- min(cap, nrow(dd))
    dd[sample.int(nrow(dd), size = m), , drop = FALSE]
  }))
}
set.seed(123)
qc_cells_ds <- downsample_per_specimen(qc_cells_all, cap = 10000L)

qc_cols <- c(keep = "cornflowerblue", remove = "brown2")
scale_qc_keep <- function(name = "") scale_color_manual(values = qc_cols, drop = FALSE, name = name)

## ---- 08-10: Histogramas globales con línea de corte -------------------

p_hist_mt <- ggplot(qc_cells_all, aes(x = percent_mt + 1)) +
  geom_histogram(binwidth = 0.25, fill = "cornflowerblue", colour = "black") +
  geom_vline(xintercept = mt_max + 1, linetype = 2, color = "red", linewidth = 0.8) +
  scale_x_log10() +
  ggtitle("Distribución global de percent_mt (todas las células)") +
  labs(x = "percent_mt (log10; +1 shift)", y = "células") +
  base_theme
save_plot_png(p_hist_mt, snakemake@output[["hist_mt"]])

p_hist_nf <- ggplot(qc_cells_all, aes(nFeature_RNA)) +
  geom_histogram(bins = 80, fill = "cornflowerblue", colour = "black") +
  geom_vline(xintercept = min_feat, linetype = 2, color = "red", linewidth = 0.8) +
  ggtitle("Distribución global de nFeature_RNA") +
  labs(x = "nFeature_RNA", y = "células") +
  base_theme
save_plot_png(p_hist_nf, snakemake@output[["hist_nfeature"]])

p_hist_nc <- ggplot(qc_cells_all, aes(nCount_RNA)) +
  geom_histogram(bins = 80, fill = "cornflowerblue", colour = "black") +
  geom_vline(xintercept = min_count, linetype = 2, color = "red", linewidth = 0.8) +
  scale_x_log10() +
  ggtitle("Distribución global de nCount_RNA (log10)") +
  labs(x = "nCount_RNA (log10)", y = "células") +
  base_theme
save_plot_png(p_hist_nc, snakemake@output[["hist_ncount"]])

## ---- 11: scatter gradiente mitocondrial (global) -----------------------

p_mito_grad <- qc_cells_ds %>%
  arrange(percent_mt) %>%
  ggplot(aes(nCount_RNA, nFeature_RNA, colour = percent_mt)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
  scale_x_log10() +
  scale_y_log10() +
  ggtitle("Nivel de expresión mitocondrial (todas las células, log10)") +
  labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)", colour = "percent_mt") +
  base_theme
save_plot_png(p_mito_grad, snakemake@output[["scatter_mito_gradient"]], width = 6.8, height = 6.2)

## ---- 12-14: scatters keep/remove (global) -------------------------------

p_keep_cf <- ggplot(qc_cells_ds, aes(nCount_RNA, nFeature_RNA, color = keep_cf)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_vline(xintercept = min_count, linetype = 2, color = "black") +
  geom_hline(yintercept = min_feat, linetype = 2, color = "black") +
  scale_x_log10() + scale_y_log10() +
  scale_qc_keep() +
  ggtitle("Keep vs remove — umbral counts + features (global)") +
  labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  base_theme
save_plot_png(p_keep_cf, snakemake@output[["scatter_keep_cf"]], width = 6.8, height = 6.2)

p_keep_mt <- ggplot(qc_cells_ds, aes(nCount_RNA, nFeature_RNA, color = keep_mt)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_x_log10() + scale_y_log10() +
  scale_qc_keep() +
  ggtitle("Keep vs remove — umbral mitocondrial (global)") +
  labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  base_theme
save_plot_png(p_keep_mt, snakemake@output[["scatter_keep_mt"]], width = 6.8, height = 6.2)

p_keep_all <- ggplot(qc_cells_ds, aes(nCount_RNA, nFeature_RNA, color = keep_all)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_vline(xintercept = min_count, linetype = 2, color = "black") +
  geom_hline(yintercept = min_feat, linetype = 2, color = "black") +
  scale_x_log10() + scale_y_log10() +
  scale_qc_keep() +
  ggtitle("Keep vs remove — combinado (global)") +
  labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  base_theme
save_plot_png(p_keep_all, snakemake@output[["scatter_keep_all"]], width = 6.8, height = 6.2)

## ---- 15-17: violines globales por keep_all ------------------------------

p_vln_mt <- ggplot(qc_cells_all, aes(x = keep_all, y = percent_mt)) +
  geom_violin(fill = "grey85", colour = "black", scale = "width") +
  geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.8) +
  geom_hline(yintercept = mt_max, linetype = 2, color = "red") +
  ggtitle("percent_mt global, keep vs remove") +
  labs(x = "", y = "percent_mt") +
  base_theme
save_plot_png(p_vln_mt, snakemake@output[["violin_mt"]], width = 6, height = 6)

p_vln_nf <- ggplot(qc_cells_all, aes(x = keep_all, y = nFeature_RNA)) +
  geom_violin(fill = "grey85", colour = "black", scale = "width") +
  geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.8) +
  geom_hline(yintercept = min_feat, linetype = 2, color = "red") +
  ggtitle("nFeature_RNA global, keep vs remove") +
  labs(x = "", y = "nFeature_RNA") +
  base_theme
save_plot_png(p_vln_nf, snakemake@output[["violin_nfeature"]], width = 6, height = 6)

p_vln_nc <- ggplot(qc_cells_all, aes(x = keep_all, y = nCount_RNA)) +
  geom_violin(fill = "grey85", colour = "black", scale = "width") +
  geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.8) +
  geom_hline(yintercept = min_count, linetype = 2, color = "red") +
  scale_y_log10() +
  ggtitle("nCount_RNA global, keep vs remove (log10)") +
  labs(x = "", y = "nCount_RNA (log10)") +
  base_theme
save_plot_png(p_vln_nc, snakemake@output[["violin_ncount"]], width = 6, height = 6)

## ---- 18-20: doublets globales -------------------------------------------

dbl_cols <- c(singlet = "cornflowerblue", doublet = "brown2")

doublet_cells_all <- doublet_cells_all %>%
  mutate(
    scDblFinder.class = factor(scDblFinder.class, levels = c("singlet", "doublet")),
    keep_dbl = factor(if_else(scDblFinder.class == "doublet", "remove", "keep"),
                       levels = c("keep", "remove"))
  )
set.seed(123)
doublet_cells_ds <- downsample_per_specimen(doublet_cells_all, cap = 10000L)

p_dbl_keep <- ggplot(doublet_cells_ds, aes(nCount_RNA, nFeature_RNA, color = keep_dbl)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_x_log10() + scale_y_log10() +
  scale_qc_keep() +
  ggtitle("Doublets — keep vs remove (global, scDblFinder)") +
  labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  base_theme
save_plot_png(p_dbl_keep, snakemake@output[["scatter_doublets_keep"]], width = 6.8, height = 6.2)

p_dbl_class <- ggplot(doublet_cells_ds, aes(nCount_RNA, nFeature_RNA, color = scDblFinder.class)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = dbl_cols, drop = FALSE, name = "") +
  ggtitle("Clasificación scDblFinder — singlet vs doublet (global)") +
  labs(x = "nCount_RNA (log10)", y = "nFeature_RNA (log10)") +
  base_theme
save_plot_png(p_dbl_class, snakemake@output[["scatter_doublets_class"]], width = 6.8, height = 6.2)

p_dbl_score <- ggplot(doublet_cells_all, aes(x = scDblFinder.class, y = scDblFinder.score)) +
  geom_violin(fill = "grey85", colour = "black", scale = "width") +
  geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.8) +
  ggtitle("scDblFinder.score por clase (global)") +
  labs(x = "", y = "scDblFinder.score") +
  base_theme
save_plot_png(p_dbl_score, snakemake@output[["violin_doublet_score"]], width = 6, height = 6)