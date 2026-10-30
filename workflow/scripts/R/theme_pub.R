# theme_pub.R
# Shared plotting style for the pipeline: colors, font, ggplot2 themes and a save helper.
# Load it from any R script with: source(snakemake@input[["theme"]])
# (declare it as an input of the rule, so changing a color re-runs the plots)

suppressPackageStartupMessages(library(ggplot2))

# Font and sizes
# "sans" maps to Helvetica/Arial when available, and always exists inside conda envs
pub_font      <- "sans"
pub_base_size <- 11
pub_dpi       <- 300

# Colors
# Phenotypes, ordered by clinical progression (cool = control, warm = disease)
pal_phenotype <- c(
  "control"    = "#3E7CB1",
  "MCI"        = "#E0A526",
  "AD-NC_ASYM" = "#D2691E",
  "AD-NC_SYM"  = "#9E2A2B"
)

# Cell types, grouped by lineage (same colors as Final_plots_cells.R)
pal_celltype <- c(
  "Excitatory Neurons" = "#1B2A5C",
  "Inhibitory Neurons" = "#8E4A8F",
  "Oligodendrocytes"   = "#C4A044",
  "OPCs"               = "#E8D07A",
  "Astrocyte"          = "#2F6577",
  "Microglia"          = "#B5413B",
  "Endothelial"        = "#6B4F40"
)

# Neutral colors for secondary elements (edges, grid, non-highlighted points)
pal_neutral <- c(
  text      = "black",
  secondary = "grey35",
  line      = "grey60",
  light     = "grey85"
)

# Accent color to highlight genes of interest (e.g. SLC1A2) in any plot
pal_highlight <- "#C2185B"

# Order in which phenotypes are shown in legends and panels
phenotype_levels <- names(pal_phenotype)

# Color of one or more phenotypes; unknown names get grey instead of failing
phenotype_color <- function(x) {
  col <- unname(pal_phenotype[x])
  col[is.na(col)] <- pal_neutral[["line"]]
  col
}

# Themes

# Plots with axes (boxplots, bars, scatter plots)
theme_pub <- function(base_size = pub_base_size, base_family = pub_font) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.line        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks       = element_line(colour = "black", linewidth = 0.4),
      axis.text        = element_text(colour = "black"),
      plot.title       = element_text(face = "bold", size = rel(1.8)),
      plot.subtitle    = element_text(colour = pal_neutral[["secondary"]], size = rel(1.5)),
      plot.caption     = element_text(colour = pal_neutral[["secondary"]], size = rel(0.8)),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      legend.key       = element_blank(),
      legend.title     = element_text(face = "bold")
    )
}

# Plots without axes (networks, UMAPs)
theme_pub_void <- function(base_size = pub_base_size, base_family = pub_font) {
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title      = element_text(face = "bold", size = rel(1.8)),
      plot.subtitle   = element_text(colour = pal_neutral[["secondary"]], size = rel(1.5)),
      plot.caption    = element_text(colour = pal_neutral[["secondary"]], size = rel(0.8)),
      plot.margin     = margin(8, 8, 8, 8),
      legend.title    = element_text(face = "bold")
    )
}

# Color scales
scale_color_phenotype <- function(...) {
  scale_color_manual(values = pal_phenotype, breaks = phenotype_levels, ...)
}
scale_fill_phenotype <- function(...) {
  scale_fill_manual(values = pal_phenotype, breaks = phenotype_levels, ...)
}
scale_color_celltype <- function(...) {
  scale_color_manual(values = pal_celltype, ...)
}
scale_fill_celltype <- function(...) {
  scale_fill_manual(values = pal_celltype, ...)
}

# Continuous scale for heatmaps of values from 0 to 1 (NMI, correlations, Jaccard...)
pal_sequential <- c("#FFFFCC", "#A1DAB4", "#41B6C4", "#2C7FB8", "#253494")
scale_fill_pub_seq <- function(limits = c(0, 1), ...) {
  scale_fill_gradientn(colours = pal_sequential, limits = limits, ...)
}

# Saving
# Writes the same plot as PNG and PDF.
# png_path / pdf_path: full output paths (as given by Snakemake)
save_pub <- function(p, png_path, pdf_path, width = 7, height = 6) {
  ggsave(png_path, p, width = width, height = height, dpi = pub_dpi, bg = "white")
  # cairo_pdf embeds fonts correctly; fall back to the default PDF device if missing
  pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"
  ggsave(pdf_path, p, width = width, height = height, device = pdf_device, bg = "white")
  invisible(c(png_path, pdf_path))
}
