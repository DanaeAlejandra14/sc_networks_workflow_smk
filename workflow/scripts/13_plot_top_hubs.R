# 13_plot_top_hubs.R
# Top PageRank hubs of the ARACNe networks of one cell type and gamma:
# one panel per phenotype (clinical order, theme_pub.R) with its top_n hubs
# as bars, from most to least central. Genes of interest are highlighted.
# PageRank comes from node_metrics.tsv (11_network_metrics.R).

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source(snakemake@input[["theme"]])

# Inputs, outputs and parameters from the Snakemake rule
metrics_files <- unlist(snakemake@input[["metrics"]])
phenotypes    <- unlist(snakemake@params[["phenotypes"]])  # same order as metrics_files
top_n         <- snakemake@params[["top_n"]]
highlight     <- unlist(snakemake@params[["highlight_genes"]])

out_png   <- snakemake@output[["png"]]
out_pdf   <- snakemake@output[["pdf"]]

stopifnot(length(metrics_files) == length(phenotypes))

# Top hubs of every network
top_hubs <- function(path, phenotype) {
  m <- read.delim(path, stringsAsFactors = FALSE)
  m <- head(m[order(m$pagerank, decreasing = TRUE), ], top_n)
  message(phenotype, ": ", paste(m$gene, collapse = ", "))
  data.frame(phenotype = phenotype, rank = seq_len(nrow(m)), gene = m$gene, pagerank = m$pagerank)
}
hubs <- do.call(rbind, Map(top_hubs, metrics_files, phenotypes))
rownames(hubs) <- NULL

# Plot
# Panels in clinical order; phenotypes without a defined color go last
ord <- phenotype_levels[phenotype_levels %in% phenotypes]
ord <- c(ord, setdiff(phenotypes, ord))

# Same x axis in every panel so bar lengths can be compared
x_max <- max(hubs$pagerank) * 1.05

plot_hubs <- function(phenotype) {
  d <- hubs[hubs$phenotype == phenotype, ]
  d$gene <- factor(d$gene, levels = rev(d$gene))  # most central on top
  d$is_hl <- as.character(d$gene) %in% highlight
  col <- phenotype_color(phenotype)
  # Highlighted genes: bar and label in the highlight color, label in bold
  label_face <- ifelse(rev(d$is_hl), "bold", "plain")
  label_col  <- ifelse(rev(d$is_hl), pal_highlight, "black")

  ggplot(d, aes(x = pagerank, y = gene)) +
    geom_col(aes(fill = is_hl), width = 0.75) +
    scale_fill_manual(values = c(`FALSE` = col, `TRUE` = pal_highlight), guide = "none") +
    scale_x_continuous(limits = c(0, x_max), labels = scales::label_scientific(digits = 2),
                       expand = expansion(mult = c(0, 0))) +
    labs(title = phenotype, x = "PageRank", y = NULL) +
    theme_pub() +
    theme(
      plot.title  = element_text(colour = col),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      axis.text.x  = element_text(size = rel(1.25)),
      axis.title.x = element_text(size = rel(1.4)),
      # Vectorized label style (one value per gene, bottom to top)
      axis.text.y  = element_text(size = rel(1.5), face = label_face, colour = label_col)
    )
}

# element_text() with one value per label works but ggplot2 warns about it
panels <- suppressWarnings(lapply(ord, plot_hubs))
fig <- wrap_plots(panels, ncol = 2)

n_rows <- ceiling(length(ord) / 2)
suppressWarnings(save_pub(fig, out_png, out_pdf, width = 11, height = 4 * n_rows))

message("\nDone.")
