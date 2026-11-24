# 20_network_jaccard.R  (Snakemake script)
#
# Jaccard similarity between the ARACNe networks built with different gammas for ONE
# phenotype x cell type. Edges are undirected pairs of genes.
#   J(A, B) = |edges in both| / |edges in A or B|
# Networks built with different gammas keep different sets of genes (the gene filter is applied
# after the metacells are built). All edges of each network are compared (optionally only the
# top N by MI, config: jaccard.top_edges).
# Outputs: jaccard matrix (tsv), summary table (tsv) and a heatmap + table figure (png/pdf).

library(ggplot2)
library(gridExtra)

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")
message("SCRIPT VERSION: blues + triangle + table")

files       <- unlist(snakemake@input[["networks"]])
mats        <- unlist(snakemake@input[["matrices"]])
gammas      <- as.character(unlist(snakemake@params[["gammas"]]))
top_edges   <- as.numeric(snakemake@params[["top_edges"]])
out_jaccard <- snakemake@output[["jaccard"]]
out_summary <- snakemake@output[["summary"]]
out_heatmap <- snakemake@output[["heatmap"]]
out_pdf <- snakemake@output[["heatmap_pdf"]]
stopifnot(length(files) == length(gammas), length(mats) == length(gammas))
for (f in c(out_jaccard, out_summary, out_heatmap, out_pdf)) {
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
}

#  Read every network: columns 1-3 = regulator, target, MI -
empty_net <- data.frame(a = character(), b = character(), mi = numeric(), key = character(),
                        stringsAsFactors = FALSE)

nets <- lapply(files, function(f) {
  message("Reading ", f)
  d <- tryCatch(read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || ncol(d) < 3 || nrow(d) == 0) {
    message("  empty network: ", f)
    return(empty_net)
  }
  d <- d[, 1:3]; names(d) <- c("a", "b", "mi")
  d$a <- as.character(d$a); d$b <- as.character(d$b)
  d <- d[order(-d$mi), ]
  d$key <- ifelse(d$a < d$b, paste(d$a, d$b, sep = "|"), paste(d$b, d$a, sep = "|"))
  d[!duplicated(d$key), ]            # an undirected edge counts once (keeps the highest MI)
})
names(nets) <- gammas

genes <- lapply(nets, function(d) unique(c(d$a, d$b)))
shared_genes <- Reduce(intersect, genes)
message("Genes per network: ", paste(gammas, lengths(genes), sep = "=", collapse = "; "),
        " | in all networks: ", length(shared_genes))

# Edge sets (optionally the top N edges by MI) 
sets <- lapply(nets, function(d) {
  if (!is.na(top_edges) && top_edges > 0) d <- head(d[order(-d$mi), ], top_edges)
  d$key
})

# Jaccard matrix (all vs all) 
jac <- function(x, y) {
  u <- length(union(x, y))
  if (u == 0) NA_real_ else length(intersect(x, y)) / u
}
J <- outer(seq_along(sets), seq_along(sets),
           Vectorize(function(i, j) jac(sets[[i]], sets[[j]])))
dimnames(J) <- list(gammas, gammas)
print(round(J, 3))

write.table(data.frame(gamma = gammas, J, check.names = FALSE),
            out_jaccard, sep = "\t", quote = FALSE, row.names = FALSE)

# Summary: sparsity of the ARACNe input matrix and size of each network 
sparsity <- vapply(mats, function(f) {
  x <- read.delim(f, row.names = 1, check.names = FALSE)
  mean(as.matrix(x) == 0)                      # fraction of zeros of the ARACNe input
}, numeric(1))

summ <- data.frame(
  gamma             = gammas,
  sparsity          = round(unname(sparsity), 3),
  n_genes           = lengths(genes),
  n_edges           = lengths(sets),
  genes_in_all_nets = length(shared_genes)
)
write.table(summ, out_summary, sep = "\t", quote = FALSE, row.names = FALSE)
#  Figure: heatmap + summary table
jac_df <- read.delim(out_jaccard, check.names = FALSE)
m <- as.matrix(jac_df[, -1])
rownames(m) <- jac_df[[1]]

df <- as.data.frame(as.table(m))
names(df) <- c("gamma_y", "gamma_x", "jaccard")
df$jaccard <- as.numeric(df$jaccard)

if (all(is.na(df$jaccard))) {
  message("All Jaccard values are NA (empty networks); skipping the figure.")
  file.create(out_heatmap)
  file.create(out_pdf)
} else {
  df$ix <- match(as.character(df$gamma_x), colnames(m))
  df$iy <- match(as.character(df$gamma_y), rownames(m))
  df <- df[df$ix <= df$iy, ]                       # lower triangle, diagonal included
  df$gamma_x <- factor(df$gamma_x, levels = colnames(m))
  df$gamma_y <- factor(df$gamma_y, levels = rev(rownames(m)))  # first gamma on top

  p <- ggplot(df, aes(gamma_x, gamma_y, fill = jaccard)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", jaccard),
                  colour = !is.na(jaccard) & jaccard > 0.6), size = 4) +
    scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
    scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(9, "Blues"),   # light blue (low) -> dark blue (high)
      limits = c(0, 1), na.value = "grey85", name = "Jaccard"
    ) +
    coord_fixed() +
    labs(
      title = paste(snakemake@wildcards[["phenotype"]], "-", snakemake@wildcards[["cell_type"]]),
      subtitle = "Edge overlap between networks (Jaccard similarity)",
      x = "SuperCell gamma", y = "SuperCell gamma"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())

  # Table with borders (same data as networks_summary.tsv)
  tab <- summ
  tab$n_genes <- format(tab$n_genes, big.mark = ",")
  tab$n_edges <- format(tab$n_edges, big.mark = ",")
  names(tab) <- c("gamma", "sparsity", "genes", "edges", "genes in\nall networks")

  tbl_theme <- gridExtra::ttheme_minimal(
    base_size = 10,
    core    = list(bg_params = list(fill = "white", col = "grey30", lwd = 1)),
    colhead = list(bg_params = list(fill = "grey90", col = "grey30", lwd = 1),
                   fg_params = list(fontface = "bold"))
  )
  tbl <- gridExtra::tableGrob(tab, rows = NULL, theme = tbl_theme)

  fig <- gridExtra::arrangeGrob(p, tbl, ncol = 1, heights = c(4, 1.4))

  ggsave(out_heatmap, fig, width = 6.5, height = 7.5, dpi = 300)
  ggsave(out_pdf,     fig, width = 6.5, height = 7.5)
}


message("Saved: ", normalizePath(out_pdf), " (", file.info(out_pdf)$size, " bytes)")
message("Saved: ", normalizePath(out_heatmap), " (", file.info(out_heatmap)$size, " bytes)")
message("Done.")
sink(type = "message"); sink()