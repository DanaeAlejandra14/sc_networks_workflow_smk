# Networks built with different gammas keep different sets of genes (the gene filter is applied
# after the metacells are built). All edges of each network are compared (optionally only the
# top N by MI, config: jaccard.top_edges).

library(ggplot2)

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

files        <- unlist(snakemake@input[["networks"]])
gammas       <- as.character(unlist(snakemake@params[["gammas"]]))
top_edges    <- as.numeric(snakemake@params[["top_edges"]])
out_jaccard  <- snakemake@output[["jaccard"]]
out_summary  <- snakemake@output[["summary"]]
stopifnot(length(files) == length(gammas))
for (f in c(out_jaccard, out_summary)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

# Read every network: 
nets <- lapply(files, function(f) {
  message("Reading ", f)
  d <- read.delim(f, stringsAsFactors = FALSE)
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

# Restrict edges 
sets <- lapply(nets, function(d) {
  if (!is.na(top_edges) && top_edges > 0) d <- head(d[order(-d$mi), ], top_edges)
  d$key
})

# Jaccard matrix 
# Jaccard similarity = |A ∩ B| / |A ∪ B|, 0 = no overlap, 1 = identical 
jac <- function(x, y) {
  u <- length(union(x, y))
  if (u == 0) NA_real_ else length(intersect(x, y)) / u
}

#The jacard is calculated in every pair of networks and stored in a matrix
J <- outer(seq_along(sets), seq_along(sets),
           Vectorize(function(i, j) jac(sets[[i]], sets[[j]])))
dimnames(J) <- list(gammas, gammas)
print(round(J, 3))

write.table(data.frame(gamma = gammas, J, check.names = FALSE),
            out_jaccard, sep = "\t", quote = FALSE, row.names = FALSE)

#  Size of each network (to interpret the Jaccard)
summ <- data.frame(
  gamma             = gammas,
  n_genes           = lengths(genes),
  n_edges_original  = vapply(nets, nrow, integer(1)),
  n_edges_compared  = lengths(sets),
  genes_in_all_nets = length(shared_genes)
)
write.table(summ, out_summary, sep = "\t", quote = FALSE, row.names = FALSE)

#  Heatmap of the Jaccard matrix 

# Re-read the matrix written above (first column = gamma labels)
jac_df <- read.delim(snakemake@output[["jaccard"]], check.names = FALSE)
m <- as.matrix(jac_df[, -1])
rownames(m) <- jac_df[[1]]

#  Heatmap of the Jaccard matrix
df <- as.data.frame(as.table(m))
names(df) <- c("gamma_y", "gamma_x", "jaccard")
df$gamma_x <- factor(df$gamma_x, levels = colnames(m))
df$gamma_y <- factor(df$gamma_y, levels = rev(rownames(m)))  # first gamma on top

p <- ggplot(df, aes(gamma_x, gamma_y, fill = jaccard)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 4) +
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, 1), name = "Jaccard") +
  coord_fixed() +
  labs(
    title = paste(snakemake@wildcards[["phenotype"]], "-", snakemake@wildcards[["cell_type"]]),
    subtitle = "Edge overlap between networks (Jaccard similarity)",
    x = "SuperCell gamma", y = "SuperCell gamma"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

ggsave(snakemake@output[["heatmap"]], p, width = 5.5, height = 5, dpi = 300)
ggsave(sub("\\.png$", ".pdf", snakemake@output[["heatmap"]]), p, width = 5.5, height = 5)


message("Done.")
sink(type = "message"); sink()