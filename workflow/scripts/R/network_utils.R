# network_utils.R
# Shared functions to read ARACNe networks as igraph objects.
# Load it from any R script with: source(snakemake@input[["utils"]])

suppressPackageStartupMessages(library(igraph))

# Reads an ARACNe network.txt and returns one row per gene pair (undirected).
# ARACNe lists many pairs twice (A->B and B->A) because every gene is used as a
# regulator; MI is symmetric, so both rows are the same edge.
read_aracne_network <- function(path) {
  net <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  # Sort each pair alphabetically so A-B and B-A get the same key
  net$gene_a <- pmin(net$Regulator, net$Target)
  net$gene_b <- pmax(net$Regulator, net$Target)
  net <- net[net$gene_a != net$gene_b, ]  # self-loops (ARACNe does not produce them)
  # Keep the highest MI and lowest p-value of each pair
  mi   <- aggregate(MI ~ gene_a + gene_b, data = net, FUN = max)
  pval <- aggregate(pvalue ~ gene_a + gene_b, data = net, FUN = min)
  merge(mi, pval, by = c("gene_a", "gene_b"))
}

# Undirected igraph object from the output of read_aracne_network(), with MI as edge weight
aracne_graph <- function(edges) {
  g <- graph_from_data_frame(
    edges[, c("gene_a", "gene_b", "MI", "pvalue")], directed = FALSE
  )
  E(g)$weight <- E(g)$MI
  g
}
