# network_utils.R
# Shared functions to read ARACNe networks as igraph objects.
# Load it from any R script with: source(snakemake@input[["utils"]])

suppressPackageStartupMessages(library(igraph))

# Reads an ARACNe network and returns one row per gene pair (undirected).
# Works for the final network.txt (Regulator, Target, MI, pvalue) and for the
# bootstrap networks (Regulator, Target, MI; no p-value).
# ARACNe lists many pairs twice (A->B and B->A) because every gene is used as a
# regulator; MI is symmetric, so both rows are the same edge.
read_aracne_network <- function(path) {
  net <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  # Sort each pair alphabetically so A-B and B-A get the same key
  net$gene_a <- pmin(net$Regulator, net$Target)
  net$gene_b <- pmax(net$Regulator, net$Target)
  net <- net[net$gene_a != net$gene_b, ]  # self-loops (ARACNe does not produce them)
  key <- paste(net$gene_a, net$gene_b, sep = "\t")
  # Lowest p-value of each pair (final network only)
  if ("pvalue" %in% names(net)) net$pvalue <- ave(net$pvalue, key, FUN = min)
  # Highest MI of each pair: sort by MI and keep the first row of each pair
  o <- order(net$MI, decreasing = TRUE)
  net <- net[o, ]
  net <- net[!duplicated(key[o]), ]
  cols <- intersect(c("gene_a", "gene_b", "MI", "pvalue"), names(net))
  net <- net[, cols]
  rownames(net) <- NULL
  net
}

# Undirected igraph object from the output of read_aracne_network(), with MI as edge weight
aracne_graph <- function(edges) {
  g <- graph_from_data_frame(edges, directed = FALSE)
  E(g)$weight <- E(g)$MI
  g
}
