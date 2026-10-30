# 11_network_metrics.R
# Centrality of every gene of one ARACNe network (one phenotype, cell type and gamma).
# The table is read by the plots and by the ego network analysis, so centrality
# is computed only once per network.
#   degree:      number of neighbors
#   strength:    sum of the MI of its edges
#   pagerank:    PageRank with MI as edge weight
#   betweenness: shortest paths through the gene, with distance = 1 / MI

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(igraph))

source(snakemake@input[["utils"]])  # read_aracne_network(), aracne_graph()

g <- aracne_graph(read_aracne_network(snakemake@input[["network"]]))
message("Network: ", vcount(g), " genes, ", ecount(g), " edges")

pr  <- page_rank(g)$vector
btw <- betweenness(g, weights = 1 / E(g)$weight)

metrics <- data.frame(
  gene             = V(g)$name,
  degree           = degree(g),
  strength         = strength(g),
  pagerank         = pr,
  pagerank_rank    = rank(-pr, ties.method = "min"),
  betweenness      = btw,
  betweenness_rank = rank(-btw, ties.method = "min"),
  row.names = NULL
)
metrics <- metrics[order(metrics$pagerank_rank), ]

write.table(metrics, snakemake@output[["metrics"]], sep = "\t", quote = FALSE, row.names = FALSE)
message("Top 5 PageRank: ", paste(head(metrics$gene, 5), collapse = ", "))
