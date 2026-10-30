# 14_ego_network.R
# Ego network of one gene in one ARACNe network (one phenotype, cell type and gamma).
# Writes tables only; the figures comparing phenotypes are made by 16_plot_ego_network.R
#   neighbors.tsv         first neighbors, with MI and bootstrap support
#   edges.tsv             edges of the ego network with first neighbors (for the figure)
#   members.tsv           genes up to analysis_order and their distance to the gene
#                         (input of the enrichment, 15_ego_enrichment.R)
#   metrics.tsv           one row: position of the gene in the network
#   bootstrap_support.tsv every partner of the gene in the bootstrap networks
# If the gene is not in the final network, the tables are written with headers only
# (metrics.tsv says in_network = FALSE), and bootstrap_support.tsv still shows its
# partners in the bootstraps, which explains why it is missing.

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(igraph))

source(snakemake@input[["utils"]])  # read_aracne_network(), aracne_graph()

# Inputs, outputs and parameters from the Snakemake rule
network_file   <- snakemake@input[["network"]]
metrics_file   <- snakemake@input[["metrics"]]
boot_files     <- unlist(snakemake@input[["bootstraps"]])
ego_gene       <- snakemake@wildcards[["gene"]]
phenotype      <- snakemake@wildcards[["phenotype"]]
analysis_order <- snakemake@params[["analysis_order"]]
out            <- snakemake@output

write_tsv <- function(df, path) write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)

g <- aracne_graph(read_aracne_network(network_file))
node_metrics <- read.delim(metrics_file, stringsAsFactors = FALSE)
in_net <- ego_gene %in% V(g)$name
message(phenotype, ": ", ego_gene, if (in_net) " is" else " is NOT", " in the network")

# Bootstrap support
# For each bootstrap, the set of partners of the gene (A-B and B-A count once).
# Support = fraction of bootstraps in which the pair appears
bootstrap_partners <- function(path) {
  lines <- readLines(path)
  hits <- lines[grepl(ego_gene, lines, fixed = TRUE)]  # fast pre-filter before splitting
  if (length(hits) == 0) return(character(0))
  f <- do.call(rbind, strsplit(hits, "\t", fixed = TRUE))
  partner <- ifelse(f[, 1] == ego_gene, f[, 2], ifelse(f[, 2] == ego_gene, f[, 1], NA))
  unique(na.omit(partner))
}
n_boot <- length(boot_files)
partners <- unlist(lapply(boot_files, bootstrap_partners))
final_neighbors <- if (in_net) neighbors(g, ego_gene)$name else character(0)

counts <- sort(table(partners), decreasing = TRUE)
support <- data.frame(
  phenotype        = rep(phenotype, length(counts)),
  partner          = names(counts),
  n_bootstraps     = as.integer(counts),
  n_total          = rep(n_boot, length(counts)),
  support          = as.integer(counts) / n_boot,
  in_final_network = names(counts) %in% final_neighbors
)
write_tsv(support, out[["support"]])
message("Bootstraps: ", n_boot, " | distinct partners: ", nrow(support),
        " | max support: ", if (nrow(support) > 0) max(support$n_bootstraps) else 0, "/", n_boot)

# First neighbors and ego network edges
if (in_net) {
  mi <- sapply(final_neighbors, function(n) E(g, P = c(ego_gene, n))$MI)
  neighbors_tab <- data.frame(phenotype = phenotype, neighbor = final_neighbors, MI = unname(mi))
  neighbors_tab <- neighbors_tab[order(-neighbors_tab$MI), ]
  neighbors_tab$n_bootstraps <- support$n_bootstraps[match(neighbors_tab$neighbor, support$partner)]
  neighbors_tab$n_total <- n_boot

  ego1 <- make_ego_graph(g, order = 1, nodes = ego_gene)[[1]]
  edges_tab <- igraph::as_data_frame(ego1, what = "edges")[, c("from", "to", "MI")]
} else {
  neighbors_tab <- data.frame(phenotype = character(0), neighbor = character(0), MI = numeric(0),
                              n_bootstraps = integer(0), n_total = integer(0))
  edges_tab <- data.frame(from = character(0), to = character(0), MI = numeric(0))
}
write_tsv(neighbors_tab, out[["neighbors"]])
write_tsv(edges_tab, out[["edges"]])

# Genes up to analysis_order
if (in_net) {
  ego_nodes <- ego(g, order = analysis_order, nodes = ego_gene)[[1]]
  d <- distances(g, v = ego_gene, to = ego_nodes)
  members <- data.frame(phenotype = phenotype, gene = colnames(d), order = as.integer(d[1, ]))
} else {
  members <- data.frame(phenotype = character(0), gene = character(0), order = integer(0))
}
write_tsv(members, out[["members"]])

# Metrics: position of the gene in the network
m <- node_metrics[node_metrics$gene == ego_gene, ]
k <- if (in_net) m$degree else NA
links_between_neighbors <- if (in_net) ecount(make_ego_graph(g, order = 1, nodes = ego_gene)[[1]]) - k else NA
metrics <- data.frame(
  phenotype               = phenotype,
  gene                    = ego_gene,
  in_network              = in_net,
  network_genes           = vcount(g),
  degree                  = k,
  pagerank                = if (in_net) m$pagerank else NA,
  pagerank_rank           = if (in_net) m$pagerank_rank else NA,
  betweenness_rank        = if (in_net) m$betweenness_rank else NA,
  # Density among neighbors: high = cohesive module, low = the gene bridges groups
  links_between_neighbors = links_between_neighbors,
  neighbor_density        = if (in_net && k > 1) links_between_neighbors / choose(k, 2) else NA,
  genes_up_to_order       = if (in_net) nrow(members) - 1 else NA,
  analysis_order          = analysis_order,
  bootstrap_partners      = nrow(support),
  max_bootstrap_support   = if (nrow(support) > 0) max(support$support) else 0
)
print(t(metrics))
write_tsv(metrics, out[["metrics"]])

message("\nDone.")
