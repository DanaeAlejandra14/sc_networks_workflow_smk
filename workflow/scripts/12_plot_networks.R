# 12_plot_networks.R
# Plots the ARACNe networks of one cell type and gamma: one figure per phenotype
# plus a panel with all phenotypes side by side.
#
# For each network:
#   1. Read network.txt and collapse A->B / B->A into one undirected edge (MI is symmetric)
#   2. Read the PageRank of every gene (node_metrics.tsv, from 11_network_metrics.R)
#   3. Keep the strongest edges (highest MI) until n_nodes genes are reached
#   4. Highlight the top n_hubs genes by PageRank among the genes shown
# Node size (PageRank) and edge width (MI) use the same scale in every network,
# so sizes can be compared across phenotypes and the panel has a single legend.

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(patchwork)
})

source(snakemake@input[["theme"]])
source(snakemake@input[["utils"]])  # read_aracne_network(), aracne_graph()

# Inputs, outputs and parameters from the Snakemake rule
network_files <- unlist(snakemake@input[["networks"]])
metrics_files <- unlist(snakemake@input[["metrics"]])      # same order as network_files
phenotypes    <- unlist(snakemake@params[["phenotypes"]])  # same order as network_files
cell_type     <- snakemake@wildcards[["cell_type"]]
gamma         <- snakemake@wildcards[["gamma"]]
n_nodes       <- snakemake@params[["n_nodes"]]
n_hubs        <- snakemake@params[["n_hubs"]]
seed          <- snakemake@params[["seed"]]

out_png   <- unlist(snakemake@output[["png"]])
out_pdf   <- unlist(snakemake@output[["pdf"]])
out_panel_png <- snakemake@output[["panel_png"]]
out_panel_pdf <- snakemake@output[["panel_pdf"]]

stopifnot(length(network_files) == length(phenotypes),
          length(metrics_files) == length(phenotypes))

# Functions

# Genes of the n_nodes strongest edges: walks edges from highest to lowest MI and
# takes genes in the order they first appear
select_top_mi_genes <- function(edges, n_nodes) {
  edges <- edges[order(edges$MI, decreasing = TRUE), ]
  genes_in_order <- unique(as.vector(rbind(edges$gene_a, edges$gene_b)))
  head(genes_in_order, n_nodes)
}

# Builds the reduced graph of one phenotype
build_network <- function(path, metrics_path, phenotype) {
  message("\n== ", phenotype, " | ", path)
  edges <- read_aracne_network(path)

  g_full <- aracne_graph(edges)
  node_metrics <- read.delim(metrics_path, stringsAsFactors = FALSE)
  pagerank <- setNames(node_metrics$pagerank, node_metrics$gene)
  message("Full network: ", vcount(g_full), " genes, ", ecount(g_full), " edges")

  # Reduced graph: selected genes plus every edge between them
  genes <- select_top_mi_genes(edges, n_nodes)
  g <- induced_subgraph(g_full, genes)
  g <- delete_vertices(g, degree(g) == 0)  # the last gene added can be left alone
  message("Shown: ", vcount(g), " genes, ", ecount(g), " edges, ",
          components(g)$no, " components")

  # Hubs: top n_hubs by PageRank of the full network, among the genes shown
  V(g)$pagerank <- pagerank[V(g)$name]
  hub_genes <- head(V(g)$name[order(V(g)$pagerank, decreasing = TRUE)], n_hubs)
  V(g)$hub <- V(g)$name %in% hub_genes
  message("Hubs: ", paste(hub_genes, collapse = ", "))

  list(graph = g, phenotype = phenotype)
}

# Draws one network. pr_limits / mi_limits are shared by all networks
draw_network <- function(net, pr_limits, mi_limits) {
  g <- net$graph
  col <- phenotype_color(net$phenotype)

  # Layout with a fixed seed so the figure is reproducible
  set.seed(seed)
  lay <- layout_with_fr(g, weights = E(g)$MI)

  ggraph(g, layout = "manual", x = lay[, 1], y = lay[, 2]) +
    geom_edge_link(aes(edge_width = MI), colour = pal_neutral[["line"]], alpha = 0.6) +
    geom_node_point(aes(size = pagerank, colour = hub), fill = col, shape = 21, stroke = 0.7) +
    geom_node_label(
      aes(filter = hub, label = name),
      repel = TRUE, size = 4.5, fontface = "bold", family = pub_font,
      label.size = 0.25, label.padding = unit(0.2, "lines"),
      min.segment.length = 0, max.overlaps = Inf
    ) +
    scale_edge_width(
      range = c(0.2, 1.6), limits = mi_limits,
      breaks = pretty(mi_limits, n = 4), name = "MI"
    ) +
    scale_size(
      range = c(2.5, 9.5), limits = pr_limits,
      breaks = signif(pretty(pr_limits, n = 4), 2), name = "PageRank",
      labels = scales::label_scientific(digits = 2)
    ) +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
    guides(
      # Grey circles: the legend is shared by all phenotypes
      size = guide_legend(order = 1, override.aes = list(fill = "grey60", colour = "grey30")),
      edge_width = guide_legend(order = 2)
    ) +
    labs(
      title    = net$phenotype,
      subtitle = paste0(cell_type, " | gamma ", gamma, " | ",
                        vcount(g), " genes, ", ecount(g), " edges")
    ) +
    theme_pub_void() +
    theme(plot.title = element_text(colour = col))
}

# Run

networks <- Map(build_network, network_files, metrics_files, phenotypes)

# Shared scales: range of PageRank and MI over all the networks shown
pr_limits <- range(unlist(lapply(networks, function(n) V(n$graph)$pagerank)))
mi_limits <- range(unlist(lapply(networks, function(n) E(n$graph)$MI)))
message("\nShared PageRank range: ", paste(signif(pr_limits, 3), collapse = " - "),
        " | MI range: ", paste(round(mi_limits, 2), collapse = " - "))

plots <- lapply(networks, draw_network, pr_limits = pr_limits, mi_limits = mi_limits)

# One figure per phenotype
for (i in seq_along(plots)) {
  save_pub(plots[[i]], out_png[i], out_pdf[i], width = 7, height = 6)
}

# Panel: all phenotypes in clinical order (theme_pub.R), 2 columns.
# The scales are identical, so guides = "collect" leaves a single legend
panel_order <- order(match(phenotypes, phenotype_levels))
panel <- wrap_plots(plots[panel_order], ncol = 2, guides = "collect") +
  plot_annotation(
    title = paste0(cell_type, " ARACNe networks (gamma ", gamma, ")"),
    theme = theme_pub_void(base_size = pub_base_size + 2)
  )
n_rows <- ceiling(length(plots) / 2)
save_pub(panel, out_panel_png, out_panel_pdf, width = 14, height = 6 * n_rows)

message("\nDone.")
