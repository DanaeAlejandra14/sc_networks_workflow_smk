# 16_plot_ego_network.R
# Figures of the ego network of one gene across phenotypes (clinical order,
# theme_pub.R), from the tables of 14_ego_network.R:
#   ego_network        the gene and its first neighbors, one panel per phenotype;
#                      every gene keeps its position in all panels
#   ego_neighbors      dot matrix of neighbors x phenotypes (size = MI)
#   bootstrap_support  partners of the gene and the fraction of bootstraps with the edge
# The GO enrichment figure is made by 17_plot_ego_enrichment.R

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(patchwork)
})

source(snakemake@input[["theme"]])

# Inputs, outputs and parameters from the Snakemake rule
# Every input is a list with one file per phenotype, in the order of params.phenotypes
phenotypes  <- unlist(snakemake@params[["phenotypes"]])
ego_gene    <- snakemake@wildcards[["gene"]]
top_support <- snakemake@params[["top_support"]]
seed        <- snakemake@params[["seed"]]
out         <- snakemake@output

read_all <- function(name) {
  files <- unlist(snakemake@input[[name]])
  stopifnot(length(files) == length(phenotypes))
  setNames(lapply(files, read.delim, stringsAsFactors = FALSE), phenotypes)
}
neighbors <- read_all("neighbors")
edges     <- read_all("edges")
metrics   <- read_all("metrics")
support   <- read_all("support")

# Phenotypes in clinical order; phenotypes without a defined color go last
ord <- phenotype_levels[phenotype_levels %in% phenotypes]
ord <- c(ord, setdiff(phenotypes, ord))
n_rows <- ceiling(length(ord) / 2)

in_net <- sapply(ord, function(p) isTRUE(as.logical(metrics[[p]]$in_network)))
first_neighbors <- lapply(setNames(ord, ord), function(p) as.character(neighbors[[p]]$neighbor))
all_neighbors <- unique(unlist(first_neighbors))
neighbor_counts <- table(unlist(first_neighbors))
shared <- names(neighbor_counts)[neighbor_counts >= 2]  # neighbors in >= 2 phenotypes

# Shared MI range for the figures that use MI
mi_values <- unlist(lapply(ord, function(p) c(neighbors[[p]]$MI, edges[[p]]$MI)))
mi_limits <- if (length(mi_values) > 0) range(mi_values) else c(0, 1)

# Placeholder when there is nothing to draw
empty_plot <- function(text, title = NULL, col = "black") {
  ggplot() + annotate("text", x = 0, y = 0, label = text, size = 5, family = pub_font) +
    labs(title = title) + theme_pub_void() + theme(plot.title = element_text(colour = col))
}

# Ego network (first neighbors)
if (length(all_neighbors) > 0) {
  # Layout computed once on the union of all ego networks
  union_nodes <- unique(c(ego_gene, all_neighbors))
  union_edges <- unique(do.call(rbind, lapply(ord, function(p) edges[[p]][, c("from", "to")])))
  g_union <- graph_from_data_frame(union_edges, directed = FALSE,
                                   vertices = data.frame(name = union_nodes))
  set.seed(seed)
  lay <- layout_with_fr(g_union)
  lay <- sweep(lay, 2, lay[match(ego_gene, union_nodes), ])  # gene at the center
  rownames(lay) <- union_nodes

  plot_ego <- function(p) {
    col <- phenotype_color(p)
    present <- c(ego_gene, first_neighbors[[p]])
    nodes <- data.frame(name = union_nodes)
    nodes$status <- ifelse(!in_net[[p]] | !nodes$name %in% present, "absent",
                    ifelse(nodes$name == ego_gene, "ego", "neighbor"))
    nodes$outline <- ifelse(nodes$name %in% shared, "Shared (>= 2 phenotypes)", "Phenotype-specific")
    g <- graph_from_data_frame(edges[[p]], directed = FALSE, vertices = nodes)

    subtitle <- if (in_net[[p]]) {
      paste0(length(first_neighbors[[p]]), " neighbors")
    } else {
      paste0(ego_gene, " not in the network")
    }

    ggraph(g, layout = "manual", x = lay[V(g)$name, 1], y = lay[V(g)$name, 2]) +
      geom_edge_link(aes(edge_width = MI), colour = pal_neutral[["line"]], alpha = 0.7) +
      # Neighbors in other phenotypes only: faded, so each gene keeps its place
      geom_node_point(aes(filter = status == "absent"), size = 4, colour = pal_neutral[["light"]]) +
      geom_node_text(aes(filter = status == "absent", label = name), size = 3.4,
                     colour = "grey70", family = pub_font, repel = TRUE) +
      # show.legend for colour only: draw both "Neighbor" keys even if a panel has no
      # shared (or no) neighbors, so all panels have the same legend and patchwork
      # merges them; the layer stays out of the MI legend
      geom_node_point(aes(filter = status == "neighbor", colour = outline), fill = col,
                      shape = 21, size = 8.5, stroke = 1,
                      show.legend = c(colour = TRUE, edge_width = FALSE)) +
      geom_node_point(aes(filter = status == "ego"), fill = pal_highlight, colour = "black",
                      shape = 21, size = 13, stroke = 1.2) +
      # Gene names in a box, as in the network figures (12_plot_networks.R)
      geom_node_label(aes(filter = status != "absent", label = name,
                          fontface = ifelse(status == "ego", "bold", "plain")),
                      size = 5, family = pub_font, repel = TRUE,
                      label.size = 0.25, label.padding = unit(0.2, "lines"),
                      point.padding = unit(0.6, "lines"), min.segment.length = 0,
                      max.overlaps = Inf) +
      scale_edge_width(range = c(0.4, 2.5), limits = mi_limits, name = "MI") +
      scale_colour_manual(values = c("Shared (>= 2 phenotypes)" = "black", "Phenotype-specific" = "white"),
                          name = "Neighbor", drop = FALSE,
                          # Both values always listed, so every panel has the same legend
                          limits = c("Shared (>= 2 phenotypes)", "Phenotype-specific")) +
      guides(colour = guide_legend(override.aes = list(fill = "grey60", size = 6))) +
      labs(title = p, subtitle = subtitle) +
      theme_pub_void() +
      theme(plot.title = element_text(colour = col))
  }
  p_net <- wrap_plots(lapply(ord, plot_ego), ncol = 2, guides = "collect") +
    plot_annotation(title = paste0(ego_gene, " ego network"),
                    theme = theme_pub_void(base_size = pub_base_size + 2))
} else {
  p_net <- empty_plot(paste(ego_gene, "has no neighbors in any network"))
}
save_pub(p_net, out[["network_png"]], out[["network_pdf"]], width = 14, height = 6.5 * n_rows)

# Neighbors x phenotypes
if (length(all_neighbors) > 0) {
  dm <- do.call(rbind, lapply(ord, function(p) neighbors[[p]][, c("phenotype", "neighbor", "MI")]))
  # Rows: neighbors shared by more phenotypes first, then by strongest MI
  n_pheno <- tapply(dm$phenotype, dm$neighbor, length)
  max_mi  <- tapply(dm$MI, dm$neighbor, max)
  row_order <- names(sort(n_pheno * 10 + max_mi / max(max_mi), decreasing = TRUE))
  dm$neighbor  <- factor(dm$neighbor, levels = rev(row_order))
  dm$phenotype <- factor(dm$phenotype, levels = ord)
  dm$fill <- phenotype_color(as.character(dm$phenotype))

  p_dm <- ggplot(dm, aes(x = phenotype, y = neighbor)) +
    geom_point(aes(size = MI, fill = fill), shape = 21, colour = "black", stroke = 0.4) +
    scale_fill_identity() +
    scale_size(range = c(2.5, 9), limits = mi_limits, name = paste0("MI with\n", ego_gene)) +
    scale_x_discrete(drop = FALSE, position = "top") +
    labs(title = paste0(ego_gene, " neighbors"), x = NULL, y = NULL) +
    theme_pub() +
    theme(
      panel.grid.major = element_line(colour = pal_neutral[["light"]], linewidth = 0.3),
      axis.line   = element_blank(),
      axis.ticks  = element_blank(),
      axis.text.x = element_text(size = rel(1.2), face = "bold"),
      axis.text.y = element_text(size = rel(1.2))
    )
  height_dm <- max(4, 1.5 + 0.3 * length(all_neighbors))
} else {
  p_dm <- empty_plot(paste(ego_gene, "has no neighbors in any network"))
  height_dm <- 4
}
save_pub(p_dm, out[["neighbors_png"]], out[["neighbors_pdf"]], width = 7, height = height_dm)

# Bootstrap support
plot_support <- function(p) {
  col <- phenotype_color(p)
  d <- support[[p]]
  n_partners <- nrow(d)
  if (n_partners == 0) return(empty_plot(paste(ego_gene, "has no edges in the bootstraps"), p, col))
  d <- head(d[order(-d$n_bootstraps, d$partner), ], top_support)
  d$partner <- factor(d$partner, levels = rev(d$partner))
  d$final <- ifelse(as.logical(d$in_final_network), "Yes", "No")
  ggplot(d, aes(x = support, y = partner, alpha = final)) +
    geom_col(fill = col, width = 0.75) +
    scale_alpha_manual(values = c(Yes = 1, No = 0.3), name = "In final\nnetwork", drop = FALSE,
                       limits = c("Yes", "No")) +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
    labs(title = p,
         subtitle = paste0(n_partners, " distinct partners in ", d$n_total[1], " bootstraps"),
         x = "Bootstraps with the edge", y = NULL) +
    theme_pub() +
    theme(
      plot.title   = element_text(colour = col),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      axis.text.y  = element_text(size = rel(1.2)),
      axis.text.x  = element_text(size = rel(1.1))
    )
}
p_sup <- wrap_plots(lapply(ord, plot_support), ncol = 2, guides = "collect") +
  plot_annotation(title = paste0(ego_gene, " partners in the bootstrap networks"),
                  theme = theme_pub_void(base_size = pub_base_size + 2))
save_pub(p_sup, out[["support_png"]], out[["support_pdf"]], width = 13,
         height = n_rows * (1.6 + 0.3 * top_support))

message("Done.")
