# 18_network_proximity.R
# Network proximity of one gene to predefined gene sets (config proximity.gene_sets),
# in the final ARACNe network of one phenotype and in each of its bootstrap networks.
#
# Proximity: mean shortest-path distance (number of edges) from the gene to the
# genes of the set present in the network. Genes the gene cannot reach count as
# (longest finite distance in the network + 1).
#
# Null model: the same distance from random genes to random sets with a similar
# number of connections (degree bins of >= bin_size genes), so a hub is not called
# "close" just because it is a hub, and each network is compared with itself
# (networks differ in size and density).
#   z < 0  -> the gene is closer to the set than expected by chance
#
# The bootstrap networks act as replicates: they give a distribution per phenotype,
# and also a value when the gene is missing from the final network.

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(igraph))

source(snakemake@input[["utils"]])  # read_aracne_network(), aracne_graph()

# Inputs, outputs and parameters from the Snakemake rule
network_file <- snakemake@input[["network"]]
boot_files   <- unlist(snakemake@input[["bootstraps"]])
ego_gene     <- snakemake@wildcards[["gene"]]
phenotype    <- snakemake@wildcards[["phenotype"]]
gene_sets    <- lapply(snakemake@params[["gene_sets"]], unlist)
n_sources    <- snakemake@params[["n_random_sources"]]
n_sets       <- snakemake@params[["n_random_sets"]]
bin_size     <- snakemake@params[["bin_size"]]
seed         <- snakemake@params[["seed"]]

# Degree bins: genes sorted by degree, cut into groups of >= bin_size genes,
# never splitting genes with the same degree
degree_bins <- function(deg) {
  levels_k <- sort(unique(deg))
  counts <- table(factor(deg, levels = levels_k))
  bin_of_k <- integer(length(levels_k))
  b <- 1
  n <- 0
  for (i in seq_along(levels_k)) {
    bin_of_k[i] <- b
    n <- n + counts[[i]]
    if (n >= bin_size) {
      b <- b + 1
      n <- 0
    }
  }
  # A last bin smaller than bin_size is merged with the previous one
  if (n > 0 && n < bin_size && b > 1) bin_of_k[bin_of_k == b] <- b - 1
  setNames(bin_of_k[match(deg, levels_k)], names(deg))
}

# Random gene with the same degree bin as each gene in `genes`
sample_matched <- function(genes, bins, pool_by_bin) {
  vapply(genes, function(x) {
    pool <- pool_by_bin[[as.character(bins[[x]])]]
    pool[sample.int(length(pool), 1)]
  }, character(1))
}

# Proximity of the gene to every gene set in one network
proximity_in_network <- function(g, network_id) {
  genes <- V(g)$name
  base <- data.frame(phenotype = phenotype, gene = ego_gene, network = network_id,
                     gene_set = names(gene_sets),
                     set_size = sapply(gene_sets, function(s) length(setdiff(s, ego_gene))),
                     set_in_network = sapply(gene_sets, function(s) sum(setdiff(s, ego_gene) %in% genes)))
  if (!ego_gene %in% genes) {
    message(network_id, ": ", ego_gene, " not in the network")
    return(cbind(base, source_in_network = FALSE, reachable = NA, distance = NA,
                 null_mean = NA, null_sd = NA, z = NA, p_value = NA))
  }

  deg <- degree(g)
  bins <- degree_bins(deg)
  pool_by_bin <- split(names(bins), bins)

  # Random sources (degree-matched to the gene); one BFS per source
  set.seed(seed)
  sources <- sample_matched(rep(ego_gene, n_sources), bins, pool_by_bin)
  # Rows: gene + random sources. weights = NA: count edges (igraph would
  # otherwise use the MI weights as lengths)
  d_all <- distances(g, v = c(ego_gene, sources), weights = NA)
  max_finite <- max(d_all[is.finite(d_all)])
  d_all[!is.finite(d_all)] <- max_finite + 1      # unreachable genes

  res <- lapply(names(gene_sets), function(s) {
    targets <- intersect(setdiff(gene_sets[[s]], ego_gene), genes)
    if (length(targets) == 0) {
      return(data.frame(source_in_network = TRUE, reachable = NA, distance = NA,
                        null_mean = NA, null_sd = NA, z = NA, p_value = NA))
    }
    observed <- mean(d_all[1, targets])
    reachable <- sum(d_all[1, targets] <= max_finite)
    # Null: each random source with n_sets degree-matched random sets
    null <- unlist(lapply(seq_len(n_sources), function(i) {
      replicate(n_sets, mean(d_all[i + 1, sample_matched(targets, bins, pool_by_bin)]))
    }))
    data.frame(
      source_in_network = TRUE,
      reachable = reachable,
      distance  = observed,
      null_mean = mean(null),
      null_sd   = sd(null),
      z         = if (sd(null) > 0) (observed - mean(null)) / sd(null) else NA,
      # One-sided: is the gene closer than chance?
      p_value   = (sum(null <= observed) + 1) / (length(null) + 1)
    )
  })
  out <- cbind(base, do.call(rbind, res))
  message(network_id, ": ", paste(sprintf("%s z=%.2f", out$gene_set, out$z), collapse = " | "))
  out
}

# Final network + every bootstrap network
networks <- c(final = network_file,
              setNames(boot_files, sub("\\.txt$", "", sub("^bootstrapNetwork_", "", basename(boot_files)))))

results <- do.call(rbind, lapply(names(networks), function(id) {
  g <- aracne_graph(read_aracne_network(networks[[id]]))
  proximity_in_network(g, id)
}))
results$network_type <- ifelse(results$network == "final", "final", "bootstrap")
rownames(results) <- NULL

write.table(results, snakemake@output[["table"]], sep = "\t", quote = FALSE, row.names = FALSE)
message("\nDone.")
