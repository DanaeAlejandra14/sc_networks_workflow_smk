# 15_ego_enrichment.R
# GO enrichment of the ego network of one gene in one ARACNe network (one
# phenotype, cell type and gamma): genes up to analysis_order (members.tsv,
# from 14_ego_network.R).
# Universe: every gene of the expression matrix given to ARACNe for this
# phenotype (centroids.txt), i.e. every gene that could have been in the network.
# Writes every tested term; the figure is made by 16_plot_ego_network.R

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

members   <- read.delim(snakemake@input[["members"]], stringsAsFactors = FALSE)
universe  <- readLines(snakemake@input[["universe"]])
phenotype <- snakemake@wildcards[["phenotype"]]
ego_gene  <- snakemake@wildcards[["gene"]]
out_table <- snakemake@output[["table"]]

# Columns kept in the output, the same whether there are results or not
cols <- c("phenotype", "ID", "Description", "GeneRatio", "BgRatio",
          "pvalue", "p.adjust", "qvalue", "Count", "geneID")
empty <- setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)

genes <- unique(members$gene)
if (length(genes) == 0) {
  message(ego_gene, " is not in the ", phenotype, " network: no enrichment")
  write.table(empty, out_table, sep = "\t", quote = FALSE, row.names = FALSE)
  quit(save = "no")
}
message(phenotype, ": ", length(genes), " ego network genes | universe: ", length(universe), " genes")

res <- enrichGO(
  gene = genes, universe = universe, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
  ont = snakemake@params[["ont"]], pAdjustMethod = "BH",
  pvalueCutoff = 1, qvalueCutoff = 1,  # keep every tested term; the cutoff is applied in the plot
  minGSSize = snakemake@params[["min_gs_size"]], maxGSSize = snakemake@params[["max_gs_size"]]
)

df <- if (is.null(res)) empty else as.data.frame(res)
if (nrow(df) > 0) {
  df$phenotype <- phenotype
  df <- df[order(df$p.adjust), cols]
}
message("Terms tested: ", nrow(df), " | p.adjust < 0.05: ", sum(df$p.adjust < 0.05))
write.table(df, out_table, sep = "\t", quote = FALSE, row.names = FALSE)

message("\nDone.")
