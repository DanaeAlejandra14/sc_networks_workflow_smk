# 17_plot_ego_enrichment.R
# GO enrichment of the ego network of one gene across phenotypes (clinical
# order, theme_pub.R): top significant terms of each phenotype, from the
# tables of 15_ego_enrichment.R.

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(ggplot2))

source(snakemake@input[["theme"]])

# Inputs, outputs and parameters from the Snakemake rule
go_files    <- unlist(snakemake@input[["go"]])
phenotypes  <- unlist(snakemake@params[["phenotypes"]])  # same order as go_files
ego_gene    <- snakemake@wildcards[["gene"]]
padj_cutoff <- snakemake@params[["padj_cutoff"]]
top_terms   <- snakemake@params[["top_terms"]]
out         <- snakemake@output

stopifnot(length(go_files) == length(phenotypes))
go <- setNames(lapply(go_files, read.delim, stringsAsFactors = FALSE), phenotypes)

# Phenotypes in clinical order; phenotypes without a defined color go last
ord <- phenotype_levels[phenotype_levels %in% phenotypes]
ord <- c(ord, setdiff(phenotypes, ord))

keep_cols <- c("phenotype", "Description", "GeneRatio", "p.adjust", "Count")
go_all <- do.call(rbind, lapply(ord, function(p) go[[p]][, keep_cols]))
sig <- go_all[!is.na(go_all$p.adjust) & go_all$p.adjust < padj_cutoff, ]

if (nrow(sig) > 0) {
  # Top terms of each phenotype; each selected term is shown wherever it is significant
  top <- unlist(lapply(split(sig, sig$phenotype), function(d) head(d$Description[order(d$p.adjust)], top_terms)))
  go_df <- sig[sig$Description %in% top, ]
  go_df$phenotype <- factor(go_df$phenotype, levels = ord)
  best <- tapply(go_df$p.adjust, go_df$Description, min)
  go_df$term <- factor(go_df$Description, levels = rev(names(sort(best))))

  p_go <- ggplot(go_df, aes(x = phenotype, y = term)) +
    geom_point(aes(size = Count, fill = -log10(p.adjust)), shape = 21, colour = "black", stroke = 0.4) +
    scale_fill_gradientn(colours = pal_sequential, name = "-log10\n(adj. p)") +
    scale_size(range = c(3, 9), name = "Genes") +
    scale_x_discrete(drop = FALSE, position = "top") +
    # Long GO term names wrapped at 50 characters
    scale_y_discrete(labels = function(x) sapply(x, function(s) paste(strwrap(s, 50), collapse = "\n"))) +
    labs(title = paste0(ego_gene, " ego network: GO enrichment"), x = NULL, y = NULL) +
    theme_pub() +
    theme(
      panel.grid.major = element_line(colour = pal_neutral[["light"]], linewidth = 0.3),
      axis.line   = element_blank(),
      axis.ticks  = element_blank(),
      axis.text.x = element_text(size = rel(1.2), face = "bold"),
      axis.text.y = element_text(size = rel(1.05))
    )
  height_go <- max(4, 1.5 + 0.45 * length(unique(go_df$term)))
} else {
  p_go <- ggplot() +
    annotate("text", x = 0, y = 0, size = 5, family = pub_font,
             label = paste0("No GO terms with adjusted p < ", padj_cutoff)) +
    theme_pub_void()
  height_go <- 4
}
save_pub(p_go, out[["png"]], out[["pdf"]], width = 10, height = height_go)

message("Done.")
