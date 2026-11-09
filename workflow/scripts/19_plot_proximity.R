# 19_plot_proximity.R
# Differential network proximity of one gene across phenotypes (clinical order,
# theme_pub.R), from the tables of 18_network_proximity.R.
#   Figure: one panel per gene set; z-score of each bootstrap network (points),
#           their distribution (box) and the final network (diamond)
#   Table:  each phenotype against the reference phenotype (config
#           proximity.reference), using the bootstrap z-scores as replicates:
#           difference of medians and Wilcoxon test (BH-adjusted)

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(ggplot2))

source(snakemake@input[["theme"]])

tables     <- unlist(snakemake@input[["tables"]])
phenotypes <- unlist(snakemake@params[["phenotypes"]])
reference  <- snakemake@params[["reference"]]
ego_gene   <- snakemake@wildcards[["gene"]]
out        <- snakemake@output

prox <- do.call(rbind, lapply(tables, read.delim, stringsAsFactors = FALSE))

# Phenotypes in clinical order; phenotypes without a defined color go last
ord <- phenotype_levels[phenotype_levels %in% phenotypes]
ord <- c(ord, setdiff(phenotypes, ord))
prox$phenotype <- factor(prox$phenotype, levels = ord)
# Panel titles from the config names: "ad_risk" -> "AD risk",
# "glutamate_glutamine" -> "Glutamate glutamine"
gene_set_label <- function(x) {
  s <- gsub("_", " ", x)
  s <- gsub("\\bad\\b", "AD", s)  # acronyms in upper case
  paste0(toupper(substr(s, 1, 1)), substring(s, 2))
}
prox$gene_set_label <- gene_set_label(prox$gene_set)

boot  <- prox[prox$network_type == "bootstrap" & !is.na(prox$z), ]
final <- prox[prox$network_type == "final" & !is.na(prox$z), ]

# Each phenotype against the reference
diff_tab <- do.call(rbind, lapply(unique(prox$gene_set), function(s) {
  ref_z <- boot$z[boot$gene_set == s & boot$phenotype == reference]
  do.call(rbind, lapply(setdiff(ord, reference), function(p) {
    z <- boot$z[boot$gene_set == s & boot$phenotype == p]
    test <- if (length(z) >= 3 && length(ref_z) >= 3) suppressWarnings(wilcox.test(z, ref_z)) else NULL
    data.frame(
      gene = ego_gene, gene_set = s, phenotype = p, reference = reference,
      n_bootstraps = length(z), n_bootstraps_reference = length(ref_z),
      median_z = if (length(z)) median(z) else NA,
      median_z_reference = if (length(ref_z)) median(ref_z) else NA,
      # Negative delta: the gene is closer to the set than in the reference
      delta_z = if (length(z) && length(ref_z)) median(z) - median(ref_z) else NA,
      p_value = if (is.null(test)) NA else test$p.value
    )
  }))
}))
diff_tab$p_adj <- p.adjust(diff_tab$p_value, method = "BH")
print(diff_tab, digits = 3)
write.table(diff_tab, out[["table"]], sep = "\t", quote = FALSE, row.names = FALSE)

# Figure
# Adjusted p-value of each phenotype against the reference, above its box
diff_tab$gene_set_label <- gene_set_label(diff_tab$gene_set)
diff_tab$phenotype <- factor(diff_tab$phenotype, levels = ord)
diff_tab$label <- ifelse(is.na(diff_tab$p_adj), "",
                         ifelse(diff_tab$p_adj < 0.001, "p<0.001", sprintf("p=%.3f", diff_tab$p_adj)))
y_top <- max(c(boot$z, final$z), na.rm = TRUE)
y_bottom <- min(c(boot$z, final$z, -2), na.rm = TRUE)

p <- ggplot(boot, aes(x = phenotype, y = z)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = pal_neutral[["secondary"]]) +
  geom_hline(yintercept = -2, linetype = "dotted", colour = pal_neutral[["line"]]) +
  geom_boxplot(aes(colour = phenotype), fill = NA, width = 0.55, outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(aes(fill = phenotype), shape = 21, colour = "white", stroke = 0.3,
              size = 2.6, width = 0.12, height = 0, alpha = 0.85) +
  geom_point(data = final, aes(fill = phenotype), shape = 23, size = 4.5,
             colour = "black", stroke = 0.9) +
  geom_text(data = diff_tab, aes(x = phenotype, y = y_top + 0.15 * (y_top - y_bottom), label = label),
            size = 3.6, family = pub_font, colour = pal_neutral[["secondary"]]) +
  scale_fill_phenotype(guide = "none") +
  scale_color_phenotype(guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  facet_wrap(~ gene_set_label, nrow = 1) +
  labs(
    title = paste0(ego_gene, " network proximity"),
    x = NULL, y = "Proximity z-score\n(negative = closer than chance)"
  ) +
  theme_pub() +
  theme(
    strip.text  = element_text(size = rel(1.2), face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = rel(1.15), face = "bold"),
    axis.text.y = element_text(size = rel(1.1)),
    axis.title.y = element_text(size = rel(1.1)),
    panel.spacing.x = unit(1.5, "lines")
  )

n_sets <- length(unique(prox$gene_set))
save_pub(p, out[["png"]], out[["pdf"]], width = 3.6 * n_sets + 1.5, height = 6)

# Vertical version: one gene set per row
p_vertical <- p + facet_wrap(~ gene_set_label, ncol = 1)
save_pub(p_vertical, out[["png_vertical"]], out[["pdf_vertical"]],
         width = 7.5, height = 3.2 * n_sets + 1.5)

message("\nDone.")
