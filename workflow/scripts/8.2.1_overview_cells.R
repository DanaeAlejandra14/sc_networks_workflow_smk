# Overview of cells per phenotype x cell type x donor at one stage of the workflow
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(scales)
  library(patchwork)
})

files      <- unlist(snakemake@input[["objects"]])
out_dir    <- snakemake@output[[1]]
phenotypes <- unlist(snakemake@params[["phenotypes"]])
cell_types <- unlist(snakemake@params[["cell_types"]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Project colors and theme
source(snakemake@input[["theme"]])

## Data: cells per phenotype x cell type x donor 

rows <- lapply(files, function(f) {
  ct <- basename(dirname(f))
  ph <- basename(dirname(dirname(f)))
  message(ph, " / ", ct)
  obj <- readRDS(f)
  donors <- as.character(obj$individualID)
  rm(obj); invisible(gc())
  tibble(phenotype = ph, cell_type = ct, individualID = donors) %>%
    count(phenotype, cell_type, individualID, name = "n_cells")
})


df <- bind_rows(rows) %>%
  mutate(
    phenotype = factor(phenotype, levels = phenotypes),
    cell_type = factor(cell_type, levels = cell_types)
  )

df_by_donor <- df %>%
  group_by(phenotype, individualID) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

df_by_celltype <- df %>%
  group_by(phenotype, cell_type) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

df_total_by_phenotype <- df %>%
  group_by(phenotype) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

## Palettes
# Cell type colors come from theme_pub.R (grey if a cell type is missing there)
cell_colors <- setNames(
  unname(ifelse(cell_types %in% names(pal_celltype), pal_celltype[cell_types], "grey60")),
  cell_types
)

donors <- sort(unique(na.omit(df_by_donor$individualID)))
donor_colors <- setNames(
  colorRampPalette(c("#2B3445", "#A3A93A", "#C8C3A8"))(length(donors)),
  donors
)
show_donor_legend <- length(donors) <= 12   # with many donors the legend is unreadable

no_x <- theme(axis.text.x  = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank())

## Panel a: cells per phenotype, stacked by donor 
p_a <- ggplot(df_by_donor, aes(x = phenotype, y = n_cells, fill = individualID)) +
  geom_col(width = 0.6, colour = NA) +
  scale_fill_manual(values = donor_colors, name = "Donor",
                    guide = if (show_donor_legend) guide_legend() else "none") +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Number of cells") +
  theme_pub() + no_x

## Panel b: cell type fraction per phenotype 
p_b <- ggplot(df_by_celltype, aes(x = phenotype, y = n_cells, fill = cell_type)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.6, colour = NA) +
  scale_fill_manual(values = cell_colors, guide = "none") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(y = "Fraction of cells") +
  theme_pub() + no_x

## Panel c: bubbles (size = number of cells) 
p_c <- ggplot(df_by_celltype,
              aes(x = phenotype, y = cell_type, size = n_cells, colour = cell_type)) +
  geom_point() +
  scale_colour_manual(values = cell_colors, guide = "none") +
  scale_size_continuous(name = "Cells", range = c(1.5, 14), trans = "log10",
                        breaks = c(1, 10, 100, 1000, 10000), labels = comma_format()) +
  scale_y_discrete(limits = rev(cell_types)) +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
        axis.line        = element_blank(),
        axis.ticks       = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1))

# Save 
save_plot <- function(p, name, w, h) {
  save_pub(p, file.path(out_dir, paste0(name, ".png")),
           file.path(out_dir, paste0(name, ".pdf")), width = w, height = h)
}

show_x <- theme(axis.text.x = element_text(), axis.ticks.x = element_line(),
                axis.title.x = element_blank())

save_plot(p_a + show_x, "01_cells_per_phenotype_total", 7, 6)
save_plot(p_b + show_x + scale_fill_manual(values = cell_colors, name = "Cell type"),
          "02_celltype_fraction_per_phenotype", 8, 6)
save_plot(p_c, "03_celltype_bubble_per_phenotype", 7, 6)

fig <- (p_a / p_b / p_c) +
  plot_layout(heights = c(1, 1.2, 1.4)) +
  plot_annotation(tag_levels = "a", title = "After subsampling")
save_plot(fig, "00_overview_combined", 8, 13)

##  Tables
write_csv(df_total_by_phenotype, file.path(out_dir, "summary_cells_per_phenotype.csv"))
write_csv(df_by_donor,           file.path(out_dir, "summary_cells_per_phenotype_donor.csv"))
write_csv(df_by_celltype %>% group_by(phenotype) %>%
            mutate(fraction = n_cells / sum(n_cells)) %>% ungroup(),
          file.path(out_dir, "summary_cells_per_phenotype_celltype.csv"))

message("Done. Outputs in ", out_dir)
sink(type = "message"); sink()