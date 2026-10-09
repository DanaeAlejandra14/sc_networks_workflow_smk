suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(scales)
})

phenotypes <- c("control", "MCI", "AD-NC_SYM", "AD-NC_ASYM")
cell_types <- c("Excitatory Neurons", "Oligodendrocytes", "Inhibitory Neurons",
                 "Astrocyte", "Endothelial", "OPCs", "Microglia")

# Paleta ggsci "npg" que pediste
npg_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)

celltype_colors  <- setNames(npg_colors[1:7], cell_types)      # 7 colores para tipo celular
phenotype_colors <- setNames(npg_colors[c(8, 2, 5, 9)], phenotypes)  # 4 colores distintos para fenotipo

## ---- 1. Cargar y agregar los datos ----------------------------------------

rows <- list()

for (ph in phenotypes) {
  for (ct in cell_types) {
    dir_path <- file.path("results/celltype_split", ph, ct)
    if (!dir.exists(dir_path)) next
    rds_files <- list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)
    message(ph, " / ", ct, ": ", length(rds_files), " specimen(s)")

    for (f in rds_files) {
      obj <- readRDS(f)
      donor <- unique(obj$individualID)
      donor <- donor[!is.na(donor)]
      donor <- if (length(donor) >= 1) donor[1] else NA_character_

      rows[[length(rows) + 1]] <- tibble(
        phenotype    = ph,
        cell_type    = ct,
        individualID = donor,
        n_cells      = ncol(obj)
      )
      rm(obj)
    }
  }
}

df <- bind_rows(rows) %>%
  group_by(phenotype, cell_type, individualID) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop") %>%
  mutate(
    phenotype = factor(phenotype, levels = phenotypes),
    cell_type = factor(cell_type, levels = cell_types)
  )

dir.create("results/plots/overview", recursive = TRUE, showWarnings = FALSE)

base_theme <- theme_minimal(base_size = 13) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line         = element_line(colour = "black", linewidth = 0.4),
    axis.ticks        = element_line(colour = "black", linewidth = 0.4),
    axis.text         = element_text(colour = "black")
  )

## ---- 2. Grafico 1: n_cells totales por fenotipo ----------------------------

df_total_by_phenotype <- df %>%
  group_by(phenotype) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

p1 <- ggplot(df_total_by_phenotype, aes(x = phenotype, y = n_cells, fill = phenotype)) +
  geom_col(width = 0.65, colour = NA) +
  geom_text(aes(label = comma(n_cells)), vjust = -0.6, size = 4, colour = "black") +
  scale_fill_manual(values = phenotype_colors, guide = "none") +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Número de células",
       title = "Número total de células por fenotipo") +
  base_theme

ggsave("results/plots/overview/01_cells_per_phenotype_total.png", p1,
       width = 7, height = 6, dpi = 150, bg = "white")

## ---- 3. Grafico 2: fraccion de cell_type por fenotipo (100% stacked) ------

df_by_celltype <- df %>%
  group_by(phenotype, cell_type) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

p2 <- ggplot(df_by_celltype, aes(x = phenotype, y = n_cells, fill = cell_type)) +
  geom_col(position = "fill", width = 0.65, colour = NA) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = celltype_colors) +
  labs(x = NULL, y = "Fracción de células", fill = "Tipo celular",
       title = "Fracción de tipos celulares por fenotipo") +
  base_theme +
  theme(legend.position = "right")

ggsave("results/plots/overview/02_celltype_fraction_per_phenotype.png", p2,
       width = 8, height = 6, dpi = 150, bg = "white")

## ---- 4. Grafico 3: bubble plot con numero gradual dentro del circulo -----

df_bubble <- df_by_celltype %>%
  mutate(
    text_size  = scales::rescale(n_cells, to = c(2.2, 6.2)),
    text_alpha = scales::rescale(n_cells, to = c(0.25, 1))
  )

p3 <- ggplot(df_by_celltype, aes(x = phenotype, y = cell_type)) +
  geom_point(aes(size = n_cells, colour = cell_type)) +
  scale_size_continuous(
    name   = "Cells",
    trans  = "log10",
    range  = c(4, 22),
    breaks = c(1, 10, 100, 1000, 10000),
    labels = comma_format()
  ) +
  scale_colour_manual(values = celltype_colors, guide = "none") +
  labs(x = NULL, y = NULL, title = "Número de células por fenotipo y tipo celular") +
  base_theme +
  theme(legend.position = "right")

ggsave("results/plots/overview/03_celltype_bubble_per_phenotype.png", p3,
       width = 9, height = 6.5, dpi = 150, bg = "white")

message("Listo. Graficas guardadas en results/plots/overview/")