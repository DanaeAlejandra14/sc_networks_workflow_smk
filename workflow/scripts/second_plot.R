suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(scales)
  library(ggsci)
  library(ggpubr)
})

phenotypes <- c("control", "MCI", "AD-NC_SYM", "AD-NC_ASYM")
cell_types <- c("Excitatory Neurons", "Oligodendrocytes", "Inhibitory Neurons",
                 "Astrocyte", "Endothelial", "OPCs", "Microglia")

##  1. Cargar y agregar los datos 

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
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line          = element_line(colour = "black", linewidth = 0.4),
    axis.ticks         = element_line(colour = "black", linewidth = 0.4),
    axis.text          = element_text(colour = "black")
  )

## 1: n_cells totales por fenotipo

df_total_by_phenotype <- df %>%
  group_by(phenotype) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

p1 <- ggplot(df_total_by_phenotype, aes(x = phenotype, y = n_cells, fill = phenotype)) +
  geom_col(width = 0.65, colour = NA) +
  geom_text(aes(label = comma(n_cells)), vjust = -0.6, size = 4, colour = "black") +
  scale_fill_npg(guide = "none") +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Número de células",
       title = "Número total de células por fenotipo") +
  base_theme

ggsave("results/plots/overview/01_cells_per_phenotype_total.png", p1,
       width = 7, height = 6, dpi = 150, bg = "white")

## -2. fraccion de cell_type por fenotipo (100% stacked) ------

df_by_celltype <- df %>%
  group_by(phenotype, cell_type) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

p2 <- ggplot(df_by_celltype, aes(x = phenotype, y = n_cells, fill = cell_type)) +
  geom_col(position = "fill", width = 0.65, colour = NA) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_npg() +
  labs(x = NULL, y = "Fracción de células", fill = "Tipo celular",
       title = "Fracción de tipos celulares por fenotipo") +
  base_theme +
  theme(legend.position = "right")

ggsave("results/plots/overview/02_celltype_fraction_per_phenotype.png", p2,
       width = 8, height = 6, dpi = 150, bg = "white")

## 4. Grafico 3: bubble plot con ggballoonplot 

p3 <- ggballoonplot(
  df_by_celltype,
  x    = "phenotype",
  y    = "cell_type",
  size = "n_cells",
  fill = "cell_type",
  size.range = c(3, 18),
  ggtheme = base_theme
) +
  scale_fill_npg(guide = "none") +
  scale_size_continuous(
    range  = c(3, 18),
    trans  = "log10",
    breaks = c(1, 10, 50, 100, 500, 1000, 5000, 10000),
    labels = comma_format()
  ) +
  labs(x = NULL, y = NULL, size = "Cells",
       title = "Número de células por fenotipo y tipo celular")

ggsave("results/plots/overview/03_celltype_bubble_per_phenotype.png", p3,
       width = 9, height = 6.5, dpi = 150, bg = "white")

message("plots guardados en  results/plots/overview/")

## ---- 5. Resumen en texto/CSV con los valores reales -----------------------

# Tabla 1: total de celulas por fenotipo
write_csv(df_total_by_phenotype, "results/plots/overview/summary_cells_per_phenotype.csv")

# Tabla 2: celulas por fenotipo x tipo celular (con fraccion incluida)
df_by_celltype_pct <- df_by_celltype %>%
  group_by(phenotype) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()

write_csv(df_by_celltype_pct, "results/plots/overview/summary_cells_per_phenotype_celltype.csv")

# Resumen en texto plano, facil de leer de un vistazo
resumen_txt <- c(
  "RESUMEN DE CELULAS POR FENOTIPO (7 tipos celulares seleccionados)",
  "====================================================================",
  "",
  "Total por fenotipo:",
  paste0("  ", df_total_by_phenotype$phenotype, ": ",
         comma(df_total_by_phenotype$n_cells), " celulas"),
  "",
  "Desglose por tipo celular dentro de cada fenotipo:",
  unlist(lapply(unique(df_by_celltype_pct$phenotype), function(ph) {
    sub <- df_by_celltype_pct %>% filter(phenotype == ph) %>% arrange(desc(n_cells))
    c(
      paste0("  ", ph, ":"),
      paste0("    ", sub$cell_type, ": ", comma(sub$n_cells),
             " celulas (", percent(sub$fraction, accuracy = 0.1), ")")
    )
  }))
)

writeLines(resumen_txt, "results/plots/overview/summary.txt")

message("Resumen guardado en:")
message("  results/plots/overview/summary_cells_per_phenotype.csv")
message("  results/plots/overview/summary_cells_per_phenotype_celltype.csv")
message("  results/plots/overview/summary.txt")