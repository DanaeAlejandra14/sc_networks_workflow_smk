suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(scales)
  library(patchwork)
})

phenotypes <- c("control", "MCI", "AD-NC_SYM", "AD-NC_ASYM")
cell_types <- c("Excitatory Neurons", "Inhibitory Neurons", "Oligodendrocytes",
                "Astrocyte", "OPCs", "Microglia", "Endothelial")

out_dir <- "results/plots/overview"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Cargar y agregar los datos (misma lógica que el script original) ----

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

df_by_donor <- df %>%
  group_by(phenotype, individualID) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

df_by_celltype <- df %>%
  group_by(phenotype, cell_type) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

df_total_by_phenotype <- df %>%
  group_by(phenotype) %>%
  summarise(n_cells = sum(n_cells), .groups = "drop")

## ---- 2. Paletas y tema (estilo de la figura de referencia) ------------------

# Un color por tipo celular, agrupados por linaje:
# neuronas (azul/morado), glía (dorado/verde azulado/amarillo/rojo), vascular (café)
cell_colors <- c(
  "Excitatory Neurons" = "#1B2A5C",
  "Inhibitory Neurons" = "#8E4A8F",
  "Oligodendrocytes"   = "#C4A044",
  "Astrocyte"          = "#2F6577",
  "OPCs"               = "#E8D07A",
  "Microglia"          = "#B5413B",
  "Endothelial"        = "#6B4F40"
)

donors <- sort(unique(na.omit(df_by_donor$individualID)))
donor_colors <- setNames(
  colorRampPalette(c("#2B3445", "#A3A93A", "#C8C3A8"))(length(donors)),
  donors
)
show_donor_legend <- length(donors) <= 12   # con muchos donadores la leyenda no se lee

base_theme <- theme_classic(base_size = 13) +
  theme(
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA),
    axis.line          = element_line(colour = "black", linewidth = 0.4),
    axis.ticks         = element_line(colour = "black", linewidth = 0.4),
    axis.text          = element_text(colour = "black"),
    legend.key.size    = unit(0.45, "cm"),
    legend.title       = element_text(face = "bold"),
    plot.tag           = element_text(face = "bold", size = 14)
  )

no_x <- theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank())

## ---- 3. Panel a: número de células por fenotipo, apilado por donador --------

p_a <- ggplot(df_by_donor, aes(x = phenotype, y = n_cells, fill = individualID)) +
  geom_col(width = 0.6, colour = NA) +
  scale_fill_manual(
    values = donor_colors, name = "Donor",
    guide  = if (show_donor_legend) guide_legend() else "none"
  ) +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Number of cells") +
  base_theme + no_x

## ---- 4. Panel b: fracción de tipos celulares por fenotipo -------------------

p_b <- ggplot(df_by_celltype, aes(x = phenotype, y = n_cells, fill = cell_type)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.6, colour = NA) +
  scale_fill_manual(values = cell_colors, guide = "none") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(y = "Fraction of cells") +
  base_theme + no_x

## ---- 5. Panel c: bubble plot (tipo celular x fenotipo) ----------------------

p_c <- ggplot(df_by_celltype,
              aes(x = phenotype, y = cell_type, size = n_cells, colour = cell_type)) +
  geom_point() +
  scale_colour_manual(values = cell_colors, guide = "none") +
  scale_size_continuous(
    name   = "Cells",
    range  = c(1.5, 14),
    trans  = "log10",
    breaks = c(1, 10, 100, 1000, 10000),
    labels = comma_format()
  ) +
  scale_y_discrete(limits = rev(cell_types)) +
  labs(x = NULL, y = NULL) +
  base_theme +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1)
  )

## ---- 6. Guardar: paneles individuales + figura combinada --------------------

save_plot <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
}

save_plot(p_a + theme(axis.text.x = element_text(), axis.ticks.x = element_line(),
                      axis.title.x = element_blank()),
          "01_cells_per_phenotype_total", 7, 6)
save_plot(p_b + theme(axis.text.x = element_text(), axis.ticks.x = element_line(),
                      axis.title.x = element_blank()) +
            scale_fill_manual(values = cell_colors, name = "Cell type"),
          "02_celltype_fraction_per_phenotype", 8, 6)
save_plot(p_c, "03_celltype_bubble_per_phenotype", 7, 6)

fig <- (p_a / p_b / p_c) +
  plot_layout(heights = c(1, 1.2, 1.4)) +
  plot_annotation(tag_levels = "a")

save_plot(fig, "00_overview_combined", 8, 13)

message("plots guardados en ", out_dir)

## ---- 7. Resumen en texto/CSV con los valores reales -------------------------

write_csv(df_total_by_phenotype, file.path(out_dir, "summary_cells_per_phenotype.csv"))

df_by_celltype_pct <- df_by_celltype %>%
  group_by(phenotype) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()

write_csv(df_by_celltype_pct, file.path(out_dir, "summary_cells_per_phenotype_celltype.csv"))
write_csv(df_by_donor, file.path(out_dir, "summary_cells_per_phenotype_donor.csv"))

resumen_txt <- c(
  "RESUMEN DE CELULAS POR FENOTIPO (7 tipos celulares seleccionados)",
  "====================================================================",
  "",
  "Total por fenotipo:",
  paste0("  ", df_total_by_phenotype$phenotype, ": ",
         comma(df_total_by_phenotype$n_cells), " celulas"),
  "",
  "Desglose por tipo celular dentro de cada fenotipo:",
  unlist(lapply(levels(df_by_celltype_pct$phenotype), function(ph) {
    sub <- df_by_celltype_pct %>% filter(phenotype == ph) %>% arrange(desc(n_cells))
    c(
      paste0("  ", ph, ":"),
      paste0("    ", sub$cell_type, ": ", comma(sub$n_cells),
             " celulas (", percent(sub$fraction, accuracy = 0.1), ")")
    )
  }))
)

writeLines(resumen_txt, file.path(out_dir, "summary.txt"))
message("Resumen guardado en ", out_dir)