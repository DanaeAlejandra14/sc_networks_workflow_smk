log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

annotated_dirs  <- snakemake@input[["annotated_dirs"]]
diagnostics_dir <- snakemake@output[["diagnostics_dir"]]
max_cells_per_specimen <- 1000

dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)

message("Doing diagnostics for ", length(annotated_dirs), " batch(es)")

obj_list  <- list()
meta_list <- list()  # Complete metadata

for (dir in annotated_dirs) {
  library_batch <- basename(dir)
  rds_files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)

  for (rds_file in rds_files) {
    specimen_id <- tools::file_path_sans_ext(basename(rds_file))
    message("  Downsampling ", library_batch, " / ", specimen_id)

    obj <- readRDS(rds_file)
    obj$specimen_id <- specimen_id

    # Complete metadata BEFORE downsampling — the counts/percentages for
    # the diagnoses must be exact, not based on a sample
    key <- paste0(library_batch, "_", specimen_id)
    meta_list[[key]] <- obj@meta.data %>%
      as_tibble() %>%
      select(individualID, libraryBatch, phenotype_label, cell_type, annotation_source)

    if (ncol(obj) > max_cells_per_specimen) {
      set.seed(42)
      cells_keep <- sample(colnames(obj), max_cells_per_specimen)
      obj <- subset(obj, cells = cells_keep)
    }

    obj <- DietSeurat(obj, layers = "counts", dimreducs = NULL, graphs = NULL)
    obj_list[[key]] <- obj
  }
}

message("Total specimens: ", length(obj_list))

## diagnostics of the cell anotation, before and after CellTypist rescue
all_meta <- bind_rows(meta_list)
n_total <- nrow(all_meta)
message("Total of cells (with out downsample): ", n_total)

write_csv(all_meta, file.path(diagnostics_dir, "annotation_diagnostics_metadata.csv"))

source_summary <- all_meta %>%
  count(annotation_source) %>%
  mutate(pct_del_total = round(100 * n / n_total, 2))

get_n <- function(df, val) {
  x <- df$n[df$annotation_source == val]
  if (length(x) == 0) 0 else x
}

n_original   <- get_n(source_summary, "original")
n_celltypist <- get_n(source_summary, "celltypist")
n_unassigned <- get_n(source_summary, "unassigned")
n_na_before_celltypist <- n_celltypist + n_unassigned

resumen_texto <- tibble(
  metrica = c(
    "Original anotation (before CellTypist)",
    "NA of the reference CSV (before CellTypist)",
    "Annotated by CellTypist (over the NA)",
    "Annotated by CellTypist (over the total of cells)",
    "NA after CellTypist (unassigned)"
  ),
  n = c(n_original, n_na_before_celltypist, n_celltypist, n_celltypist, n_unassigned),
  pct = c(
    round(100 * n_original / n_total, 2),
    round(100 * n_na_before_celltypist / n_total, 2),
    round(100 * n_celltypist / n_na_before_celltypist, 2),
    round(100 * n_celltypist / n_total, 2),
    round(100 * n_unassigned / n_total, 2)
  )
)
write_csv(resumen_texto, file.path(diagnostics_dir, "01_summary_origin_annotation.csv.csv"))

p1 <- source_summary %>%
  ggplot(aes(x = annotation_source, y = n, fill = annotation_source)) +
  geom_col() +
  geom_text(aes(label = paste0(n, " (", pct_del_total, "%)")), vjust = -0.3) +
  labs(x = "annotation_source", y = "Number of cells", title = "Origin of Cell Annotation") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(diagnostics_dir, "01_origen_anotacion.png"), p1, width = 6, height = 5, dpi = 150)

celltype_summary <- all_meta %>%
  filter(!is.na(cell_type)) %>%
  count(cell_type) %>%
  mutate(pct = round(100 * n / sum(n), 2)) %>%
  arrange(desc(n))
write_csv(celltype_summary, file.path(diagnostics_dir, "02_cell_type_summary.csv"))

p2 <- celltype_summary %>%
  ggplot(aes(x = reorder(cell_type, n), y = n)) +
  geom_col() +
  geom_text(aes(label = paste0(pct, "%")), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(x = "cell_type", y = "Number of cells", title = "Distribution of cell_type (all individuals)") +
  theme_minimal()
ggsave(file.path(diagnostics_dir, "02_cell_type_global.png"), p2, width = 8, height = 6, dpi = 150)

celltype_by_specimen <- all_meta %>%
  filter(!is.na(cell_type)) %>%
  count(individualID, phenotype_label, cell_type) %>%
  group_by(individualID) %>%
  mutate(pct_specimen = round(100 * n / sum(n), 2)) %>%
  ungroup()
write_csv(celltype_by_specimen, file.path(diagnostics_dir, "03_cell_type_por_specimen.csv"))

p3 <- celltype_by_specimen %>%
  ggplot(aes(x = individualID, y = pct_specimen, fill = cell_type)) +
  geom_col() +
  facet_grid(~ phenotype_label, scales = "free_x", space = "free_x") +
  coord_flip() +
  labs(x = "individualID", y = "% of cells",
       title = "Composition of cell_type by individual, grouped by phenotype") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(diagnostics_dir, "03_cell_type_por_specimen_y_phenotype.png"), p3, width = 10, height = 8, dpi = 150)

celltype_by_phenotype <- all_meta %>%
  filter(!is.na(cell_type), !is.na(phenotype_label)) %>%
  count(phenotype_label, cell_type) %>%
  group_by(phenotype_label) %>%
  mutate(pct_phenotype = round(100 * n / sum(n), 2)) %>%
  ungroup()
write_csv(celltype_by_phenotype, file.path(diagnostics_dir, "04_cell_type_por_phenotype.csv"))

p4 <- celltype_by_phenotype %>%
  ggplot(aes(x = phenotype_label, y = pct_phenotype, fill = cell_type)) +
  geom_col(position = "stack") +
  labs(x = "phenotype", y = "% of cells", fill = "cell_type",
       title = "Composition of cell_type by phenotype (aggregated)") +
  theme_minimal()
ggsave(file.path(diagnostics_dir, "04_cell_type_por_phenotype.png"), p4, width = 8, height = 6, dpi = 150)

message("Diagnostics of cell annotation saved in: ", diagnostics_dir)

## Unambiguous UMAP of all specimens together, colored by cell_type, phenotype and library_batch

gc()
merged <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
merged <- JoinLayers(merged)
rm(obj_list)
gc()
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:30, verbose = FALSE)

p_celltype <- DimPlot(merged, group.by = "cell_type", label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("All the specimens - cell_type")

p_source <- DimPlot(merged, group.by = "annotation_source") +
  ggtitle("All the specimens - annotation_source")

p_batch <- DimPlot(merged, group.by = "libraryBatch") +
  ggtitle("All the specimens - library_batch (checking batch effect)")

ggsave(snakemake@output[["umap_celltype"]], p_celltype, width = 8, height = 7, dpi = 150)
ggsave(snakemake@output[["umap_source"]], p_source, width = 8, height = 7, dpi = 150)
ggsave(snakemake@output[["umap_batch"]], p_batch, width = 8, height = 7, dpi = 150)

message("Combined UMAPs saved in results/annotation_plots/")
sink(type = "message"); sink()