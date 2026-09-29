log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

batch_name <- snakemake@wildcards[["library_batch"]]
annotated_dir <- snakemake@input[["annotated_dir"]]
predictions_path <- snakemake@input[["predictions"]]
final_dir <- snakemake@output[["final_dir"]]

dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

message("Merging CellTypist predictions back for batch: ", batch_name)

predictions <- readr::read_csv(predictions_path, show_col_types = FALSE) %>%
  mutate(across(c(cell, celltypist_label, celltype_simple), as.character))

lab_map  <- setNames(predictions$celltypist_label, predictions$cell)
sim_map  <- setNames(predictions$celltype_simple, predictions$cell)
prob_map <- setNames(predictions$max_prob, predictions$cell)
mar_map  <- setNames(predictions$margin_top1_top2, predictions$cell)

rds_files <- list.files(annotated_dir, pattern = "\\.rds$", full.names = TRUE)
message("Found ", length(rds_files), " specimen(s) in this batch")

for (rds_file in rds_files) {

  specimen_id <- tools::file_path_sans_ext(basename(rds_file))
  message("Processing specimen: ", specimen_id)

  obj <- readRDS(rds_file)

  obj$celltypist_label   <- NA_character_
  obj$celltypist_simple  <- NA_character_
  obj$celltypist_max_prob <- NA_real_
  obj$celltypist_margin  <- NA_real_

  # Reconstruimos el MISMO prefijo que se usó en export_unassigned.R
  # para poder cruzar contra los IDs que trae el CSV de predicciones
  prefixed_cells <- paste0(specimen_id, "_", colnames(obj))
  in_predictions <- prefixed_cells %in% predictions$cell

  obj$celltypist_label[in_predictions]    <- lab_map[prefixed_cells[in_predictions]]
  obj$celltypist_simple[in_predictions]   <- sim_map[prefixed_cells[in_predictions]]
  obj$celltypist_max_prob[in_predictions] <- prob_map[prefixed_cells[in_predictions]]
  obj$celltypist_margin[in_predictions]   <- mar_map[prefixed_cells[in_predictions]]

  na_cells <- is.na(obj$cell_type)
  fill_cells <- na_cells & in_predictions &
    !is.na(obj$celltypist_simple) & obj$celltypist_simple != "Unassigned"

  obj$cell_type[fill_cells] <- obj$celltypist_simple[fill_cells]

  obj$annotation_source <- "original"
  obj$annotation_source[fill_cells] <- "celltypist"
  obj$annotation_source[is.na(obj$cell_type)] <- "unassigned"

  message("  original: ", sum(obj$annotation_source == "original"),
          " | rescued: ", sum(obj$annotation_source == "celltypist"),
          " | still unassigned: ", sum(obj$annotation_source == "unassigned"))

  saveRDS(obj, file.path(final_dir, paste0(specimen_id, ".rds")), compress = "gzip")
}

message("Merge complete for batch: ", batch_name)
sink(type = "message"); sink()