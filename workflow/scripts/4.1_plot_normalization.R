# Write all stdout and messages/warnings/errors to the Snakemake log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)
library(Seurat)

batch_name <- snakemake@wildcards[["library_batch"]]
norm_dir   <- snakemake@input[["norm_dir"]]
out_png    <- snakemake@output[[1]]

dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

message("Normalization plot for batch: ", batch_name)

rds_files <- list.files(norm_dir, pattern = "\\.rds$", full.names = TRUE)
message("Found ", length(rds_files), " specimen(s) in this batch")

# Pick one specimen at random from this batch
set.seed(123)  # for reproducibility (always the same "random" choice)
rds_file <- sample(rds_files, 1)
specimen_id <- tools::file_path_sans_ext(basename(rds_file))

message("Randomly selected specimen: ", specimen_id)

obj <- readRDS(rds_file)

counts_mat <- GetAssayData(obj, assay = "RNA", slot = "counts")
data_mat   <- GetAssayData(obj, assay = "RNA", slot = "data")

# Non-zero values stored in the sparse matrices
counts_values <- counts_mat@x
data_values   <- data_mat@x

message("Counts values (non-zero): ", length(counts_values))
message("Data values (non-zero): ", length(data_values))

df_plot <- bind_rows(
  tibble(valor = counts_values, etapa = "Before (counts, raw)"),
  tibble(valor = data_values, etapa = "After (data, log-normalized)")
) %>%
  mutate(etapa = factor(etapa, levels = c("Before (counts, raw)", "After (data, log-normalized)")))

p <- df_plot %>%
  ggplot(aes(x = valor, fill = etapa)) +
  geom_histogram(bins = 80, alpha = 0.8, colour = "black") +
  facet_wrap(~etapa, scales = "free", ncol = 2) +
  labs(x = "expression value", y = "frequency",
       title = paste0("Before vs after LogNormalize\n",
                       "Batch: ", batch_name, " | Specimen: ", specimen_id)) +
  theme_minimal() +
  theme(
    legend.position   = "none",
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    strip.background  = element_rect(fill = "grey90", colour = NA),
    strip.text        = element_text(face = "bold")
  )

ggsave(out_png, p, width = 11, height = 6, dpi = 150, bg = "white")

message("Plot saved to: ", out_png)
sink(type = "message"); sink()