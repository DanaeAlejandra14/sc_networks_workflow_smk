log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))

csv_files <- unlist(snakemake@input)
message("aggregate", length(csv_files), " files in summary...")

df <- map_dfr(csv_files, read_csv, show_col_types = FALSE)

write_csv(df, snakemake@output[[1]])
message("Summary save in: ", snakemake@output[[1]])
message("Total de rows: ", nrow(df))

sink(type = "message"); sink()