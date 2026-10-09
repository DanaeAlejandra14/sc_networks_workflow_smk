log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))

csv_files <- unlist(snakemake@input)
message("Juntando ", length(csv_files), " archivos de resumen...")

df <- map_dfr(csv_files, read_csv, show_col_types = FALSE)

write_csv(df, snakemake@output[[1]])
message("Resumen maestro guardado en: ", snakemake@output[[1]])
message("Total de filas: ", nrow(df))

sink(type = "message"); sink()