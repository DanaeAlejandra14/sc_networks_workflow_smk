log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)

input_files <- unlist(snakemake@input)
output_path <- snakemake@output[[1]]

message("Aggregating ", length(input_files), " QC summary files")

qc_all <- map_dfr(input_files, read_csv, show_col_types = FALSE)

write_csv(qc_all, output_path)

message("Consolidated QC summary written to: ", output_path)

sink(type = "message"); sink()