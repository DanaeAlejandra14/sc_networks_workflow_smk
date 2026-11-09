#!/usr/bin/env Rscript
# trace_library_map.R
#
# Read-only check: which pipeline libraries (libraryBatch) have NO batch with the
# same name in the authors' atlas, and which atlas batch (MAP<number>) corresponds
# to each of them, matching by individualID.
#
# It does not change anything in the pipeline. It only prints a report and writes
# a CSV table (library_batch -> atlas_batch) that can become the alias table later.
#
# Usage (from the project root):
#   Rscript trace_library_map.R <atlas_csv> <demux_csv> [out_csv]
# Example:
#   Rscript trace_library_map.R resources/cell-annotation.full-atlas.csv \
#       resources/ROSMAPsnRNAseqdemultiplexedIDmapping.csv \
#       results/qc/library_alias_check.csv

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript trace_library_map.R <atlas_csv> <demux_csv> [out_csv]")
}
atlas_csv <- args[1]
demux_csv <- args[2]
out_csv   <- if (length(args) >= 3) args[3] else "results/qc/library_alias_check.csv"

# Atlas batches that need an alias are named like MAP36492755
map_pattern <- "^MAP[0-9]+$"
n_probe <- 2000  # rows used to detect which column is which

# ---------------------------------------------------------------- 1. Atlas ----
message("Reading atlas: ", atlas_csv)
atlas <- readr::read_csv(atlas_csv, show_col_types = FALSE,
                         col_types = cols(.default = col_character()))
message("  rows: ", nrow(atlas), " | columns: ", paste(names(atlas), collapse = ", "))

# Detect columns by their content, so the script does not depend on header names:
#  - cell column: values like <batch>_<barcode>-1
#  - batch column: its value + "_" is the start of the cell column
#  - individual column: values like R1234567
probe <- head(atlas, n_probe)

cell_col <- names(probe)[sapply(probe, function(x)
  all(grepl("_[ACGT]+-[0-9]+$", x[!is.na(x)])) && any(!is.na(x)))]
if (length(cell_col) != 1) {
  stop("Could not identify the cell column. Candidates: ",
       paste(cell_col, collapse = ", "), "\nColumns: ", paste(names(atlas), collapse = ", "))
}

batch_col <- names(probe)[sapply(probe, function(x)
  all(!is.na(x)) && all(startsWith(probe[[cell_col]], paste0(x, "_"))))]
batch_col <- setdiff(batch_col, cell_col)
if (length(batch_col) != 1) {
  stop("Could not identify the batch column. Candidates: ",
       paste(batch_col, collapse = ", "), "\nColumns: ", paste(names(atlas), collapse = ", "))
}

ind_col <- names(probe)[sapply(probe, function(x)
  all(grepl("^R[0-9]+$", x[!is.na(x)])) && any(!is.na(x)))]
if (length(ind_col) != 1) {
  stop("Could not identify the individualID column. Candidates: ",
       paste(ind_col, collapse = ", "), "\nColumns: ", paste(names(atlas), collapse = ", "))
}
message("  detected -> cell: ", cell_col, " | batch: ", batch_col, " | individual: ", ind_col)

atlas_batches <- atlas %>%
  transmute(batch = .data[[batch_col]], individualID = .data[[ind_col]]) %>%
  group_by(batch) %>%
  summarise(n_cells = n(),
            n_individuals = n_distinct(individualID),
            individualID = first(individualID),
            .groups = "drop")

map_batches <- atlas_batches %>% filter(str_detect(batch, map_pattern))

message("\nAtlas: ", nrow(atlas_batches), " batches in total, ",
        nrow(map_batches), " named MAP<number>")
bad_map <- map_batches %>% filter(n_individuals != 1)
message("MAP batches with more than one individual (expected 0): ", nrow(bad_map))

# ------------------------------------------------------------- 2. Demux file --
message("\nReading demultiplexing file: ", demux_csv)
demux <- readr::read_csv(demux_csv, show_col_types = FALSE)
message("  columns: ", paste(names(demux), collapse = ", "))

find_col <- function(df, wanted) {
  hit <- names(df)[tolower(names(df)) == tolower(wanted)]
  if (length(hit) != 1) {
    stop("Column '", wanted, "' not found in the demux file. Columns: ",
         paste(names(df), collapse = ", "))
  }
  hit
}
lib_col  <- find_col(demux, "libraryBatch")
dind_col <- find_col(demux, "individualID")

lib_ind <- demux %>%
  transmute(library_batch = .data[[lib_col]], individualID = .data[[dind_col]]) %>%
  filter(!is.na(library_batch), !is.na(individualID)) %>%
  distinct()

lib_summary <- lib_ind %>%
  group_by(library_batch) %>%
  summarise(n_ind_lib = n_distinct(individualID),
            individualID = first(individualID),
            .groups = "drop")

direct  <- lib_summary %>% filter(library_batch %in% atlas_batches$batch)
targets <- lib_summary %>% filter(!library_batch %in% atlas_batches$batch)
message("Libraries in demux: ", nrow(lib_summary),
        " | found directly in the atlas: ", nrow(direct),
        " | NOT found (need an alias): ", nrow(targets))

# ------------------------------------------------------------------ 3. Match ---
cand <- targets %>%
  left_join(map_batches %>% select(atlas_batch = batch, individualID, n_cells),
            by = "individualID")

lib_counts <- cand %>% filter(!is.na(atlas_batch)) %>% count(library_batch, name = "n_map_for_lib")
map_counts <- cand %>% filter(!is.na(atlas_batch)) %>% count(atlas_batch, name = "n_lib_for_map")

res <- cand %>%
  left_join(lib_counts, by = "library_batch") %>%
  left_join(map_counts, by = "atlas_batch") %>%
  mutate(n_map_for_lib = replace_na(n_map_for_lib, 0L),
         id_in_name = str_extract(library_batch, "R[0-9]+"),
         name_matches_id = if_else(is.na(id_in_name), NA, id_in_name == individualID),
         status = case_when(
           n_map_for_lib == 0                           ~ "no_atlas_match",
           n_ind_lib != 1                               ~ "library_has_several_individuals",
           n_map_for_lib > 1 | n_lib_for_map > 1        ~ "ambiguous",
           isFALSE(name_matches_id)                     ~ "name_mismatch",
           TRUE                                         ~ "ok"
         )) %>%
  select(library_batch, atlas_batch, individualID, n_cells_atlas = n_cells,
         id_in_name, name_matches_id, n_ind_lib, n_map_for_lib, n_lib_for_map, status) %>%
  arrange(status, library_batch)

# ------------------------------------------------------------------ 4. Report --
message("\n=== Library -> atlas batch table ===")
print(res %>% select(library_batch, atlas_batch, individualID, n_cells_atlas, status),
      n = Inf, width = Inf)

message("\n=== Status summary ===")
print(count(res, status))

orphans <- map_batches %>% filter(!batch %in% res$atlas_batch[res$status == "ok"])
message("\nMAP batches in the atlas not matched to any library of this demux file: ",
        nrow(orphans))

message("\n=== Checks ===")
chk <- function(label, ok) message(sprintf("  [%s] %s", if (isTRUE(ok)) "PASS" else "FAIL", label))
chk("Every MAP batch has exactly one individual", nrow(bad_map) == 0)
chk("Every library without a direct match found a MAP batch", all(res$n_map_for_lib >= 1))
chk("Matches are one-to-one (no library/MAP repeated)",
    all(res$n_map_for_lib <= 1, na.rm = TRUE) && all(res$n_lib_for_map <= 1, na.rm = TRUE))
chk("Individual in the library name equals individualID", !any(res$name_matches_id %in% FALSE))
chk("All rows have status 'ok'", all(res$status == "ok"))

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(res, out_csv)
message("\nTable written to: ", out_csv)