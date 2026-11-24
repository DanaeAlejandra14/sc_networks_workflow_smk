log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# One gene set per cell type, shared by every phenotype and gamma, so that the
# networks are comparable. Computed on single cells (not metacells), with the
# subsampled cells of all phenotypes pooled:
#   1. keep genes detected (count > 0) in >= min_frac_cells of the cells
#   2. drop genes matching exclude_patterns (MT-, ribosomal)
#   3. keep the n_variable most variable of those (Seurat vst)
in_rds           <- unlist(snakemake@input[["subsampled"]])
out_genes        <- snakemake@output[["genes"]]
out_report       <- snakemake@output[["report"]]
min_frac_cells   <- snakemake@params[["min_frac_cells"]]
n_variable       <- snakemake@params[["n_variable"]]
exclude_patterns <- unlist(snakemake@params[["exclude_patterns"]])

message("Cell type: ", snakemake@wildcards[["cell_type"]])
message("min_frac_cells: ", min_frac_cells, " | n_variable: ", n_variable,
        " | exclude: ", paste(exclude_patterns, collapse = ", "))

# Raw counts of every phenotype, one column block per phenotype
counts <- lapply(in_rds, function(f) {
  obj <- readRDS(f)
  # Merged objects keep one counts layer per specimen; join them into one matrix
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  m <- LayerData(obj, assay = "RNA", layer = "counts")
  message("  ", f, ": ", nrow(m), " genes x ", ncol(m), " cells")
  m
})
stopifnot("Gene order differs between phenotypes" =
            all(vapply(counts, function(m) identical(rownames(m), rownames(counts[[1]])), logical(1))))
counts <- do.call(cbind, counts)
message("Pooled: ", nrow(counts), " genes x ", ncol(counts), " cells")

# 1. Detection
frac_cells <- Matrix::rowMeans(counts > 0)
detected <- frac_cells >= min_frac_cells
message("Detected in >= ", min_frac_cells, " of cells: ", sum(detected))

# 2. Exclusion
excluded <- Reduce(`|`, lapply(exclude_patterns, function(p) grepl(p, rownames(counts))))
message("Excluded by pattern (among detected): ", sum(detected & excluded))
candidates <- rownames(counts)[detected & !excluded]
message("Candidates for variable gene selection: ", length(candidates))
stopifnot("Fewer candidate genes than n_variable" = length(candidates) >= n_variable)

# 3. Most variable genes (vst on the raw counts of the candidates)
obj <- CreateSeuratObject(counts = counts[candidates, ])
rm(counts); invisible(gc())
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_variable)
genes <- VariableFeatures(obj)
message("Selected variable genes: ", length(genes))

hvf <- HVFInfo(obj, method = "vst")
col_of <- function(pattern) grep(pattern, colnames(hvf), value = TRUE)[1]
report <- data.frame(
  gene = candidates,
  frac_cells = round(frac_cells[candidates], 4),
  mean = hvf[candidates, col_of("mean$")],
  variance_standardized = hvf[candidates, col_of("variance.standardized$")],
  selected = candidates %in% genes
)
report <- report[order(-report$variance_standardized), ]

dir.create(dirname(out_genes), recursive = TRUE, showWarnings = FALSE)
writeLines(genes, out_genes)
write.table(report, out_report, sep = "\t", quote = FALSE, row.names = FALSE)
message("Written: ", out_genes, " and ", out_report)

sink(type = "message"); sink()
