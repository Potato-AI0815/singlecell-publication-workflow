#!/usr/bin/env Rscript
# Conservative preparation template reconstructed from the validation notes.
# It is not the original byte-for-byte locally validated script.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Usage: Rscript 00_prepare_GSE132465_template.R raw_umi.txt.gz annotation.txt.gz output_dir [SMC01,SMC02,...]")
raw_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
ann_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[[3]], winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
patients <- if (length(args) >= 4L) strsplit(args[[4]], ",", fixed = TRUE)[[1]] else sprintf("SMC%02d", 1:6)
if (!requireNamespace("data.table", quietly = TRUE) || !requireNamespace("Seurat", quietly = TRUE)) stop("data.table and Seurat are required")
ann <- data.table::fread(ann_path, data.table = FALSE)
# Inspect and edit these aliases when the upstream annotation file changes.
find_col <- function(candidates) {
  hit <- candidates[candidates %in% names(ann)]
  if (!length(hit)) stop("Missing annotation column; tried: ", paste(candidates, collapse = ", "))
  hit[[1]]
}
cell_col <- find_col(c("Index", "cell", "cell_id", "barcode"))
patient_col <- find_col(c("Patient", "patient", "patient_id"))
condition_col <- find_col(c("Class", "condition", "Condition", "sample_type"))
author_type_col <- find_col(c("Cell_type", "cell_type", "CellType", "Cluster"))
ann <- ann[as.character(ann[[patient_col]]) %in% patients, , drop = FALSE]
keep_cells <- as.character(ann[[cell_col]])
# The matrix is assumed to have gene IDs in the first column and cells in subsequent columns.
hdr <- names(data.table::fread(raw_path, nrows = 0L, data.table = FALSE, check.names = FALSE))
gene_col <- hdr[[1]]
selected <- c(gene_col, intersect(keep_cells, hdr))
if (length(selected) < 2L) stop("No selected annotation barcodes matched matrix columns.")
mat_df <- data.table::fread(raw_path, select = selected, data.table = FALSE, check.names = FALSE)
genes <- make.unique(as.character(mat_df[[1]]))
mat <- as.matrix(mat_df[, -1, drop = FALSE]); storage.mode(mat) <- "numeric"; rownames(mat) <- genes
obj <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(mat, sparse = TRUE), project = "GSE132465_6pair", min.cells = 0, min.features = 0)
mm <- ann[match(colnames(obj), as.character(ann[[cell_col]])), , drop = FALSE]
obj$sample_id <- paste0(as.character(mm[[patient_col]]), "-", substr(as.character(mm[[condition_col]]), 1, 1))
obj$patient_id <- as.character(mm[[patient_col]])
obj$condition <- as.character(mm[[condition_col]])
obj$author_cell_type <- as.character(mm[[author_type_col]])
obj$tumor_epithelial_flag <- obj$condition == "Tumor" & grepl("epithelial", obj$author_cell_type, ignore.case = TRUE)
saveRDS(obj, file.path(out_dir, "GSE132465_6pair_raw_counts.rds"), compress = FALSE)
meta <- unique(obj[[]][, c("sample_id", "patient_id", "condition"), drop = FALSE])
data.table::fwrite(meta, file.path(out_dir, "sample_metadata.csv"))
cat("OUTPUT=", normalizePath(out_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
