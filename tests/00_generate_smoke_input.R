#!/usr/bin/env Rscript
# Creates a small synthetic four-sample dataset that exercises QC, integration,
# provisional annotation, composition, marker figures, continuous embeddings and split violins.
args <- commandArgs(trailingOnly = TRUE)
out_root <- if (length(args)) args[[1]] else file.path("tests", "tmp_smoke_project")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
set.seed(20260804)

genes <- unique(c(
  paste0("G", 1:1200),
  "CD3D", "CD3E", "TRAC", "IL7R", "LTB", "NKG7", "GNLY", "KLRD1",
  "MS4A1", "CD79A", "CD37", "CD74", "LST1", "S100A8", "S100A9", "CTSS",
  "FCGR3A", "LYZ", "COL1A1", "DCN", "EPCAM", "KRT8", "KRT18",
  "MKI67", "TOP2A", "MT-CO1", "MT-CO2"
))

make_sample <- function(sample_id, n = 320L, treatment = FALSE) {
  cell_type <- sample(c("T", "NK", "B", "Myeloid"), n, replace = TRUE, prob = c(0.38, 0.17, 0.18, 0.27))
  m <- matrix(rpois(length(genes) * n, lambda = 0.12), nrow = length(genes), dimnames = list(genes, paste0(sample_id, "_C", seq_len(n))))
  add_program <- function(idx, gs, lambda) {
    if (!length(idx)) return()
    m[gs, idx] <<- m[gs, idx, drop = FALSE] + matrix(rpois(length(gs) * length(idx), lambda), nrow = length(gs))
  }
  add_program(which(cell_type == "T"), c("CD3D", "CD3E", "TRAC", "IL7R", "LTB"), if (treatment) 3.4 else 3.0)
  add_program(which(cell_type == "NK"), c("NKG7", "GNLY", "KLRD1"), 4.0)
  add_program(which(cell_type == "B"), c("MS4A1", "CD79A", "CD37", "CD74"), 3.5)
  add_program(which(cell_type == "Myeloid"), c("LST1", "S100A8", "S100A9", "CTSS", "FCGR3A", "LYZ"), if (treatment) 3.6 else 3.0)
  Matrix::Matrix(m, sparse = TRUE)
}

sample_ids <- c("C1", "C2", "T1", "T2")
conditions <- c("Control", "Control", "Treatment", "Treatment")
counts <- setNames(lapply(seq_along(sample_ids), function(i) make_sample(sample_ids[i], treatment = conditions[i] == "Treatment")), sample_ids)
obj <- Seurat::CreateSeuratObject(counts = counts, min.features = 0)
obj$sample_id <- sub("_C.*$", "", colnames(obj))
obj$condition <- conditions[match(obj$sample_id, sample_ids)]
obj$patient_id <- paste0("P", match(obj$sample_id, sample_ids))
obj$batch <- ifelse(obj$sample_id %in% c("C1", "T1"), "B1", "B2")
obj$tissue <- "PBMC"
input_path <- file.path(out_root, "synthetic_multisample.rds")
saveRDS(obj, input_path, compress = FALSE)

cat("Synthetic input written to:", normalizePath(input_path, winslash = "/", mustWork = TRUE), "\n")
cat("Use templates/config.full.yml with input.type=seurat_rds, input.path set to this file, and analysis.mode=fast or full.\n")
