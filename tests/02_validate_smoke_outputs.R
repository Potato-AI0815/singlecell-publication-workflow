#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript tests/02_validate_smoke_outputs.R /path/to/result_root")
root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
req <- c(
  "01_objects/06_annotated.rds",
  "14_qa/figure_generation_status.csv",
  "14_qa/final_QA_checks.csv",
  "14_qa/FINAL_STATUS.txt",
  "15_manifests/figure_export_manifest.csv"
)
missing <- req[!file.exists(file.path(root, req))]
if (length(missing)) stop("Missing smoke outputs: ", paste(missing, collapse = ", "))
status_path <- file.path(root, "14_qa/FINAL_STATUS.txt")
status <- trimws(paste(readLines(status_path, warn = FALSE), collapse = "\n"))
if (!status %in% c("PASS", "PASS_WITH_WARNINGS")) {
  stop("FINAL_STATUS must be PASS or PASS_WITH_WARNINGS; got: ", status)
}
m <- utils::read.csv(file.path(root, "15_manifests/figure_export_manifest.csv"), check.names = FALSE)
needed <- c("figure_id", "pdf", "tiff", "png", "proof_png", "plot_rds", "parameters", "visual_qa", "visual_qa_status")
if (!all(needed %in% names(m))) stop("Manifest missing columns: ", paste(setdiff(needed, names(m)), collapse = ", "))
if (!nrow(m)) stop("No figure was exported")
for (nm in c("pdf", "tiff", "png", "proof_png", "plot_rds", "parameters", "visual_qa")) {
  if (!all(file.exists(m[[nm]]))) stop("Missing files in manifest column: ", nm)
}
if (any(m$visual_qa_status == "FAIL")) stop("At least one figure failed visual QA")
cat("SMOKE_OUTPUT_VALIDATION PASS\n")
