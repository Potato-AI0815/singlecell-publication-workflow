#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript 03_validate_GSE132465_outputs.R result_root [core|cnv|cellchat|nmf|hdwgcna|all]")
root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
scope <- if (length(args) >= 2L) tolower(args[[2]]) else "all"
checks <- data.frame(check = character(), status = character(), evidence = character(), stringsAsFactors = FALSE)
add <- function(id, ok, evidence) checks <<- rbind(checks, data.frame(check=id, status=if(ok)"PASS" else "FAIL", evidence=as.character(evidence), stringsAsFactors=FALSE))
need <- function(rel) { p <- file.path(root, rel); add(paste0("file:", rel), file.exists(p) && file.info(p)$size > 0, p); p }
if (scope %in% c("core","all")) {
  need("01_objects/06_annotated.rds"); need("02_qc/sample_qc_summary.csv"); need("04_annotation/annotation_evidence.csv")
  need("14_qa/FINAL_STATUS.txt"); need("15_manifests/figure_export_manifest.csv")
}
if (scope %in% c("cnv","all")) { need("17_advanced_modules/cnv_cell_scores.csv"); need("17_advanced_modules/cnv_parameters.yml") }
if (scope %in% c("cellchat","all")) { need("17_advanced_modules/cellchat_by_sample.rds"); need("17_advanced_modules/cellchat_interactions_by_sample.csv") }
if (scope %in% c("nmf","all")) { need("17_advanced_modules/nmf_meta_program_genes.csv"); need("17_advanced_modules/nmf_scored_object.rds") }
if (scope %in% c("hdwgcna","all")) { need("17_advanced_modules/hdwgcna_parameters.yml") }
utils::write.csv(checks, file.path(root, paste0("14_qa/GSE132465_validation_", scope, ".csv")), row.names = FALSE)
print(checks)
quit(status = if (any(checks$status == "FAIL")) 1 else 0)
