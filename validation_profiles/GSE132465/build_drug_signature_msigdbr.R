#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else "msigdb_c2_cgp_directional_signatures.csv.gz"
if (!requireNamespace("msigdbr", quietly = TRUE) || !requireNamespace("data.table", quietly = TRUE)) stop("msigdbr and data.table are required")
x <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
name_col <- if ("gs_name" %in% names(x)) "gs_name" else stop("msigdbr result lacks gs_name")
gene_col <- if ("gene_symbol" %in% names(x)) "gene_symbol" else stop("msigdbr result lacks gene_symbol")
name <- as.character(x[[name_col]])
direction <- ifelse(grepl("_UP$", name), 1L, ifelse(grepl("_DN$", name), -1L, NA_integer_))
base <- sub("_(UP|DN)$", "", name)
pairs <- intersect(unique(base[direction == 1L]), unique(base[direction == -1L]))
keep <- !is.na(direction) & base %in% pairs
out_tab <- data.frame(drug = base[keep], gene = as.character(x[[gene_col]][keep]), effect = direction[keep],
                      source = "MSigDB C2:CGP directional pair", original_set = name[keep], stringsAsFactors = FALSE)
out_tab <- unique(out_tab[nzchar(out_tab$gene), ])
data.table::fwrite(out_tab, out)
cat("SIGNATURES=", length(unique(out_tab$drug)), "\nROWS=", nrow(out_tab), "\nOUTPUT=", normalizePath(out, winslash = "/", mustWork = FALSE), "\n", sep = "")
