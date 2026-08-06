# Gated Slingshot + optional tradeSeq module.
# Root, analysis subset and topology are explicit. The module never infers a biological direction from the visual embedding alone.

write_tradeseq_result <- function(x, path) {
  if (is.null(x)) return(invisible(FALSE))
  z <- as.data.frame(x, check.names = FALSE)
  z$gene <- rownames(z)
  if ("pvalue" %in% names(z)) z$FDR <- stats::p.adjust(z$pvalue, method = "BH")
  write_csv_safe(z, path)
  invisible(TRUE)
}

run_tradeseq_layer <- function(ctx, sub, pseudotime, cell_weights, cfg) {
  if (!as_flag(cfg$run_tradeSeq, TRUE)) return(status_table("NOT_RUN", "tradeSeq disabled"))
  if (!requireNamespace("tradeSeq", quietly = TRUE)) return(status_table("NOT_RUN", "tradeSeq unavailable"))
  sub <- join_layers_safe(sub, "RNA")
  counts <- get_layer_safe(sub, "RNA", "counts")
  counts <- counts[, rownames(pseudotime), drop = FALSE]
  expressed <- Matrix::rowSums(counts > 0) >= (cfg$tradeSeq_minimum_expressing_cells %||% 20L)
  genes <- rownames(counts)[expressed]
  maximum_genes <- cfg$tradeSeq_maximum_genes %||% 12000L
  if (length(genes) > maximum_genes) {
    mu <- Matrix::rowMeans(counts[genes, , drop = FALSE])
    genes <- names(sort(mu, decreasing = TRUE))[seq_len(maximum_genes)]
  }
  if (length(genes) < (cfg$tradeSeq_minimum_genes %||% 200L)) return(status_table("NOT_EVALUABLE", "too few expressed genes for tradeSeq"))

  knots <- as.integer(unlist(cfg$tradeSeq_knots %||% 3:7))
  selected_knots <- as.integer(cfg$tradeSeq_selected_knots %||% 6L)
  if (as_flag(cfg$tradeSeq_evaluate_knots, TRUE)) {
    set.seed(ctx$cfg$runtime$seed %||% 20260804)
    eval_res <- tryCatch(
      tradeSeq::evaluateK(
        counts = counts[genes, , drop = FALSE], pseudotime = pseudotime,
        cellWeights = cell_weights, k = knots,
        nGenes = min(cfg$tradeSeq_evaluate_n_genes %||% 500L, length(genes)),
        plot = FALSE, verbose = FALSE,
        parallel = as_flag(cfg$tradeSeq_parallel, FALSE)
      ), error = function(e) e
    )
    if (!inherits(eval_res, "error")) {
      saveRDS(eval_res, file.path(ctx$dirs$advanced, "tradeSeq_knot_evaluation.rds"))
      aic <- if (is.list(eval_res) && length(eval_res) >= 2L) eval_res[[2L]] else if (is.matrix(eval_res)) eval_res else NULL
      if (!is.null(aic)) {
        aic <- as.matrix(aic)
        knot_names <- suppressWarnings(as.integer(gsub("[^0-9]", "", colnames(aic))))
        if (all(is.na(knot_names))) knot_names <- knots[seq_len(min(length(knots), ncol(aic)))]
        med <- apply(aic, 2, stats::median, na.rm = TRUE)
        diag <- data.frame(knots = knot_names[seq_along(med)], median_AIC = med, stringsAsFactors = FALSE)
        write_csv_safe(diag, file.path(ctx$dirs$advanced, "tradeSeq_knot_diagnostics.csv"))
        if (any(is.finite(med))) selected_knots <- diag$knots[which.min(diag$median_AIC)]
      }
    } else append_issue(ctx, "TRAJECTORY", "WARNING", paste("tradeSeq evaluateK failed:", conditionMessage(eval_res)))
  }
  if (!is.finite(selected_knots) || selected_knots < 3L) selected_knots <- 6L

  U <- NULL
  fixed <- unlist(cfg$tradeSeq_fixed_effects %||% list())
  fixed <- fixed[fixed %in% colnames(sub[[]])]
  if (length(fixed)) {
    mm <- stats::model.matrix(stats::as.formula(paste("~ 0 +", paste(fixed, collapse = " + "))), data = sub[[]][rownames(pseudotime), , drop = FALSE])
    if (ncol(mm) && qr(mm)$rank == ncol(mm)) U <- mm else append_issue(ctx, "TRAJECTORY", "WARNING", "tradeSeq fixed-effect matrix was rank deficient and was omitted")
  }

  fit <- tryCatch(
    tradeSeq::fitGAM(
      counts = counts[genes, , drop = FALSE], pseudotime = pseudotime,
      cellWeights = cell_weights, U = U, nknots = selected_knots,
      verbose = FALSE, parallel = as_flag(cfg$tradeSeq_parallel, FALSE), sce = TRUE
    ), error = function(e) e
  )
  if (inherits(fit, "error")) return(status_table("FAILED", paste("tradeSeq fitGAM failed:", conditionMessage(fit))))
  saveRDS(fit, file.path(ctx$dirs$advanced, "tradeSeq_fit.rds"), compress = FALSE)

  assoc <- tryCatch(tradeSeq::associationTest(fit, global = TRUE, lineages = TRUE, l2fc = cfg$tradeSeq_l2fc %||% 0), error = function(e) NULL)
  start_end <- tryCatch(tradeSeq::startVsEndTest(fit, global = TRUE, lineages = TRUE, l2fc = cfg$tradeSeq_l2fc %||% 0), error = function(e) NULL)
  write_tradeseq_result(assoc, file.path(ctx$dirs$advanced, "tradeSeq_association_test.csv"))
  write_tradeseq_result(start_end, file.path(ctx$dirs$advanced, "tradeSeq_start_vs_end_test.csv"))
  if (ncol(pseudotime) > 1L) {
    pattern <- tryCatch(tradeSeq::patternTest(fit, global = TRUE, pairwise = TRUE, l2fc = cfg$tradeSeq_l2fc %||% 0), error = function(e) NULL)
    endpoints <- tryCatch(tradeSeq::diffEndTest(fit, global = TRUE, pairwise = TRUE, l2fc = cfg$tradeSeq_l2fc %||% 0), error = function(e) NULL)
    write_tradeseq_result(pattern, file.path(ctx$dirs$advanced, "tradeSeq_pattern_test.csv"))
    write_tradeseq_result(endpoints, file.path(ctx$dirs$advanced, "tradeSeq_differential_end_test.csv"))
  }
  write_yaml(list(selected_knots = selected_knots, genes_fitted = length(genes), fixed_effects = fixed,
                  interpretation = "exploratory NB-GAM tests along a prespecified fitted trajectory"),
             file.path(ctx$dirs$advanced, "tradeSeq_parameters.yml"))
  status_table("COMPLETED", paste("tradeSeq fitted", length(genes), "genes with", selected_knots, "knots"), length(genes))
}

run_trajectory_module <- function(ctx, object_path = file.path(ctx$dirs$objects, "06_annotated.rds")) {
  cfg <- ctx$cfg$advanced_modules$trajectory
  if (!as_flag(cfg$enabled, FALSE)) return(status_table("NOT_RUN", "disabled"))
  if (!requireNamespace("slingshot", quietly = TRUE) || !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    return(status_table("NOT_RUN", "slingshot/SingleCellExperiment unavailable"))
  }
  if (is.null(cfg$subset_column) || !length(cfg$subset_values) || is.null(cfg$root_cluster)) {
    return(status_table("NOT_EVALUABLE", "explicit subset_column, subset_values and root_cluster are required"))
  }
  obj <- readRDS(object_path)
  if (!cfg$subset_column %in% colnames(obj[[]])) return(status_table("NOT_EVALUABLE", "trajectory subset column is absent"))
  keep <- obj[[cfg$subset_column]][, 1] %in% unlist(cfg$subset_values)
  sub <- subset(obj, cells = colnames(obj)[keep])
  if (ncol(sub) < (cfg$minimum_cells %||% 200L)) return(status_table("NOT_EVALUABLE", "too few cells in the declared lineage subset"))
  cluster_col <- cfg$cluster_column %||% "cluster_raw"
  if (!cluster_col %in% colnames(sub[[]])) return(status_table("NOT_EVALUABLE", "trajectory cluster column is absent"))
  if (!as.character(cfg$root_cluster) %in% unique(as.character(sub[[cluster_col]][, 1]))) return(status_table("NOT_EVALUABLE", "root cluster is absent from the subset"))
  reduction <- cfg$reduction %||% "umap"
  if (!reduction %in% names(sub@reductions)) return(status_table("NOT_EVALUABLE", paste(reduction, "reduction is required")))

  sce <- Seurat::as.SingleCellExperiment(sub, assay = "RNA")
  SingleCellExperiment::reducedDim(sce, toupper(reduction)) <- Seurat::Embeddings(sub, reduction)
  sce <- slingshot::slingshot(
    sce, clusterLabels = as.character(sub[[cluster_col]][, 1]), reducedDim = toupper(reduction),
    start.clus = as.character(cfg$root_cluster), end.clus = unlist(cfg$end_clusters %||% list()),
    allow.breaks = isTRUE(cfg$allow_breaks %||% FALSE), shrink = cfg$shrink %||% TRUE
  )
  saveRDS(sce, file.path(ctx$dirs$advanced, "slingshot_sce.rds"), compress = FALSE)

  coords <- SingleCellExperiment::reducedDim(sce, toupper(reduction))[, 1:2, drop = FALSE]
  pt <- slingshot::slingPseudotime(sce, na = TRUE)
  wt <- slingshot::slingCurveWeights(sce)
  if (is.null(dim(pt))) pt <- matrix(pt, ncol = 1L)
  if (is.null(dim(wt))) wt <- matrix(wt, ncol = 1L)
  rownames(pt) <- rownames(coords); rownames(wt) <- rownames(coords)
  colnames(pt) <- paste0("lineage_", seq_len(ncol(pt))); colnames(wt) <- paste0("lineage_", seq_len(ncol(wt)))
  weighted_pt <- rowSums(pt * wt, na.rm = TRUE) / rowSums(wt, na.rm = TRUE)
  weighted_pt[!is.finite(weighted_pt)] <- NA_real_

  md <- sub[[]]
  trend_columns <- unlist(cfg$trend_columns %||% list())
  trend_columns <- trend_columns[trend_columns %in% colnames(md)]
  cell_data <- data.frame(cell_id = rownames(coords), embedding_1 = coords[, 1], embedding_2 = coords[, 2],
                          cluster_raw = as.character(sub[[cluster_col]][, 1]), pseudotime = weighted_pt, stringsAsFactors = FALSE)
  for (nm in intersect(c("sample_id", "patient_id", "condition", "cell_type_l1", "cell_type_l2"), colnames(md))) cell_data[[nm]] <- md[rownames(coords), nm]
  for (nm in trend_columns) cell_data[[nm]] <- md[rownames(coords), nm]
  for (j in seq_len(ncol(pt))) {cell_data[[paste0("pseudotime_", j)]] <- pt[, j]; cell_data[[paste0("weight_", j)]] <- wt[, j]}
  write_csv_safe(cell_data, file.path(ctx$dirs$advanced, "trajectory_cell_data.csv"))

  curves <- slingshot::slingCurves(sce)
  curve_rows <- list()
  for (j in seq_along(curves)) {
    crv <- curves[[j]]; s <- crv$s
    if (is.null(s) || !nrow(s)) next
    ord <- crv$ord %||% seq_len(nrow(s)); ord <- ord[ord >= 1L & ord <= nrow(s)]
    z <- s[ord, 1:2, drop = FALSE]
    curve_rows[[j]] <- data.frame(lineage = paste0("lineage_", j), curve_order = seq_len(nrow(z)), embedding_1 = z[, 1], embedding_2 = z[, 2], stringsAsFactors = FALSE)
  }
  if (!length(curve_rows)) return(status_table("NOT_EVALUABLE", "Slingshot fitted no drawable lineage curve"))
  write_csv_safe(do.call(rbind, curve_rows), file.path(ctx$dirs$advanced, "trajectory_curve_data.csv"))

  node_rows <- lapply(split(seq_len(nrow(cell_data)), cell_data$cluster_raw), function(ii) data.frame(
    cluster_raw = cell_data$cluster_raw[ii[1]], embedding_1 = stats::median(cell_data$embedding_1[ii], na.rm = TRUE),
    embedding_2 = stats::median(cell_data$embedding_2[ii], na.rm = TRUE), n_cells = length(ii), stringsAsFactors = FALSE))
  write_csv_safe(do.call(rbind, node_rows), file.path(ctx$dirs$advanced, "trajectory_node_data.csv"))

  ts_status <- run_tradeseq_layer(ctx, sub, pt, wt, cfg)
  write_csv_safe(ts_status, file.path(ctx$dirs$advanced, "tradeSeq_status.csv"))
  write_yaml(list(root_cluster = cfg$root_cluster, end_clusters = unlist(cfg$end_clusters %||% list()), subset_column = cfg$subset_column,
                  subset_values = unlist(cfg$subset_values), reduction = reduction, allow_breaks = cfg$allow_breaks %||% FALSE,
                  shrink = cfg$shrink %||% TRUE, n_cells = nrow(cell_data), n_lineages = length(curves),
                  interpretation = "exploratory trajectory with prespecified root; topology and direction require biological validation"),
             file.path(ctx$dirs$advanced, "trajectory_parameters.yml"))
  status_table("COMPLETED", paste("Slingshot exported", length(curves), "lineage(s); tradeSeq:", ts_status$status[1]), nrow(cell_data))
}
