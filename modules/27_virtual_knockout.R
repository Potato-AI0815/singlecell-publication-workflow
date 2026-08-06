# Exploratory virtual-knockout module.
# Primary backend: scTenifoldKnk, using raw RNA counts from an already QC-filtered Seurat object.
# The module estimates network perturbation magnitude. scTenifoldKnk distances are unsigned and
# must not be described as predicted up- or down-regulation.

vk_clean_scalar <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[[1]]) || !nzchar(trimws(as.character(x[[1]])))) return(NULL)
  trimws(as.character(x[[1]]))
}

vk_read_targets <- function(cfg) {
  direct <- unique(trimws(as.character(unlist(cfg$target_genes %||% list()))))
  direct <- direct[!is.na(direct) & nzchar(direct)]
  tab <- data.frame(target_gene = direct, target_label = direct, stringsAsFactors = FALSE)
  path <- vk_clean_scalar(cfg$targets_path)
  if (!is.null(path)) {
    if (!file.exists(path)) stop("virtual-knockout targets_path not found: ", path)
    x <- read_csv_safe(path)
    nms <- tolower(names(x))
    gcol <- names(x)[match(TRUE, nms %in% c("target_gene", "gene", "target"), nomatch = 0L)]
    if (!length(gcol)) stop("targets_path requires a target_gene or gene column")
    out <- data.frame(target_gene = trimws(as.character(x[[gcol]])), stringsAsFactors = FALSE)
    out$target_label <- if ("target_label" %in% names(x)) as.character(x$target_label) else out$target_gene
    for (nm in intersect(c("subset_value", "priority", "rationale"), names(x))) out[[nm]] <- x[[nm]]
    all_cols <- union(names(tab), names(out))
    for (nm in setdiff(all_cols, names(tab))) tab[[nm]] <- NA
    for (nm in setdiff(all_cols, names(out))) out[[nm]] <- NA
    tab <- rbind(tab[, all_cols, drop = FALSE], out[, all_cols, drop = FALSE])
  }
  tab <- tab[!is.na(tab$target_gene) & nzchar(tab$target_gene), , drop = FALSE]
  tab <- tab[!duplicated(toupper(tab$target_gene)), , drop = FALSE]
  rownames(tab) <- NULL
  tab
}

vk_resolve_symbols <- function(symbols, available) {
  lookup <- split(available, toupper(available))
  rows <- lapply(symbols, function(g) {
    hit <- lookup[[toupper(g)]]
    if (is.null(hit)) return(data.frame(requested = g, resolved = NA_character_, status = "ABSENT", stringsAsFactors = FALSE))
    if (length(hit) > 1L) return(data.frame(requested = g, resolved = NA_character_, status = "AMBIGUOUS_CASE_MATCH", stringsAsFactors = FALSE))
    data.frame(requested = g, resolved = hit[[1]], status = if (identical(g, hit[[1]])) "EXACT" else "CASE_MATCH", stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

vk_exclude_gene_mask <- function(genes, patterns) {
  keep <- rep(TRUE, length(genes))
  for (p in patterns) keep <- keep & !grepl(p, genes, ignore.case = TRUE, perl = TRUE)
  keep
}

vk_sparse_logcpm_variance <- function(counts, maximum_cells = 5000L, seed = 1L) {
  if (ncol(counts) > maximum_cells) {
    set.seed(seed)
    counts <- counts[, sample.int(ncol(counts), maximum_cells), drop = FALSE]
  }
  lib <- Matrix::colSums(counts)
  lib[!is.finite(lib) | lib <= 0] <- 1
  norm <- counts %*% Matrix::Diagonal(x = 1e4 / lib)
  norm@x <- log1p(norm@x)
  mu <- Matrix::rowMeans(norm)
  mu2 <- Matrix::rowMeans(norm ^ 2)
  vv <- pmax(mu2 - mu ^ 2, 0)
  data.frame(gene = rownames(counts), mean_logcpm = as.numeric(mu), variance_logcpm = as.numeric(vv), stringsAsFactors = FALSE)
}

vk_build_gene_universe <- function(counts, md, target_genes, cfg, seed) {
  genes <- rownames(counts)
  patterns <- unlist(cfg$exclude_gene_patterns %||% c("^MT-", "^RPL", "^RPS", "^HB[AB]"))
  keep <- vk_exclude_gene_mask(genes, patterns)
  counts0 <- counts[keep, , drop = FALSE]
  genes0 <- rownames(counts0)
  detect_global <- Matrix::rowMeans(counts0 > 0)
  sample_ids <- unique(as.character(md$sample_id))
  det_by_sample <- lapply(sample_ids, function(s) Matrix::rowMeans(counts0[, rownames(md)[md$sample_id == s], drop = FALSE] > 0))
  if (length(det_by_sample)) {
    det_mat <- do.call(cbind, det_by_sample)
    sample_support <- rowMeans(det_mat >= (cfg$minimum_detection_fraction %||% 0.01))
  } else sample_support <- rep(0, nrow(counts0))
  candidates <- genes0[detect_global >= (cfg$minimum_global_detection_fraction %||% 0.01) & sample_support >= (cfg$minimum_sample_fraction_for_gene %||% 0.5)]
  configured_path <- vk_clean_scalar(cfg$gene_universe_path)
  method <- tolower(cfg$gene_universe_method %||% "variable_features")
  if (!is.null(configured_path)) {
    if (!file.exists(configured_path)) stop("virtual-knockout gene_universe_path not found: ", configured_path)
    gx <- read_csv_safe(configured_path)
    gcol <- intersect(c("gene", "gene_symbol", "feature"), names(gx))[1]
    if (is.na(gcol)) stop("gene_universe_path requires a gene column")
    chosen <- intersect(as.character(gx[[gcol]]), genes0)
    method <- "configured"
  } else {
    metrics <- vk_sparse_logcpm_variance(counts0[candidates, , drop = FALSE], cfg$variance_estimation_maximum_cells %||% 5000L, seed)
    metrics <- metrics[order(metrics$variance_logcpm, decreasing = TRUE), , drop = FALSE]
    maxg <- cfg$maximum_genes %||% 2500L
    chosen <- head(metrics$gene, maxg)
    if (identical(method, "expressed")) chosen <- head(candidates, maxg)
  }
  chosen <- unique(c(target_genes, chosen))
  chosen <- intersect(chosen, genes0)
  min_genes <- cfg$minimum_genes %||% 500L
  if (length(chosen) < min_genes) {
    metrics2 <- vk_sparse_logcpm_variance(counts0, cfg$variance_estimation_maximum_cells %||% 5000L, seed)
    extra <- metrics2$gene[order(metrics2$variance_logcpm, decreasing = TRUE)]
    chosen <- unique(c(chosen, extra))
    chosen <- head(chosen, max(cfg$maximum_genes %||% 2500L, min_genes))
  }
  audit <- data.frame(
    gene = genes0,
    global_detection_fraction = as.numeric(detect_global),
    sample_support_fraction = as.numeric(sample_support),
    selected = genes0 %in% chosen,
    target_gene = genes0 %in% target_genes,
    selection_method = method,
    stringsAsFactors = FALSE
  )
  list(genes = chosen, audit = audit)
}

vk_target_expression_audit <- function(counts, md, targets, subset_col, subset_values, baseline_conditions, condition_col) {
  rows <- list(); z <- 0L
  for (sv in subset_values) {
    keep0 <- as.character(md[[subset_col]]) == sv
    if (length(baseline_conditions) && condition_col %in% names(md)) keep0 <- keep0 & as.character(md[[condition_col]]) %in% baseline_conditions
    keep0[is.na(keep0)] <- FALSE
    for (s in unique(as.character(md$sample_id[keep0]))) {
      cells <- rownames(md)[keep0 & as.character(md$sample_id) == s]
      for (g in targets) {
        z <- z + 1L
        vals <- if (length(cells)) as.numeric(counts[g, cells, drop = TRUE]) else numeric()
        rows[[z]] <- data.frame(
          target_gene = g,
          subset_value = sv,
          sample_id = s,
          n_cells = length(cells),
          n_target_expressing_cells = sum(vals > 0),
          target_expressing_fraction = if (length(vals)) mean(vals > 0) else NA_real_,
          target_mean_raw_count = if (length(vals)) mean(vals) else NA_real_,
          target_total_raw_count = if (length(vals)) sum(vals) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

vk_select_cells <- function(cells, maximum_cells, seed, fraction = 1, minimum_cells = 1L) {
  cells <- unique(as.character(cells))
  if (!length(cells)) return(character())
  fraction <- suppressWarnings(as.numeric(fraction %||% 1))
  if (!is.finite(fraction) || fraction <= 0 || fraction > 1) fraction <- 1
  target_n <- min(as.integer(maximum_cells), length(cells))
  if (fraction < 1) {
    target_n <- min(target_n, max(as.integer(minimum_cells), floor(length(cells) * fraction)))
  }
  target_n <- min(max(1L, target_n), length(cells))
  if (length(cells) <= target_n) return(cells)
  set.seed(seed)
  sample(cells, target_n, replace = FALSE)
}

vk_extract_manifold <- function(ma, target_gene, subset_value, sample_id, repeat_id, analysis_unit) {
  if (is.null(ma) || is.null(rownames(ma)) || !nrow(ma) || !ncol(ma)) return(data.frame())
  rn <- rownames(ma)
  # Current and historical scTenifoldNet builds have used X_/Y_ or WT_/KO_ row prefixes.
  xrows <- grep("^(X|WT)_", rn, value = TRUE)
  yrows <- grep("^(Y|KO)_", rn, value = TRUE)
  if (!length(xrows) || length(xrows) != length(yrows)) return(data.frame())
  gx <- sub("^(X|WT)_", "", xrows)
  gy <- sub("^(Y|KO)_", "", yrows)
  if (!setequal(gx, gy)) return(data.frame())
  yrows <- yrows[match(gx, gy)]
  dims <- min(2L, ncol(ma))
  make <- function(rows, state) {
    out <- data.frame(
      target_gene = target_gene,
      subset_value = subset_value,
      sample_id = sample_id,
      repeat_id = repeat_id,
      analysis_unit = analysis_unit,
      gene = gx,
      state = state,
      dim1 = as.numeric(ma[rows, 1]),
      dim2 = if (dims >= 2L) as.numeric(ma[rows, 2]) else 0,
      stringsAsFactors = FALSE
    )
    out
  }
  rbind(make(xrows, "WT"), make(yrows, "KO"))
}

vk_get_tensor_network <- function(ans, state = c("WT", "KO")) {
  state <- match.arg(state)
  nets <- ans$tensorNetworks
  if (is.null(nets)) return(NULL)
  candidates <- if (state == "WT") c("WT", "X") else c("KO", "Y")
  for (nm in candidates) if (!is.null(nets[[nm]])) return(nets[[nm]])
  NULL
}

vk_extract_network_edges <- function(network, dr, target_gene, subset_value, sample_id, repeat_id, analysis_unit, top_genes = 50L, edge_quantile = 0.97) {
  if (is.null(network) || !target_gene %in% rownames(network)) return(data.frame())
  d <- dr[order(dr$p.adj, -dr$distance), , drop = FALSE]
  genes <- unique(c(target_gene, head(as.character(d$gene), top_genes)))
  genes <- intersect(genes, rownames(network))
  if (length(genes) < 2L) return(data.frame())
  m <- as.matrix(network[genes, genes, drop = FALSE])
  diag(m) <- 0
  nz <- abs(m[is.finite(m) & m != 0])
  if (!length(nz)) return(data.frame())
  cut <- as.numeric(stats::quantile(nz, probs = edge_quantile, na.rm = TRUE, names = FALSE))
  idx <- which(is.finite(m) & m != 0 & abs(m) >= cut, arr.ind = TRUE)
  if (!nrow(idx)) return(data.frame())
  data.frame(
    target_gene = target_gene,
    subset_value = subset_value,
    sample_id = sample_id,
    repeat_id = repeat_id,
    analysis_unit = analysis_unit,
    from = rownames(m)[idx[, 1]],
    to = colnames(m)[idx[, 2]],
    weight = as.numeric(m[idx]),
    absolute_weight = abs(as.numeric(m[idx])),
    edge_sign = ifelse(as.numeric(m[idx]) >= 0, "positive", "negative"),
    target_outgoing = rownames(m)[idx[, 1]] == target_gene,
    stringsAsFactors = FALSE
  )
}

vk_run_one <- function(counts, target_gene, cfg, seed, prior_network = NULL) {
  set.seed(seed)
  scTenifoldKnk::scTenifoldKnk(
    countMatrix = counts,
    gKO = target_gene,
    qc = FALSE,
    nc_lambda = cfg$directionality_lambda %||% 0,
    nc_nNet = cfg$n_networks %||% 10L,
    nc_nCells = min(cfg$cells_per_network %||% 500L, ncol(counts)),
    nc_nComp = min(cfg$principal_components %||% 3L, max(3L, nrow(counts) - 1L)),
    nc_scaleScores = as_flag(cfg$scale_scores, TRUE),
    nc_symmetric = as_flag(cfg$symmetric, FALSE),
    nc_q = cfg$network_quantile %||% 0.9,
    nc_priorNetwork = prior_network,
    td_K = cfg$tensor_rank %||% 3L,
    td_maxIter = cfg$tensor_max_iter %||% 1000L,
    td_maxError = cfg$tensor_max_error %||% 1e-5,
    td_nDecimal = cfg$tensor_decimal_places %||% 3L,
    ma_nDim = cfg$manifold_dimensions %||% 2L,
    nCores = cfg$workers %||% 1L
  )
}

vk_prepare_backend_result <- function(ans) {
  if (!is.list(ans)) stop("scTenifoldKnk backend did not return a list")
  if (is.null(ans$diffRegulation)) stop("scTenifoldKnk result is missing diffRegulation")
  dr <- as.data.frame(ans$diffRegulation)
  required <- c("gene", "distance", "Z", "FC", "p.value", "p.adj")
  missing <- setdiff(required, names(dr))
  if (length(missing)) stop("scTenifoldKnk diffRegulation is missing columns: ", paste(missing, collapse = ", "))
  if (!nrow(dr)) stop("scTenifoldKnk diffRegulation is empty")
  dr$gene <- as.character(dr$gene)
  for (nm in setdiff(required, "gene")) dr[[nm]] <- suppressWarnings(as.numeric(dr[[nm]]))
  dr
}

vk_fisher_p <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  stats::pchisq(-2 * sum(log(pmax(p, .Machine$double.xmin))), df = 2 * length(p), lower.tail = FALSE)
}

vk_aggregate_results <- function(all_dr, fdr_threshold, minimum_repeat_recurrence, minimum_sample_recurrence) {
  if (!nrow(all_dr)) return(list(sample = data.frame(), consensus = data.frame()))
  keys <- interaction(all_dr$target_gene, all_dr$subset_value, all_dr$sample_id, all_dr$gene, drop = TRUE, lex.order = TRUE)
  sample_rows <- lapply(split(all_dr, keys), function(x) {
    data.frame(
      target_gene = x$target_gene[1],
      subset_value = x$subset_value[1],
      sample_id = x$sample_id[1],
      gene = x$gene[1],
      n_repeats = length(unique(x$repeat_id)),
      repeat_hit_count = sum(x$p.adj < fdr_threshold, na.rm = TRUE),
      repeat_recurrence_fraction = mean(x$p.adj < fdr_threshold, na.rm = TRUE),
      median_distance = stats::median(x$distance, na.rm = TRUE),
      median_Z = stats::median(x$Z, na.rm = TRUE),
      median_FC = stats::median(x$FC, na.rm = TRUE),
      median_p_value = stats::median(x$p.value, na.rm = TRUE),
      minimum_p_adj = min(x$p.adj, na.rm = TRUE),
      sample_supported = mean(x$p.adj < fdr_threshold, na.rm = TRUE) >= minimum_repeat_recurrence,
      stringsAsFactors = FALSE
    )
  })
  sample_tab <- do.call(rbind, sample_rows)
  keys2 <- interaction(sample_tab$target_gene, sample_tab$subset_value, sample_tab$gene, drop = TRUE, lex.order = TRUE)
  con_rows <- lapply(split(sample_tab, keys2), function(x) {
    p_comb <- vk_fisher_p(x$median_p_value)
    data.frame(
      target_gene = x$target_gene[1],
      subset_value = x$subset_value[1],
      gene = x$gene[1],
      n_samples = length(unique(x$sample_id)),
      n_supported_samples = sum(x$sample_supported, na.rm = TRUE),
      sample_recurrence_fraction = mean(x$sample_supported, na.rm = TRUE),
      median_sample_distance = stats::median(x$median_distance, na.rm = TRUE),
      median_sample_Z = stats::median(x$median_Z, na.rm = TRUE),
      median_sample_FC = stats::median(x$median_FC, na.rm = TRUE),
      fisher_p_value = p_comb,
      stringsAsFactors = FALSE
    )
  })
  con <- do.call(rbind, con_rows)
  con$fisher_fdr <- ave(con$fisher_p_value, interaction(con$target_gene, con$subset_value, drop = TRUE), FUN = function(z) stats::p.adjust(z, method = "BH"))
  con$consensus_impact_score <- con$sample_recurrence_fraction * pmax(con$median_sample_Z, 0) * (-log10(pmax(con$fisher_fdr, .Machine$double.xmin)))
  con$evidence_class <- ifelse(
    con$sample_recurrence_fraction >= minimum_sample_recurrence & con$fisher_fdr < fdr_threshold,
    "Replicate-consistent",
    ifelse(con$n_supported_samples > 0, "Sample-restricted", "Not-supported")
  )
  con <- con[order(con$target_gene, con$subset_value, -con$consensus_impact_score, con$fisher_fdr), , drop = FALSE]
  list(sample = sample_tab, consensus = con)
}

vk_aggregate_edges <- function(edges, minimum_recurrence = 0.5) {
  if (!nrow(edges)) return(data.frame())
  group_run_totals <- tapply(
    interaction(edges$sample_id, edges$repeat_id, edges$analysis_unit, drop = TRUE),
    interaction(edges$target_gene, edges$subset_value, drop = TRUE, lex.order = TRUE),
    function(z) length(unique(z))
  )
  keys <- interaction(edges$target_gene, edges$subset_value, edges$from, edges$to, drop = TRUE, lex.order = TRUE)
  out <- lapply(split(edges, keys), function(x) {
    group_key <- as.character(interaction(x$target_gene[1], x$subset_value[1], drop = TRUE, lex.order = TRUE))
    total_runs <- as.integer(group_run_totals[[group_key]])
    n_runs <- length(unique(interaction(x$sample_id, x$repeat_id, x$analysis_unit, drop = TRUE)))
    data.frame(
    target_gene = x$target_gene[1], subset_value = x$subset_value[1], from = x$from[1], to = x$to[1],
    n_runs = n_runs,
    total_runs = total_runs,
    edge_recurrence_fraction = if (is.finite(total_runs) && total_runs > 0L) n_runs / total_runs else NA_real_,
    median_weight = stats::median(x$weight, na.rm = TRUE),
    median_absolute_weight = stats::median(x$absolute_weight, na.rm = TRUE),
    target_outgoing = any(x$target_outgoing),
    stringsAsFactors = FALSE
  )})
  ans <- do.call(rbind, out)
  ans$edge_sign <- ifelse(ans$median_weight >= 0, "positive", "negative")
  ans <- ans[ans$edge_recurrence_fraction >= minimum_recurrence | ans$target_outgoing, , drop = FALSE]
  ans[order(ans$target_gene, ans$subset_value, -ans$target_outgoing, -ans$edge_recurrence_fraction, -ans$median_absolute_weight), , drop = FALSE]
}

vk_run_similarity <- function(all_dr, top_n = 100L) {
  if (!nrow(all_dr)) return(data.frame())
  rows <- list(); z <- 0L
  groups <- split(all_dr, interaction(all_dr$target_gene, all_dr$subset_value, drop = TRUE))
  for (x in groups) {
    x$run_id <- paste(x$sample_id, paste0("R", x$repeat_id), x$analysis_unit, sep = "|")
    runs <- unique(x$run_id)
    top <- lapply(runs, function(r) head(x$gene[x$run_id == r][order(x$p.adj[x$run_id == r], -x$distance[x$run_id == r])], top_n))
    names(top) <- runs
    for (i in seq_along(runs)) for (j in seq_along(runs)) {
      a <- unique(top[[i]]); b <- unique(top[[j]])
      z <- z + 1L
      rows[[z]] <- data.frame(
        target_gene = x$target_gene[1], subset_value = x$subset_value[1],
        run_1 = runs[i], run_2 = runs[j],
        jaccard_top_genes = if (length(union(a, b))) length(intersect(a, b)) / length(union(a, b)) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

vk_load_pathways <- function(ctx, cfg) {
  gmt <- vk_clean_scalar(cfg$gmt_path %||% ctx$cfg$pathways$gmt_path)
  if (!is.null(gmt) && file.exists(gmt) && requireNamespace("fgsea", quietly = TRUE)) return(fgsea::gmtPathways(gmt))
  if (!requireNamespace("msigdbr", quietly = TRUE)) return(NULL)
  species <- if (tolower(ctx$cfg$input$species %||% "human") == "mouse") "Mus musculus" else "Homo sapiens"
  fm <- names(formals(msigdbr::msigdbr))
  args <- list(species = species)
  collection <- cfg$pathway_collection %||% "H"
  if ("collection" %in% fm) args$collection <- collection else args$category <- collection
  ms <- tryCatch(do.call(msigdbr::msigdbr, args), error = function(e) NULL)
  if (is.null(ms) || !nrow(ms)) return(NULL)
  split(ms$gene_symbol, ms$gs_name)
}

vk_pathway_ora <- function(consensus, universe, pathways, cfg) {
  if (is.null(pathways) || !length(pathways) || !nrow(consensus)) return(data.frame())
  rows <- list(); z <- 0L
  fdr <- cfg$fdr_threshold %||% 0.05
  min_query <- cfg$minimum_pathway_query_genes %||% 10L
  min_overlap <- cfg$minimum_pathway_overlap %||% 3L
  groups <- split(consensus, interaction(consensus$target_gene, consensus$subset_value, drop = TRUE))
  for (x in groups) {
    query <- unique(x$gene[x$evidence_class == "Replicate-consistent" & x$fisher_fdr < fdr])
    query <- setdiff(intersect(query, universe), x$target_gene[1])
    if (length(query) < min_query) next
    for (nm in names(pathways)) {
      gs <- unique(intersect(pathways[[nm]], universe))
      ov <- intersect(query, gs)
      if (length(ov) < min_overlap) next
      p <- stats::phyper(length(ov) - 1L, length(gs), length(universe) - length(gs), length(query), lower.tail = FALSE)
      z <- z + 1L
      rows[[z]] <- data.frame(
        target_gene = x$target_gene[1], subset_value = x$subset_value[1], pathway = nm,
        overlap = length(ov), pathway_size = length(gs), query_size = length(query), universe_size = length(universe),
        p_value = p, genes = paste(sort(ov), collapse = ";"), stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(data.frame())
  ans <- do.call(rbind, rows)
  ans$FDR <- ave(ans$p_value, interaction(ans$target_gene, ans$subset_value, drop = TRUE), FUN = function(z) stats::p.adjust(z, "BH"))
  ans$enrichment_ratio <- ans$overlap / pmax(ans$pathway_size, 1)
  ans[order(ans$target_gene, ans$subset_value, ans$FDR, -ans$overlap), , drop = FALSE]
}

run_virtual_knockout_module <- function(ctx, object_path = file.path(ctx$dirs$objects, "06_annotated.rds")) {
  cfg <- ctx$cfg$advanced_modules$virtual_knockout
  if (!as_flag(cfg$enabled, FALSE)) return(status_table("NOT_RUN", "disabled"))
  backend <- tolower(cfg$backend %||% "sctenifoldknk")
  if (!identical(backend, "sctenifoldknk")) return(status_table("NOT_EVALUABLE", paste("unsupported virtual-knockout backend in this R engine:", backend)))
  if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) return(status_table("NOT_RUN", "scTenifoldKnk package unavailable"))
  if (!file.exists(object_path)) return(status_table("NOT_EVALUABLE", "annotated Seurat object is required"))

  targets_tab <- vk_read_targets(cfg)
  if (!nrow(targets_tab)) return(status_table("NOT_EVALUABLE", "predeclared target_genes or targets_path is required"))
  max_targets <- cfg$maximum_targets %||% 20L
  if (nrow(targets_tab) > max_targets && !as_flag(cfg$allow_large_screen, FALSE)) return(status_table("NOT_EVALUABLE", paste("target count exceeds maximum_targets", max_targets)))

  obj <- readRDS(object_path)
  assay <- cfg$assay %||% "RNA"
  obj <- join_layers_safe(obj, assay)
  counts <- get_layer_safe(obj, assay, "counts")
  md <- obj[[]]
  if (!"sample_id" %in% names(md)) return(status_table("NOT_EVALUABLE", "sample_id metadata is required"))
  bad_sample <- is.na(md$sample_id) | !nzchar(trimws(as.character(md$sample_id)))
  if (any(bad_sample)) return(status_table("NOT_EVALUABLE", paste(sum(bad_sample), "cells have missing sample_id")))
  subset_col <- cfg$subset_column %||% "cell_type_l1"
  if (!subset_col %in% names(md)) return(status_table("NOT_EVALUABLE", paste("subset_column missing:", subset_col)))
  subset_values <- unique(as.character(unlist(cfg$subset_values %||% list())))
  subset_values <- subset_values[!is.na(subset_values) & nzchar(subset_values)]
  if (!length(subset_values)) {
    if (!as_flag(cfg$allow_all_cells, FALSE)) return(status_table("NOT_EVALUABLE", "explicit subset_values are required unless allow_all_cells is true"))
    md$.__virtual_knockout_all__ <- "All_cells"
    subset_col <- ".__virtual_knockout_all__"
    subset_values <- "All_cells"
  }
  missing_subset <- setdiff(subset_values, unique(as.character(md[[subset_col]])))
  if (length(missing_subset)) append_issue(ctx, "ADVANCED", "WARNING", "virtual-knockout subset values absent", paste(missing_subset, collapse = ","))
  subset_values <- intersect(subset_values, unique(as.character(md[[subset_col]])))
  if (!length(subset_values)) return(status_table("NOT_EVALUABLE", "none of the configured subset_values are present"))

  condition_col <- cfg$condition_column %||% "condition"
  baseline_conditions <- unique(as.character(unlist(cfg$baseline_conditions %||% list())))
  baseline_conditions <- baseline_conditions[!is.na(baseline_conditions) & nzchar(baseline_conditions)]
  if (condition_col %in% names(md)) {
    observed_conditions <- unique(stats::na.omit(as.character(md[[condition_col]])))
    if (length(observed_conditions) > 1L && !length(baseline_conditions) && as_flag(cfg$require_baseline_condition_when_multiple, TRUE)) {
      return(status_table("NOT_EVALUABLE", "multiple conditions are present; baseline_conditions must be predeclared"))
    }
    if (length(baseline_conditions)) {
      absent_cond <- setdiff(baseline_conditions, observed_conditions)
      if (length(absent_cond)) return(status_table("NOT_EVALUABLE", paste("baseline condition absent:", paste(absent_cond, collapse = ","))))
    }
  } else if (length(baseline_conditions)) return(status_table("NOT_EVALUABLE", paste("condition column missing:", condition_col)))

  resolved <- vk_resolve_symbols(targets_tab$target_gene, rownames(counts))
  write_csv_safe(resolved, file.path(ctx$dirs$advanced, "virtual_knockout_target_symbol_audit.csv"))
  resolved_targets <- resolved$resolved[resolved$status %in% c("EXACT", "CASE_MATCH")]
  if (!length(resolved_targets)) return(status_table("NOT_EVALUABLE", "none of the predeclared target genes are present"))
  targets_tab$resolved_target_gene <- resolved$resolved[match(targets_tab$target_gene, resolved$requested)]
  targets_tab <- targets_tab[!is.na(targets_tab$resolved_target_gene), , drop = FALSE]
  if (!"subset_value" %in% names(targets_tab)) targets_tab$subset_value <- NA_character_
  targets_tab$subset_value <- trimws(as.character(targets_tab$subset_value))
  scoped_subset_values <- unique(targets_tab$subset_value[!is.na(targets_tab$subset_value) & nzchar(targets_tab$subset_value)])
  absent_scoped <- setdiff(scoped_subset_values, subset_values)
  if (length(absent_scoped)) return(status_table("NOT_EVALUABLE", paste("targets_path subset_value absent from configured/present subsets:", paste(absent_scoped, collapse = ","))))

  eligible0 <- as.character(md[[subset_col]]) %in% subset_values
  if (length(baseline_conditions)) eligible0 <- eligible0 & as.character(md[[condition_col]]) %in% baseline_conditions
  allowed_conf <- unique(as.character(unlist(cfg$allowed_annotation_confidence %||% c("High", "Medium"))))
  if ("annotation_confidence" %in% names(md) && length(allowed_conf)) eligible0 <- eligible0 & as.character(md$annotation_confidence) %in% allowed_conf
  eligible0[is.na(eligible0)] <- FALSE
  eligible_cells <- rownames(md)[eligible0]
  if (!length(eligible_cells)) return(status_table("NOT_EVALUABLE", "no cells pass subset, baseline and annotation-confidence gates"))
  counts_eligible <- counts[, eligible_cells, drop = FALSE]
  md_eligible <- md[eligible_cells, , drop = FALSE]

  gene_universe <- vk_build_gene_universe(counts_eligible, md_eligible, resolved_targets, cfg, ctx$cfg$runtime$seed %||% 20260804L)
  write_csv_safe(gene_universe$audit, file.path(ctx$dirs$advanced, "virtual_knockout_gene_universe.csv"))
  if (length(gene_universe$genes) < (cfg$minimum_genes %||% 500L)) return(status_table("NOT_EVALUABLE", "gene universe is smaller than minimum_genes"))

  target_audit <- vk_target_expression_audit(counts_eligible, md_eligible, resolved_targets, subset_col, subset_values, character(), condition_col)
  if (nrow(target_audit)) {
    scope_map <- setNames(targets_tab$subset_value, targets_tab$resolved_target_gene)
    intended <- vapply(seq_len(nrow(target_audit)), function(i) {
      sv0 <- scope_map[[target_audit$target_gene[i]]]
      is.null(sv0) || is.na(sv0) || !nzchar(sv0) || identical(as.character(sv0), as.character(target_audit$subset_value[i]))
    }, logical(1))
    target_audit <- target_audit[intended, , drop = FALSE]
  }
  write_csv_safe(target_audit, file.path(ctx$dirs$advanced, "virtual_knockout_target_expression_audit.csv"))
  if (!nrow(target_audit)) return(status_table("NOT_EVALUABLE", "target-expression audit produced no sample-level rows"))
  target_audit$eligible_by_cells <- target_audit$n_cells >= (cfg$minimum_cells_per_sample %||% 500L)
  target_audit$eligible_by_expression <- target_audit$n_target_expressing_cells >= (cfg$minimum_target_expressing_cells %||% 20L) &
    target_audit$target_expressing_fraction >= (cfg$minimum_target_expressing_fraction %||% 0.05)
  target_audit$eligible_for_knockout <- target_audit$eligible_by_cells & target_audit$eligible_by_expression
  write_csv_safe(target_audit, file.path(ctx$dirs$advanced, "virtual_knockout_target_expression_audit.csv"))
  target_subset_groups <- split(target_audit, interaction(target_audit$target_gene, target_audit$subset_value, drop = TRUE, lex.order = TRUE))
  target_subset_eval <- do.call(rbind, lapply(target_subset_groups, function(z) data.frame(
    target_gene = z$target_gene[1],
    subset_value = z$subset_value[1],
    n_eligible_samples = sum(z$eligible_for_knockout, na.rm = TRUE),
    n_total_samples = length(unique(z$sample_id)),
    stringsAsFactors = FALSE
  )))
  rownames(target_subset_eval) <- NULL
  target_subset_eval$replicate_consensus_evaluable <- target_subset_eval$n_eligible_samples >= (cfg$minimum_eligible_samples %||% 2L)
  write_csv_safe(target_subset_eval, file.path(ctx$dirs$advanced, "virtual_knockout_target_subset_evaluability.csv"))

  prior_network <- NULL
  prior_path <- vk_clean_scalar(cfg$prior_network_path)
  if (!is.null(prior_path)) {
    if (!file.exists(prior_path)) return(status_table("NOT_EVALUABLE", "prior_network_path does not exist"))
    prior_network <- read_csv_safe(prior_path)
    if (!all(c("regulators", "targets") %in% names(prior_network))) return(status_table("NOT_EVALUABLE", "prior network requires regulators and targets columns"))
  }

  analysis_unit <- tolower(cfg$analysis_unit %||% "sample_stratified_consensus")
  if (!analysis_unit %in% c("sample_stratified_consensus", "balanced_pooled", "both")) return(status_table("NOT_EVALUABLE", paste("unsupported analysis_unit:", analysis_unit)))
  repeats <- max(1L, cfg$subsampling_repeats %||% 3L)
  subsampling_fraction <- suppressWarnings(as.numeric(cfg$subsampling_fraction %||% 0.8))
  if (!is.finite(subsampling_fraction) || subsampling_fraction <= 0 || subsampling_fraction > 1) return(status_table("NOT_EVALUABLE", "subsampling_fraction must be in (0,1]"))
  minimum_cells <- cfg$minimum_cells_per_sample %||% 500L
  maximum_cells <- cfg$maximum_cells_per_sample %||% 5000L
  minimum_target_fraction <- cfg$minimum_target_expressing_fraction %||% 0.05
  minimum_target_cells <- cfg$minimum_target_expressing_cells %||% 20L
  seed0 <- ctx$cfg$runtime$seed %||% 20260804L
  save_objects <- as_flag(cfg$save_full_backend_objects, FALSE)
  run_dir <- safe_dir(file.path(ctx$dirs$advanced, "virtual_knockout_runs"))

  all_dr <- list(); all_ma <- list(); all_edges <- list(); manifest <- list(); k <- 0L; km <- 0L; ke <- 0L; kr <- 0L
  run_counter <- 0L
  for (sv in subset_values) {
    base_cells <- rownames(md)[as.character(md[[subset_col]]) == sv & eligible0]
    sample_split <- split(base_cells, as.character(md[base_cells, "sample_id", drop = TRUE]))
    for (tg in resolved_targets) {
      target_scope <- targets_tab$subset_value[match(tg, targets_tab$resolved_target_gene)]
      if (length(target_scope) && !is.na(target_scope) && nzchar(target_scope) && !identical(as.character(target_scope), as.character(sv))) next
      ta <- target_audit[target_audit$target_gene == tg & target_audit$subset_value == sv, , drop = FALSE]
      eligible_samples <- ta$sample_id[ta$eligible_for_knockout]
      eligible_samples <- intersect(eligible_samples, names(sample_split))
      if (!length(eligible_samples)) {
        kr <- kr + 1L
        manifest[[kr]] <- data.frame(target_gene = tg, subset_value = sv, sample_id = NA_character_, repeat_id = NA_integer_, analysis_unit = analysis_unit, n_cells = 0L, n_genes = length(gene_universe$genes), status = "NOT_EVALUABLE", reason = "no sample passed cell-count and target-expression gates", elapsed_seconds = 0, stringsAsFactors = FALSE)
        next
      }
      if (analysis_unit %in% c("sample_stratified_consensus", "both")) {
        for (s in eligible_samples) {
          cells0 <- sample_split[[s]]
          for (r in seq_len(repeats)) {
            run_counter <- run_counter + 1L
            seed <- seed0 + run_counter * 1009L
            cells <- vk_select_cells(cells0, maximum_cells, seed, fraction = subsampling_fraction, minimum_cells = minimum_cells)
            mat <- counts[gene_universe$genes, cells, drop = FALSE]
            run_name <- paste(sanitize_stem(tg), sanitize_stem(sv), sanitize_stem(s), paste0("R", r), sep = "__")
            started <- Sys.time()
            ans <- tryCatch(vk_run_one(mat, tg, cfg, seed, prior_network), error = function(e) e)
            elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
            if (inherits(ans, "error")) {
              kr <- kr + 1L
              manifest[[kr]] <- data.frame(target_gene = tg, subset_value = sv, sample_id = s, repeat_id = r, analysis_unit = "sample", n_cells = ncol(mat), n_genes = nrow(mat), status = "FAILED", reason = conditionMessage(ans), elapsed_seconds = elapsed, stringsAsFactors = FALSE)
              next
            }
            dr <- tryCatch(vk_prepare_backend_result(ans), error = function(e) e)
            if (inherits(dr, "error")) {
              kr <- kr + 1L
              manifest[[kr]] <- data.frame(target_gene = tg, subset_value = sv, sample_id = s, repeat_id = r, analysis_unit = "sample", n_cells = ncol(mat), n_genes = nrow(mat), status = "FAILED", reason = conditionMessage(dr), elapsed_seconds = elapsed, stringsAsFactors = FALSE)
              next
            }
            dr$target_gene <- tg; dr$subset_value <- sv; dr$sample_id <- s; dr$repeat_id <- r; dr$analysis_unit <- "sample"
            k <- k + 1L; all_dr[[k]] <- dr
            ma <- vk_extract_manifold(ans$manifoldAlignment, tg, sv, s, r, "sample")
            if (nrow(ma)) { km <- km + 1L; all_ma[[km]] <- ma }
            ed <- vk_extract_network_edges(vk_get_tensor_network(ans, "WT"), dr, tg, sv, s, r, "sample", cfg$network_top_genes %||% 50L, cfg$network_edge_quantile %||% 0.97)
            if (nrow(ed)) { ke <- ke + 1L; all_edges[[ke]] <- ed }
            if (save_objects) saveRDS(ans, file.path(run_dir, paste0(run_name, ".rds")), compress = FALSE)
            kr <- kr + 1L
            manifest[[kr]] <- data.frame(target_gene = tg, subset_value = sv, sample_id = s, repeat_id = r, analysis_unit = "sample", n_cells = ncol(mat), n_genes = nrow(mat), status = "COMPLETED", reason = NA_character_, elapsed_seconds = elapsed, stringsAsFactors = FALSE)
          }
        }
      }
      if (analysis_unit %in% c("balanced_pooled", "both") && length(eligible_samples)) {
        balanced_per_sample <- min(vapply(sample_split[eligible_samples], length, integer(1)), floor((cfg$maximum_pooled_cells %||% 10000L) / length(eligible_samples)))
        if (is.finite(balanced_per_sample) && balanced_per_sample >= (cfg$minimum_balanced_cells_per_sample %||% 100L)) {
          for (r in seq_len(repeats)) {
            run_counter <- run_counter + 1L
            seed <- seed0 + run_counter * 1009L
            selected <- unlist(lapply(seq_along(eligible_samples), function(i) vk_select_cells(sample_split[[eligible_samples[i]]], balanced_per_sample, seed + i * 53L, fraction = 1, minimum_cells = balanced_per_sample)), use.names = FALSE)
            mat <- counts[gene_universe$genes, selected, drop = FALSE]
            started <- Sys.time()
            ans <- tryCatch(vk_run_one(mat, tg, cfg, seed, prior_network), error = function(e) e)
            elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
            if (inherits(ans, "error")) {
              kr <- kr + 1L
              manifest[[kr]] <- data.frame(target_gene = tg, subset_value = sv, sample_id = "BALANCED_POOLED", repeat_id = r, analysis_unit = "balanced_pooled", n_cells = ncol(mat), n_genes = nrow(mat), status = "FAILED", reason = conditionMessage(ans), elapsed_seconds = elapsed, stringsAsFactors = FALSE)
              next
            }
            dr <- tryCatch(vk_prepare_backend_result(ans), error = function(e) e)
            if (inherits(dr, "error")) {
              kr <- kr + 1L
              manifest[[kr]] <- data.frame(target_gene = tg, subset_value = sv, sample_id = "BALANCED_POOLED", repeat_id = r, analysis_unit = "balanced_pooled", n_cells = ncol(mat), n_genes = nrow(mat), status = "FAILED", reason = conditionMessage(dr), elapsed_seconds = elapsed, stringsAsFactors = FALSE)
              next
            }
            dr$target_gene <- tg; dr$subset_value <- sv; dr$sample_id <- "BALANCED_POOLED"; dr$repeat_id <- r; dr$analysis_unit <- "balanced_pooled"
            k <- k + 1L; all_dr[[k]] <- dr
            ma <- vk_extract_manifold(ans$manifoldAlignment, tg, sv, "BALANCED_POOLED", r, "balanced_pooled")
            if (nrow(ma)) { km <- km + 1L; all_ma[[km]] <- ma }
            ed <- vk_extract_network_edges(vk_get_tensor_network(ans, "WT"), dr, tg, sv, "BALANCED_POOLED", r, "balanced_pooled", cfg$network_top_genes %||% 50L, cfg$network_edge_quantile %||% 0.97)
            if (nrow(ed)) { ke <- ke + 1L; all_edges[[ke]] <- ed }
            if (save_objects) saveRDS(ans, file.path(run_dir, paste0(sanitize_stem(tg), "__", sanitize_stem(sv), "__BALANCED_POOLED__R", r, ".rds")), compress = FALSE)
            kr <- kr + 1L
            manifest[[kr]] <- data.frame(target_gene = tg, subset_value = sv, sample_id = "BALANCED_POOLED", repeat_id = r, analysis_unit = "balanced_pooled", n_cells = ncol(mat), n_genes = nrow(mat), status = "COMPLETED", reason = NA_character_, elapsed_seconds = elapsed, stringsAsFactors = FALSE)
          }
        }
      }
    }
  }

  manifest_tab <- if (length(manifest)) do.call(rbind, manifest) else data.frame()
  write_csv_safe(manifest_tab, file.path(ctx$dirs$advanced, "virtual_knockout_run_manifest.csv"))
  if (!length(all_dr)) return(status_table("NOT_EVALUABLE", "no virtual-knockout run completed after scientific and expression gates"))
  dr_all <- do.call(rbind, all_dr)
  write_source_table(dr_all, file.path(ctx$dirs$advanced, "virtual_knockout_diffregulation_all.csv.gz"))
  ma_all <- if (length(all_ma)) do.call(rbind, all_ma) else data.frame()
  if (nrow(ma_all)) write_source_table(ma_all, file.path(ctx$dirs$advanced, "virtual_knockout_manifold_coordinates.csv.gz"))
  edge_all <- if (length(all_edges)) do.call(rbind, all_edges) else data.frame()
  if (nrow(edge_all)) write_source_table(edge_all, file.path(ctx$dirs$advanced, "virtual_knockout_network_edges_all.csv.gz"))

  sample_dr <- dr_all[dr_all$analysis_unit == "sample", , drop = FALSE]
  if (nrow(sample_dr)) {
    agg <- vk_aggregate_results(sample_dr, cfg$fdr_threshold %||% 0.05, cfg$minimum_repeat_recurrence %||% 0.5, cfg$minimum_sample_recurrence %||% 0.5)
    write_csv_safe(agg$sample, file.path(ctx$dirs$advanced, "virtual_knockout_sample_consensus.csv"))
    write_csv_safe(agg$consensus, file.path(ctx$dirs$advanced, "virtual_knockout_consensus_genes.csv"))
    consensus <- agg$consensus
    minimum_samples_required <- cfg$minimum_eligible_samples %||% 2L
    consensus$replicate_consensus_evaluable <- consensus$n_samples >= minimum_samples_required
    consensus$evidence_class[!consensus$replicate_consensus_evaluable & consensus$evidence_class == "Replicate-consistent"] <- "Single-sample-descriptive"
    write_csv_safe(consensus, file.path(ctx$dirs$advanced, "virtual_knockout_consensus_genes.csv"))
  } else {
    pooled <- dr_all[dr_all$analysis_unit == "balanced_pooled", , drop = FALSE]
    pooled$sample_id <- "BALANCED_POOLED"
    agg <- vk_aggregate_results(pooled, cfg$fdr_threshold %||% 0.05, cfg$minimum_repeat_recurrence %||% 0.5, 1)
    agg$consensus$evidence_class <- ifelse(agg$consensus$n_supported_samples > 0, "Pooled-descriptive", "Not-supported")
    agg$consensus$replicate_consensus_evaluable <- FALSE
    write_csv_safe(agg$sample, file.path(ctx$dirs$advanced, "virtual_knockout_sample_consensus.csv"))
    write_csv_safe(agg$consensus, file.path(ctx$dirs$advanced, "virtual_knockout_consensus_genes.csv"))
    consensus <- agg$consensus
  }

  if (nrow(edge_all)) {
    edge_cons <- vk_aggregate_edges(edge_all, cfg$minimum_edge_recurrence %||% 0.5)
    write_csv_safe(edge_cons, file.path(ctx$dirs$advanced, "virtual_knockout_consensus_network_edges.csv"))
  }
  sim <- vk_run_similarity(dr_all, cfg$similarity_top_genes %||% 100L)
  if (nrow(sim)) write_csv_safe(sim, file.path(ctx$dirs$advanced, "virtual_knockout_run_similarity.csv"))
  pathways <- vk_load_pathways(ctx, cfg)
  ora <- vk_pathway_ora(consensus, gene_universe$genes, pathways, cfg)
  if (nrow(ora)) write_csv_safe(ora, file.path(ctx$dirs$advanced, "virtual_knockout_pathway_enrichment.csv"))

  write_yaml(list(
    backend = backend,
    interpretation = "Exploratory GRN perturbation magnitude. scTenifoldKnk differential-regulation distances are unsigned and do not predict expression direction; WT network edges are displayed for topology context only.",
    outer_subsampling_fraction = subsampling_fraction,
    outer_subsampling_reason = "Independent cell subsets are drawn before each backend run so repeated stability checks are not dependent on backend-internal random-seed behavior.",
    analysis_unit = analysis_unit,
    target_genes = as.list(resolved_targets),
    subset_column = subset_col,
    subset_values = as.list(subset_values),
    baseline_conditions = as.list(baseline_conditions),
    gene_universe_size = length(gene_universe$genes),
    fdr_threshold = cfg$fdr_threshold %||% 0.05,
    minimum_repeat_recurrence = cfg$minimum_repeat_recurrence %||% 0.5,
    minimum_sample_recurrence = cfg$minimum_sample_recurrence %||% 0.5,
    required_validation = c("independent biological samples", "real perturbation data when available", "orthogonal molecular validation", "functional experiment")
  ), file.path(ctx$dirs$advanced, "virtual_knockout_parameters.yml"))

  n_completed <- sum(manifest_tab$status == "COMPLETED", na.rm = TRUE)
  n_failed <- sum(manifest_tab$status == "FAILED", na.rm = TRUE)
  reason <- paste("scTenifoldKnk completed", n_completed, "runs;", n_failed, "runs failed; outputs are exploratory and unsigned")
  out <- status_table("COMPLETED", reason, n_completed)
  out$backend <- backend
  out$evidence_level <- if (nrow(sample_dr) && length(unique(sample_dr$sample_id)) >= (cfg$minimum_eligible_samples %||% 2L)) "replicate_consensus_exploratory" else "pooled_or_single_sample_descriptive"
  out
}
