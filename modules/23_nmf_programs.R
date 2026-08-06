# Sample-aware NMF program discovery.
# Preferred backend: GeneNMF multiNMF/getMetaPrograms. A deterministic NMF-package fallback preserves the same output contract.

cosine_matrix <- function(x) {
  x <- as.matrix(x); den <- sqrt(colSums(x^2)); den[den == 0] <- NA_real_
  z <- crossprod(x) / outer(den, den)
  z[!is.finite(z)] <- 0
  z
}

score_gene_programs <- function(obj, program_genes, assay = "RNA") {
  obj <- join_layers_safe(obj, assay)
  mat <- get_layer_safe(obj, assay, "data")
  out <- matrix(NA_real_, nrow = ncol(mat), ncol = length(program_genes), dimnames = list(colnames(mat), names(program_genes)))
  for (j in seq_along(program_genes)) {
    g <- intersect(program_genes[[j]], rownames(mat))
    if (length(g) >= 3L) out[, j] <- as.numeric(Matrix::colMeans(mat[g, , drop = FALSE]))
  }
  scale(out)
}

run_genenmf_backend <- function(obj_list, ranks, cfg) {
  if (!requireNamespace("GeneNMF", quietly = TRUE)) stop("GeneNMF unavailable")
  nmf_res <- GeneNMF::multiNMF(obj_list, k = ranks,
                               nfeatures = cfg$maximum_features %||% 2000L,
                               min.cells = cfg$minimum_cells_per_sample %||% 100L,
                               seed = cfg$seed %||% 20260804)
  mp <- GeneNMF::getMetaPrograms(
    nmf_res, nMP = cfg$n_meta_programs %||% 6L,
    specificity.weight = cfg$specificity_weight %||% 5,
    weight.explained = cfg$weight_explained %||% 0.5,
    max.genes = cfg$maximum_genes_per_program %||% 100L,
    metric = cfg$similarity_metric %||% "cosine",
    min.confidence = cfg$minimum_gene_confidence %||% 0.5
  )
  list(raw = nmf_res, meta = mp, backend = "GeneNMF")
}

extract_genenmf_contract <- function(ans) {
  mp <- ans$meta
  genes <- mp$metaprograms.genes %||% mp$meta.programs.genes %||% NULL
  if (is.null(genes)) stop("GeneNMF result did not expose metaprogram genes")
  if (is.data.frame(genes)) {
    candidates <- names(genes)[grepl("gene", names(genes), ignore.case = TRUE)]
    program_col <- names(genes)[grepl("program|meta", names(genes), ignore.case = TRUE)][1]
    gene_col <- candidates[1]
    if (is.na(program_col) || is.na(gene_col)) stop("Could not standardize GeneNMF metaprogram gene table")
    gene_list <- split(as.character(genes[[gene_col]]), as.character(genes[[program_col]]))
  } else if (is.list(genes)) gene_list <- lapply(genes, as.character) else stop("Unsupported GeneNMF metaprogram gene representation")
  names(gene_list) <- paste0("NMF_MP", seq_along(gene_list))
  similarity <- mp$programs.similarity %||% NULL
  metrics <- as.data.frame(mp$metaprograms.metrics %||% data.frame(), stringsAsFactors = FALSE)
  clusters <- mp$programs.clusters %||% NULL
  list(program_genes = gene_list, similarity = similarity, metrics = metrics, program_clusters = clusters)
}

run_nmf_fallback <- function(obj_list, ranks, cfg, features) {
  if (!requireNamespace("NMF", quietly = TRUE)) stop("Neither GeneNMF nor NMF is available")
  programs <- list(); diagnostics <- list(); fits <- list(); z <- 0L
  for (sn in names(obj_list)) {
    x <- get_layer_safe(join_layers_safe(obj_list[[sn]], "RNA"), "RNA", "data")
    g <- intersect(features, rownames(x)); x <- as.matrix(x[g, , drop = FALSE]); x[x < 0] <- 0
    for (k in ranks) {
      set.seed((cfg$seed %||% 20260804) + k)
      fit <- NMF::nmf(x, rank = k, nrun = cfg$n_runs %||% 10L, .options = "p")
      W <- NMF::basis(fit); H <- NMF::coef(fit)
      recon <- sum((x - W %*% H)^2) / pmax(sum(x^2), .Machine$double.eps)
      diagnostics[[paste(sn, k)]] <- data.frame(sample_id = sn, rank = k, metric = "relative_reconstruction_error", metric_value = recon, stringsAsFactors = FALSE)
      for (j in seq_len(ncol(W))) {z <- z + 1L; programs[[z]] <- W[, j] / pmax(sqrt(sum(W[, j]^2)), .Machine$double.eps); names(programs)[z] <- paste(sn, paste0("k", k), paste0("p", j), sep = "__")}
      fits[[paste(sn, k, sep = "__")]] <- fit
    }
  }
  all_genes <- unique(unlist(lapply(programs, names)))
  load <- matrix(0, nrow = length(all_genes), ncol = length(programs), dimnames = list(all_genes, names(programs)))
  for (j in seq_along(programs)) load[names(programs[[j]]), j] <- programs[[j]]
  sim <- cosine_matrix(load)
  nmp <- min(cfg$n_meta_programs %||% 6L, ncol(load))
  hc <- stats::hclust(stats::as.dist(1 - pmin(pmax(sim, 0), 1)), method = "ward.D2")
  cl <- stats::cutree(hc, k = nmp)
  meta_weights <- sapply(split(seq_len(ncol(load)), cl), function(ii) rowMeans(load[, ii, drop = FALSE]))
  colnames(meta_weights) <- paste0("NMF_MP", seq_len(ncol(meta_weights)))
  maxg <- cfg$maximum_genes_per_program %||% 100L
  gene_list <- lapply(seq_len(ncol(meta_weights)), function(j) names(sort(meta_weights[, j], decreasing = TRUE))[seq_len(min(maxg, sum(meta_weights[, j] > 0)))])
  names(gene_list) <- colnames(meta_weights)
  metrics <- do.call(rbind, lapply(seq_len(ncol(meta_weights)), function(j) {
    ii <- which(cl == j); ss <- unique(sub("__.*", "", colnames(load)[ii]))
    data.frame(program = colnames(meta_weights)[j], n_component_programs = length(ii), n_samples = length(ss),
               sample_fraction = length(ss) / length(obj_list), mean_similarity = if (length(ii) > 1L) mean(sim[ii, ii][upper.tri(sim[ii, ii])]) else 1,
               stringsAsFactors = FALSE)
  }))
  list(raw = fits, contract = list(program_genes = gene_list, similarity = sim, metrics = metrics, program_clusters = cl,
                                   loadings = load, meta_weights = meta_weights), diagnostics = do.call(rbind, diagnostics), backend = "NMF_fallback")
}

run_nmf_module <- function(ctx, object_path = file.path(ctx$dirs$objects, "06_annotated.rds")) {
  cfg <- ctx$cfg$advanced_modules$nmf
  if (!as_flag(cfg$enabled, FALSE)) return(status_table("NOT_RUN", "disabled"))
  if (!file.exists(object_path)) return(status_table("NOT_EVALUABLE", "annotated object absent"))
  obj <- join_layers_safe(readRDS(object_path), "RNA")
  md <- obj[[]]
  subset_col <- cfg$subset_column %||% NULL
  subset_values <- unlist(cfg$subset_values %||% list())
  if (!is.null(subset_col)) {
    if (!subset_col %in% names(md) || !length(subset_values)) return(status_table("NOT_EVALUABLE", "NMF subset_column requires explicit subset_values"))
    cells <- rownames(md)[as.character(md[[subset_col]]) %in% subset_values]
    obj <- subset(obj, cells = cells); md <- obj[[]]
  }
  if (!"sample_id" %in% names(md) || length(unique(md$sample_id)) < (cfg$minimum_samples %||% 2L)) return(status_table("NOT_EVALUABLE", "NMF requires multiple biological samples"))
  min_cells <- cfg$minimum_cells_per_sample %||% 100L
  valid_samples <- names(table(md$sample_id))[table(md$sample_id) >= min_cells]
  if (length(valid_samples) < (cfg$minimum_samples %||% 2L)) return(status_table("NOT_EVALUABLE", "too few samples pass the minimum cell gate"))
  obj <- subset(obj, cells = rownames(md)[md$sample_id %in% valid_samples]); md <- obj[[]]
  maxc <- cfg$maximum_cells_per_sample %||% 5000L
  set.seed(cfg$seed %||% ctx$cfg$runtime$seed %||% 20260804)
  keep <- unlist(lapply(split(rownames(md), md$sample_id), function(z) if (length(z) > maxc) sample(z, maxc) else z), use.names = FALSE)
  discovery <- subset(obj, cells = keep)
  obj_list <- Seurat::SplitObject(discovery, split.by = "sample_id")
  ranks <- as.integer(unlist(cfg$ranks %||% 4:9)); ranks <- ranks[ranks >= 2L]
  if (!length(ranks)) return(status_table("NOT_EVALUABLE", "no valid NMF ranks"))

  features <- Seurat::VariableFeatures(obj)
  if (length(features) < 500L) {
    tmp <- tryCatch(Seurat::FindVariableFeatures(obj, nfeatures = cfg$maximum_features %||% 2000L, verbose = FALSE), error = function(e) obj)
    features <- Seurat::VariableFeatures(tmp)
  }
  features <- intersect(features, rownames(obj))
  features <- head(features, cfg$maximum_features %||% 2000L)
  if (length(features) < (cfg$minimum_features %||% 500L)) return(status_table("NOT_EVALUABLE", "fewer than the minimum variable features"))

  backend <- tolower(cfg$backend %||% "auto")
  ans <- NULL
  if (backend %in% c("auto", "genenmf") && requireNamespace("GeneNMF", quietly = TRUE)) {
    ans <- tryCatch(run_genenmf_backend(obj_list, ranks, cfg), error = function(e) {append_issue(ctx, "NMF", "WARNING", paste("GeneNMF backend failed:", conditionMessage(e))); NULL})
  }
  if (!is.null(ans)) {
    contract <- extract_genenmf_contract(ans)
    raw <- ans$raw; diagnostics <- data.frame(); backend_used <- ans$backend
  } else {
    fb <- run_nmf_fallback(obj_list, ranks, cfg, features)
    contract <- fb$contract; raw <- fb$raw; diagnostics <- fb$diagnostics; backend_used <- fb$backend
  }
  saveRDS(raw, file.path(ctx$dirs$advanced, "nmf_backend_objects.rds"), compress = FALSE)
  saveRDS(contract, file.path(ctx$dirs$advanced, "nmf_program_contract.rds"))
  if (nrow(diagnostics)) write_csv_safe(diagnostics, file.path(ctx$dirs$advanced, "nmf_rank_diagnostics.csv"))

  program_genes <- contract$program_genes
  gene_rows <- do.call(rbind, lapply(names(program_genes), function(p) data.frame(program = p, rank = seq_along(program_genes[[p]]), gene = program_genes[[p]], stringsAsFactors = FALSE)))
  write_csv_safe(gene_rows, file.path(ctx$dirs$advanced, "nmf_meta_program_genes.csv"))
  if (!is.null(contract$metrics) && nrow(contract$metrics)) write_csv_safe(contract$metrics, file.path(ctx$dirs$advanced, "nmf_meta_program_metrics.csv"))
  if (!is.null(contract$similarity)) {
    sm <- as.matrix(contract$similarity)
    sim_long <- data.frame(program_1 = rep(rownames(sm), times = ncol(sm)), program_2 = rep(colnames(sm), each = nrow(sm)), similarity = as.vector(sm), stringsAsFactors = FALSE)
    write_source_table(sim_long, file.path(ctx$dirs$advanced, "nmf_program_similarity.csv.gz"))
  }
  if (!is.null(contract$meta_weights)) {
    mw <- contract$meta_weights
    weight_long <- data.frame(gene = rep(rownames(mw), times = ncol(mw)), program = rep(colnames(mw), each = nrow(mw)), weight = as.vector(mw), stringsAsFactors = FALSE)
    write_source_table(weight_long, file.path(ctx$dirs$advanced, "nmf_meta_program_gene_weights.csv.gz"))
  }

  scores <- score_gene_programs(obj, program_genes)
  score_names <- colnames(scores)
  for (p in score_names) obj[[p]] <- scores[colnames(obj), p]
  saveRDS(obj, file.path(ctx$dirs$advanced, "nmf_scored_object.rds"), compress = FALSE)
  cell_tab <- data.frame(cell_id = rownames(scores), scores, check.names = FALSE, stringsAsFactors = FALSE)
  for (nm in intersect(c("sample_id", "patient_id", "condition", "cell_type_l1", "cell_type_l2"), names(obj[[]]))) cell_tab[[nm]] <- obj[[]][rownames(scores), nm]
  write_source_table(cell_tab, file.path(ctx$dirs$advanced, "nmf_cell_program_scores.csv.gz"))
  group_col <- cfg$summary_group_column %||% "cell_type_l1"
  if (!group_col %in% names(cell_tab)) group_col <- "sample_id"
  summaries <- list()
  for (p in score_names) {
    sp <- split(cell_tab[[p]], cell_tab[[group_col]])
    summaries[[p]] <- do.call(rbind, lapply(names(sp), function(g) {
      v <- sp[[g]]; cut <- stats::quantile(cell_tab[[p]], 0.75, na.rm = TRUE)
      data.frame(program = p, group = g, mean_score = mean(v, na.rm = TRUE), median_score = stats::median(v, na.rm = TRUE), percent_high = 100 * mean(v >= cut, na.rm = TRUE), n_cells = sum(is.finite(v)), stringsAsFactors = FALSE)
    }))
  }
  write_csv_safe(do.call(rbind, summaries), file.path(ctx$dirs$advanced, "nmf_program_group_summary.csv"))
  write_yaml(list(backend = backend_used, ranks = ranks, n_meta_programs = length(program_genes), samples = valid_samples,
                  subset_column = subset_col, subset_values = subset_values,
                  interpretation = "recurrent transcriptional programs; names require marker/pathway review and external validation"),
             file.path(ctx$dirs$advanced, "nmf_parameters.yml"))
  status_table("COMPLETED", paste(backend_used, "identified", length(program_genes), "meta-programs"), length(program_genes))
}
