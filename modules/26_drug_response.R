# Exploratory drug-response module.
# Backends: transparent perturbation-signature reversal, gene-set score matrices, or oncoPredict on sample x cell-type pseudobulk expression.
# Outputs are hypothesis-generating and never constitute treatment recommendations.

read_drug_signature <- function(path) {
  if (!file.exists(path)) stop("drug signature file not found: ", path)
  x <- read_csv_safe(path)
  names(x) <- tolower(names(x))
  if (!all(c("drug", "gene") %in% names(x))) stop("drug signature requires drug and gene columns")
  x$drug <- as.character(x$drug); x$gene <- as.character(x$gene)
  x
}

signature_reversal_engine <- function(de, sig, min_overlap = 20L) {
  effect_col <- intersect(c("effect", "logfc", "z", "score", "drug_effect"), names(sig))[1]
  if (is.na(effect_col)) stop("signature reversal requires numeric effect/drug_effect column")
  sig$drug_effect <- suppressWarnings(as.numeric(sig[[effect_col]]))
  de$disease_effect <- suppressWarnings(as.numeric(de$logFC))
  de_groups <- split(de, as.character(de$cell_type))
  sig_groups <- split(sig, as.character(sig$drug))
  ranks <- list(); contrib <- list(); z <- 0L
  for (ct in names(de_groups)) for (drug in names(sig_groups)) {
    m <- merge(de_groups[[ct]][, c("gene", "disease_effect", "FDR"), drop = FALSE], sig_groups[[drug]][, c("gene", "drug_effect"), drop = FALSE], by = "gene")
    m <- m[is.finite(m$disease_effect) & is.finite(m$drug_effect), , drop = FALSE]
    if (nrow(m) < min_overlap) next
    ctst <- tryCatch(stats::cor.test(m$disease_effect, m$drug_effect, method = "spearman", exact = FALSE), error = function(e) NULL)
    if (is.null(ctst)) next
    z <- z + 1L
    ranks[[z]] <- data.frame(cell_type = ct, drug = drug, n_overlap = nrow(m), correlation = unname(ctst$estimate), reversal_score = -unname(ctst$estimate), p_value = ctst$p.value, stringsAsFactors = FALSE)
    m$cell_type <- ct; m$drug <- drug; m$reversal_contribution <- -m$disease_effect * m$drug_effect
    contrib[[z]] <- m
  }
  if (!length(ranks)) return(list(ranking = data.frame(), contribution = data.frame()))
  r <- do.call(rbind, ranks); r$FDR <- stats::p.adjust(r$p_value, method = "BH")
  list(ranking = r, contribution = do.call(rbind, contrib))
}

cell_signature_score_engine <- function(obj, sig, cfg) {
  dir_col <- intersect(c("direction", "effect", "sign"), names(sig))[1]
  if (is.na(dir_col)) stop("cell signature scoring requires direction/effect/sign column")
  sig$direction_num <- if (is.numeric(sig[[dir_col]])) sign(sig[[dir_col]]) else ifelse(tolower(as.character(sig[[dir_col]])) %in% c("up", "positive", "+", "1"), 1, -1)
  mat <- get_layer_safe(join_layers_safe(obj, "RNA"), "RNA", "data")
  drugs <- unique(sig$drug); scores <- matrix(NA_real_, nrow = ncol(mat), ncol = length(drugs), dimnames = list(colnames(mat), drugs))
  for (j in seq_along(drugs)) {
    ss <- sig[sig$drug == drugs[j], , drop = FALSE]
    up <- intersect(ss$gene[ss$direction_num > 0], rownames(mat)); dn <- intersect(ss$gene[ss$direction_num < 0], rownames(mat))
    if (length(c(up, dn)) < (cfg$minimum_signature_genes %||% 5L)) next
    scores[, j] <- as.numeric((if (length(up)) Matrix::colMeans(mat[up, , drop = FALSE]) else 0) - (if (length(dn)) Matrix::colMeans(mat[dn, , drop = FALSE]) else 0))
  }
  scores <- scale(scores)
  data.frame(cell_id = rownames(scores), scores, check.names = FALSE, stringsAsFactors = FALSE)
}

aggregate_drug_scores <- function(cell_scores, md) {
  drug_cols <- setdiff(names(cell_scores), "cell_id")
  common <- intersect(cell_scores$cell_id, rownames(md)); x <- cell_scores[match(common, cell_scores$cell_id), , drop = FALSE]; md <- md[common, , drop = FALSE]
  rows <- list(); z <- 0L
  key <- paste(md$sample_id, md$cell_type_l1, sep = "|||")
  for (k in unique(key)) {
    ii <- which(key == k); parts <- strsplit(k, "\\|\\|\\|")[[1]]
    for (d in drug_cols) {z <- z + 1L; rows[[z]] <- data.frame(sample_id = parts[1], cell_type = parts[2], condition = md$condition[ii[1]], drug = d,
                                                               score = mean(x[[d]][ii], na.rm = TRUE), n_cells = length(ii), stringsAsFactors = FALSE)}
  }
  do.call(rbind, rows)
}

run_oncopredict_engine <- function(ctx, cfg) {
  if (!requireNamespace("oncoPredict", quietly = TRUE)) return(list(status = "NOT_RUN", reason = "oncoPredict unavailable"))
  req <- c(cfg$training_expression_path, cfg$training_response_path)
  if (any(vapply(req, is.null, logical(1))) || any(!file.exists(req))) return(list(status = "NOT_EVALUABLE", reason = "oncoPredict training expression/response files are required"))
  pb_path <- file.path(ctx$dirs$de, "pseudobulk_raw_counts.rds"); meta_path <- file.path(ctx$dirs$de, "pseudobulk_sample_metadata.csv")
  if (!file.exists(pb_path) || !file.exists(meta_path)) return(list(status = "NOT_EVALUABLE", reason = "pseudobulk counts and metadata are required"))
  pb <- readRDS(pb_path); pm <- read_csv_safe(meta_path)
  test <- log2(as.matrix(pb) + 1)
  train_expr <- as.matrix(utils::read.csv(cfg$training_expression_path, row.names = 1, check.names = FALSE))
  train_resp <- as.matrix(utils::read.csv(cfg$training_response_path, row.names = 1, check.names = FALSE))
  pred <- oncoPredict::calcPhenotype(trainingExprData = train_expr, trainingPtype = train_resp, testExprData = test,
                                     batchCorrect = cfg$batch_correct %||% "eb", powerTransformPhenotype = as_flag(cfg$power_transform, TRUE),
                                     removeLowVaryingGenes = cfg$remove_low_varying_genes %||% 0.2,
                                     minNumSamples = cfg$minimum_training_samples %||% 10L, printOutput = FALSE)
  pred <- as.matrix(pred)
  long <- data.frame(pb_id = rep(rownames(pred), times = ncol(pred)), drug = rep(colnames(pred), each = nrow(pred)), predicted_sensitivity = as.vector(pred), stringsAsFactors = FALSE)
  long <- merge(long, pm, by = "pb_id", all.x = TRUE)
  list(status = "COMPLETED", predictions = long)
}

run_drug_response_module <- function(ctx, object_path = file.path(ctx$dirs$objects, "06_annotated.rds")) {
  cfg <- ctx$cfg$advanced_modules$drug_response
  if (!as_flag(cfg$enabled, FALSE)) return(status_table("NOT_RUN", "disabled"))
  backend <- tolower(cfg$backend %||% "signature_reversal")
  if (backend == "signature_reversal") {
    de_path <- cfg$differential_expression_path %||% file.path(ctx$dirs$de, "pseudobulk_edgeR_all_results.csv")
    if (!file.exists(de_path) || is.null(cfg$signature_path)) return(status_table("NOT_EVALUABLE", "pseudobulk DE and drug signature_path are required"))
    de <- read_csv_safe(de_path); sig <- read_drug_signature(cfg$signature_path)
    if (!all(c("gene", "logFC", "cell_type") %in% names(de))) return(status_table("NOT_EVALUABLE", "DE table lacks gene/logFC/cell_type"))
    ans <- signature_reversal_engine(de, sig, cfg$minimum_gene_overlap %||% 20L)
    if (!nrow(ans$ranking)) return(status_table("NOT_EVALUABLE", "no drug-cell-type pair passed the overlap gate"))
    write_csv_safe(ans$ranking, file.path(ctx$dirs$advanced, "drug_reversal_ranking.csv"))
    write_source_table(ans$contribution, file.path(ctx$dirs$advanced, "drug_reversal_gene_contributions.csv.gz"))
    n <- nrow(ans$ranking)
  } else if (backend == "signature_score") {
    if (!file.exists(object_path) || is.null(cfg$signature_path)) return(status_table("NOT_EVALUABLE", "annotated object and signature_path are required"))
    obj <- readRDS(object_path); sig <- read_drug_signature(cfg$signature_path)
    scores <- cell_signature_score_engine(obj, sig, cfg)
    write_source_table(scores, file.path(ctx$dirs$advanced, "drug_cell_scores.csv.gz"))
    ag <- aggregate_drug_scores(scores, obj[[]]); write_csv_safe(ag, file.path(ctx$dirs$advanced, "drug_sample_celltype_scores.csv"))
    top_drugs <- names(sort(colMeans(abs(as.matrix(scores[, setdiff(names(scores), "cell_id"), drop = FALSE])), na.rm = TRUE), decreasing = TRUE))
    top_drugs <- head(top_drugs, cfg$therapeutic_cluster_drugs %||% 30L)
    if (length(top_drugs) >= 3L && nrow(scores) >= (cfg$minimum_cells_for_therapeutic_clusters %||% 100L)) {
      x <- as.matrix(scores[, top_drugs, drop = FALSE]); x[!is.finite(x)] <- 0
      k <- min(cfg$therapeutic_clusters %||% 5L, max(2L, floor(sqrt(nrow(x) / 50))))
      set.seed(ctx$cfg$runtime$seed %||% 20260804); km <- stats::kmeans(x, centers = k, nstart = 20)
      tc <- data.frame(cell_id = scores$cell_id, therapeutic_cluster = paste0("TC", km$cluster), stringsAsFactors = FALSE)
      write_csv_safe(tc, file.path(ctx$dirs$advanced, "drug_therapeutic_clusters.csv"))
    }
    n <- nrow(ag)
  } else if (backend == "oncopredict") {
    ans <- tryCatch(run_oncopredict_engine(ctx, cfg), error = function(e) list(status = "FAILED", reason = conditionMessage(e)))
    if (!identical(ans$status, "COMPLETED")) return(status_table(ans$status, ans$reason))
    write_csv_safe(ans$predictions, file.path(ctx$dirs$advanced, "oncoPredict_pseudobulk_predictions.csv")); n <- nrow(ans$predictions)
  } else return(status_table("NOT_EVALUABLE", paste("unsupported drug-response backend", backend)))

  write_yaml(list(backend = backend, signature_path = cfg$signature_path %||% NULL,
                  interpretation = "exploratory transcriptome-based hypothesis generation; not clinical efficacy or treatment recommendation",
                  required_validation = c("independent dataset", "orthogonal pharmacology", "clinical context and exposure review")),
             file.path(ctx$dirs$advanced, "drug_response_parameters.yml"))
  status_table("COMPLETED", paste("exploratory drug-response backend completed:", backend), n)
}
