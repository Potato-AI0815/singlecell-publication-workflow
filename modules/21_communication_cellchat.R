# Gated, sample-aware CellChat workflow. Inference is descriptive and not causal.

run_cellchat_one <- function(obj, group_col, species, min_cells, cfg) {
  md <- obj[[]]
  data.input <- get_layer_safe(join_layers_safe(obj, "RNA"), "RNA", "data")
  meta <- md[, unique(c(group_col, "sample_id", "condition", "patient_id")), drop = FALSE]
  meta[[group_col]] <- droplevels(factor(as.character(meta[[group_col]])))
  cc <- CellChat::createCellChat(object = data.input, meta = meta, group.by = group_col)
  cc@DB <- if (tolower(species) == "mouse") CellChat::CellChatDB.mouse else CellChat::CellChatDB.human
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc, type = cfg$probability_type %||% "triMean",
                                    population.size = as_flag(cfg$population_size, FALSE),
                                    raw.use = TRUE)
  cc <- CellChat::filterCommunication(cc, min.cells = min_cells)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)
  cc <- tryCatch(CellChat::netAnalysis_computeCentrality(cc, slot.name = "netP"), error = function(e) cc)
  cc
}

standardize_cellchat_table <- function(tab, sample_id, condition, patient_id = NA_character_) {
  if (is.null(tab) || !nrow(tab)) return(data.frame())
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  aliases <- list(source = c("source", "source.cell"), target = c("target", "target.cell"),
                  ligand = c("ligand"), receptor = c("receptor"), interaction = c("interaction_name", "interaction"),
                  pathway = c("pathway_name", "pathway"), probability = c("prob", "probability"), p_value = c("pval", "p.value"))
  out <- data.frame(sample_id = rep(sample_id, nrow(tab)), patient_id = rep(patient_id, nrow(tab)), condition = rep(condition, nrow(tab)), stringsAsFactors = FALSE)
  for (nm in names(aliases)) {
    hit <- aliases[[nm]][aliases[[nm]] %in% names(tab)][1]
    out[[nm]] <- if (length(hit) && !is.na(hit)) tab[[hit]] else NA
  }
  out$probability <- suppressWarnings(as.numeric(out$probability))
  out$p_value <- suppressWarnings(as.numeric(out$p_value))
  out
}

communication_sample_summary <- function(interactions) {
  if (!nrow(interactions)) return(data.frame())
  interactions$pair <- paste(interactions$source, interactions$target, interactions$interaction, interactions$pathway, sep = "||")
  interactions$patient_key <- as.character(interactions$patient_id)
  interactions$patient_key[is.na(interactions$patient_key) | !nzchar(interactions$patient_key)] <- as.character(interactions$sample_id[is.na(interactions$patient_key) | !nzchar(interactions$patient_key)])
  agg <- aggregate(probability ~ sample_id + patient_key + condition + source + target + interaction + pathway,
                   interactions, sum, na.rm = TRUE)
  names(agg)[names(agg) == "patient_key"] <- "patient_id"
  agg
}

communication_condition_contrast <- function(sample_tab, reference, target, minimum_replicates = 3L, paired = FALSE) {
  if (!nrow(sample_tab) || any(!c(reference, target) %in% unique(sample_tab$condition))) return(data.frame())
  key_cols <- c("source", "target", "interaction", "pathway")
  keys <- unique(sample_tab[, key_cols, drop = FALSE])
  sample_cols <- intersect(c("sample_id", "patient_id", "condition"), names(sample_tab))
  samples <- unique(sample_tab[, sample_cols, drop = FALSE])
  grid <- merge(samples, keys, all = TRUE)
  full <- merge(grid, sample_tab, by = c(sample_cols, key_cols), all.x = TRUE)
  full$probability[is.na(full$probability)] <- 0
  spl <- split(full, interaction(full$source, full$target, full$interaction, full$pathway, drop = TRUE))
  rows <- lapply(spl, function(z) {
    a <- z[z$condition == reference, , drop = FALSE]
    b <- z[z$condition == target, , drop = FALSE]
    p <- NA_real_; delta <- mean(b$probability) - mean(a$probability); n_pairs <- NA_integer_
    if (isTRUE(paired) && "patient_id" %in% names(z) && any(!is.na(z$patient_id))) {
      aa <- aggregate(probability ~ patient_id, a, mean, na.rm = TRUE)
      bb <- aggregate(probability ~ patient_id, b, mean, na.rm = TRUE)
      names(aa)[2] <- "reference_probability"; names(bb)[2] <- "target_probability"
      mm <- merge(aa, bb, by = "patient_id")
      n_pairs <- nrow(mm)
      if (n_pairs >= minimum_replicates) p <- tryCatch(stats::wilcox.test(mm$target_probability, mm$reference_probability, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_)
      if (n_pairs) delta <- mean(mm$target_probability - mm$reference_probability, na.rm = TRUE)
    } else if (nrow(a) >= minimum_replicates && nrow(b) >= minimum_replicates) {
      p <- tryCatch(stats::wilcox.test(b$probability, a$probability, exact = FALSE)$p.value, error = function(e) NA_real_)
    }
    data.frame(source = z$source[1], target = z$target[1], interaction = z$interaction[1], pathway = z$pathway[1],
               reference = reference, target_condition = target, n_reference = nrow(a), n_target = nrow(b), n_pairs = n_pairs,
               mean_reference = mean(a$probability), mean_target = mean(b$probability), delta_probability = delta,
               p_value = p, paired = isTRUE(paired), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$FDR <- stats::p.adjust(out$p_value, method = "BH")
  out
}

run_communication_module <- function(ctx, object_path = file.path(ctx$dirs$objects, "06_annotated.rds")) {
  cfg <- ctx$cfg$advanced_modules$communication
  if (!as_flag(cfg$enabled, FALSE)) return(status_table("NOT_RUN", "disabled"))
  if (!requireNamespace("CellChat", quietly = TRUE)) return(status_table("NOT_RUN", "CellChat unavailable"))
  if (!file.exists(object_path)) return(status_table("NOT_EVALUABLE", "annotated object absent"))
  obj <- join_layers_safe(readRDS(object_path), "RNA")
  md <- obj[[]]
  group_col <- cfg$grouping_column %||% "cell_type_l1"
  if (!all(c("sample_id", group_col) %in% names(md))) return(status_table("NOT_EVALUABLE", "sample_id or grouping column absent"))
  excluded <- unlist(cfg$exclude_labels %||% c("Unresolved", "Ambiguous", "Doublet-like"))
  keep <- !as.character(md[[group_col]]) %in% excluded & !is.na(md[[group_col]])
  obj <- subset(obj, cells = rownames(md)[keep]); md <- obj[[]]
  min_cells <- cfg$minimum_cells_per_group_sample %||% 20L
  samples <- unique(as.character(md$sample_id))
  objects <- list(); tables <- list(); status <- list()
  for (s in samples) {
    cells <- rownames(md)[as.character(md$sample_id) == s]
    ms <- md[cells, , drop = FALSE]
    tab <- table(ms[[group_col]])
    valid_groups <- names(tab)[tab >= min_cells]
    if (length(valid_groups) < 2L) {
      status[[s]] <- data.frame(sample_id = s, status = "NOT_EVALUABLE", reason = "fewer than two supported cell groups", stringsAsFactors = FALSE); next
    }
    cells <- rownames(ms)[as.character(ms[[group_col]]) %in% valid_groups]
    sub <- subset(obj, cells = cells)
    cc <- tryCatch(run_cellchat_one(sub, group_col, ctx$cfg$input$species, min_cells, cfg),
                   error = function(e) {append_issue(ctx, "COMMUNICATION", "WARNING", conditionMessage(e), s); NULL})
    if (is.null(cc)) {status[[s]] <- data.frame(sample_id = s, status = "FAILED", reason = "CellChat error", stringsAsFactors = FALSE); next}
    objects[[s]] <- cc
    cond <- unique(stats::na.omit(as.character(ms$condition))); cond <- if (length(cond) == 1L) cond else "Mixed_or_unknown"
    tab2 <- tryCatch(CellChat::subsetCommunication(cc), error = function(e) data.frame())
    pid <- unique(stats::na.omit(as.character(ms$patient_id))); pid <- if (length(pid) == 1L) pid else NA_character_
    tables[[s]] <- standardize_cellchat_table(tab2, s, cond, pid)
    status[[s]] <- data.frame(sample_id = s, status = "COMPLETED", reason = paste(length(valid_groups), "supported groups"), stringsAsFactors = FALSE)
  }
  status_df <- do.call(rbind, status); write_csv_safe(status_df, file.path(ctx$dirs$advanced, "cellchat_sample_status.csv"))
  if (!length(objects)) return(status_table("NOT_EVALUABLE", "no sample passed CellChat gates"))
  saveRDS(objects, file.path(ctx$dirs$advanced, "cellchat_by_sample.rds"))
  nonempty <- Filter(function(x) !is.null(x) && nrow(x) > 0L, tables)
  interactions <- if (length(nonempty)) do.call(rbind, nonempty) else data.frame()
  if (!nrow(interactions)) return(status_table("NOT_EVALUABLE", "CellChat produced no supported interactions"))
  write_csv_safe(interactions, file.path(ctx$dirs$advanced, "cellchat_interactions_by_sample.csv"))
  sample_summary <- communication_sample_summary(interactions)
  write_csv_safe(sample_summary, file.path(ctx$dirs$advanced, "cellchat_sample_interaction_summary.csv"))

  ref <- ctx$cfg$comparisons$reference %||% NULL; target <- ctx$cfg$comparisons$target %||% NULL
  contrast <- if (!is.null(ref) && !is.null(target)) communication_condition_contrast(sample_summary, ref, target,
      ctx$cfg$comparisons$minimum_replicates_per_condition %||% 3L, as_flag(ctx$cfg$metadata$paired, FALSE)) else data.frame()
  if (nrow(contrast)) write_csv_safe(contrast, file.path(ctx$dirs$advanced, "cellchat_condition_contrast.csv"))

  edges <- aggregate(probability ~ condition + source + target, sample_summary, mean, na.rm = TRUE)
  names(edges)[names(edges) == "probability"] <- "weight"
  write_csv_safe(edges, file.path(ctx$dirs$advanced, "cellchat_network_edges.csv"))
  pathways <- aggregate(probability ~ condition + pathway + source + target, sample_summary, mean, na.rm = TRUE)
  write_csv_safe(pathways, file.path(ctx$dirs$advanced, "cellchat_pathway_edges.csv"))
  write_yaml(list(grouping_column = group_col, minimum_cells = min_cells, excluded_labels = excluded,
                  inference = "sample-wise CellChat; condition summaries use biological samples and are descriptive unless replicate gate passes"),
             file.path(ctx$dirs$advanced, "cellchat_parameters.yml"))
  status_table("COMPLETED", "sample-wise CellChat completed; probabilities are model-derived and non-causal", length(objects))
}
