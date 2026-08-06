select_publication_marker_genes <- function(ctx, obj, maximum = 36L) {
  configured <- unlist(ctx$cfg$figures$marker_genes %||% list())
  configured <- configured[!is.na(configured) & nzchar(configured)]
  if (length(configured) && !identical(configured[[1]], "auto")) return(head(unique(configured), maximum))
  labels <- unique(as.character(obj$cell_type_l1))
  sp <- tolower(ctx$cfg$input$species %||% "human")
  dict_path <- if (sp == "mouse") file.path(ctx$skill_root, "resources", "marker_dictionary_mouse_l1.csv") else file.path(ctx$skill_root, "resources", "marker_dictionary_human_l1.csv")
  genes <- character()
  if (file.exists(dict_path)) {
    d <- read_csv_safe(dict_path)
    d <- d[d$label %in% labels, , drop = FALSE]
    if (nrow(d)) genes <- unique(unlist(lapply(d$positive, function(x) head(split_genes(x), 3L))))
  }
  top_path <- file.path(ctx$dirs$source, "top_cluster_markers.csv")
  if (length(genes) < 8L && file.exists(top_path)) {
    top <- read_csv_safe(top_path)
    if ("gene" %in% names(top)) genes <- unique(c(genes, top$gene))
  }
  head(intersect(genes, rownames(obj)), maximum)
}

figure_status_row <- function(id, type, status, reason = NA_character_) {
  data.frame(figure_id = id, type = type, status = status, reason = as.character(reason), stringsAsFactors = FALSE)
}

try_export_figure <- function(ctx, id, type, expr, collection = "main", width_mm = NULL, height_mm = NULL) {
  ans <- tryCatch(force(expr), error = function(e) e)
  if (inherits(ans, "error")) {
    append_issue(ctx, "11", "WARNING", paste0("Figure ", id, " not generated: ", conditionMessage(ans)), type)
    return(figure_status_row(id, type, "NOT_EVALUABLE", conditionMessage(ans)))
  }
  ans$parameters <- ans$parameters %||% list()
  ans$parameters$template_type <- ans$parameters$type %||% NA_character_
  ans$parameters$type <- type
  exp <- tryCatch(export_publication_figure(ctx, ans, id, collection, width_mm, height_mm), error = function(e) e)
  if (inherits(exp, "error")) {
    append_issue(ctx, "11", "WARNING", paste0("Figure ", id, " export failed: ", conditionMessage(exp)), type)
    return(figure_status_row(id, type, "FAILED", conditionMessage(exp)))
  }
  export_status <- if ("status" %in% names(exp)) as.character(exp$status[1]) else "EXPORTED"
  export_reason <- if ("reason" %in% names(exp)) as.character(exp$reason[1]) else NA_character_
  figure_status_row(id, type, export_status, export_reason)
}

stage_main <- function(ctx) {
  if (!as_flag(ctx$cfg$figures$export_individual_figures, TRUE)) stop("Publication Figure Engine alpha requires figures.export_individual_figures: true")
  if (!as_flag(ctx$cfg$figures$composite_figures_forbidden, TRUE)) stop("Alpha release requires figures.composite_figures_forbidden: true")
  if (as_flag(ctx$cfg$figures$build_contact_sheet, FALSE)) stop("Alpha release forbids contact sheets; set figures.build_contact_sheet: false")
  objpath <- file.path(ctx$dirs$objects, "06_annotated.rds")
  if (!file.exists(objpath)) stop("Annotated object missing for figure engine.")
  obj <- readRDS(objpath)
  reduction <- ctx$cfg$figures$reduction %||% "auto"
  statuses <- list(); ks <- 0L
  add_status <- function(x) { ks <<- ks + 1L; statuses[[ks]] <<- x; invisible(x) }

  # 1. Embeddings: every plot is exported independently with its own Source Data and parameters.
  add_status(try_export_figure(ctx, "Embedding_01_clusters", "embedding_discrete",
    plot_embedding_discrete(ctx, obj, "cluster_raw", reduction, "Unsupervised clusters",
      labels = TRUE, hulls = as_flag(ctx$cfg$figures$cluster_hulls, TRUE)), "main", 183, 155))
  add_status(try_export_figure(ctx, "Embedding_02_cell_types", "embedding_discrete",
    plot_embedding_discrete(ctx, obj, "cell_type_l1", reduction, "Provisional broad cell types",
      labels = TRUE, hulls = as_flag(ctx$cfg$figures$celltype_hulls, TRUE)), "main", 183, 155))

  if ("condition" %in% colnames(obj[[]]) && length(unique(stats::na.omit(obj$condition))) > 1L) {
    add_status(try_export_figure(ctx, "Embedding_03_condition", "embedding_discrete",
      plot_embedding_discrete(ctx, obj, "condition", reduction, "Cells by biological condition", labels = FALSE, hulls = FALSE), "main", 183, 155))
  } else add_status(figure_status_row("Embedding_03_condition", "embedding_discrete", "NOT_EVALUABLE", "condition has fewer than two observed levels"))

  if ("sample_id" %in% colnames(obj[[]]) && length(unique(obj$sample_id)) <= (ctx$cfg$figures$max_sample_colours %||% 16)) {
    add_status(try_export_figure(ctx, "Embedding_04_samples", "embedding_discrete",
      plot_embedding_discrete(ctx, obj, "sample_id", reduction, "Cells by sample", labels = FALSE, hulls = FALSE), "extended", 183, 155))
  } else add_status(figure_status_row("Embedding_04_samples", "embedding_discrete", "NOT_EVALUABLE", "too many samples for a legible categorical embedding"))

  if ("annotation_confidence" %in% colnames(obj[[]])) {
    conf_pal <- c(High = "#2A9D8F", Medium = "#73BFB8", Low = "#E9C46A", Ambiguous = "#F4A261", `Doublet-like` = "#9C6644", Unresolved = "#B8B8B8")
    add_status(try_export_figure(ctx, "Embedding_05_annotation_confidence", "embedding_discrete",
      plot_embedding_discrete(ctx, obj, "annotation_confidence", reduction, "Annotation confidence", labels = FALSE, hulls = FALSE, palette = conf_pal), "extended", 183, 155))
  }

  # 2. Marker evidence.
  marker_genes <- select_publication_marker_genes(ctx, obj, maximum = ctx$cfg$figures$maximum_marker_genes %||% 36L)
  if (length(marker_genes)) {
    dot_h <- max(100, min(250, 5.0 * length(marker_genes)))
    add_status(try_export_figure(ctx, "Markers_01_dotplot", "marker_dotplot",
      plot_marker_dotplot(ctx, obj, marker_genes, "cell_type_l1", "Canonical marker evidence"), "main", 210, dot_h))
    add_status(try_export_figure(ctx, "Markers_02_heatmap", "marker_heatmap",
      plot_marker_heatmap(ctx, obj, marker_genes, "cell_type_l1", "Average marker expression"), "extended", 210, dot_h))
  } else {
    add_status(figure_status_row("Markers_01_dotplot", "marker_dotplot", "NOT_EVALUABLE", "no configured or automatically selected marker genes were present"))
    add_status(figure_status_row("Markers_02_heatmap", "marker_heatmap", "NOT_EVALUABLE", "no configured or automatically selected marker genes were present"))
  }

  # 3. Continuous feature embeddings and split violins.
  continuous <- unlist(ctx$cfg$figures$continuous_features %||% list())
  if (length(continuous) && identical(continuous[[1]], "auto")) continuous <- head(marker_genes, ctx$cfg$figures$maximum_continuous_features %||% 4L)
  continuous <- unique(continuous[!is.na(continuous) & nzchar(continuous)])
  for (feature in continuous) {
    id <- paste0("Feature_", sanitize_stem(feature), "_embedding")
    add_status(try_export_figure(ctx, id, "embedding_continuous",
      plot_embedding_continuous(ctx, obj, feature, reduction, paste0(feature, " expression")), "extended", 150, 135))
  }

  split_features <- unlist(ctx$cfg$figures$split_violin_features %||% list())
  if (length(split_features) && identical(split_features[[1]], "auto")) split_features <- head(marker_genes, ctx$cfg$figures$maximum_split_violin_features %||% 4L)
  split_features <- unique(split_features[!is.na(split_features) & nzchar(split_features)])
  cond_levels <- if ("condition" %in% colnames(obj[[]])) unique(stats::na.omit(as.character(obj$condition))) else character()
  if (length(split_features) && length(cond_levels) == 2L) {
    md <- obj[[]]
    for (feature in split_features) {
      vals <- tryCatch(get_feature_vector(obj, feature), error = function(e) NULL)
      if (is.null(vals)) {
        add_status(figure_status_row(paste0("Violin_", sanitize_stem(feature)), "split_violin_box", "NOT_EVALUABLE", "feature absent")); next
      }
      dd <- data.frame(cell_type_l1 = md$cell_type_l1, condition = md$condition, value = vals, stringsAsFactors = FALSE)
      # Restrict to configured confidence levels when available, but preserve all source values in the exported table.
      keep_types <- names(sort(table(dd$cell_type_l1), decreasing = TRUE))[seq_len(min(length(unique(dd$cell_type_l1)), ctx$cfg$figures$maximum_violin_celltypes %||% 12L))]
      dd <- dd[dd$cell_type_l1 %in% keep_types, , drop = FALSE]
      add_status(try_export_figure(ctx, paste0("Violin_", sanitize_stem(feature)), "split_violin_box",
        plot_split_violin_box(ctx, dd, "cell_type_l1", "condition", "value", paste0(feature, " by cell type and condition"), feature), "extended", 210, 135))
    }
  } else if (length(split_features)) {
    for (feature in split_features) add_status(figure_status_row(paste0("Violin_", sanitize_stem(feature)), "split_violin_box", "NOT_EVALUABLE", "exactly two condition levels are required"))
  }

  # 4. Composition figures.
  comp_path <- file.path(ctx$dirs$source, "sample_celltype_composition.csv")
  if (file.exists(comp_path)) {
    comp <- read_csv_safe(comp_path)
    add_status(try_export_figure(ctx, "Composition_01_stacked_samples", "composition_stacked",
      plot_composition_stacked(ctx, comp), "main", 210, 135))
    if ("condition" %in% names(comp) && length(unique(stats::na.omit(comp$condition))) > 1L) {
      add_status(try_export_figure(ctx, "Composition_02_by_condition", "composition_by_condition",
        plot_composition_by_condition(ctx, comp), "extended", 297, 210))
    }
    if (as_flag(ctx$cfg$metadata$paired, FALSE) && "patient_id" %in% names(comp)) {
      add_status(try_export_figure(ctx, "Composition_03_paired", "paired_composition",
        plot_paired_composition(ctx, comp), "extended", 297, 210))
    }
  } else add_status(figure_status_row("Composition_01_stacked_samples", "composition_stacked", "NOT_EVALUABLE", "composition Source Data absent"))

  # 5. Pseudobulk DE: one independent publication figure per evaluable cell type.
  de_path <- file.path(ctx$dirs$de, "pseudobulk_edgeR_all_results.csv")
  if (file.exists(de_path)) {
    de <- read_csv_safe(de_path)
    best <- aggregate(FDR ~ cell_type, de, min, na.rm = TRUE)
    best <- best[order(best$FDR), , drop = FALSE]
    celltypes <- head(best$cell_type, ctx$cfg$figures$maximum_de_celltypes %||% 12L)
    for (ct in celltypes) {
      add_status(try_export_figure(ctx, paste0("DE_volcano_", sanitize_stem(ct)), "volcano",
        plot_volcano_single(ctx, de, ct, label_n = ctx$cfg$figures$volcano_label_n %||% 10L), "extended", 150, 135))
    }
  }

  # 6. Enrichment: one independent figure per cell type.
  fg_path <- file.path(ctx$dirs$pathways, "fgsea_all_results.csv")
  if (file.exists(fg_path)) {
    fg <- read_csv_safe(fg_path)
    best <- aggregate(padj ~ cell_type, fg, min, na.rm = TRUE)
    best <- best[order(best$padj), , drop = FALSE]
    celltypes <- head(best$cell_type, ctx$cfg$figures$maximum_enrichment_celltypes %||% 12L)
    for (ct in celltypes) {
      add_status(try_export_figure(ctx, paste0("Enrichment_", sanitize_stem(ct)), "enrichment_dotplot",
        plot_enrichment_single(ctx, fg, ct, top_n = ctx$cfg$figures$enrichment_top_n %||% 15L), "extended", 183, 150))
    }
  }

  # 7. Fitted trajectory and pseudotime trends. No hand-drawn arrows are allowed.
  traj_cells <- file.path(ctx$dirs$advanced, "trajectory_cell_data.csv")
  traj_curve <- file.path(ctx$dirs$advanced, "trajectory_curve_data.csv")
  traj_nodes <- file.path(ctx$dirs$advanced, "trajectory_node_data.csv")
  if (file.exists(traj_cells) && file.exists(traj_curve)) {
    tc <- read_csv_safe(traj_cells); tr <- read_csv_safe(traj_curve); tn <- if (file.exists(traj_nodes)) read_csv_safe(traj_nodes) else NULL
    add_status(try_export_figure(ctx, "Trajectory_01_embedding", "trajectory_embedding",
      plot_trajectory_embedding(ctx, tc, tr, tn, "Fitted lineage trajectory"), "main", 183, 155))
    trend_cols <- unlist(ctx$cfg$figures$pseudotime_trend_columns %||% list())
    trend_cols <- trend_cols[trend_cols %in% names(tc)]
    group_col <- ctx$cfg$figures$pseudotime_group_column %||% NULL
    if (!is.null(group_col) && !group_col %in% names(tc)) group_col <- NULL
    for (v in trend_cols) {
      add_status(try_export_figure(ctx, paste0("Trajectory_trend_", sanitize_stem(v)), "pseudotime_trend",
        plot_pseudotime_trend(ctx, tc, "pseudotime", v, group_col, paste0(v, " along fitted pseudotime")), "extended", 170, 125))
    }
  } else {
    add_status(figure_status_row("Trajectory_01_embedding", "trajectory_embedding", "NOT_EVALUABLE", "validated trajectory Source Data were not produced"))
  }
  knot_path <- file.path(ctx$dirs$advanced, "tradeSeq_knot_diagnostics.csv")
  if (file.exists(knot_path)) {
    kd <- read_csv_safe(knot_path)
    add_status(try_export_figure(ctx, "Trajectory_02_tradeSeq_knots", "tradeseq_knot_diagnostics",
      plot_metric_curve(ctx, kd, "knots", "median_AIC", "tradeSeq knot selection",
                        x_label = "Number of knots", y_label = "Median AIC", lower_is_better = TRUE),
      "extended", 150, 125))
  }
  trade_tests <- list(
    association = c("tradeSeq_association_test.csv", "Genes associated with fitted pseudotime"),
    start_vs_end = c("tradeSeq_start_vs_end_test.csv", "Genes differing between trajectory start and end"),
    pattern = c("tradeSeq_pattern_test.csv", "Genes with lineage-dependent expression patterns"),
    differential_end = c("tradeSeq_differential_end_test.csv", "Genes differing between lineage endpoints")
  )
  for (nm in names(trade_tests)) {
    fp <- file.path(ctx$dirs$advanced, trade_tests[[nm]][1])
    if (file.exists(fp)) {
      tt <- read_csv_safe(fp)
      if (all(c("gene", "FDR") %in% names(tt))) add_status(try_export_figure(ctx,
        paste0("Trajectory_tradeSeq_", nm), "tradeseq_gene_ranking",
        plot_gene_test_ranking(ctx, tt, trade_tests[[nm]][2], top_n = ctx$cfg$figures$tradeSeq_top_genes %||% 30L),
        "extended", 160, 155))
    }
  }

  # 8. Expression-inferred CNA. Every sample receives an independent heatmap.
  cnv_cells_path <- file.path(ctx$dirs$advanced, "cnv_cell_scores.csv")
  cnv_heat_path <- file.path(ctx$dirs$advanced, "cnv_heatmap_long.csv.gz")
  if (file.exists(cnv_cells_path)) {
    cc <- read_csv_safe(cnv_cells_path)
    if (all(c("cna_correlation", "cna_signal", "cnv_class") %in% names(cc))) {
      add_status(try_export_figure(ctx, "CNV_01_signal_scatter", "cnv_scatter",
        plot_cnv_scatter(ctx, cc), "extended", 150, 135))
    }
  }
  if (file.exists(cnv_heat_path)) {
    ch <- read_csv_safe(cnv_heat_path)
    samples <- head(unique(as.character(ch$sample_id)), ctx$cfg$figures$maximum_cnv_samples %||% 12L)
    for (sn in samples) {
      dd <- ch[ch$sample_id == sn, , drop = FALSE]
      add_status(try_export_figure(ctx, paste0("CNV_heatmap_", sanitize_stem(sn)), "cnv_heatmap",
        plot_cnv_heatmap(ctx, dd, title = paste("Expression-inferred CNA:", sn)), "extended", 183, 150))
    }
  }

  # 9. Sample-aware CellChat outputs. Pooled cells are not treated as biological replicates.
  comm_path <- file.path(ctx$dirs$advanced, "cellchat_sample_interaction_summary.csv")
  edge_path <- file.path(ctx$dirs$advanced, "cellchat_network_edges.csv")
  if (file.exists(comm_path)) {
    cm <- read_csv_safe(comm_path)
    for (cond in head(unique(as.character(cm$condition)), ctx$cfg$figures$maximum_communication_conditions %||% 6L)) {
      dd <- cm[cm$condition == cond, , drop = FALSE]
      add_status(try_export_figure(ctx, paste0("Communication_bubble_", sanitize_stem(cond)), "communication_bubble",
        plot_communication_bubble(ctx, dd, title = paste("Cell-cell communication:", cond),
                                  top_n = ctx$cfg$figures$communication_top_interactions %||% 40L), "extended", 210, 160))
    }
  }
  if (file.exists(edge_path)) {
    ee <- read_csv_safe(edge_path)
    for (cond in head(unique(as.character(ee$condition)), ctx$cfg$figures$maximum_communication_conditions %||% 6L)) {
      dd <- ee[ee$condition == cond, , drop = FALSE]
      add_status(try_export_figure(ctx, paste0("Communication_network_", sanitize_stem(cond)), "communication_network",
        plot_network_edges(ctx, dd, title = paste("Communication network:", cond),
                           top_n = ctx$cfg$figures$communication_top_edges %||% 80L), "extended", 170, 155))
    }
  }
  contrast_path <- file.path(ctx$dirs$advanced, "cellchat_condition_contrast.csv")
  if (file.exists(contrast_path)) {
    cx <- read_csv_safe(contrast_path)
    cx$label <- paste(cx$source, "→", cx$target, cx$interaction, sep = " | ")
    add_status(try_export_figure(ctx, "Communication_condition_changes", "communication_contrast",
      plot_lollipop_ranking(ctx, cx, "label", "delta_probability", "Sample-level communication changes",
                            top_n = ctx$cfg$figures$communication_contrast_top_n %||% 30L, x_label = "Mean probability difference"),
      "extended", 183, 170))
  }

  # 10. Recurrent NMF meta-programs.
  nmf_diag <- file.path(ctx$dirs$advanced, "nmf_rank_diagnostics.csv")
  if (file.exists(nmf_diag)) add_status(try_export_figure(ctx, "NMF_01_rank_diagnostics", "nmf_rank_diagnostics",
    plot_rank_diagnostics(ctx, read_csv_safe(nmf_diag), title = "NMF rank diagnostics"), "extended", 160, 125))
  nmf_sim <- file.path(ctx$dirs$advanced, "nmf_program_similarity.csv.gz")
  if (file.exists(nmf_sim)) add_status(try_export_figure(ctx, "NMF_02_program_similarity", "nmf_similarity_heatmap",
    plot_advanced_heatmap(ctx, read_csv_safe(nmf_sim), "program_1", "program_2", "similarity",
                          title = "Similarity among recurrent NMF programs", midpoint = 0.5, diverging = FALSE,
                          cluster_x = TRUE, cluster_y = TRUE, fill_label = "Cosine similarity"), "extended", 183, 170))
  nmf_genes <- file.path(ctx$dirs$advanced, "nmf_meta_program_genes.csv")
  if (file.exists(nmf_genes)) {
    ng <- read_csv_safe(nmf_genes); ng$gene_weight_display <- 1 / pmax(ng$rank, 1)
    ng <- do.call(rbind, lapply(split(ng, ng$program), head, n = ctx$cfg$figures$nmf_top_genes_per_program %||% 20L))
    add_status(try_export_figure(ctx, "NMF_03_meta_program_genes", "nmf_gene_heatmap",
      plot_advanced_heatmap(ctx, ng, "program", "gene", "gene_weight_display", title = "Top genes defining NMF meta-programs",
                            midpoint = 0, diverging = FALSE, cluster_x = FALSE, cluster_y = FALSE, fill_label = "Reciprocal rank"),
      "extended", 183, max(130, min(250, nrow(ng) * 1.8))))
  }
  nmf_group <- file.path(ctx$dirs$advanced, "nmf_program_group_summary.csv")
  if (file.exists(nmf_group)) add_status(try_export_figure(ctx, "NMF_04_program_activity", "nmf_program_dotplot",
    plot_program_prevalence_dotplot(ctx, read_csv_safe(nmf_group), title = "NMF program activity across groups"), "extended", 183, 145))
  nmf_obj_path <- file.path(ctx$dirs$advanced, "nmf_scored_object.rds")
  if (file.exists(nmf_obj_path)) {
    no <- readRDS(nmf_obj_path); progs <- grep("^NMF_MP", names(no[[]]), value = TRUE)
    for (pr in head(progs, ctx$cfg$figures$maximum_nmf_program_embeddings %||% 12L)) {
      add_status(try_export_figure(ctx, paste0("NMF_embedding_", sanitize_stem(pr)), "nmf_program_embedding",
        plot_embedding_continuous(ctx, no, pr, reduction, paste(pr, "activity")), "extended", 150, 135))
    }
  }

  # 11. hdWGCNA co-expression modules.
  hp <- file.path(ctx$dirs$advanced, "hdWGCNA_soft_power_table.csv")
  if (file.exists(hp)) {
    hd <- read_csv_safe(hp)
    pcol <- names(hd)[tolower(names(hd)) %in% c("power", "softpower", "soft_power")][1]
    rcol <- names(hd)[grepl("sft.*r|signed.*r|scale.*fit", names(hd), ignore.case = TRUE)][1]
    if (!is.na(pcol) && !is.na(rcol)) add_status(try_export_figure(ctx, "hdWGCNA_01_soft_power", "hdwgcna_soft_power",
      plot_soft_power_metric(ctx, hd, pcol, rcol, ctx$cfg$advanced_modules$hdWGCNA$scale_free_threshold %||% 0.8), "extended", 150, 125))
    kcol <- names(hd)[grepl("mean.*connect|connectivity", names(hd), ignore.case = TRUE)][1]
    if (!is.na(pcol) && !is.na(kcol)) add_status(try_export_figure(ctx, "hdWGCNA_02_connectivity", "hdwgcna_connectivity",
      plot_soft_power_metric(ctx, hd, pcol, kcol, threshold = 0, title = "Mean network connectivity"), "extended", 150, 125))
  }
  hh <- file.path(ctx$dirs$advanced, "hdWGCNA_hub_genes.csv")
  if (file.exists(hh)) {
    hubs <- read_csv_safe(hh); kcol <- names(hubs)[grepl("^kme$|kme", names(hubs), ignore.case = TRUE)][1]
    if (!is.na(kcol)) add_status(try_export_figure(ctx, "hdWGCNA_03_hub_genes", "hdwgcna_hub_genes",
      plot_hub_genes(ctx, hubs, value_col = kcol, top_n_per_module = ctx$cfg$figures$hdwgcna_hubs_per_module %||% 8L),
      "extended", 183, max(140, min(260, nrow(hubs) * 2))))
  }
  ht <- file.path(ctx$dirs$advanced, "hdWGCNA_module_trait_correlations.csv")
  if (file.exists(ht)) add_status(try_export_figure(ctx, "hdWGCNA_04_module_trait", "hdwgcna_trait_heatmap",
    plot_advanced_heatmap(ctx, read_csv_safe(ht), "trait", "module", "correlation", title = "Module-trait correlations",
                          midpoint = 0, diverging = TRUE, limits = c(-1, 1), text_col = NULL, fill_label = "Spearman r"),
    "extended", 160, 145))
  hdme <- file.path(ctx$dirs$advanced, "hdWGCNA_differential_module_eigengenes.csv")
  if (file.exists(hdme)) add_status(try_export_figure(ctx, "hdWGCNA_05_differential_modules", "hdwgcna_dme",
    plot_lollipop_ranking(ctx, read_csv_safe(hdme), "module", "delta_ME", "Differential module eigengenes",
                          top_n = ctx$cfg$figures$maximum_hdwgcna_modules %||% 20L, x_label = "Mean eigengene difference"),
    "extended", 160, 145))
  hss <- file.path(ctx$dirs$advanced, "hdWGCNA_sample_module_summary.csv")
  if (file.exists(hss)) add_status(try_export_figure(ctx, "hdWGCNA_06_sample_module_activity", "hdwgcna_sample_heatmap",
    plot_advanced_heatmap(ctx, read_csv_safe(hss), "sample_id", "module", "mean_ME",
                          title = "Sample-level module eigengene activity", midpoint = 0, diverging = TRUE,
                          cluster_x = TRUE, cluster_y = TRUE, fill_label = "Mean eigengene"),
    "extended", 183, 160))
  hd_dendro <- file.path(ctx$dirs$advanced, "hdWGCNA_dendrogram_plot.rds")
  if (file.exists(hd_dendro)) {
    dp <- tryCatch(readRDS(hd_dendro), error = function(e) NULL)
    if (!is.null(dp) && inherits(dp, c("gg", "ggplot")) && !inherits(dp, "patchwork")) {
      add_status(try_export_figure(ctx, "hdWGCNA_07_module_dendrogram", "hdwgcna_dendrogram",
        list(plot = dp,
             source_data = list(modules = if (file.exists(file.path(ctx$dirs$advanced,"hdWGCNA_modules.csv"))) read_csv_safe(file.path(ctx$dirs$advanced,"hdWGCNA_modules.csv")) else data.frame()),
             parameters = list(type = "hdwgcna_dendrogram")),
        "extended", 183, 135))
    }
  }
  hdobj <- file.path(ctx$dirs$advanced, "hdWGCNA_object.rds")
  hdcells <- file.path(ctx$dirs$advanced, "hdWGCNA_cell_module_eigengenes.csv.gz")
  if (file.exists(hdobj) && file.exists(hdcells)) {
    ho <- readRDS(hdobj); hm <- read_csv_safe(hdcells); modules <- setdiff(names(hm), "cell_id")
    for (m in head(modules, ctx$cfg$figures$maximum_hdwgcna_module_embeddings %||% 12L)) {
      ho[[paste0("V4_", m)]] <- hm[[m]][match(colnames(ho), hm$cell_id)]
      add_status(try_export_figure(ctx, paste0("hdWGCNA_embedding_", sanitize_stem(m)), "hdwgcna_module_embedding",
        plot_embedding_continuous(ctx, ho, paste0("V4_", m), reduction, paste(m, "eigengene")), "extended", 150, 135))
    }
  }

  # 12. Spatial transcriptomics: each slice, cell-type score, spatial gene and neighborhood result is independent.
  scoord <- file.path(ctx$dirs$advanced, "spatial_coordinates_metadata.csv.gz")
  if (file.exists(scoord)) {
    sd <- read_csv_safe(scoord)
    slice_col <- if ("spatial_sample_id" %in% names(sd)) "spatial_sample_id" else "image"
    slices <- head(unique(as.character(sd[[slice_col]])), ctx$cfg$figures$maximum_spatial_slices %||% 12L)
    for (sl in slices) {
      dd <- sd[as.character(sd[[slice_col]]) == sl, , drop = FALSE]
      if ("spatial_cluster" %in% names(dd)) {
        names(dd)[names(dd) == "spatial_cluster"] <- "group"
        add_status(try_export_figure(ctx, paste0("Spatial_cluster_", sanitize_stem(sl)), "spatial_discrete",
          plot_spatial_discrete(ctx, dd, title = paste("Spatial clusters:", sl)), "main", 160, 145))
        names(dd)[names(dd) == "group"] <- "spatial_cluster"
      }
      if ("predicted_cell_type" %in% names(dd)) {
        names(dd)[names(dd) == "predicted_cell_type"] <- "group"
        add_status(try_export_figure(ctx, paste0("Spatial_predicted_type_", sanitize_stem(sl)), "spatial_discrete",
          plot_spatial_discrete(ctx, dd, title = paste("Transferred cell types:", sl)), "main", 160, 145))
      }
    }
  }
  slabel <- file.path(ctx$dirs$advanced, "spatial_label_transfer_predictions.csv")
  if (file.exists(slabel) && file.exists(scoord)) {
    lp <- read_csv_safe(slabel); sd <- read_csv_safe(scoord)
    dd_all <- merge(sd[, intersect(c("cell_id","x","y","spatial_sample_id","image"), names(sd)), drop=FALSE], lp, by="cell_id", all.x=TRUE)
    names(dd_all)[names(dd_all)=="prediction_score_max"] <- "value"
    slice_col <- if ("spatial_sample_id" %in% names(dd_all)) "spatial_sample_id" else if ("image" %in% names(dd_all)) "image" else NULL
    if ("value" %in% names(dd_all)) {
      split_conf <- if (!is.null(slice_col)) split(dd_all, as.character(dd_all[[slice_col]]), drop=TRUE) else list(slice1=dd_all)
      split_conf <- head(split_conf, ctx$cfg$figures$maximum_spatial_slices %||% 12L)
      for (sl in names(split_conf)) add_status(try_export_figure(ctx,
        paste0("Spatial_label_transfer_confidence_", sanitize_stem(sl)), "spatial_continuous",
        plot_spatial_continuous(ctx, split_conf[[sl]], title = paste("Label-transfer confidence:", sl)),
        "extended", 150, 135))
    }
  }
  sscores <- file.path(ctx$dirs$advanced, "spatial_label_transfer_scores.csv.gz")
  if (file.exists(sscores) && file.exists(scoord)) {
    sc <- read_csv_safe(sscores); sd <- read_csv_safe(scoord)
    keep <- head(names(sort(tapply(sc$score, sc$cell_type, mean, na.rm=TRUE), decreasing=TRUE)),
                 ctx$cfg$figures$maximum_spatial_celltype_maps %||% 12L)
    base_cols <- intersect(c("cell_id","x","y","spatial_sample_id","image"), names(sd))
    for (ct in keep) {
      all_ct <- merge(sd[, base_cols, drop=FALSE], sc[sc$cell_type==ct,c("cell_id","score")], by="cell_id", all.x=TRUE)
      names(all_ct)[names(all_ct)=="score"] <- "value"
      slice_col <- if ("spatial_sample_id" %in% names(all_ct)) "spatial_sample_id" else if ("image" %in% names(all_ct)) "image" else NULL
      split_ct <- if (!is.null(slice_col)) split(all_ct, as.character(all_ct[[slice_col]]), drop=TRUE) else list(slice1=all_ct)
      split_ct <- head(split_ct, ctx$cfg$figures$maximum_spatial_slices %||% 12L)
      for (sl in names(split_ct)) add_status(try_export_figure(ctx,
        paste0("Spatial_transfer_score_", sanitize_stem(ct), "_", sanitize_stem(sl)), "spatial_continuous",
        plot_spatial_continuous(ctx, split_ct[[sl]], title = paste("Transfer score:", ct, "|", sl)),
        "extended", 150, 135))
    }
  }
  srctd <- file.path(ctx$dirs$advanced, "spatial_RCTD_weights.csv.gz")
  if (file.exists(srctd) && file.exists(scoord)) {
    rw <- read_csv_safe(srctd); sd <- read_csv_safe(scoord)
    means <- aggregate(weight ~ cell_type, rw, mean, na.rm = TRUE)
    ctypes <- head(means$cell_type[order(means$weight, decreasing = TRUE)], ctx$cfg$figures$maximum_spatial_celltype_maps %||% 12L)
    base_cols <- intersect(c("cell_id","x","y","spatial_sample_id","image"), names(sd))
    for (ct in ctypes) {
      all_ct <- merge(sd[, base_cols, drop=FALSE], rw[rw$cell_type == ct, c("cell_id", "weight")], by = "cell_id", all.x = TRUE)
      names(all_ct)[names(all_ct) == "weight"] <- "value"
      slice_col <- if ("spatial_sample_id" %in% names(all_ct)) "spatial_sample_id" else if ("image" %in% names(all_ct)) "image" else NULL
      split_ct <- if (!is.null(slice_col)) split(all_ct, as.character(all_ct[[slice_col]]), drop=TRUE) else list(slice1=all_ct)
      split_ct <- head(split_ct, ctx$cfg$figures$maximum_spatial_slices %||% 12L)
      for (sl in names(split_ct)) add_status(try_export_figure(ctx,
        paste0("Spatial_RCTD_", sanitize_stem(ct), "_", sanitize_stem(sl)), "spatial_continuous",
        plot_spatial_continuous(ctx, split_ct[[sl]], title = paste("RCTD weight:", ct, "|", sl)),
        "extended", 150, 135))
    }
  }
  sobj <- file.path(ctx$dirs$advanced, "spatial_processed_object.rds")
  svgf <- file.path(ctx$dirs$advanced, "spatial_variable_features.csv")
  if (file.exists(sobj) && file.exists(scoord) && file.exists(svgf)) {
    so <- readRDS(sobj); sd <- read_csv_safe(scoord); sf <- read_csv_safe(svgf)
    base_cols <- intersect(c("cell_id","x","y","spatial_sample_id","image"), names(sd))
    for (g in head(sf$gene, ctx$cfg$figures$maximum_spatial_gene_maps %||% 8L)) {
      vals <- tryCatch(get_feature_vector(so, g), error = function(e) NULL)
      if (is.null(vals)) next
      all_g <- merge(sd[, base_cols, drop=FALSE], data.frame(cell_id = names(vals), value = as.numeric(vals)), by = "cell_id")
      slice_col <- if ("spatial_sample_id" %in% names(all_g)) "spatial_sample_id" else if ("image" %in% names(all_g)) "image" else NULL
      split_g <- if (!is.null(slice_col)) split(all_g, as.character(all_g[[slice_col]]), drop=TRUE) else list(slice1=all_g)
      split_g <- head(split_g, ctx$cfg$figures$maximum_spatial_slices %||% 12L)
      for (sl in names(split_g)) add_status(try_export_figure(ctx,
        paste0("Spatial_gene_", sanitize_stem(g), "_", sanitize_stem(sl)), "spatial_continuous",
        plot_spatial_continuous(ctx, split_g[[sl]], title = paste(g, "spatial expression |", sl)),
        "extended", 150, 135))
    }
  }
  if (file.exists(svgf)) {
    sf <- read_csv_safe(svgf)
    if (all(c("gene","rank") %in% names(sf))) {
      sf$priority <- 1 / pmax(sf$rank, 1)
      add_status(try_export_figure(ctx, "Spatial_variable_gene_ranking", "spatial_variable_gene_ranking",
        plot_lollipop_ranking(ctx, sf, "gene", "priority", "Top spatially variable genes",
                              top_n = ctx$cfg$figures$maximum_spatial_variable_genes %||% 30L,
                              x_label = "Reciprocal rank"), "extended", 155, 155))
    }
  }
  sne <- file.path(ctx$dirs$advanced, "spatial_neighborhood_enrichment.csv")
  if (file.exists(sne)) {
    ne <- read_csv_safe(sne)
    split_ne <- if ("spatial_sample_id" %in% names(ne)) split(ne, as.character(ne$spatial_sample_id), drop=TRUE) else list(slice1=ne)
    split_ne <- head(split_ne, ctx$cfg$figures$maximum_spatial_slices %||% 12L)
    for (sl in names(split_ne)) add_status(try_export_figure(ctx,
      paste0("Spatial_neighborhood_enrichment_", sanitize_stem(sl)), "spatial_neighborhood_heatmap",
      plot_advanced_heatmap(ctx, split_ne[[sl]], "source", "target", "z_score",
                            title = paste("Spatial neighborhood enrichment:", sl), midpoint = 0,
                            diverging = TRUE, fill_label = "Permutation z"),
      "extended", 170, 155))
  }

  # 13. Exploratory drug-response hypotheses. No plot is labelled as a clinical recommendation.
  dr <- file.path(ctx$dirs$advanced, "drug_reversal_ranking.csv")
  dc <- file.path(ctx$dirs$advanced, "drug_reversal_gene_contributions.csv.gz")
  if (file.exists(dr)) {
    rr <- read_csv_safe(dr)
    for (ct in head(unique(as.character(rr$cell_type)), ctx$cfg$figures$maximum_drug_celltypes %||% 12L)) {
      dd <- rr[rr$cell_type == ct, , drop = FALSE]
      add_status(try_export_figure(ctx, paste0("Drug_reversal_rank_", sanitize_stem(ct)), "drug_reversal_ranking",
        plot_lollipop_ranking(ctx, dd, "drug", "reversal_score", paste("Exploratory signature reversal:", ct),
                              top_n = ctx$cfg$figures$drug_top_n %||% 25L, x_label = "Reversal score"), "extended", 160, 150))
      if (file.exists(dc) && nrow(dd)) {
        contrib <- read_csv_safe(dc); best <- dd$drug[which.max(dd$reversal_score)]
        cc2 <- contrib[contrib$cell_type == ct & contrib$drug == best, , drop = FALSE]
        add_status(try_export_figure(ctx, paste0("Drug_reversal_genes_", sanitize_stem(ct), "_", sanitize_stem(best)), "drug_reversal_scatter",
          plot_drug_reversal_scatter(ctx, cc2, title = paste("Signature reversal:", best, "in", ct)), "extended", 150, 135))
      }
    }
  }
  ds <- file.path(ctx$dirs$advanced, "drug_sample_celltype_scores.csv")
  if (file.exists(ds)) {
    dscore <- read_csv_safe(ds); dscore$group <- paste(dscore$sample_id, dscore$cell_type, sep = " | ")
    add_status(try_export_figure(ctx, "Drug_signature_score_heatmap", "drug_score_heatmap",
      plot_drug_score_heatmap(ctx, dscore, group_col = "group", title = "Exploratory drug-signature scores",
                              top_n = ctx$cfg$figures$drug_top_n %||% 30L), "extended", 210, 180))
  }
  tc_path <- file.path(ctx$dirs$advanced, "drug_therapeutic_clusters.csv")
  if (file.exists(tc_path)) {
    tc <- read_csv_safe(tc_path); dobj <- obj
    dobj$therapeutic_cluster <- tc$therapeutic_cluster[match(colnames(dobj), tc$cell_id)]
    add_status(try_export_figure(ctx, "Drug_therapeutic_clusters_embedding", "drug_therapeutic_cluster_embedding",
      plot_embedding_discrete(ctx, dobj, "therapeutic_cluster", reduction,
                              "Exploratory transcriptomic therapeutic clusters", labels=TRUE, hulls=TRUE),
      "extended", 183, 155))
  }
  if (file.exists(ds)) {
    dscore <- read_csv_safe(ds)
    top_drugs <- head(names(sort(tapply(abs(dscore$score), dscore$drug, mean, na.rm=TRUE), decreasing=TRUE)),
                      ctx$cfg$figures$maximum_drug_distribution_plots %||% 8L)
    for (drug in top_drugs) {
      dd <- dscore[dscore$drug==drug, , drop=FALSE]
      split <- if ("condition" %in% names(dd) && length(unique(stats::na.omit(dd$condition))) > 1L) "condition" else NULL
      add_status(try_export_figure(ctx, paste0("Drug_score_distribution_", sanitize_stem(drug)), "drug_score_distribution",
        plot_group_score_distribution(ctx, dd, "cell_type", "score", split,
                                      paste("Exploratory signature score:", drug), "Standardized score"),
        "extended", 183, 145))
    }
  }
  op <- file.path(ctx$dirs$advanced, "oncoPredict_pseudobulk_predictions.csv")
  if (file.exists(op)) {
    od <- read_csv_safe(op); od$group <- paste(od$sample_id, od$cell_type, sep = " | "); names(od)[names(od) == "predicted_sensitivity"] <- "score"
    add_status(try_export_figure(ctx, "oncoPredict_pseudobulk_heatmap", "oncopredict_heatmap",
      plot_drug_score_heatmap(ctx, od, group_col = "group", title = "Exploratory pseudobulk drug-sensitivity prediction",
                              top_n = ctx$cfg$figures$drug_top_n %||% 30L), "extended", 210, 180))
  }

  # 14. Exploratory virtual knockout. Every target/subset/readout is an independent figure.
  vk_consensus_path <- file.path(ctx$dirs$advanced, "virtual_knockout_consensus_genes.csv")
  vk_sample_path <- file.path(ctx$dirs$advanced, "virtual_knockout_sample_consensus.csv")
  vk_audit_path <- file.path(ctx$dirs$advanced, "virtual_knockout_target_expression_audit.csv")
  vk_edges_path <- file.path(ctx$dirs$advanced, "virtual_knockout_consensus_network_edges.csv")
  vk_similarity_path <- file.path(ctx$dirs$advanced, "virtual_knockout_run_similarity.csv")
  vk_pathways_path <- file.path(ctx$dirs$advanced, "virtual_knockout_pathway_enrichment.csv")
  vk_manifold_path <- file.path(ctx$dirs$advanced, "virtual_knockout_manifold_coordinates.csv.gz")
  vk_diff_path <- file.path(ctx$dirs$advanced, "virtual_knockout_diffregulation_all.csv.gz")
  if (file.exists(vk_consensus_path)) {
    vk_consensus <- read_csv_safe(vk_consensus_path)
    vk_sample <- if (file.exists(vk_sample_path)) read_csv_safe(vk_sample_path) else data.frame()
    vk_audit <- if (file.exists(vk_audit_path)) read_csv_safe(vk_audit_path) else data.frame()
    vk_edges <- if (file.exists(vk_edges_path)) read_csv_safe(vk_edges_path) else data.frame()
    vk_similarity <- if (file.exists(vk_similarity_path)) read_csv_safe(vk_similarity_path) else data.frame()
    vk_pathways <- if (file.exists(vk_pathways_path)) read_csv_safe(vk_pathways_path) else data.frame()
    vk_manifold <- if (file.exists(vk_manifold_path)) read_csv_safe(vk_manifold_path) else data.frame()
    vk_diff <- if (file.exists(vk_diff_path)) read_csv_safe(vk_diff_path) else data.frame()
    vk_groups <- split(vk_consensus, interaction(vk_consensus$target_gene, vk_consensus$subset_value, drop = TRUE, lex.order = TRUE))
    for (gname in names(vk_groups)) {
      cc <- vk_groups[[gname]]
      tg <- as.character(cc$target_gene[1]); sv <- as.character(cc$subset_value[1])
      stem <- paste(sanitize_stem(tg), sanitize_stem(sv), sep = "__")
      add_status(try_export_figure(ctx, paste0("VirtualKO_rank__", stem), "virtual_knockout_ranking",
        plot_virtual_knockout_ranking(ctx, cc,
          title = paste0("Virtual knockout of ", tg, " | ", sv, ": regulatory impact"),
          top_n = ctx$cfg$figures$virtual_knockout_top_genes %||% 30L,
          fdr_threshold = ctx$cfg$advanced_modules$virtual_knockout$fdr_threshold %||% 0.05),
        "extended", 175, 170))

      if (nrow(vk_sample)) {
        ss <- vk_sample[vk_sample$target_gene == tg & vk_sample$subset_value == sv, , drop = FALSE]
        if (nrow(ss)) add_status(try_export_figure(ctx, paste0("VirtualKO_sample_support__", stem), "virtual_knockout_sample_support",
          plot_virtual_knockout_sample_support(ctx, ss, cc,
            title = paste0("Virtual knockout of ", tg, " | ", sv, ": sample support"),
            top_n = ctx$cfg$figures$virtual_knockout_top_genes %||% 30L),
          "extended", 190, 175))
      }

      if (nrow(vk_audit)) {
        aa <- vk_audit[vk_audit$target_gene == tg & vk_audit$subset_value == sv, , drop = FALSE]
        if (nrow(aa)) add_status(try_export_figure(ctx, paste0("VirtualKO_target_expression__", stem), "virtual_knockout_target_expression",
          plot_virtual_knockout_target_expression(ctx, aa,
            title = paste0(tg, " expression eligibility | ", sv)),
          "extended", 155, max(100, min(190, 12 * nrow(aa) + 65))))
      }

      if (nrow(vk_edges)) {
        ee <- vk_edges[vk_edges$target_gene == tg & vk_edges$subset_value == sv, , drop = FALSE]
        if (nrow(ee)) add_status(try_export_figure(ctx, paste0("VirtualKO_network__", stem), "virtual_knockout_network",
          plot_virtual_knockout_network(ctx, ee, tg,
            title = paste0(tg, "-centred recurrent WT-network context | ", sv),
            top_n = ctx$cfg$figures$virtual_knockout_network_edges %||% 80L),
          "extended", 175, 165))
      }

      if (nrow(vk_pathways)) {
        pp <- vk_pathways[vk_pathways$target_gene == tg & vk_pathways$subset_value == sv, , drop = FALSE]
        if (nrow(pp)) add_status(try_export_figure(ctx, paste0("VirtualKO_pathways__", stem), "virtual_knockout_pathways",
          plot_virtual_knockout_pathways(ctx, pp,
            title = paste0("Pathways among ", tg, " virtual-KO impacted genes | ", sv),
            top_n = ctx$cfg$figures$virtual_knockout_top_pathways %||% 20L),
          "extended", 175, 160))
      }

      if (nrow(vk_similarity)) {
        rr <- vk_similarity[vk_similarity$target_gene == tg & vk_similarity$subset_value == sv, , drop = FALSE]
        if (nrow(rr)) add_status(try_export_figure(ctx, paste0("VirtualKO_run_similarity__", stem), "virtual_knockout_run_similarity",
          plot_advanced_heatmap(ctx, rr, "run_1", "run_2", "jaccard_top_genes",
            title = paste0(tg, " virtual-KO run stability | ", sv), midpoint = 0.5,
            diverging = FALSE, limits = c(0, 1), fill_label = "Top-gene Jaccard",
            cluster_x = TRUE, cluster_y = TRUE),
          "extended", 175, 165))
      }

      if (nrow(vk_manifold) && nrow(vk_diff)) {
        mm0 <- vk_manifold[vk_manifold$target_gene == tg & vk_manifold$subset_value == sv, , drop = FALSE]
        dd0 <- vk_diff[vk_diff$target_gene == tg & vk_diff$subset_value == sv, , drop = FALSE]
        if (nrow(mm0) && nrow(dd0)) {
          candidates <- unique(mm0[, c("sample_id", "repeat_id", "analysis_unit"), drop = FALSE])
          candidates$priority <- ifelse(candidates$analysis_unit == "balanced_pooled", 0L, 1L)
          candidates <- candidates[order(candidates$priority, candidates$sample_id, candidates$repeat_id), , drop = FALSE]
          pick <- candidates[1, , drop = FALSE]
          same_run <- function(z) as.character(z$sample_id) == as.character(pick$sample_id) &
            as.character(z$repeat_id) == as.character(pick$repeat_id) &
            as.character(z$analysis_unit) == as.character(pick$analysis_unit)
          mm <- mm0[same_run(mm0), , drop = FALSE]
          dd <- dd0[same_run(dd0), , drop = FALSE]
          if (nrow(mm) && nrow(dd)) add_status(try_export_figure(ctx, paste0("VirtualKO_manifold__", stem), "virtual_knockout_manifold",
            plot_virtual_knockout_manifold(ctx, mm, dd,
              title = paste0(tg, " WT–KO network-state displacement | ", sv, " | ", pick$sample_id, " R", pick$repeat_id),
              top_n = ctx$cfg$figures$virtual_knockout_manifold_genes %||% 20L,
              fdr_threshold = ctx$cfg$advanced_modules$virtual_knockout$fdr_threshold %||% 0.05),
            "extended", 175, 160))
        }
      }
    }
  }

  # Composite or contact-sheet outputs are intentionally forbidden in V5.
  status <- do.call(rbind, statuses)
  write_csv_safe(status, file.path(ctx$dirs$qa, "figure_generation_status.csv"))
  if (!any(status$status %in% c("EXPORTED", "SKIPPED_EXISTING"))) stop("The figure engine did not export or reuse any figure.")
}
