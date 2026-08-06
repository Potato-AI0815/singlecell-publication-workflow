# Publication Figure Engine v4: advanced single-figure templates.
# Every function returns list(plot, source_data, parameters) and is exported separately.

advanced_diverging_scale <- function(ctx, midpoint = 0, limits = NULL, name = NULL) {
  ggplot2::scale_fill_gradient2(
    low = figure_cfg(ctx, "continuous_low", "#E67861"), mid = "#F7F7F7",
    high = figure_cfg(ctx, "continuous_high", "#557FB5"), midpoint = midpoint,
    limits = limits, oob = scales::squish, na.value = "#D0D0D0", name = name
  )
}

advanced_sequential_scale <- function(ctx, limits = NULL, name = NULL, option = "D") {
  ggplot2::scale_fill_viridis_c(option = option, limits = limits, oob = scales::squish, na.value = "#D0D0D0", name = name)
}

plot_advanced_heatmap <- function(ctx, df, x_col, y_col, value_col, title = NULL,
                                  midpoint = 0, diverging = TRUE, limits = NULL,
                                  x_label = NULL, y_label = NULL, fill_label = NULL,
                                  text_col = NULL, cluster_x = FALSE, cluster_y = FALSE) {
  req <- c(x_col, y_col, value_col)
  if (!all(req %in% names(df))) stop("Heatmap data missing: ", paste(setdiff(req, names(df)), collapse = ", "))
  raw <- df[complete.cases(df[, req, drop = FALSE]) & is.finite(df[[value_col]]), , drop = FALSE]
  if (!nrow(raw)) stop("Heatmap has no finite rows.")
  key <- interaction(as.character(raw[[x_col]]), as.character(raw[[y_col]]), drop = TRUE, lex.order = TRUE)
  spl <- split(seq_len(nrow(raw)), key)
  dd <- do.call(rbind, lapply(spl, function(ii) {
    z <- raw[ii, , drop = FALSE]
    out <- z[1, , drop = FALSE]
    out[[value_col]] <- mean(z[[value_col]], na.rm = TRUE)
    if (!is.null(text_col) && text_col %in% names(z)) out[[text_col]] <- as.character(z[[text_col]][1])
    out
  }))
  rownames(dd) <- NULL
  xlev <- unique(as.character(dd[[x_col]])); ylev <- unique(as.character(dd[[y_col]]))
  safe_cluster_order <- function(mat, labels, by_columns = TRUE) {
    z <- if (by_columns) t(mat) else mat
    z <- as.matrix(z); z[!is.finite(z)] <- 0
    if (nrow(z) <= 2L) return(labels)
    d <- tryCatch(stats::dist(z), error = function(e) NULL)
    if (is.null(d) || !length(d) || all(as.numeric(d) == 0)) return(labels)
    labels[stats::hclust(d, method = "ward.D2")$order]
  }
  if ((isTRUE(cluster_x) && length(xlev) > 2L) || (isTRUE(cluster_y) && length(ylev) > 2L)) {
    wide <- xtabs(stats::as.formula(paste(value_col, "~", y_col, "+", x_col)), dd)
    if (isTRUE(cluster_x) && length(xlev) > 2L) xlev <- safe_cluster_order(wide, colnames(wide), TRUE)
    if (isTRUE(cluster_y) && length(ylev) > 2L) ylev <- safe_cluster_order(wide, rownames(wide), FALSE)
  }
  dd[[x_col]] <- factor(as.character(dd[[x_col]]), levels = xlev)
  dd[[y_col]] <- factor(as.character(dd[[y_col]]), levels = rev(ylev))
  p <- ggplot2::ggplot(dd, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[value_col]])) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.18) +
    pub_theme(ctx) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank()) +
    ggplot2::labs(title = title, x = x_label %||% gsub("_", " ", x_col),
                  y = y_label %||% gsub("_", " ", y_col), fill = fill_label %||% value_col)
  p <- p + if (isTRUE(diverging)) advanced_diverging_scale(ctx, midpoint, limits, fill_label %||% value_col) else advanced_sequential_scale(ctx, limits, fill_label %||% value_col)
  if (!is.null(text_col) && text_col %in% names(dd)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = .data[[text_col]]), size = 2.4, colour = "black")
  }
  list(plot = p, source_data = list(raw = raw, display = dd),
       parameters = list(type = "advanced_heatmap", x_col = x_col, y_col = y_col, value_col = value_col,
                         midpoint = midpoint, diverging = diverging, limits = limits, cluster_x = cluster_x, cluster_y = cluster_y,
                         duplicate_cells_aggregated_by = "mean"))
}

plot_rank_diagnostics <- function(ctx, df, rank_col = "rank", metric_col = "metric_value",
                                  metric_name_col = "metric", title = "NMF rank diagnostics") {
  req <- c(rank_col, metric_col, metric_name_col)
  if (!all(req %in% names(df))) stop("Rank diagnostics missing required columns.")
  dd <- df[is.finite(df[[metric_col]]), , drop = FALSE]
  p <- ggplot2::ggplot(dd, ggplot2::aes(.data[[rank_col]], .data[[metric_col]], colour = .data[[metric_name_col]], group = .data[[metric_name_col]])) +
    ggplot2::geom_line(linewidth = 0.75) + ggplot2::geom_point(size = 2.1) +
    pub_theme(ctx) + ggplot2::labs(title = title, x = "NMF rank", y = "Diagnostic value", colour = "Metric")
  list(plot = p, source_data = dd, parameters = list(type = "rank_diagnostics"))
}

plot_program_prevalence_dotplot <- function(ctx, df, program_col = "program", group_col = "group",
                                            mean_col = "mean_score", pct_col = "percent_high",
                                            title = "Program activity across groups") {
  req <- c(program_col, group_col, mean_col, pct_col)
  if (!all(req %in% names(df))) stop("Program prevalence table is incomplete.")
  dd <- df[complete.cases(df[, req, drop = FALSE]), , drop = FALSE]
  p <- ggplot2::ggplot(dd, ggplot2::aes(.data[[program_col]], .data[[group_col]], size = .data[[pct_col]], colour = .data[[mean_col]])) +
    ggplot2::geom_point(alpha = 0.92) +
    ggplot2::scale_size_area(max_size = 8, limits = c(0, max(100, max(dd[[pct_col]], na.rm = TRUE)))) +
    ggplot2::scale_colour_gradient2(low = "#E67861", mid = "#F3F0F7", high = "#54278F", midpoint = 0, oob = scales::squish) +
    pub_theme(ctx) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = title, x = "Program", y = gsub("_", " ", group_col), size = "% high", colour = "Mean score")
  list(plot = p, source_data = dd, parameters = list(type = "program_prevalence_dotplot"))
}

plot_lollipop_ranking <- function(ctx, df, label_col, value_col, title = NULL, top_n = 25,
                                  colour_col = NULL, x_label = NULL) {
  req <- c(label_col, value_col)
  if (!all(req %in% names(df))) stop("Ranking table missing required columns.")
  dd <- df[is.finite(df[[value_col]]) & !is.na(df[[label_col]]), , drop = FALSE]
  dd <- dd[order(dd[[value_col]], decreasing = TRUE), , drop = FALSE]
  dd <- head(dd, top_n)
  dd[[label_col]] <- factor(as.character(dd[[label_col]]), levels = rev(as.character(dd[[label_col]])))
  aes_col <- if (!is.null(colour_col) && colour_col %in% names(dd)) ggplot2::aes(colour = .data[[colour_col]]) else NULL
  p <- ggplot2::ggplot(dd, ggplot2::aes(y = .data[[label_col]], x = .data[[value_col]])) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data[[value_col]], yend = .data[[label_col]]), linewidth = 0.55, colour = "#B7B7B7") +
    ggplot2::geom_point(mapping = aes_col, size = 2.8, colour = if (is.null(aes_col)) "#4C78A8" else NULL) +
    pub_theme(ctx) + ggplot2::labs(title = title, x = x_label %||% value_col, y = NULL, colour = colour_col)
  if (!is.null(aes_col)) p <- p + ggplot2::scale_colour_viridis_c(option = "D")
  list(plot = p, source_data = dd, parameters = list(type = "lollipop_ranking", top_n = top_n))
}

plot_cnv_heatmap <- function(ctx, df, cell_col = "cell_id", bin_col = "genomic_bin", value_col = "cna_value",
                             cell_group_col = "cnv_class", title = "Expression-inferred copy-number landscape",
                             display_limits = c(-1, 1)) {
  req <- c(cell_col, bin_col, value_col)
  if (!all(req %in% names(df))) stop("CNV heatmap table is incomplete.")
  dd <- df[is.finite(df[[value_col]]), , drop = FALSE]
  if (!nrow(dd)) stop("No finite CNV values.")
  if (!"cell_order" %in% names(dd)) dd$cell_order <- match(dd[[cell_col]], unique(dd[[cell_col]]))
  if (!"bin_order" %in% names(dd)) dd$bin_order <- match(dd[[bin_col]], unique(dd[[bin_col]]))
  dd$value_display <- pmin(pmax(dd[[value_col]], display_limits[1]), display_limits[2])
  p <- ggplot2::ggplot(dd, ggplot2::aes(bin_order, cell_order, fill = value_display)) +
    ggplot2::geom_raster() +
    advanced_diverging_scale(ctx, midpoint = 0, limits = display_limits, name = "Inferred CNA") +
    ggplot2::scale_x_continuous(expand = c(0, 0)) + ggplot2::scale_y_continuous(expand = c(0, 0)) +
    pub_theme(ctx) + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                                    axis.line = ggplot2::element_blank()) +
    ggplot2::labs(title = title, x = "Ordered genomic position", y = "Cells")
  list(plot = p, source_data = dd,
       parameters = list(type = "cnv_heatmap", display_limits = display_limits, cell_group_col = cell_group_col))
}

plot_cnv_scatter <- function(ctx, df, x_col = "cna_correlation", y_col = "cna_signal",
                             class_col = "cnv_class", title = "CNA signal and profile correlation") {
  req <- c(x_col, y_col, class_col)
  if (!all(req %in% names(df))) stop("CNV scatter table is incomplete.")
  dd <- df[is.finite(df[[x_col]]) & is.finite(df[[y_col]]), , drop = FALSE]
  lev <- sort(unique(as.character(dd[[class_col]])))
  pal <- setNames(categorical_palette(length(lev)), lev)
  p <- ggplot2::ggplot(dd, ggplot2::aes(.data[[x_col]], .data[[y_col]], colour = .data[[class_col]])) +
    ggplot2::geom_point(size = 0.65, alpha = 0.65) + ggplot2::scale_colour_manual(values = pal, na.value = "#BDBDBD") +
    pub_theme(ctx) + ggplot2::labs(title = title, x = "CNA-profile correlation", y = "CNA signal", colour = "Class")
  list(plot = p, source_data = dd, parameters = list(type = "cnv_scatter", palette = as.list(pal)))
}

plot_communication_bubble <- function(ctx, df, source_col = "source", target_col = "target",
                                      interaction_col = "interaction", probability_col = "probability",
                                      pvalue_col = "p_value", title = "Selected cell-cell communications", top_n = 40) {
  req <- c(source_col, target_col, interaction_col, probability_col)
  if (!all(req %in% names(df))) stop("Communication table is incomplete.")
  dd <- df[is.finite(df[[probability_col]]), , drop = FALSE]
  if (pvalue_col %in% names(dd)) dd <- dd[is.na(dd[[pvalue_col]]) | dd[[pvalue_col]] <= 0.05, , drop = FALSE]
  dd <- dd[order(dd[[probability_col]], decreasing = TRUE), , drop = FALSE]
  dd <- head(dd, top_n)
  dd$pair <- paste(dd[[source_col]], "→", dd[[target_col]])
  dd$interaction <- factor(as.character(dd[[interaction_col]]), levels = rev(unique(as.character(dd[[interaction_col]]))))
  p <- ggplot2::ggplot(dd, ggplot2::aes(pair, interaction, size = .data[[probability_col]], colour = .data[[probability_col]])) +
    ggplot2::geom_point(alpha = 0.9) + ggplot2::scale_colour_viridis_c(option = "D") + ggplot2::scale_size_area(max_size = 7) +
    pub_theme(ctx) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = title, x = "Sender → receiver", y = "Ligand-receptor interaction", size = "Probability", colour = "Probability")
  list(plot = p, source_data = dd, parameters = list(type = "communication_bubble", top_n = top_n))
}

plot_network_edges <- function(ctx, edges, source_col = "source", target_col = "target", weight_col = "weight",
                               title = "Interaction network", node_size = NULL, directed = TRUE, top_n = 80) {
  req <- c(source_col, target_col, weight_col)
  if (!all(req %in% names(edges))) stop("Network edge table is incomplete.")
  ee <- edges[is.finite(edges[[weight_col]]) & edges[[weight_col]] > 0, , drop = FALSE]
  ee <- head(ee[order(ee[[weight_col]], decreasing = TRUE), , drop = FALSE], top_n)
  nodes <- sort(unique(c(as.character(ee[[source_col]]), as.character(ee[[target_col]]))))
  theta <- seq(0, 2 * pi, length.out = length(nodes) + 1L)[-1L]
  nd <- data.frame(node = nodes, x = cos(theta), y = sin(theta), stringsAsFactors = FALSE)
  if (!is.null(node_size) && is.data.frame(node_size) && all(c("node", "size") %in% names(node_size))) nd <- merge(nd, node_size, by = "node", all.x = TRUE) else nd$size <- 1
  ee$x <- nd$x[match(ee[[source_col]], nd$node)]; ee$y <- nd$y[match(ee[[source_col]], nd$node)]
  ee$xend <- nd$x[match(ee[[target_col]], nd$node)]; ee$yend <- nd$y[match(ee[[target_col]], nd$node)]
  self_loop <- as.character(ee[[source_col]]) == as.character(ee[[target_col]])
  n_self_loops <- sum(self_loop, na.rm = TRUE)
  if (n_self_loops > 0L) {
    # geom_curve rejects identical endpoints. Offset self-loop endpoints tangentially
    # around the node while preserving the original source/target labels in Source Data.
    angle <- atan2(ee$y[self_loop], ee$x[self_loop])
    dx <- 0.16 * (-sin(angle)); dy <- 0.16 * cos(angle)
    ee$x[self_loop] <- ee$x[self_loop] - dx
    ee$y[self_loop] <- ee$y[self_loop] - dy
    ee$xend[self_loop] <- ee$xend[self_loop] + dx
    ee$yend[self_loop] <- ee$yend[self_loop] + dy
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_curve(data = ee, ggplot2::aes(x, y, xend = xend, yend = yend, linewidth = .data[[weight_col]]),
                        curvature = 0.18, alpha = 0.42, colour = "#6B6B6B",
                        arrow = if (directed) grid::arrow(length = grid::unit(1.5, "mm"), type = "closed") else NULL) +
    ggplot2::geom_point(data = nd, ggplot2::aes(x, y, size = size), colour = "#F2F2F2", fill = "#BDBDBD", shape = 21, stroke = 0.35) +
    ggplot2::geom_text(data = nd, ggplot2::aes(x, y, label = node), size = 2.7, vjust = -0.9) +
    ggplot2::scale_linewidth(range = c(0.2, 2.1)) + ggplot2::scale_size(range = c(3, 8), guide = "none") +
    ggplot2::coord_equal(clip = "off") + embedding_theme(ctx) + ggplot2::theme(legend.position = "right") +
    ggplot2::labs(title = title, linewidth = "Weight")
  list(plot = p, source_data = list(edges = ee, nodes = nd), parameters = list(type = "network_edges", top_n = top_n, directed = directed, n_self_loops = n_self_loops))
}

plot_soft_power_metric <- function(ctx, df, power_col = "Power", metric_col = "SFT.R.sq",
                                   threshold = 0.8, title = "Scale-free topology fit") {
  req <- c(power_col, metric_col)
  if (!all(req %in% names(df))) stop("Soft-power table is incomplete.")
  dd <- df[is.finite(df[[power_col]]) & is.finite(df[[metric_col]]), , drop = FALSE]
  p <- ggplot2::ggplot(dd, ggplot2::aes(.data[[power_col]], .data[[metric_col]])) +
    ggplot2::geom_line(linewidth = 0.7, colour = "#4C78A8") + ggplot2::geom_point(size = 2.0, colour = "#4C78A8") +
    ggplot2::geom_hline(yintercept = threshold, linetype = 2, colour = "#E45756", linewidth = 0.5) +
    pub_theme(ctx) + ggplot2::labs(title = title, x = "Soft-thresholding power", y = metric_col)
  list(plot = p, source_data = dd, parameters = list(type = "soft_power_metric", threshold = threshold, metric = metric_col))
}

plot_hub_genes <- function(ctx, df, module_col = "module", gene_col = "gene_name", value_col = "kME",
                           title = "Top intramodular hub genes", top_n_per_module = 8) {
  req <- c(module_col, gene_col, value_col)
  if (!all(req %in% names(df))) stop("Hub-gene table is incomplete.")
  spl <- split(df, as.character(df[[module_col]]), drop = TRUE)
  dd <- do.call(rbind, lapply(spl, function(z) head(z[order(z[[value_col]], decreasing = TRUE), , drop = FALSE], top_n_per_module)))
  dd$label <- paste(dd[[module_col]], dd[[gene_col]], sep = " | ")
  dd$label <- factor(dd$label, levels = rev(dd$label))
  p <- ggplot2::ggplot(dd, ggplot2::aes(.data[[value_col]], label, fill = .data[[module_col]])) +
    ggplot2::geom_col(width = 0.72, show.legend = FALSE) + pub_theme(ctx) +
    ggplot2::labs(title = title, x = "kME", y = NULL)
  list(plot = p, source_data = dd, parameters = list(type = "hub_gene_bar", top_n_per_module = top_n_per_module))
}

plot_spatial_discrete <- function(ctx, df, x_col = "x", y_col = "y", group_col = "group",
                                  title = "Spatial domains", palette = NULL, point_size = 1.6, flip_y = TRUE) {
  req <- c(x_col, y_col, group_col)
  if (!all(req %in% names(df))) stop("Spatial coordinate table is incomplete.")
  dd <- df[is.finite(df[[x_col]]) & is.finite(df[[y_col]]) & !is.na(df[[group_col]]), , drop = FALSE]
  if (!nrow(dd)) stop("Spatial discrete plot has no finite labelled locations.")
  lev <- sort(unique(as.character(dd[[group_col]]))); palette <- palette %||% setNames(categorical_palette(length(lev)), lev)
  p <- ggplot2::ggplot(dd, ggplot2::aes(.data[[x_col]], .data[[y_col]], colour = .data[[group_col]])) +
    ggplot2::geom_point(size = point_size, alpha = 0.92, stroke = 0) + ggplot2::scale_colour_manual(values = palette, na.value = "#BDBDBD") +
    ggplot2::coord_fixed() + embedding_theme(ctx) + ggplot2::labs(title = title, colour = gsub("_", " ", group_col))
  if (isTRUE(flip_y)) p <- p + ggplot2::scale_y_reverse()
  list(plot = p, source_data = dd, parameters = list(type = "spatial_discrete", flip_y = flip_y, palette = as.list(palette)))
}

plot_spatial_continuous <- function(ctx, df, x_col = "x", y_col = "y", value_col = "value",
                                    title = NULL, point_size = 1.6, quantile_limits = c(0.01, 0.99), flip_y = TRUE) {
  req <- c(x_col, y_col, value_col)
  if (!all(req %in% names(df))) stop("Spatial continuous table is incomplete.")
  dd <- df[is.finite(df[[x_col]]) & is.finite(df[[y_col]]), , drop = FALSE]
  finite_values <- dd[[value_col]][is.finite(dd[[value_col]])]
  if (!nrow(dd) || !length(finite_values)) stop("Spatial continuous plot has no finite values.")
  lim <- stats::quantile(finite_values, quantile_limits, na.rm = TRUE, names = FALSE)
  if (length(lim) != 2L || any(!is.finite(lim)) || lim[1] == lim[2]) lim <- range(finite_values, finite = TRUE)
  if (lim[1] == lim[2]) lim <- lim + c(-0.5, 0.5) * max(abs(lim[1]), 1) * 0.02
  dd$value_display <- pmin(pmax(dd[[value_col]], lim[1]), lim[2])
  p <- ggplot2::ggplot(dd, ggplot2::aes(.data[[x_col]], .data[[y_col]], colour = value_display)) +
    ggplot2::geom_point(size = point_size, alpha = 0.95, stroke = 0) +
    ggplot2::scale_colour_viridis_c(option = "D", limits = lim, oob = scales::squish, na.value = "#D0D0D0") +
    ggplot2::coord_fixed() + embedding_theme(ctx) + ggplot2::labs(title = title %||% value_col, colour = value_col)
  if (isTRUE(flip_y)) p <- p + ggplot2::scale_y_reverse()
  list(plot = p, source_data = dd, parameters = list(type = "spatial_continuous", limits = lim, quantile_limits = quantile_limits, flip_y = flip_y))
}

plot_drug_score_heatmap <- function(ctx, df, drug_col = "drug", group_col = "group", score_col = "score",
                                    title = "Exploratory drug-response scores", top_n = 30) {
  req <- c(drug_col, group_col, score_col)
  if (!all(req %in% names(df))) stop("Drug score table is incomplete.")
  dd <- df[is.finite(df[[score_col]]), , drop = FALSE]
  ranks <- aggregate(abs(dd[[score_col]]), list(drug = dd[[drug_col]]), mean, na.rm = TRUE)
  keep <- head(ranks$drug[order(ranks$x, decreasing = TRUE)], top_n)
  dd <- dd[dd[[drug_col]] %in% keep, , drop = FALSE]
  names(dd)[names(dd) == drug_col] <- "drug_plot"; names(dd)[names(dd) == group_col] <- "group_plot"; names(dd)[names(dd) == score_col] <- "score_plot"
  plot_advanced_heatmap(ctx, dd, "group_plot", "drug_plot", "score_plot", title = title,
                        midpoint = 0, diverging = TRUE, fill_label = "Score", cluster_x = FALSE, cluster_y = TRUE)
}

plot_drug_reversal_scatter <- function(ctx, df, disease_col = "disease_effect", drug_col = "drug_effect",
                                       gene_col = "gene", title = "Disease–drug signature reversal", label_n = 12) {
  req <- c(disease_col, drug_col, gene_col)
  if (!all(req %in% names(df))) stop("Drug reversal table is incomplete.")
  dd <- df[is.finite(df[[disease_col]]) & is.finite(df[[drug_col]]), , drop = FALSE]
  dd$reversal_contribution <- -dd[[disease_col]] * dd[[drug_col]]
  lab <- head(dd[order(abs(dd$reversal_contribution), decreasing = TRUE), , drop = FALSE], label_n)
  p <- ggplot2::ggplot(dd, ggplot2::aes(.data[[disease_col]], .data[[drug_col]], colour = reversal_contribution)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "#888888") + ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "#888888") +
    ggplot2::geom_point(size = 0.75, alpha = 0.65) + ggplot2::scale_colour_gradient2(low = "#E67861", mid = "#F7F7F7", high = "#557FB5", midpoint = 0) +
    pub_theme(ctx) + ggplot2::labs(title = title, x = "Disease effect", y = "Drug perturbation effect", colour = "Reversal\ncontribution")
  if (nrow(lab) && requireNamespace("ggrepel", quietly = TRUE)) p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = .data[[gene_col]]), size = 2.4, max.overlaps = Inf, seed = ctx$cfg$runtime$seed %||% 1L)
  list(plot = p, source_data = dd, parameters = list(type = "drug_reversal_scatter", label_n = label_n))
}

plot_metric_curve <- function(ctx, df, x_col, y_col, title = NULL, x_label = NULL, y_label = NULL,
                              selected_x = NULL, lower_is_better = FALSE) {
  req <- c(x_col, y_col)
  if (!all(req %in% names(df))) stop("Metric curve is missing required columns.")
  dd <- df[is.finite(df[[x_col]]) & is.finite(df[[y_col]]), , drop = FALSE]
  if (!nrow(dd)) stop("Metric curve has no finite values.")
  dd <- dd[order(dd[[x_col]]), , drop = FALSE]
  if (is.null(selected_x)) selected_x <- dd[[x_col]][if (isTRUE(lower_is_better)) which.min(dd[[y_col]]) else which.max(dd[[y_col]])]
  p <- ggplot2::ggplot(dd, ggplot2::aes(.data[[x_col]], .data[[y_col]])) +
    ggplot2::geom_line(linewidth = 0.75, colour = "#4C78A8") +
    ggplot2::geom_point(size = 2.0, colour = "#4C78A8") +
    ggplot2::geom_vline(xintercept = selected_x, linetype = 2, linewidth = 0.5, colour = "#E45756") +
    pub_theme(ctx) +
    ggplot2::labs(title = title, x = x_label %||% x_col, y = y_label %||% y_col)
  list(plot = p, source_data = dd,
       parameters = list(type = "metric_curve", x_col = x_col, y_col = y_col,
                         selected_x = selected_x, lower_is_better = lower_is_better))
}

plot_gene_test_ranking <- function(ctx, df, title = NULL, top_n = 30L, gene_col = "gene",
                                   fdr_col = "FDR", effect_col = NULL) {
  if (!all(c(gene_col, fdr_col) %in% names(df))) stop("Gene-test table requires gene and FDR columns.")
  dd <- df[is.finite(df[[fdr_col]]) & !is.na(df[[gene_col]]), , drop = FALSE]
  if (!nrow(dd)) stop("Gene-test table has no finite adjusted P values.")
  dd$minus_log10_FDR <- -log10(pmax(dd[[fdr_col]], .Machine$double.xmin))
  dd <- head(dd[order(dd$minus_log10_FDR, decreasing = TRUE), , drop = FALSE], top_n)
  dd[[gene_col]] <- factor(as.character(dd[[gene_col]]), levels = rev(as.character(dd[[gene_col]])))
  p <- ggplot2::ggplot(dd, ggplot2::aes(minus_log10_FDR, .data[[gene_col]])) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = minus_log10_FDR, yend = .data[[gene_col]]),
                          linewidth = 0.5, colour = "#B7B7B7") +
    ggplot2::geom_point(size = 2.5, colour = "#4C78A8") +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.45, colour = "#E45756") +
    pub_theme(ctx) + ggplot2::labs(title = title, x = expression(-log[10](FDR)), y = NULL)
  list(plot = p, source_data = dd,
       parameters = list(type = "gene_test_ranking", top_n = top_n, gene_col = gene_col,
                         fdr_col = fdr_col, effect_col = effect_col))
}

plot_group_score_distribution <- function(ctx, df, group_col, score_col, split_col = NULL,
                                          title = NULL, y_label = NULL) {
  req <- c(group_col, score_col, split_col)
  req <- req[!is.null(req) & !is.na(req) & nzchar(req)]
  if (!all(req %in% names(df))) stop("Score-distribution table is incomplete.")
  dd <- df[is.finite(df[[score_col]]) & !is.na(df[[group_col]]), , drop = FALSE]
  if (!nrow(dd)) stop("Score-distribution table has no finite values.")
  aes0 <- if (!is.null(split_col) && split_col %in% names(dd)) {
    ggplot2::aes(.data[[group_col]], .data[[score_col]], fill = .data[[split_col]])
  } else ggplot2::aes(.data[[group_col]], .data[[score_col]], fill = .data[[group_col]])
  p <- ggplot2::ggplot(dd, aes0) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.55) +
    ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.35, alpha = 0.8) +
    pub_theme(ctx) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = title, x = NULL, y = y_label %||% score_col, fill = split_col)
  if (is.null(split_col)) p <- p + ggplot2::guides(fill = "none")
  list(plot = p, source_data = dd,
       parameters = list(type = "group_score_distribution", group_col = group_col,
                         score_col = score_col, split_col = split_col))
}

# Publication Figure Engine v5: exploratory virtual-knockout single-figure templates.
# scTenifoldKnk differential-regulation distance is unsigned; none of these plots
# imply predicted expression up- or down-regulation.

plot_virtual_knockout_ranking <- function(ctx, df, title = "Virtual-knockout regulatory impact",
                                          top_n = 30L, fdr_threshold = 0.05) {
  req <- c("gene", "consensus_impact_score", "sample_recurrence_fraction", "fisher_fdr", "evidence_class")
  if (!all(req %in% names(df))) stop("Virtual-knockout consensus table is incomplete.")
  dd <- df[is.finite(df$consensus_impact_score) & !is.na(df$gene), , drop = FALSE]
  if ("target_gene" %in% names(dd)) dd <- dd[as.character(dd$gene) != as.character(dd$target_gene), , drop = FALSE]
  dd <- head(dd[order(dd$consensus_impact_score, decreasing = TRUE), , drop = FALSE], top_n)
  if (!nrow(dd)) stop("No virtual-knockout genes were available for ranking.")
  dd$gene <- factor(as.character(dd$gene), levels = rev(as.character(dd$gene)))
  dd$significance <- ifelse(is.finite(dd$fisher_fdr) & dd$fisher_fdr < fdr_threshold, "FDR-supported", "Not FDR-supported")
  evidence_levels <- c("Replicate-consistent", "Sample-restricted", "Single-sample-descriptive", "Pooled-descriptive", "Not-supported")
  evidence_palette <- c(
    `Replicate-consistent` = "#2A9D8F", `Sample-restricted` = "#E9C46A",
    `Single-sample-descriptive` = "#F4A261", `Pooled-descriptive` = "#8D99AE", `Not-supported` = "#B8B8B8"
  )
  dd$evidence_class <- factor(as.character(dd$evidence_class), levels = evidence_levels)
  p <- ggplot2::ggplot(dd, ggplot2::aes(x = consensus_impact_score, y = gene)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = consensus_impact_score, yend = gene), colour = "#C7C7C7", linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(size = sample_recurrence_fraction, colour = evidence_class, shape = significance), alpha = 0.95) +
    ggplot2::scale_colour_manual(values = evidence_palette, drop = FALSE, na.value = "#B8B8B8") +
    ggplot2::scale_shape_manual(values = c(`FDR-supported` = 16, `Not FDR-supported` = 1)) +
    ggplot2::scale_size_continuous(range = c(2.2, 6.5), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1)) +
    pub_theme(ctx) +
    ggplot2::labs(title = title, x = "Consensus perturbation-impact score (unsigned)", y = NULL,
                  colour = "Evidence", size = "Sample recurrence", shape = "Statistical support")
  list(plot = p, source_data = dd,
       parameters = list(type = "virtual_knockout_ranking", top_n = top_n, fdr_threshold = fdr_threshold,
                         directionality = "unsigned"))
}

plot_virtual_knockout_sample_support <- function(ctx, sample_df, consensus_df,
                                                 title = "Virtual-knockout support across samples",
                                                 top_n = 30L) {
  req_s <- c("sample_id", "gene", "repeat_recurrence_fraction")
  req_c <- c("gene", "consensus_impact_score")
  if (!all(req_s %in% names(sample_df)) || !all(req_c %in% names(consensus_df))) stop("Virtual-knockout sample-support tables are incomplete.")
  top <- head(as.character(consensus_df$gene[order(consensus_df$consensus_impact_score, decreasing = TRUE)]), top_n)
  dd <- sample_df[sample_df$gene %in% top & is.finite(sample_df$repeat_recurrence_fraction), , drop = FALSE]
  if (!nrow(dd)) stop("No sample-level virtual-knockout support rows were available.")
  dd$gene <- factor(as.character(dd$gene), levels = rev(top[top %in% unique(as.character(dd$gene))]))
  sample_order <- unique(as.character(dd$sample_id))
  dd$sample_id <- factor(as.character(dd$sample_id), levels = sample_order)
  p <- ggplot2::ggplot(dd, ggplot2::aes(x = sample_id, y = gene, fill = repeat_recurrence_fraction)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.18) +
    ggplot2::scale_fill_viridis_c(option = "D", limits = c(0, 1), oob = scales::squish, na.value = "#D0D0D0") +
    pub_theme(ctx) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank()) +
    ggplot2::labs(title = title, x = "Biological sample", y = NULL, fill = "Repeat recurrence")
  list(plot = p, source_data = dd,
       parameters = list(type = "virtual_knockout_sample_support", top_n = top_n, scale_limits = c(0, 1)))
}

plot_virtual_knockout_manifold <- function(ctx, manifold_df, diff_df,
                                           title = "WT–KO manifold displacement",
                                           top_n = 20L, fdr_threshold = 0.05) {
  req_m <- c("gene", "state", "dim1", "dim2")
  req_d <- c("gene", "distance", "p.adj")
  if (!all(req_m %in% names(manifold_df)) || !all(req_d %in% names(diff_df))) stop("Virtual-knockout manifold inputs are incomplete.")
  md <- manifold_df[is.finite(manifold_df$dim1) & is.finite(manifold_df$dim2) & manifold_df$state %in% c("WT", "KO"), , drop = FALSE]
  if (!nrow(md)) stop("No finite virtual-knockout manifold coordinates were available.")
  wt <- md[md$state == "WT", c("gene", "dim1", "dim2"), drop = FALSE]
  ko <- md[md$state == "KO", c("gene", "dim1", "dim2"), drop = FALSE]
  names(wt)[2:3] <- c("wt_dim1", "wt_dim2"); names(ko)[2:3] <- c("ko_dim1", "ko_dim2")
  seg <- merge(wt, ko, by = "gene")
  seg <- merge(seg, diff_df[, intersect(c("gene", "distance", "Z", "FC", "p.value", "p.adj"), names(diff_df)), drop = FALSE], by = "gene", all.x = TRUE)
  seg <- seg[is.finite(seg$distance), , drop = FALSE]
  top <- head(seg$gene[order(seg$p.adj, -seg$distance, na.last = TRUE)], top_n)
  seg$highlight <- seg$gene %in% top
  seg$significant <- is.finite(seg$p.adj) & seg$p.adj < fdr_threshold
  display_seg <- seg[seg$highlight, , drop = FALSE]
  if (!nrow(display_seg)) stop("No manifold displacement rows passed display selection.")
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = md, ggplot2::aes(dim1, dim2, colour = state), size = 0.65, alpha = 0.18) +
    ggplot2::geom_segment(data = display_seg,
      ggplot2::aes(x = wt_dim1, y = wt_dim2, xend = ko_dim1, yend = ko_dim2, linewidth = distance, alpha = significant),
      colour = "#4A4A4A", arrow = grid::arrow(length = grid::unit(1.6, "mm"), type = "closed")) +
    ggplot2::geom_point(data = display_seg, ggplot2::aes(wt_dim1, wt_dim2), shape = 21, fill = "#557FB5", colour = "white", stroke = 0.25, size = 2.2) +
    ggplot2::geom_point(data = display_seg, ggplot2::aes(ko_dim1, ko_dim2), shape = 21, fill = "#E67861", colour = "white", stroke = 0.25, size = 2.2) +
    ggplot2::geom_text(data = head(display_seg[order(display_seg$p.adj, -display_seg$distance), , drop = FALSE], min(12L, nrow(display_seg))),
                       ggplot2::aes(ko_dim1, ko_dim2, label = gene), size = 2.5, check_overlap = TRUE, vjust = -0.7) +
    ggplot2::scale_colour_manual(values = c(WT = "#557FB5", KO = "#E67861")) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 0.85, `FALSE` = 0.38), guide = "none") +
    ggplot2::scale_linewidth_continuous(range = c(0.25, 1.6)) +
    ggplot2::coord_equal() + embedding_theme(ctx) +
    ggplot2::labs(title = title, x = "Manifold dimension 1", y = "Manifold dimension 2", colour = "Network state", linewidth = "Distance")
  list(plot = p, source_data = list(coordinates = md, displacement = seg, displayed = display_seg),
       parameters = list(type = "virtual_knockout_manifold", top_n = top_n, fdr_threshold = fdr_threshold,
                         arrow_meaning = "WT-to-KO network-state displacement; not expression direction"))
}

plot_virtual_knockout_network <- function(ctx, edges, target_gene,
                                          title = "Target-centred WT network context", top_n = 80L) {
  req <- c("from", "to", "median_weight", "median_absolute_weight", "edge_recurrence_fraction")
  if (!all(req %in% names(edges))) stop("Virtual-knockout consensus-edge table is incomplete.")
  ee <- edges[is.finite(edges$median_absolute_weight) & edges$median_absolute_weight > 0, , drop = FALSE]
  ee <- head(ee[order(ee$target_outgoing, ee$edge_recurrence_fraction, ee$median_absolute_weight, decreasing = TRUE), , drop = FALSE], top_n)
  if (!nrow(ee)) stop("No virtual-knockout network edges were available.")
  nodes <- sort(unique(c(as.character(ee$from), as.character(ee$to), target_gene)))
  theta <- seq(0, 2 * pi, length.out = length(nodes) + 1L)[-1L]
  nd <- data.frame(node = nodes, x = cos(theta), y = sin(theta), stringsAsFactors = FALSE)
  if (target_gene %in% nd$node) { nd$x[nd$node == target_gene] <- 0; nd$y[nd$node == target_gene] <- 0 }
  nd$is_target <- nd$node == target_gene
  nd$degree <- vapply(nd$node, function(z) sum(ee$from == z | ee$to == z), numeric(1))
  ee$x <- nd$x[match(ee$from, nd$node)]; ee$y <- nd$y[match(ee$from, nd$node)]
  ee$xend <- nd$x[match(ee$to, nd$node)]; ee$yend <- nd$y[match(ee$to, nd$node)]
  ee$edge_sign <- ifelse(ee$median_weight >= 0, "Positive WT edge", "Negative WT edge")
  p <- ggplot2::ggplot() +
    ggplot2::geom_curve(data = ee,
      ggplot2::aes(x, y, xend = xend, yend = yend, linewidth = edge_recurrence_fraction, colour = edge_sign, alpha = median_absolute_weight),
      curvature = 0.14, arrow = grid::arrow(length = grid::unit(1.4, "mm"), type = "closed")) +
    ggplot2::geom_point(data = nd, ggplot2::aes(x, y, size = degree, fill = is_target), shape = 21, colour = "white", stroke = 0.45) +
    ggplot2::geom_text(data = nd, ggplot2::aes(x, y, label = node), size = 2.5, vjust = -0.85, check_overlap = TRUE) +
    ggplot2::scale_colour_manual(values = c(`Positive WT edge` = "#557FB5", `Negative WT edge` = "#E67861")) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#54278F", `FALSE` = "#BDBDBD"), guide = "none") +
    ggplot2::scale_linewidth_continuous(range = c(0.25, 1.8), limits = c(0, 1)) +
    ggplot2::scale_alpha_continuous(range = c(0.25, 0.85), guide = "none") +
    ggplot2::scale_size_continuous(range = c(3.2, 8), guide = "none") +
    ggplot2::coord_equal(clip = "off") + embedding_theme(ctx) + ggplot2::theme(legend.position = "right") +
    ggplot2::labs(title = title, colour = "WT edge sign", linewidth = "Run recurrence")
  list(plot = p, source_data = list(edges = ee, nodes = nd),
       parameters = list(type = "virtual_knockout_network", target_gene = target_gene, top_n = top_n,
                         network_state = "WT", causal_claim = FALSE))
}

plot_virtual_knockout_pathways <- function(ctx, df, title = "Pathways enriched among impacted genes", top_n = 20L) {
  req <- c("pathway", "FDR", "overlap", "enrichment_ratio")
  if (!all(req %in% names(df))) stop("Virtual-knockout pathway table is incomplete.")
  dd <- df[is.finite(df$FDR) & is.finite(df$enrichment_ratio), , drop = FALSE]
  dd <- head(dd[order(dd$FDR, -dd$overlap), , drop = FALSE], top_n)
  if (!nrow(dd)) stop("No virtual-knockout pathways were available.")
  dd$minus_log10_fdr <- -log10(pmax(dd$FDR, .Machine$double.xmin))
  dd$pathway <- factor(as.character(dd$pathway), levels = rev(as.character(dd$pathway)))
  p <- ggplot2::ggplot(dd, ggplot2::aes(enrichment_ratio, pathway, size = overlap, colour = minus_log10_fdr)) +
    ggplot2::geom_point(alpha = 0.92) + ggplot2::scale_colour_viridis_c(option = "D") +
    ggplot2::scale_size_area(max_size = 8) + pub_theme(ctx) +
    ggplot2::labs(title = title, x = "Overlap / pathway size", y = NULL, size = "Overlap genes", colour = expression(-log[10](FDR)))
  list(plot = p, source_data = dd,
       parameters = list(type = "virtual_knockout_pathways", top_n = top_n, enrichment = "over-representation analysis"))
}

plot_virtual_knockout_target_expression <- function(ctx, df,
                                                     title = "Target-expression eligibility by sample") {
  req <- c("sample_id", "target_expressing_fraction", "n_cells", "n_target_expressing_cells")
  if (!all(req %in% names(df))) stop("Virtual-knockout target-expression audit is incomplete.")
  dd <- df[!is.na(df$sample_id) & is.finite(df$target_expressing_fraction), , drop = FALSE]
  if (!nrow(dd)) stop("No target-expression audit rows were available.")
  if (!"eligible_for_knockout" %in% names(dd)) dd$eligible_for_knockout <- NA
  dd$sample_id <- factor(as.character(dd$sample_id), levels = rev(unique(as.character(dd$sample_id))))
  p <- ggplot2::ggplot(dd, ggplot2::aes(target_expressing_fraction, sample_id, size = n_target_expressing_cells, fill = eligible_for_knockout)) +
    ggplot2::geom_point(shape = 21, colour = "#4A4A4A", stroke = 0.35, alpha = 0.9) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, max(1, max(dd$target_expressing_fraction, na.rm = TRUE)))) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#2A9D8F", `FALSE` = "#D0D0D0"), na.value = "#D0D0D0") +
    ggplot2::scale_size_area(max_size = 8) + pub_theme(ctx) +
    ggplot2::labs(title = title, x = "Target-expressing cell fraction", y = "Biological sample",
                  size = "Expressing cells", fill = "Eligible")
  list(plot = p, source_data = dd,
       parameters = list(type = "virtual_knockout_target_expression", eligibility_is_predeclared = TRUE))
}
