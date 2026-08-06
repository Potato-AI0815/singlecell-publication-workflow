# Publication Figure Engine v4 — authoritative independent figures only
# Deterministic plotting helpers for single-cell figures. This file is sourced by stage 11.

sanitize_stem <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  sub("^_+|_+$", "", x)
}

figure_cfg <- function(ctx, name, default = NULL) {
  x <- ctx$cfg$figures[[name]]
  if (is.null(x)) default else x
}

resolve_reduction <- function(obj, requested = "auto") {
  available <- names(obj@reductions)
  if (!length(available)) stop("No dimensional reduction is available in the object.")
  if (!is.null(requested) && length(requested) == 1L && !identical(requested, "auto")) {
    if (!requested %in% available) stop("Requested reduction is absent: ", requested)
    return(requested)
  }
  preferred <- c("umap", "umap.integrated", "integrated.umap", "harmony.umap", "tsne", "tSNE", "pca")
  hit <- preferred[preferred %in% available]
  if (length(hit)) hit[[1]] else available[[1]]
}

embedding_source <- function(obj, reduction, metadata_cols = character()) {
  emb <- Seurat::Embeddings(obj, reduction = reduction)[, 1:2, drop = FALSE]
  out <- data.frame(cell_id = rownames(emb), embedding_1 = emb[, 1], embedding_2 = emb[, 2], stringsAsFactors = FALSE)
  md <- obj[[]]
  metadata_cols <- intersect(unique(metadata_cols), colnames(md))
  if (length(metadata_cols)) {
    m <- md[rownames(emb), metadata_cols, drop = FALSE]
    rownames(m) <- NULL
    out <- cbind(out, m)
  }
  out
}

get_feature_vector <- function(obj, feature, assay = "RNA", layer = "data") {
  md <- obj[[]]
  if (feature %in% colnames(md)) {
    out <- as.numeric(md[[feature]])
    names(out) <- rownames(md)
    return(out)
  }
  obj2 <- join_layers_safe(obj, assay)
  layers <- SeuratObject::Layers(obj2[[assay]])
  if (!layer %in% layers && layer == "data") {
    obj2 <- Seurat::NormalizeData(obj2, assay = assay, verbose = FALSE)
  }
  mat <- get_layer_safe(obj2, assay = assay, layer = layer)
  if (!feature %in% rownames(mat)) stop("Feature or metadata column is absent: ", feature)
  out <- as.numeric(mat[feature, colnames(obj2), drop = TRUE])
  names(out) <- colnames(obj2)
  out
}

stratified_downsample <- function(df, strata = NULL, max_points = Inf, seed = 1L) {
  if (!is.finite(max_points) || nrow(df) <= max_points) return(df)
  set.seed(seed)
  if (is.null(strata) || !strata %in% names(df)) return(df[sample.int(nrow(df), max_points), , drop = FALSE])
  spl <- split(seq_len(nrow(df)), as.character(df[[strata]]), drop = TRUE)
  quota <- max(1L, floor(max_points / max(1L, length(spl))))
  idx <- unlist(lapply(spl, function(z) if (length(z) <= quota) z else sample(z, quota)), use.names = FALSE)
  if (length(idx) < max_points) {
    remainder <- setdiff(seq_len(nrow(df)), idx)
    if (length(remainder)) idx <- c(idx, sample(remainder, min(length(remainder), max_points - length(idx))))
  }
  df[idx, , drop = FALSE]
}

publication_point_layer <- function(mapping, data, size = 0.25, alpha = 0.75, raster = TRUE, raster_dpi = 300) {
  if (isTRUE(raster) && requireNamespace("ggrastr", quietly = TRUE)) {
    return(ggrastr::geom_point_rast(mapping = mapping, data = data, size = size, alpha = alpha, raster.dpi = raster_dpi, stroke = 0))
  }
  ggplot2::geom_point(mapping = mapping, data = data, size = size, alpha = alpha, stroke = 0)
}

axis_arrow_layers <- function(df, x_label = "UMAP 1", y_label = "UMAP 2", fraction = 0.14, inset = 0.04, linewidth = 0.45, text_size = 3.0) {
  xr <- range(df$embedding_1, finite = TRUE)
  yr <- range(df$embedding_2, finite = TRUE)
  dx <- diff(xr); dy <- diff(yr)
  if (!is.finite(dx) || dx == 0) dx <- 1
  if (!is.finite(dy) || dy == 0) dy <- 1
  x0 <- xr[1] + inset * dx
  y0 <- yr[1] + inset * dy
  list(
    ggplot2::annotate("segment", x = x0, xend = x0 + fraction * dx, y = y0, yend = y0,
                      linewidth = linewidth, colour = "black", arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed")),
    ggplot2::annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + fraction * dy,
                      linewidth = linewidth, colour = "black", arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed")),
    ggplot2::annotate("text", x = x0 + fraction * dx / 2, y = y0 - 0.025 * dy, label = x_label, size = text_size, vjust = 1),
    ggplot2::annotate("text", x = x0 - 0.018 * dx, y = y0 + fraction * dy / 2, label = y_label, size = text_size, angle = 90, vjust = 1)
  )
}

label_positions <- function(df, group_col, minimum_cells = 20L) {
  groups <- split(df, as.character(df[[group_col]]), drop = TRUE)
  rows <- lapply(names(groups), function(g) {
    z <- groups[[g]]
    if (nrow(z) < minimum_cells) return(NULL)
    data.frame(label = g, embedding_1 = stats::median(z$embedding_1, na.rm = TRUE),
               embedding_2 = stats::median(z$embedding_2, na.rm = TRUE), n = nrow(z), stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

hull_components <- function(z, method = "auto", minimum_component_cells = 30L) {
  n <- nrow(z)
  if (n < minimum_component_cells) return(rep(NA_integer_, n))
  if (method %in% c("auto", "hdbscan") && requireNamespace("dbscan", quietly = TRUE) && n >= 100L) {
    xy <- scale(as.matrix(z[, c("embedding_1", "embedding_2")]))
    min_pts <- max(10L, min(50L, floor(n * 0.03)))
    cl <- tryCatch(dbscan::hdbscan(xy, minPts = min_pts)$cluster, error = function(e) rep(1L, n))
    cl[cl == 0L] <- NA_integer_
    if (sum(!is.na(cl)) >= minimum_component_cells) return(cl)
  }
  rep(1L, n)
}

make_embedding_hulls <- function(df, group_col, minimum_cells = 50L, component_method = "auto", concavity = 2) {
  groups <- split(df, as.character(df[[group_col]]), drop = TRUE)
  out <- list(); k <- 0L
  for (g in names(groups)) {
    z <- groups[[g]]
    if (nrow(z) < minimum_cells) next
    z$.component <- hull_components(z, component_method, minimum_component_cells = max(20L, floor(minimum_cells / 2)))
    comps <- split(z[!is.na(z$.component), , drop = FALSE], z$.component[!is.na(z$.component)], drop = TRUE)
    for (cc in names(comps)) {
      w <- comps[[cc]]
      if (nrow(w) < max(20L, floor(minimum_cells / 2))) next
      pts <- unique(w[, c("embedding_1", "embedding_2"), drop = FALSE])
      if (nrow(pts) < 3L) next
      if (nrow(pts) > 5000L) pts <- pts[unique(round(seq(1, nrow(pts), length.out = 5000))), , drop = FALSE]
      poly <- NULL
      if (requireNamespace("concaveman", quietly = TRUE) && nrow(pts) >= 4L) {
        poly <- tryCatch(as.data.frame(concaveman::concaveman(as.matrix(pts), concavity = concavity, length_threshold = 0)), error = function(e) NULL)
        if (!is.null(poly) && ncol(poly) >= 2L) names(poly)[1:2] <- c("embedding_1", "embedding_2")
      }
      if (is.null(poly) || nrow(poly) < 3L) {
        idx <- grDevices::chull(pts$embedding_1, pts$embedding_2)
        poly <- pts[c(idx, idx[1]), , drop = FALSE]
      }
      k <- k + 1L
      poly$group_value <- g
      poly$component <- paste(g, cc, sep = "__")
      out[[k]] <- poly
    }
  }
  if (!length(out)) data.frame() else do.call(rbind, out)
}

embedding_theme <- function(ctx) {
  fam <- choose_font(ctx$cfg$figures$font_family %||% "Arial")
  bs <- ctx$cfg$figures$base_size_pt %||% 10
  ggplot2::theme_void(base_family = fam, base_size = bs) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.title = ggplot2::element_text(face = "bold", size = max(8, bs - 1)),
      legend.text = ggplot2::element_text(size = max(8, bs - 2)),
      plot.title = ggplot2::element_text(face = "bold", size = bs + 2, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = max(8, bs - 1), colour = "#444444"),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
}

plot_embedding_discrete <- function(ctx, obj, group_by, reduction = "auto", title = NULL,
                                    labels = TRUE, hulls = FALSE, palette = NULL,
                                    max_points = NULL, point_size = NULL, point_alpha = NULL) {
  if (!group_by %in% colnames(obj[[]])) stop("Embedding group column is absent: ", group_by)
  red <- resolve_reduction(obj, reduction)
  df_all <- embedding_source(obj, red, c(group_by, "sample_id", "condition", "cluster_raw", "cell_type_l1"))
  df_all$group_value <- as.character(df_all[[group_by]])
  df_all$group_value[is.na(df_all$group_value) | !nzchar(df_all$group_value)] <- "NA"
  lev <- sort(unique(df_all$group_value))
  df_all$group_value <- factor(df_all$group_value, levels = lev)
  palette <- palette %||% setNames(categorical_palette(length(lev)), lev)
  palette <- palette[lev]
  max_points <- max_points %||% figure_cfg(ctx, "max_embedding_points", 80000)
  point_size <- point_size %||% figure_cfg(ctx, "embedding_point_size", if (nrow(df_all) > 50000) 0.12 else 0.22)
  point_alpha <- point_alpha %||% figure_cfg(ctx, "embedding_point_alpha", 0.78)
  plot_df <- stratified_downsample(df_all, "group_value", max_points, ctx$cfg$runtime$seed %||% 1L)
  hull_df <- if (isTRUE(hulls)) make_embedding_hulls(df_all, "group_value",
                                                     minimum_cells = figure_cfg(ctx, "hull_minimum_cells", 80),
                                                     component_method = figure_cfg(ctx, "hull_component_method", "auto"),
                                                     concavity = figure_cfg(ctx, "hull_concavity", 2)) else data.frame()
  p <- ggplot2::ggplot()
  if (nrow(hull_df)) {
    p <- p + ggplot2::geom_polygon(data = hull_df,
      ggplot2::aes(x = embedding_1, y = embedding_2, group = component, fill = group_value),
      alpha = figure_cfg(ctx, "hull_alpha", 0.14), colour = NA, show.legend = FALSE)
  }
  p <- p + publication_point_layer(ggplot2::aes(x = embedding_1, y = embedding_2, colour = group_value),
                                   plot_df, size = point_size, alpha = point_alpha,
                                   raster = nrow(plot_df) >= figure_cfg(ctx, "raster_threshold", 30000),
                                   raster_dpi = figure_cfg(ctx, "raster_dpi", 300)) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE, na.value = "#BDBDBD") +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE, na.value = "#BDBDBD") +
    ggplot2::coord_fixed(clip = "off") + embedding_theme(ctx) +
    ggplot2::labs(title = title %||% paste("Embedding by", group_by), colour = gsub("_", " ", group_by))
  if (isTRUE(labels)) {
    lab <- label_positions(df_all, "group_value", minimum_cells = figure_cfg(ctx, "label_minimum_cells", 30))
    if (nrow(lab)) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_label_repel(data = lab, ggplot2::aes(embedding_1, embedding_2, label = label),
          inherit.aes = FALSE, size = 2.7, label.size = 0.15, label.padding = grid::unit(0.12, "lines"),
          box.padding = 0.25, point.padding = 0.1, min.segment.length = 0, seed = ctx$cfg$runtime$seed %||% 1L,
          colour = "black", fill = grDevices::adjustcolor("white", alpha.f = 0.82), max.overlaps = Inf)
      } else {
        p <- p + ggplot2::geom_label(data = lab, ggplot2::aes(embedding_1, embedding_2, label = label),
          inherit.aes = FALSE, size = 2.7, label.size = 0.15, fill = grDevices::adjustcolor("white", alpha.f = 0.82))
      }
    }
  }
  axis_names <- if (grepl("tsne", red, ignore.case = TRUE)) c("tSNE 1", "tSNE 2") else if (grepl("umap", red, ignore.case = TRUE)) c("UMAP 1", "UMAP 2") else c(paste0(red, " 1"), paste0(red, " 2"))
  if (identical(figure_cfg(ctx, "embedding_axis_style", "arrows"), "arrows")) p <- p + axis_arrow_layers(df_all, axis_names[1], axis_names[2])
  list(plot = p, source_data = list(cells = df_all, hulls = hull_df),
       parameters = list(type = "embedding_discrete", reduction = red, group_by = group_by,
                         labels = labels, hulls = hulls, point_size = point_size, point_alpha = point_alpha,
                         max_points_plot = max_points, n_cells_total = nrow(df_all), n_cells_plotted = nrow(plot_df), palette = as.list(palette)))
}

plot_embedding_continuous <- function(ctx, obj, feature, reduction = "auto", title = NULL,
                                      midpoint = 0, quantile_limits = c(0.01, 0.99), max_points = NULL) {
  red <- resolve_reduction(obj, reduction)
  df_all <- embedding_source(obj, red, c("sample_id", "condition", "cluster_raw", "cell_type_l1"))
  names_vec <- colnames(obj)
  values <- get_feature_vector(obj, feature, assay = "RNA", layer = "data")
  if (length(values) != length(names_vec)) stop("Feature vector length differs from object cells: ", feature)
  names(values) <- names_vec
  df_all$value_raw <- values[df_all$cell_id]
  lim <- stats::quantile(df_all$value_raw[is.finite(df_all$value_raw)], probs = quantile_limits, na.rm = TRUE, names = FALSE)
  if (length(lim) != 2L || any(!is.finite(lim)) || lim[1] == lim[2]) lim <- range(df_all$value_raw, finite = TRUE)
  df_all$value_display <- pmin(pmax(df_all$value_raw, lim[1]), lim[2])
  max_points <- max_points %||% figure_cfg(ctx, "max_embedding_points", 80000)
  plot_df <- stratified_downsample(df_all, NULL, max_points, ctx$cfg$runtime$seed %||% 1L)
  symmetric <- is.finite(midpoint) && lim[1] < midpoint && lim[2] > midpoint
  maxabs <- max(abs(lim - midpoint))
  scale_obj <- if (symmetric) {
    ggplot2::scale_colour_gradient2(low = figure_cfg(ctx, "continuous_low", "#E67861"), mid = "#F7F7F7",
      high = figure_cfg(ctx, "continuous_high", "#557FB5"), midpoint = midpoint,
      limits = c(midpoint - maxabs, midpoint + maxabs), oob = scales::squish, na.value = "#D0D0D0")
  } else {
    ggplot2::scale_colour_viridis_c(option = "D", limits = lim, oob = scales::squish, na.value = "#D0D0D0")
  }
  p <- ggplot2::ggplot() +
    publication_point_layer(ggplot2::aes(embedding_1, embedding_2, colour = value_display), plot_df,
      size = figure_cfg(ctx, "embedding_point_size", if (nrow(df_all) > 50000) 0.12 else 0.22),
      alpha = figure_cfg(ctx, "embedding_point_alpha", 0.8),
      raster = nrow(plot_df) >= figure_cfg(ctx, "raster_threshold", 30000),
      raster_dpi = figure_cfg(ctx, "raster_dpi", 300)) + scale_obj +
    ggplot2::coord_fixed(clip = "off") + embedding_theme(ctx) +
    ggplot2::labs(title = title %||% feature, colour = feature)
  axis_names <- if (grepl("tsne", red, ignore.case = TRUE)) c("tSNE 1", "tSNE 2") else if (grepl("umap", red, ignore.case = TRUE)) c("UMAP 1", "UMAP 2") else c(paste0(red, " 1"), paste0(red, " 2"))
  if (identical(figure_cfg(ctx, "embedding_axis_style", "arrows"), "arrows")) p <- p + axis_arrow_layers(df_all, axis_names[1], axis_names[2])
  list(plot = p, source_data = df_all,
       parameters = list(type = "embedding_continuous", feature = feature, reduction = red, midpoint = midpoint,
                         quantile_limits = quantile_limits, display_limits = lim, n_cells_total = nrow(df_all), n_cells_plotted = nrow(plot_df)))
}

split_violin_polygon_data <- function(df, x_col, split_col, value_col, width = 0.9, n_density = 256L, trim = TRUE) {
  xlev <- unique(as.character(df[[x_col]]))
  slev <- unique(as.character(df[[split_col]]))
  slev <- slev[!is.na(slev)]
  if (length(slev) != 2L) stop("Split violin requires exactly two non-missing split levels; found ", length(slev))
  rows <- list(); k <- 0L
  for (ix in seq_along(xlev)) {
    for (is in seq_along(slev)) {
      vals <- as.numeric(df[df[[x_col]] == xlev[ix] & df[[split_col]] == slev[is], value_col])
      vals <- vals[is.finite(vals)]
      if (length(vals) < 3L || length(unique(vals)) < 2L) next
      den <- stats::density(vals, n = n_density, from = if (trim) min(vals) else NULL, to = if (trim) max(vals) else NULL)
      scale <- if (max(den$y) > 0) den$y / max(den$y) * width / 2 else rep(0, length(den$y))
      center <- ix
      edge <- if (is == 1L) center - scale else center + scale
      poly <- data.frame(x = c(rep(center, length(den$x)), rev(edge)),
                         y = c(den$x, rev(den$x)), x_level = xlev[ix], split_level = slev[is],
                         polygon_id = paste(ix, is, sep = "__"), stringsAsFactors = FALSE)
      k <- k + 1L; rows[[k]] <- poly
    }
  }
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

plot_split_violin_box <- function(ctx, df, x_col, split_col, value_col, title = NULL, y_label = NULL, palette = NULL) {
  req <- c(x_col, split_col, value_col)
  if (!all(req %in% names(df))) stop("Split violin data are missing: ", paste(setdiff(req, names(df)), collapse = ", "))
  dd <- df[complete.cases(df[, req, drop = FALSE]), req, drop = FALSE]
  dd[[x_col]] <- factor(as.character(dd[[x_col]]), levels = unique(as.character(dd[[x_col]])))
  dd[[split_col]] <- factor(as.character(dd[[split_col]]), levels = unique(as.character(dd[[split_col]])))
  lev <- levels(dd[[split_col]])
  if (length(lev) != 2L) stop("Split violin requires exactly two split levels.")
  palette <- palette %||% setNames(c("#6FA3D2", "#EF7D71"), lev)
  poly <- split_violin_polygon_data(dd, x_col, split_col, value_col,
                                    width = figure_cfg(ctx, "split_violin_width", 0.88),
                                    n_density = figure_cfg(ctx, "split_violin_density_points", 256))
  box_df <- dd
  box_df$.x_num <- as.numeric(box_df[[x_col]]) + ifelse(box_df[[split_col]] == lev[1], -0.10, 0.10)
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = poly, ggplot2::aes(x, y, group = polygon_id, fill = split_level),
                          colour = NA, alpha = 0.72) +
    ggplot2::geom_boxplot(data = box_df, ggplot2::aes(x = .x_num, y = .data[[value_col]], group = interaction(.data[[x_col]], .data[[split_col]]), colour = .data[[split_col]]),
                          width = 0.12, outlier.shape = NA, linewidth = 0.35, alpha = 0.55, fill = "white") +
    ggplot2::scale_x_continuous(breaks = seq_along(levels(dd[[x_col]])), labels = levels(dd[[x_col]]), expand = ggplot2::expansion(mult = c(0.04, 0.06))) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE, guide = "none") +
    pub_theme(ctx) +
    ggplot2::labs(title = title %||% paste(value_col, "by", x_col), x = NULL, y = y_label %||% value_col, fill = gsub("_", " ", split_col)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = if (length(levels(dd[[x_col]])) > 8) 45 else 0, hjust = if (length(levels(dd[[x_col]])) > 8) 1 else 0.5))
  list(plot = p, source_data = list(values = dd, density_polygons = poly),
       parameters = list(type = "split_violin_box", x_col = x_col, split_col = split_col, value_col = value_col, split_levels = lev, palette = as.list(palette)))
}

aggregate_expression_long <- function(obj, genes, group_by, assay = "RNA", layer = "data") {
  if (!group_by %in% colnames(obj[[]])) stop("Grouping column absent: ", group_by)
  obj2 <- join_layers_safe(obj, assay)
  if (!layer %in% SeuratObject::Layers(obj2[[assay]]) && layer == "data") obj2 <- Seurat::NormalizeData(obj2, assay = assay, verbose = FALSE)
  mat <- get_layer_safe(obj2, assay, layer)
  genes <- unique(intersect(genes, rownames(mat)))
  if (!length(genes)) stop("None of the requested marker genes are present.")
  groups <- as.character(obj2[[group_by]][, 1])
  lev <- unique(groups)
  rows <- list(); k <- 0L
  for (g in lev) {
    idx <- which(groups == g)
    if (!length(idx)) next
    sub <- mat[genes, idx, drop = FALSE]
    avg <- Matrix::rowMeans(sub)
    pct <- Matrix::rowMeans(sub > 0) * 100
    k <- k + 1L
    rows[[k]] <- data.frame(gene = genes, group = g, average_expression = as.numeric(avg), percent_expressed = as.numeric(pct), n_cells = length(idx), stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$scaled_average <- ave(out$average_expression, out$gene, FUN = function(x) {
    z <- as.numeric(scale(x)); z[!is.finite(z)] <- 0; pmax(-2, pmin(2, z))
  })
  out
}

plot_marker_dotplot <- function(ctx, obj, genes, group_by = "cell_type_l1", title = "Canonical marker evidence") {
  dat <- aggregate_expression_long(obj, genes, group_by)
  dat$gene <- factor(dat$gene, levels = rev(unique(genes[genes %in% dat$gene])))
  dat$group <- factor(dat$group, levels = unique(dat$group))
  p <- ggplot2::ggplot(dat, ggplot2::aes(group, gene)) +
    ggplot2::geom_point(ggplot2::aes(size = percent_expressed, colour = scaled_average)) +
    ggplot2::scale_size_continuous(range = c(0.5, 5.5), limits = c(0, 100), breaks = c(25, 50, 75, 100)) +
    ggplot2::scale_colour_gradient(low = "#E3D7F2", high = "#5A189A", limits = c(-2, 2), oob = scales::squish) +
    pub_theme(ctx) + ggplot2::labs(title = title, x = NULL, y = NULL, size = "Percent\nexpressed", colour = "Average\nexpression\nscaled") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid.major = ggplot2::element_line(colour = "#EEEEEE", linewidth = 0.25))
  list(plot = p, source_data = dat, parameters = list(type = "marker_dotplot", genes = genes, group_by = group_by, colour_limits = c(-2, 2), size_limits = c(0, 100)))
}

plot_marker_heatmap <- function(ctx, obj, genes, group_by = "cell_type_l1", title = "Average marker expression") {
  dat <- aggregate_expression_long(obj, genes, group_by)
  dat$gene <- factor(dat$gene, levels = rev(unique(genes[genes %in% dat$gene])))
  dat$group <- factor(dat$group, levels = unique(dat$group))
  p <- ggplot2::ggplot(dat, ggplot2::aes(group, gene, fill = scaled_average)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.12) +
    ggplot2::scale_fill_gradient2(low = "#557FB5", mid = "white", high = "#E67861", midpoint = 0, limits = c(-2, 2), oob = scales::squish) +
    pub_theme(ctx) + ggplot2::labs(title = title, x = NULL, y = NULL, fill = "Z score") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  list(plot = p, source_data = dat, parameters = list(type = "marker_heatmap", genes = genes, group_by = group_by, fill_limits = c(-2, 2)))
}

plot_composition_stacked <- function(ctx, dat, title = "Recovered-cell composition") {
  req <- c("sample_id", "cell_type", "percentage")
  if (!all(req %in% names(dat))) stop("Composition table is incomplete.")
  lev <- sort(unique(dat$cell_type)); pal <- setNames(categorical_palette(length(lev)), lev)
  dat$cell_type <- factor(dat$cell_type, levels = lev)
  p <- ggplot2::ggplot(dat, ggplot2::aes(sample_id, percentage, fill = cell_type)) +
    ggplot2::geom_col(width = 0.86) + ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
    pub_theme(ctx) + ggplot2::labs(title = title, x = "Sample", y = "Cell composition (%)", fill = "Cell type") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  list(plot = p, source_data = dat, parameters = list(type = "composition_stacked", palette = as.list(pal)))
}

plot_composition_by_condition <- function(ctx, dat, title = "Sample-level cell-type proportions") {
  req <- c("sample_id", "cell_type", "condition", "percentage")
  if (!all(req %in% names(dat))) stop("Composition table is missing condition metadata.")
  p <- ggplot2::ggplot(dat, ggplot2::aes(condition, percentage, group = sample_id)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.55, fill = "white", linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(fill = condition), shape = 21, size = 1.8, alpha = 0.85, position = ggplot2::position_jitter(width = 0.08, height = 0)) +
    ggplot2::facet_wrap(~cell_type, scales = "free_y") + pub_theme(ctx, 8) +
    ggplot2::labs(title = title, x = NULL, y = "Cell composition (%)", fill = "Condition")
  list(plot = p, source_data = dat, parameters = list(type = "composition_by_condition"))
}

plot_paired_composition <- function(ctx, dat, title = "Paired cell-type composition") {
  req <- c("patient_id", "cell_type", "condition", "percentage")
  if (!all(req %in% names(dat))) stop("Paired composition requires patient_id, condition and percentage.")
  d <- dat[!is.na(dat$patient_id) & nzchar(as.character(dat$patient_id)), , drop = FALSE]
  p <- ggplot2::ggplot(d, ggplot2::aes(condition, percentage, group = patient_id)) +
    ggplot2::geom_line(colour = "#777777", alpha = 0.65, linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(fill = condition), shape = 21, size = 1.9) +
    ggplot2::facet_wrap(~cell_type, scales = "free_y") + pub_theme(ctx, 8) +
    ggplot2::labs(title = title, x = NULL, y = "Cell composition (%)", fill = "Condition")
  list(plot = p, source_data = d, parameters = list(type = "paired_composition"))
}

plot_volcano_single <- function(ctx, dat, cell_type, title = NULL, label_n = 10L) {
  d <- dat[dat$cell_type == cell_type, , drop = FALSE]
  if (!nrow(d)) stop("No DE rows for cell type: ", cell_type)
  fdr_thr <- ctx$cfg$comparisons$fdr_threshold %||% 0.05
  lfc_thr <- ctx$cfg$comparisons$logfc_display_threshold %||% 0.5
  d$neglog10FDR <- -log10(pmax(d$FDR, 1e-300))
  d$status <- ifelse(d$FDR < fdr_thr & d$logFC >= lfc_thr, "Up", ifelse(d$FDR < fdr_thr & d$logFC <= -lfc_thr, "Down", "NS"))
  ord <- order(d$FDR, -abs(d$logFC), na.last = NA)
  label_idx <- head(ord[d$status[ord] != "NS"], label_n)
  labs <- d[label_idx, , drop = FALSE]
  pal <- c(Down = "#557FB5", NS = "#BDBDBD", Up = "#E67861")
  p <- ggplot2::ggplot(d, ggplot2::aes(logFC, neglog10FDR, colour = status)) +
    ggplot2::geom_point(size = 0.65, alpha = 0.72) +
    ggplot2::geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2, linewidth = 0.3, colour = "#777777") +
    ggplot2::geom_hline(yintercept = -log10(fdr_thr), linetype = 2, linewidth = 0.3, colour = "#777777") +
    ggplot2::scale_colour_manual(values = pal, breaks = c("Up", "Down", "NS")) + pub_theme(ctx) +
    ggplot2::labs(title = title %||% paste0(cell_type, ": pseudobulk differential expression"), x = "log2 fold change", y = "-log10 FDR", colour = NULL)
  if (nrow(labs) && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = labs, ggplot2::aes(label = gene), size = 2.7, max.overlaps = Inf,
      min.segment.length = 0, box.padding = 0.25, seed = ctx$cfg$runtime$seed %||% 1L, show.legend = FALSE)
  }
  list(plot = p, source_data = d, parameters = list(type = "volcano", cell_type = cell_type, fdr_threshold = fdr_thr, logfc_threshold = lfc_thr, label_n = label_n))
}

plot_enrichment_single <- function(ctx, dat, cell_type, title = NULL, top_n = 15L) {
  d <- dat[dat$cell_type == cell_type & is.finite(dat$NES), , drop = FALSE]
  if (!nrow(d)) stop("No enrichment rows for cell type: ", cell_type)
  d <- d[order(d$padj, -abs(d$NES)), , drop = FALSE]
  d <- head(d, top_n)
  d$pathway_display <- factor(d$pathway, levels = rev(d$pathway[order(d$NES)]))
  p <- ggplot2::ggplot(d, ggplot2::aes(NES, pathway_display)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#888888", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(size = size, colour = -log10(pmax(padj, 1e-300)))) +
    ggplot2::scale_colour_viridis_c(option = "C") + ggplot2::scale_size_continuous(range = c(2, 6)) +
    pub_theme(ctx) + ggplot2::labs(title = title %||% paste0(cell_type, ": pathway enrichment"), x = "Normalized enrichment score", y = NULL, size = "Gene-set size", colour = "-log10 FDR")
  list(plot = p, source_data = d, parameters = list(type = "enrichment_dotplot", cell_type = cell_type, top_n = top_n))
}

plot_trajectory_embedding <- function(ctx, cell_data, curve_data, node_data = NULL, title = "Trajectory on embedding") {
  reqc <- c("cell_id", "embedding_1", "embedding_2", "pseudotime")
  reqp <- c("embedding_1", "embedding_2", "lineage", "curve_order")
  if (!all(reqc %in% names(cell_data))) stop("Trajectory cell table is incomplete.")
  if (!all(reqp %in% names(curve_data))) stop("Trajectory curve table is incomplete.")
  p <- ggplot2::ggplot() +
    publication_point_layer(ggplot2::aes(embedding_1, embedding_2, colour = pseudotime), cell_data,
      size = figure_cfg(ctx, "embedding_point_size", 0.18), alpha = 0.68,
      raster = nrow(cell_data) >= figure_cfg(ctx, "raster_threshold", 30000), raster_dpi = figure_cfg(ctx, "raster_dpi", 300)) +
    ggplot2::geom_path(data = curve_data[order(curve_data$lineage, curve_data$curve_order), , drop = FALSE],
      ggplot2::aes(embedding_1, embedding_2, group = lineage), colour = "black", linewidth = 0.75,
      arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")) +
    ggplot2::scale_colour_viridis_c(option = "C", na.value = "#D0D0D0") +
    ggplot2::coord_fixed(clip = "off") + embedding_theme(ctx) +
    ggplot2::labs(title = title, colour = "Pseudotime")
  if (!is.null(node_data) && nrow(node_data) && all(c("embedding_1", "embedding_2") %in% names(node_data))) {
    p <- p + ggplot2::geom_point(data = node_data, ggplot2::aes(embedding_1, embedding_2), inherit.aes = FALSE,
      shape = 21, size = 3.2, stroke = 0.6, fill = "#F2F2F2", colour = "black")
  }
  p <- p + axis_arrow_layers(cell_data, "UMAP 1", "UMAP 2")
  list(plot = p, source_data = list(cells = cell_data, curves = curve_data, nodes = node_data %||% data.frame()),
       parameters = list(type = "trajectory_embedding", n_cells = nrow(cell_data), n_curve_points = nrow(curve_data), lineages = unique(curve_data$lineage)))
}

pseudotime_trend_data <- function(df, pseudotime_col, value_col, group_col = NULL, bins = 40L, span = 0.45) {
  req <- c(pseudotime_col, value_col, group_col)
  req <- req[!is.null(req) & nzchar(req)]
  d <- df[complete.cases(df[, req, drop = FALSE]), , drop = FALSE]
  if (!nrow(d)) return(list(binned = data.frame(), smooth = data.frame()))
  d$.group <- if (is.null(group_col)) "All" else as.character(d[[group_col]])
  rows_bin <- list(); rows_smooth <- list(); kb <- 0L; ks <- 0L
  for (g in unique(d$.group)) {
    z <- d[d$.group == g, , drop = FALSE]
    pt <- as.numeric(z[[pseudotime_col]]); val <- as.numeric(z[[value_col]])
    if (length(unique(pt)) < 6L) next
    br <- unique(stats::quantile(pt, probs = seq(0, 1, length.out = bins + 1L), na.rm = TRUE, names = FALSE))
    if (length(br) < 4L) next
    b <- cut(pt, breaks = br, include.lowest = TRUE, labels = FALSE)
    split_idx <- split(seq_along(pt), b)
    bb <- do.call(rbind, lapply(split_idx, function(ii) {
      data.frame(group = g, pseudotime = stats::median(pt[ii], na.rm = TRUE), mean = mean(val[ii], na.rm = TRUE),
                 se = stats::sd(val[ii], na.rm = TRUE) / sqrt(sum(is.finite(val[ii]))), n = sum(is.finite(val[ii])), stringsAsFactors = FALSE)
    }))
    bb <- bb[is.finite(bb$pseudotime) & is.finite(bb$mean), , drop = FALSE]
    if (!nrow(bb)) next
    kb <- kb + 1L; rows_bin[[kb]] <- bb
    grid <- seq(min(bb$pseudotime), max(bb$pseudotime), length.out = 240L)
    fit <- tryCatch(stats::loess(mean ~ pseudotime, data = bb, weights = pmax(bb$n, 1), span = span, degree = 1, control = stats::loess.control(surface = "direct")), error = function(e) NULL)
    if (!is.null(fit)) {
      pr <- tryCatch(stats::predict(fit, newdata = data.frame(pseudotime = grid), se = TRUE), error = function(e) NULL)
      if (!is.null(pr) && is.list(pr)) {
        ss <- data.frame(group = g, pseudotime = grid, mean_smooth = as.numeric(pr$fit),
                         low = as.numeric(pr$fit - 1.96 * pr$se.fit), high = as.numeric(pr$fit + 1.96 * pr$se.fit), stringsAsFactors = FALSE)
      } else {
        yy <- as.numeric(stats::predict(fit, newdata = data.frame(pseudotime = grid)))
        se_interp <- stats::approx(bb$pseudotime, bb$se, xout = grid, rule = 2)$y
        ss <- data.frame(group = g, pseudotime = grid, mean_smooth = yy, low = yy - 1.96 * se_interp, high = yy + 1.96 * se_interp, stringsAsFactors = FALSE)
      }
    } else {
      yy <- stats::approx(bb$pseudotime, bb$mean, xout = grid, rule = 2)$y
      se_interp <- stats::approx(bb$pseudotime, bb$se, xout = grid, rule = 2)$y
      ss <- data.frame(group = g, pseudotime = grid, mean_smooth = yy, low = yy - 1.96 * se_interp, high = yy + 1.96 * se_interp, stringsAsFactors = FALSE)
    }
    ss <- ss[is.finite(ss$mean_smooth), , drop = FALSE]
    ks <- ks + 1L; rows_smooth[[ks]] <- ss
  }
  list(binned = if (length(rows_bin)) do.call(rbind, rows_bin) else data.frame(),
       smooth = if (length(rows_smooth)) do.call(rbind, rows_smooth) else data.frame())
}

plot_pseudotime_trend <- function(ctx, df, pseudotime_col, value_col, group_col = NULL, title = NULL) {
  trend <- pseudotime_trend_data(df, pseudotime_col, value_col, group_col,
                                 bins = figure_cfg(ctx, "pseudotime_bins", 40),
                                 span = figure_cfg(ctx, "pseudotime_span", 0.45))
  if (!nrow(trend$smooth)) stop("Pseudotime trend was not estimable for ", value_col)
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = trend$smooth, ggplot2::aes(pseudotime, ymin = low, ymax = high, fill = group), alpha = 0.14, colour = NA) +
    ggplot2::geom_line(data = trend$smooth, ggplot2::aes(pseudotime, mean_smooth, colour = group), linewidth = 0.85) +
    ggplot2::geom_point(data = trend$binned, ggplot2::aes(pseudotime, mean, colour = group), size = 0.7, alpha = 0.35) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "#888888", linewidth = 0.3) +
    pub_theme(ctx) + ggplot2::labs(title = title %||% paste(value_col, "along pseudotime"), x = "Pseudotime", y = value_col, colour = group_col %||% NULL, fill = group_col %||% NULL)
  list(plot = p, source_data = list(binned = trend$binned, smooth = trend$smooth),
       parameters = list(type = "pseudotime_trend", pseudotime_col = pseudotime_col, value_col = value_col, group_col = group_col,
                         bins = figure_cfg(ctx, "pseudotime_bins", 40), span = figure_cfg(ctx, "pseudotime_span", 0.45)))
}

write_source_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.gz$", path, ignore.case = TRUE) && !requireNamespace("data.table", quietly = TRUE)) {
    con <- gzfile(path, open = "wt"); on.exit(close(con), add = TRUE)
    utils::write.csv(as.data.frame(x), con, row.names = FALSE, na = "NA")
  } else {
    write_csv_safe(as.data.frame(x), path)
  }
}

collect_font_sizes <- function(grob) {
  out <- numeric()
  if (inherits(grob, "text") && !is.null(grob$gp$fontsize)) out <- c(out, as.numeric(grob$gp$fontsize))
  if (!is.null(grob$grobs)) for (x in grob$grobs) out <- c(out, collect_font_sizes(x))
  if (!is.null(grob$children)) for (x in grob$children) out <- c(out, collect_font_sizes(x))
  out[is.finite(out)]
}

run_plot_qa <- function(ctx, plot_obj, paths, figure_dir, source_paths, width_mm, height_mm, dpi) {
  checks <- data.frame(check = character(), status = character(), evidence = character(), stringsAsFactors = FALSE)
  add <- function(id, status, evidence) checks <<- rbind(checks, data.frame(check = id, status = status, evidence = as.character(evidence), stringsAsFactors = FALSE))
  b <- tryCatch(ggplot2::ggplot_build(plot_obj), error = function(e) e)
  if (inherits(b, "error")) add("ggplot_build", "FAIL", conditionMessage(b)) else {
    nrows <- vapply(b$data, nrow, integer(1))
    add("ggplot_build", if (length(nrows) && any(nrows > 0)) "PASS" else "FAIL", paste(nrows, collapse = ","))
    finite_xy <- any(vapply(b$data, function(z) {
      xok <- if ("x" %in% names(z)) any(is.finite(z$x)) else TRUE
      yok <- if ("y" %in% names(z)) any(is.finite(z$y)) else TRUE
      xok && yok
    }, logical(1)))
    add("finite_geometry", if (finite_xy) "PASS" else "FAIL", finite_xy)
  }
  req <- c("pdf", "tiff", "png", "proof_png")
  for (nm in req) {
    p <- paths[[nm]]
    ok <- !is.null(p) && file.exists(p) && file.info(p)$size > 1000
    add(paste0("file_", nm), if (ok) "PASS" else "FAIL", p %||% "NA")
  }
  svg_ok <- !is.null(paths$svg) && !is.na(paths$svg) && file.exists(paths$svg) && file.info(paths$svg)$size > 500
  add("file_svg", if (svg_ok) "PASS" else "WARNING", paths$svg %||% "NA")
  source_ok <- length(source_paths) > 0L && all(file.exists(source_paths))
  add("source_data", if (source_ok) "PASS" else "WARNING", paste(source_paths, collapse = ";"))
  gt <- tryCatch(ggplot2::ggplotGrob(plot_obj), error = function(e) NULL)
  fs <- if (is.null(gt)) numeric() else collect_font_sizes(gt)
  min_fs <- if (length(fs)) min(fs, na.rm = TRUE) else NA_real_
  add("minimum_text_size", if (!is.finite(min_fs)) "WARNING" else if (min_fs >= 7) "PASS" else "WARNING", min_fs)
  aspect <- width_mm / height_mm
  add("canvas_geometry", if (is.finite(aspect) && aspect > 0.35 && aspect < 4) "PASS" else "WARNING", paste(width_mm, height_mm, dpi, sep = "/"))
  overall <- if (any(checks$status == "FAIL")) "FAIL" else if (any(checks$status == "WARNING")) "PASS_WITH_WARNINGS" else "PASS"
  checks$overall <- overall
  write_csv_safe(checks, file.path(figure_dir, "visual_QA.csv"))
  writeLines(overall, file.path(figure_dir, "VISUAL_QA_STATUS.txt"))
  list(overall = overall, checks = checks)
}


.figure_export_registry <- new.env(parent = emptyenv())

make_renderable_plot_proxy <- function(plot_obj) {
  grob <- ggplot2::ggplotGrob(plot_obj)
  proxy <- ggplot2::ggplot() +
    ggplot2::annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void()
  attr(proxy, "singlecell_plot_proxy") <- TRUE
  attr(proxy, "proxy_created_at") <- now_iso()
  proxy
}

existing_figure_is_complete <- function(row) {
  req <- intersect(c("pdf", "tiff", "png", "proof_png", "plot_rds", "parameters", "visual_qa"), names(row))
  length(req) > 0L && all(vapply(req, function(nm) {
    x <- as.character(row[[nm]][1]); nzchar(x) && file.exists(x) && file.info(x)$size > 100
  }, logical(1)))
}
export_publication_figure <- function(ctx, figure, figure_id, collection = "main", width_mm = NULL, height_mm = NULL) {
  if (is.null(figure$plot) || !inherits(figure$plot, c("gg", "ggplot"))) stop("Figure object must contain exactly one ggplot.")
  if (inherits(figure$plot, "patchwork") || inherits(figure$plot, "ggarrange")) stop("Composite/montage plot objects are forbidden; export each plot separately.")
  source_data <- figure$source_data %||% list()
  if (is.data.frame(source_data) || is.matrix(source_data)) source_data <- list(data = source_data)
  nonempty_source <- length(source_data) && any(vapply(source_data, function(x) {
    z <- tryCatch(as.data.frame(x), error = function(e) NULL)
    !is.null(z) && nrow(z) > 0L
  }, logical(1)))
  if (!nonempty_source) stop("Every figure requires non-empty figure-specific Source Data.")
  stem <- sanitize_stem(figure_id)
  mpath <- file.path(ctx$dirs$manifests, "figure_export_manifest.csv")
  if (exists(stem, envir = .figure_export_registry, inherits = FALSE)) {
    stop("Duplicate figure_id requested in the current R process: ", stem)
  }
  assign(stem, TRUE, envir = .figure_export_registry)
  if (file.exists(mpath)) {
    prior_manifest <- read_csv_safe(mpath)
    hit <- which("figure_id" %in% names(prior_manifest) & as.character(prior_manifest$figure_id) == stem)
    if (length(hit)) {
      prior <- prior_manifest[tail(hit, 1L), , drop = FALSE]
      if (existing_figure_is_complete(prior)) {
        prior$status <- "SKIPPED_EXISTING"
        prior$reason <- "complete figure assets already exist from a previous run/resume cycle"
        return(invisible(prior))
      }
      stop("Figure manifest contains an incomplete existing figure_id: ", stem, ". Repair or start a new result directory.")
    }
  }
  out_root <- if (identical(collection, "extended")) ctx$dirs$extended else ctx$dirs$main
  figure_path <- file.path(out_root, stem)
  if (dir.exists(figure_path) && length(list.files(figure_path, all.files = TRUE, no.. = TRUE))) {
    stop("A non-empty figure directory exists without a complete manifest row: ", figure_path)
  }
  figure_dir <- safe_dir(figure_path)
  cfg <- ctx$cfg$figures
  width_mm <- width_mm %||% cfg$single_width_mm %||% cfg$width_mm %||% 183
  height_mm <- height_mm %||% cfg$single_height_mm %||% 135
  dpi <- cfg$dpi %||% 600
  proof_width <- cfg$proof_width_mm %||% 183
  proof_height <- height_mm * proof_width / width_mm
  thumb_width <- cfg$thumbnail_width_mm %||% 90
  thumb_height <- height_mm * thumb_width / width_mm
  paths <- list(
    pdf = file.path(figure_dir, paste0(stem, ".pdf")),
    svg = file.path(figure_dir, paste0(stem, ".svg")),
    tiff = file.path(figure_dir, paste0(stem, "_600dpi.tiff")),
    png = file.path(figure_dir, paste0(stem, "_600dpi.png")),
    proof_png = file.path(figure_dir, paste0(stem, "_proof183mm.png")),
    thumbnail = file.path(figure_dir, paste0(stem, "_thumbnail.png")),
    plot_rds = file.path(figure_dir, paste0(stem, "_plot.rds")),
    parameters = file.path(figure_dir, paste0(stem, "_parameters.yml"))
  )
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  ggplot2::ggsave(paths$pdf, figure$plot, width = width_mm, height = height_mm, units = "mm", device = pdf_device, limitsize = FALSE, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(paths$svg, figure$plot, width = width_mm, height = height_mm, units = "mm", device = svglite::svglite, limitsize = FALSE, bg = "white")
  } else {
    paths$svg <- NA_character_
    append_issue(ctx, "FIGURE", "WARNING", "svglite unavailable; SVG not exported", stem)
  }
  ggplot2::ggsave(paths$tiff, figure$plot, width = width_mm, height = height_mm, units = "mm", device = "tiff", dpi = dpi, compression = "lzw", limitsize = FALSE, bg = "white")
  ggplot2::ggsave(paths$png, figure$plot, width = width_mm, height = height_mm, units = "mm", device = "png", dpi = dpi, limitsize = FALSE, bg = "white")
  ggplot2::ggsave(paths$proof_png, figure$plot, width = proof_width, height = proof_height, units = "mm", device = "png", dpi = cfg$proof_dpi %||% 300, limitsize = FALSE, bg = "white")
  ggplot2::ggsave(paths$thumbnail, figure$plot, width = thumb_width, height = thumb_height, units = "mm", device = "png", dpi = 150, limitsize = FALSE, bg = "white")
  proxy <- tryCatch(make_renderable_plot_proxy(figure$plot), error = function(e) NULL)
  if (is.null(proxy)) {
    append_issue(ctx, "FIGURE", "WARNING", "Could not create lightweight plot proxy; saving original ggplot", stem)
    saveRDS(figure$plot, paths$plot_rds, compress = "xz")
  } else {
    saveRDS(proxy, paths$plot_rds, compress = "xz")
  }
  params <- c(list(figure_id = stem, collection = collection, width_mm = width_mm, height_mm = height_mm,
                   dpi = dpi, proof_width_mm = proof_width, font_family = choose_font(cfg$font_family %||% "Arial"),
                   created_at = now_iso()), figure$parameters %||% list())
  write_yaml(params, paths$parameters)
  source_paths <- character()
  if (length(source_data)) {
    for (nm in names(source_data)) {
      if (is.null(source_data[[nm]])) next
      tab <- tryCatch(as.data.frame(source_data[[nm]]), error = function(e) NULL)
      if (is.null(tab)) next
      # Empty optional tables (for example, hulls below the minimum cell
      # threshold) must not be registered as traceable Source Data files.
      if (nrow(tab) == 0L) next
      ext <- if (nrow(tab) > 250000L) ".csv.gz" else ".csv"
      sp <- file.path(figure_dir, paste0(stem, "_source_", sanitize_stem(nm), ext))
      write_source_table(tab, sp); source_paths <- c(source_paths, sp)
    }
  }
  qa <- run_plot_qa(ctx, figure$plot, paths, figure_dir, source_paths, width_mm, height_mm, dpi)
  manifest <- data.frame(
    figure_id = stem, collection = collection, width_mm = width_mm, height_mm = height_mm,
    proof_width_mm = proof_width, dpi = dpi, font = choose_font(cfg$font_family %||% "Arial"),
    pdf = paths$pdf, svg = paths$svg, tiff = paths$tiff, png = paths$png, proof_png = paths$proof_png,
    thumbnail = paths$thumbnail, plot_rds = paths$plot_rds, parameters = paths$parameters,
    source_data = paste(source_paths, collapse = ";"), visual_qa = file.path(figure_dir, "visual_QA.csv"),
    visual_qa_status = qa$overall, status = "EXPORTED", reason = NA_character_, stringsAsFactors = FALSE)
  old <- if (file.exists(mpath)) read_csv_safe(mpath) else manifest[0, ]
  cols <- union(names(old), names(manifest))
  for (nm in setdiff(cols, names(old))) old[[nm]] <- NA
  for (nm in setdiff(cols, names(manifest))) manifest[[nm]] <- NA
  write_csv_safe(rbind(old[, cols, drop = FALSE], manifest[, cols, drop = FALSE]), mpath)
  invisible(manifest)
}

# Backward-compatible wrapper used by stages 06-10 and project profiles.
export_plot_set <- function(ctx, plot_obj, stem, subdir = "main", width_mm = NULL, height_mm = NULL,
                            source_data = NULL, parameters = NULL) {
  export_publication_figure(ctx,
    list(plot = plot_obj, source_data = source_data %||% list(), parameters = parameters %||% list(type = "legacy_plot_call")),
    figure_id = stem, collection = if (identical(subdir, "extended")) "extended" else "main",
    width_mm = width_mm, height_mm = height_mm)
}

