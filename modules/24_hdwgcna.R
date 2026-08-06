# Gated hdWGCNA workflow for a prespecified cell population.
# Biological samples are retained during metacell construction; modules remain exploratory until preservation/validation.

call_supported <- function(fun, args) {
  f <- if (is.character(fun)) getExportedValue("hdWGCNA", fun) else fun
  fm <- names(formals(f))
  if (!"..." %in% fm) args <- args[names(args) %in% fm]
  do.call(f, args)
}

hdwgcna_power_table <- function(obj, name) {
  candidates <- list(
    tryCatch(hdWGCNA::GetPowerTable(obj, wgcna_name = name), error = function(e) NULL),
    tryCatch(obj@misc[[name]]$powerTable, error = function(e) NULL),
    tryCatch(obj@misc[[name]]$power_table, error = function(e) NULL)
  )
  for (x in candidates) if (!is.null(x) && nrow(as.data.frame(x))) return(as.data.frame(x, stringsAsFactors = FALSE))
  data.frame()
}

select_hdwgcna_power <- function(tab, threshold = 0.8, fallback = 6L) {
  if (!nrow(tab)) return(as.integer(fallback))
  pcol <- names(tab)[tolower(names(tab)) %in% c("power", "softpower", "soft_power")][1]
  rcol <- names(tab)[grepl("sft.*r|signed.*r|scale.*fit", names(tab), ignore.case = TRUE)][1]
  if (is.na(pcol) || is.na(rcol)) return(as.integer(fallback))
  ok <- is.finite(tab[[pcol]]) & is.finite(tab[[rcol]]) & tab[[rcol]] >= threshold
  if (any(ok)) as.integer(min(tab[[pcol]][ok])) else as.integer(tab[[pcol]][which.max(tab[[rcol]])])
}

module_sample_summary <- function(MEs, md, trait_columns = character()) {
  common <- intersect(rownames(MEs), rownames(md)); MEs <- MEs[common, , drop = FALSE]; md <- md[common, , drop = FALSE]
  trait_columns <- unique(c("patient_id", "condition", trait_columns))
  trait_columns <- trait_columns[trait_columns %in% names(md)]
  rows <- list(); z <- 0L
  for (sn in unique(as.character(md$sample_id))) {
    ii <- which(as.character(md$sample_id) == sn)
    traits <- list()
    for (tr in trait_columns) {
      vals <- unique(md[[tr]][ii][!is.na(md[[tr]][ii])])
      traits[[tr]] <- if (length(vals) == 1L) vals[[1]] else NA
    }
    for (m in colnames(MEs)) {
      z <- z + 1L
      row <- data.frame(sample_id = sn, module = m,
                        mean_ME = mean(MEs[ii, m], na.rm = TRUE),
                        median_ME = stats::median(MEs[ii, m], na.rm = TRUE),
                        n_cells = length(ii), stringsAsFactors = FALSE)
      for (tr in names(traits)) row[[tr]] <- traits[[tr]]
      rows[[z]] <- row
    }
  }
  do.call(rbind, rows)
}

module_condition_test <- function(tab, reference, target, minimum_replicates = 3L, paired = FALSE) {
  if (!nrow(tab) || is.null(reference) || is.null(target)) return(data.frame())
  rows <- lapply(split(tab, tab$module), function(z) {
    a <- z[z$condition == reference, , drop = FALSE]; b <- z[z$condition == target, , drop = FALSE]
    if (nrow(a) < minimum_replicates || nrow(b) < minimum_replicates) {
      return(data.frame(module = z$module[1], n_reference = nrow(a), n_target = nrow(b), delta_ME = mean(b$mean_ME) - mean(a$mean_ME), p_value = NA_real_, stringsAsFactors = FALSE))
    }
    if (paired && all(!is.na(z$patient_id))) {
      aa <- a[, c("patient_id", "mean_ME")]; bb <- b[, c("patient_id", "mean_ME")]
      names(aa)[2] <- "a"; names(bb)[2] <- "b"; mm <- merge(aa, bb, by = "patient_id")
      p <- if (nrow(mm) >= minimum_replicates) tryCatch(stats::wilcox.test(mm$b, mm$a, paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_
      delta <- if (nrow(mm)) mean(mm$b - mm$a, na.rm = TRUE) else NA_real_
    } else {
      p <- tryCatch(stats::wilcox.test(b$mean_ME, a$mean_ME, exact = FALSE)$p.value, error = function(e) NA_real_)
      delta <- mean(b$mean_ME, na.rm = TRUE) - mean(a$mean_ME, na.rm = TRUE)
    }
    data.frame(module = z$module[1], n_reference = nrow(a), n_target = nrow(b), delta_ME = delta, p_value = p, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); out$FDR <- stats::p.adjust(out$p_value, method = "BH"); out
}

module_trait_correlations <- function(sample_tab, traits) {
  traits <- traits[traits %in% names(sample_tab)]
  if (!length(traits)) return(data.frame())
  wide <- reshape(sample_tab[, c("sample_id", "module", "mean_ME")], idvar = "sample_id", timevar = "module", direction = "wide")
  md <- unique(sample_tab[, c("sample_id", traits), drop = FALSE]); x <- merge(wide, md, by = "sample_id")
  module_cols <- grep("^mean_ME\\.", names(x), value = TRUE); rows <- list(); z <- 0L
  for (m in module_cols) for (tr in traits) {
    v <- x[[tr]]
    if (!is.numeric(v)) {
      lev <- unique(stats::na.omit(as.character(v)))
      if (length(lev) != 2L) next
      v <- as.numeric(factor(v, levels = lev)) - 1
    }
    ok <- is.finite(x[[m]]) & is.finite(v)
    if (sum(ok) < 4L || stats::sd(v[ok]) == 0) next
    ct <- tryCatch(stats::cor.test(x[[m]][ok], v[ok], method = "spearman", exact = FALSE), error = function(e) NULL)
    if (is.null(ct)) next
    z <- z + 1L; rows[[z]] <- data.frame(module = sub("^mean_ME\\.", "", m), trait = tr,
                                          correlation = unname(ct$estimate), p_value = ct$p.value, n_samples = sum(ok), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows); out$FDR <- stats::p.adjust(out$p_value, method = "BH"); out
}

run_hdwgcna_module <- function(ctx, object_path = file.path(ctx$dirs$objects, "06_annotated.rds")) {
  cfg <- ctx$cfg$advanced_modules$hdWGCNA
  if (!as_flag(cfg$enabled, FALSE)) return(status_table("NOT_RUN", "disabled"))
  if (!requireNamespace("hdWGCNA", quietly = TRUE) || !requireNamespace("WGCNA", quietly = TRUE)) return(status_table("NOT_RUN", "hdWGCNA/WGCNA unavailable"))
  if (!file.exists(object_path)) return(status_table("NOT_EVALUABLE", "annotated object absent"))
  obj <- join_layers_safe(readRDS(object_path), "RNA"); md <- obj[[]]
  group_col <- cfg$grouping_column %||% "cell_type_l1"; group_values <- unlist(cfg$group_values %||% list())
  if (!all(c("sample_id", group_col) %in% names(md)) || !length(group_values)) return(status_table("NOT_EVALUABLE", "explicit grouping_column and group_values are required"))
  cells <- rownames(md)[as.character(md[[group_col]]) %in% group_values]
  sub <- subset(obj, cells = cells); md <- sub[[]]
  valid <- names(table(md$sample_id))[table(md$sample_id) >= (cfg$minimum_cells_per_sample %||% 100L)]
  if (length(valid) < (cfg$minimum_samples %||% 3L)) return(status_table("NOT_EVALUABLE", "too few biological samples pass the hdWGCNA cell gate"))
  sub <- subset(sub, cells = rownames(md)[md$sample_id %in% valid]); md <- sub[[]]
  reduction <- cfg$reduction %||% if ("harmony" %in% names(sub@reductions)) "harmony" else "pca"
  if (!reduction %in% names(sub@reductions)) return(status_table("NOT_EVALUABLE", paste("missing reduction", reduction)))
  wname <- cfg$wgcna_name %||% paste0("hdWGCNA_", sanitize_stem(paste(group_values, collapse = "_")))

  sub <- call_supported("SetupForWGCNA", list(seurat_obj = sub, gene_select = cfg$gene_select %||% "fraction",
                                               fraction = cfg$gene_fraction %||% 0.05, wgcna_name = wname))
  group_by <- unique(c(group_col, "sample_id"))
  sub <- call_supported("MetacellsByGroups", list(seurat_obj = sub, group.by = group_by, ident.group = group_col,
                                                   k = cfg$metacell_k %||% 25L, reduction = reduction,
                                                   dims = seq_len(cfg$dimensions %||% 30L), assay = "RNA", layer = "counts",
                                                   mode = cfg$metacell_mode %||% "average", min_cells = cfg$minimum_cells_per_group %||% 100L,
                                                   max_shared = cfg$maximum_shared_cells %||% 10L, target_metacells = cfg$target_metacells %||% 1000L,
                                                   wgcna_name = wname, verbose = FALSE))
  sub <- call_supported("NormalizeMetacells", list(seurat_obj = sub, wgcna_name = wname))
  sub <- call_supported("SetDatExpr", list(seurat_obj = sub, group_name = group_values, group.by = group_col,
                                            assay = "RNA", layer = "data", wgcna_name = wname))
  powers <- as.integer(unlist(cfg$soft_powers %||% c(1:10, seq(12, 30, 2))))
  sub <- call_supported("TestSoftPowers", list(seurat_obj = sub, powers = powers, networkType = cfg$network_type %||% "signed",
                                                corFnc = cfg$correlation_function %||% "bicor", wgcna_name = wname))
  ptab <- hdwgcna_power_table(sub, wname)
  if (nrow(ptab)) write_csv_safe(ptab, file.path(ctx$dirs$advanced, "hdWGCNA_soft_power_table.csv"))
  soft_power <- as.integer(cfg$soft_power %||% select_hdwgcna_power(ptab, cfg$scale_free_threshold %||% 0.8, 6L))
  tom_dir <- safe_dir(file.path(ctx$dirs$advanced, "hdWGCNA_TOM"))
  sub <- call_supported("ConstructNetwork", list(seurat_obj = sub, soft_power = soft_power, tom_outdir = tom_dir,
                                                  tom_name = wname, networkType = cfg$network_type %||% "signed",
                                                  minModuleSize = cfg$minimum_module_size %||% 30L,
                                                  mergeCutHeight = cfg$merge_cut_height %||% 0.25,
                                                  deepSplit = cfg$deep_split %||% 2L,
                                                  store_tom_in_seurat = as_flag(cfg$store_tom_in_object, FALSE), wgcna_name = wname))
  harmonize <- unlist(cfg$harmonize_module_eigengenes_by %||% list())
  harmonize <- harmonize[harmonize %in% names(sub[[]])]
  sub <- call_supported("ModuleEigengenes", list(seurat_obj = sub, group.by.vars = if (length(harmonize)) harmonize else NULL,
                                                  assay = "RNA", wgcna_name = wname, verbose = FALSE))
  sub <- call_supported("ModuleConnectivity", list(seurat_obj = sub, group.by = group_col, group_name = group_values,
                                                    assay = "RNA", layer = "data", wgcna_name = wname))
  score_method <- cfg$module_score_method %||% if (requireNamespace("UCell", quietly = TRUE)) "UCell" else "Seurat"
  sub <- tryCatch(call_supported("ModuleExprScore", list(seurat_obj = sub, n_genes = cfg$module_score_genes %||% 25L,
                                                          method = score_method, wgcna_name = wname)), error = function(e) sub)
  saveRDS(sub, file.path(ctx$dirs$advanced, "hdWGCNA_object.rds"), compress = FALSE)

  modules <- as.data.frame(hdWGCNA::GetModules(sub, wgcna_name = wname), stringsAsFactors = FALSE)
  modules <- modules[as.character(modules$module) != "grey", , drop = FALSE]
  hubs <- as.data.frame(hdWGCNA::GetHubGenes(sub, n_hubs = cfg$hub_genes_per_module %||% 20L, wgcna_name = wname), stringsAsFactors = FALSE)
  MEs <- as.matrix(hdWGCNA::GetMEs(sub, harmonized = length(harmonize) > 0L, wgcna_name = wname))
  write_csv_safe(modules, file.path(ctx$dirs$advanced, "hdWGCNA_modules.csv"))
  write_csv_safe(hubs, file.path(ctx$dirs$advanced, "hdWGCNA_hub_genes.csv"))
  me_cells <- data.frame(cell_id = rownames(MEs), MEs, check.names = FALSE, stringsAsFactors = FALSE)
  write_source_table(me_cells, file.path(ctx$dirs$advanced, "hdWGCNA_cell_module_eigengenes.csv.gz"))
  trait_cfg <- unlist(cfg$trait_columns %||% c("condition"))
  sample_tab <- module_sample_summary(MEs, sub[[]], trait_cfg)
  write_csv_safe(sample_tab, file.path(ctx$dirs$advanced, "hdWGCNA_sample_module_summary.csv"))
  dme <- module_condition_test(sample_tab, ctx$cfg$comparisons$reference %||% NULL, ctx$cfg$comparisons$target %||% NULL,
                               ctx$cfg$comparisons$minimum_replicates_per_condition %||% 3L, as_flag(ctx$cfg$metadata$paired, FALSE))
  if (nrow(dme)) write_csv_safe(dme, file.path(ctx$dirs$advanced, "hdWGCNA_differential_module_eigengenes.csv"))
  traits <- trait_cfg
  traits <- traits[traits %in% names(sub[[]])]
  trait_tab <- module_trait_correlations(sample_tab, traits)
  if (nrow(trait_tab)) write_csv_safe(trait_tab, file.path(ctx$dirs$advanced, "hdWGCNA_module_trait_correlations.csv"))
  dendro <- tryCatch(hdWGCNA::PlotDendrogram(sub, wgcna_name = wname, main = "hdWGCNA module dendrogram"), error = function(e) NULL)
  if (!is.null(dendro)) saveRDS(dendro, file.path(ctx$dirs$advanced, "hdWGCNA_dendrogram_plot.rds"))
  write_yaml(list(wgcna_name = wname, group_column = group_col, group_values = group_values, samples = valid,
                  soft_power = soft_power, network_type = cfg$network_type %||% "signed", score_method = score_method,
                  interpretation = "exploratory co-expression modules; require preservation and orthogonal validation"),
             file.path(ctx$dirs$advanced, "hdWGCNA_parameters.yml"))
  status_table("COMPLETED", paste("hdWGCNA identified", length(unique(modules$module)), "non-grey modules"), length(unique(modules$module)))
}
