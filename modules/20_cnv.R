# Gated expression-inferred CNA module. CopyKAT is preferred; infercna is supported.
# inferCNV is legacy-only because upstream support has been discontinued.

cnv_numeric_matrix <- function(x) {
  if (inherits(x, "Matrix")) return(x)
  if (is.matrix(x)) return(x)
  if (is.data.frame(x)) {
    keep <- vapply(x, is.numeric, logical(1))
    return(as.matrix(x[, keep, drop = FALSE]))
  }
  stop("Could not coerce CNA output to a numeric matrix.")
}

run_copykat_one <- function(counts, sample_name, cfg, normal_cells = character(), species = "human") {
  if (!requireNamespace("copykat", quietly = TRUE)) stop("copykat unavailable")
  max_cells <- cfg$maximum_cells_per_sample %||% 12000L
  if (ncol(counts) > max_cells) {
    set.seed(cfg$seed %||% 20260804)
    keep <- sample(colnames(counts), max_cells)
    counts <- counts[, keep, drop = FALSE]
    normal_cells <- intersect(normal_cells, keep)
  }
  genome <- if (tolower(species) == "mouse") "mm10" else (cfg$copykat_genome %||% "hg20")
  args <- list(
    rawmat = as.matrix(counts), id.type = cfg$id_type %||% "S",
    ngene.chr = cfg$minimum_genes_per_chromosome %||% 5,
    win.size = cfg$window_size %||% 25,
    KS.cut = cfg$ks_cut %||% 0.1,
    sam.name = sanitize_stem(sample_name), distance = cfg$distance %||% "euclidean",
    norm.cell.names = if (length(normal_cells)) normal_cells else "",
    genome = genome,
    n.cores = if (.Platform$OS.type == "windows") 1L else as.integer(cfg$workers %||% 1L)
  )
  fit <- do.call(copykat::copykat, args)
  list(fit = fit, counts_cells = colnames(counts), backend = "copykat")
}

parse_copykat <- function(res, sample_name) {
  fit <- res$fit
  pred <- fit$prediction %||% fit$pred %||% NULL
  if (is.null(pred)) pred <- data.frame(cell_id = res$counts_cells, copykat_pred = NA_character_, stringsAsFactors = FALSE)
  pred <- as.data.frame(pred)
  cell_col <- names(pred)[grepl("cell", names(pred), ignore.case = TRUE)][1] %||% names(pred)[1]
  class_col <- names(pred)[grepl("pred|copykat", names(pred), ignore.case = TRUE)][1] %||% names(pred)[min(2, ncol(pred))]
  cell_tab <- data.frame(cell_id = as.character(pred[[cell_col]]), cnv_class = as.character(pred[[class_col]]), sample_id = sample_name, backend = "copykat", stringsAsFactors = FALSE)
  cell_tab$cnv_class <- sub("aneuploid", "aneuploid_like", tolower(cell_tab$cnv_class), fixed = TRUE)
  cell_tab$cnv_class <- sub("diploid", "diploid_like", cell_tab$cnv_class, fixed = TRUE)

  cm <- fit$CNAmat %||% fit$CNA.mat %||% fit$cnv %||% NULL
  if (is.null(cm)) return(list(cells = cell_tab, heatmap = data.frame()))
  cm_df <- as.data.frame(cm, check.names = FALSE)
  cell_cols <- intersect(res$counts_cells, names(cm_df))
  if (!length(cell_cols) && is.matrix(cm)) cell_cols <- intersect(res$counts_cells, colnames(cm))
  if (!length(cell_cols)) return(list(cells = cell_tab, heatmap = data.frame()))
  coord_cols <- setdiff(names(cm_df), cell_cols)
  chr_col <- coord_cols[grepl("chr|chrom", coord_cols, ignore.case = TRUE)][1]
  pos_col <- coord_cols[grepl("abs|pos|start", coord_cols, ignore.case = TRUE)][1]
  chromosome <- if (!is.na(chr_col) && length(chr_col)) as.character(cm_df[[chr_col]]) else "unknown"
  position <- if (!is.na(pos_col) && length(pos_col)) suppressWarnings(as.numeric(cm_df[[pos_col]])) else seq_len(nrow(cm_df))
  mat <- as.matrix(cm_df[, cell_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  max_bins <- 600L
  idx <- unique(round(seq(1, nrow(mat), length.out = min(max_bins, nrow(mat)))))
  mat <- mat[idx, , drop = FALSE]; chromosome <- chromosome[idx]; position <- position[idx]
  max_hm_cells <- 2000L
  if (ncol(mat) > max_hm_cells) {
    set.seed(20260804); keep <- sample(colnames(mat), max_hm_cells); mat <- mat[, keep, drop = FALSE]
  }
  heat <- data.frame(
    sample_id = sample_name,
    genomic_bin = paste0(chromosome, ":", position),
    chromosome = chromosome, position = position,
    bin_order = rep(seq_len(nrow(mat)), times = ncol(mat)),
    cell_id = rep(colnames(mat), each = nrow(mat)),
    cell_order = rep(seq_len(ncol(mat)), each = nrow(mat)),
    cna_value = as.vector(mat), backend = "copykat", stringsAsFactors = FALSE
  )
  list(cells = cell_tab, heatmap = heat)
}

run_infercna_one <- function(data_mat, sample_name, cfg, reference_cells, species = "human") {
  if (!requireNamespace("infercna", quietly = TRUE)) stop("infercna unavailable")
  genome <- cfg$infercna_genome %||% if (tolower(species) == "mouse") "mm10" else "hg38"
  infercna::useGenome(genome)
  ref_split <- split(reference_cells, rep("reference", length(reference_cells)))
  if (length(reference_cells) < (cfg$minimum_reference_cells %||% 50L)) ref_split <- NULL
  mat <- as.matrix(data_mat)
  max_cells <- cfg$maximum_cells_per_sample %||% 12000L
  if (ncol(mat) > max_cells) {
    set.seed(cfg$seed %||% 20260804); keep <- sample(colnames(mat), max_cells); mat <- mat[, keep, drop = FALSE]
    if (!is.null(ref_split)) ref_split <- lapply(ref_split, intersect, y = keep)
  }
  cna <- infercna::infercna(mat, refCells = ref_split, window = cfg$window_size %||% 100,
                            noise = cfg$noise %||% 0.1, isLog = TRUE, verbose = FALSE)
  sig <- infercna::cnaSignal(cna)
  cor <- infercna::cnaCor(cna, bySample = FALSE)
  cls <- tryCatch(infercna::findMalignant(cna), error = function(e) NULL)
  class_vec <- rep("unresolved", ncol(cna)); names(class_vec) <- colnames(cna)
  if (!is.null(cls)) {
    if (is.logical(cls) && length(cls) == ncol(cna)) class_vec <- ifelse(cls, "aneuploid_like", "diploid_like")
    if (is.character(cls)) class_vec[names(cls)] <- cls
  }
  cell_tab <- data.frame(cell_id = colnames(cna), sample_id = sample_name, cnv_class = class_vec[colnames(cna)],
                         cna_signal = as.numeric(sig[colnames(cna)]), cna_correlation = as.numeric(cor[colnames(cna)]), backend = "infercna", stringsAsFactors = FALSE)
  max_bins <- 600L; idx <- unique(round(seq(1, nrow(cna), length.out = min(max_bins, nrow(cna)))))
  c2 <- cna[idx, , drop = FALSE]
  if (ncol(c2) > 2000L) {set.seed(20260804); c2 <- c2[, sample(colnames(c2), 2000L), drop = FALSE]}
  heat <- data.frame(sample_id = sample_name, genomic_bin = rownames(c2), chromosome = NA_character_, position = idx,
                     bin_order = rep(seq_len(nrow(c2)), times = ncol(c2)), cell_id = rep(colnames(c2), each = nrow(c2)),
                     cell_order = rep(seq_len(ncol(c2)), each = nrow(c2)), cna_value = as.vector(c2), backend = "infercna", stringsAsFactors = FALSE)
  list(cells = cell_tab, heatmap = heat, fit = cna)
}

run_cnv_module <- function(ctx, object_path = file.path(ctx$dirs$objects, "06_annotated.rds")) {
  cfg <- ctx$cfg$advanced_modules$cnv
  if (!as_flag(cfg$enabled, FALSE)) return(status_table("NOT_RUN", "disabled"))
  if (!file.exists(object_path)) return(status_table("NOT_EVALUABLE", "annotated Seurat object absent"))
  method <- tolower(cfg$method %||% "copykat")
  if (method == "infercnv") append_issue(ctx, "CNV", "WARNING", "inferCNV is a legacy backend; prefer copykat or infercna")
  if (method == "copykat" && !requireNamespace("copykat", quietly = TRUE)) return(status_table("NOT_RUN", "copykat unavailable"))
  if (method == "infercna" && !requireNamespace("infercna", quietly = TRUE)) return(status_table("NOT_RUN", "infercna unavailable"))
  if (method == "infercnv") return(status_table("NOT_EVALUABLE", "legacy inferCNV execution was intentionally removed; select copykat or infercna"))

  obj <- join_layers_safe(readRDS(object_path), "RNA")
  md <- obj[[]]
  label_col <- cfg$label_column %||% "cell_type_l1"
  if (!all(c("sample_id", label_col) %in% names(md))) return(status_table("NOT_EVALUABLE", "sample_id or configured label column absent"))
  refs <- unlist(cfg$reference_labels %||% list())
  candidates <- unlist(cfg$candidate_labels %||% list())
  if (!length(refs)) return(status_table("NOT_EVALUABLE", "explicit normal reference_labels are required"))
  ref_cells_all <- rownames(md)[as.character(md[[label_col]]) %in% refs]
  if (length(ref_cells_all) < (cfg$minimum_reference_cells %||% 50L)) return(status_table("NOT_EVALUABLE", "too few explicit normal reference cells"))
  if (length(candidates)) {
    keep <- as.character(md[[label_col]]) %in% c(refs, candidates)
    obj <- subset(obj, cells = rownames(md)[keep]); md <- obj[[]]
  }
  counts <- get_layer_safe(obj, "RNA", "counts")
  data_mat <- get_layer_safe(obj, "RNA", "data")
  samples <- unique(as.character(md$sample_id))
  cell_rows <- list(); heat_rows <- list(); fits <- list()
  for (s in samples) {
    cells <- rownames(md)[as.character(md$sample_id) == s]
    ref_cells <- intersect(cells, rownames(md)[as.character(md[[label_col]]) %in% refs])
    if (length(cells) < (cfg$minimum_cells_per_sample %||% 200L) || length(ref_cells) < (cfg$minimum_reference_cells_per_sample %||% 20L)) {
      append_issue(ctx, "CNV", "WARNING", "sample skipped for insufficient cells/reference", s); next
    }
    res <- tryCatch({
      if (method == "copykat") parse_copykat(run_copykat_one(counts[, cells, drop = FALSE], s, cfg, ref_cells, ctx$cfg$input$species), s)
      else run_infercna_one(data_mat[, cells, drop = FALSE], s, cfg, ref_cells, ctx$cfg$input$species)
    }, error = function(e) {append_issue(ctx, "CNV", "WARNING", conditionMessage(e), s); NULL})
    if (is.null(res)) next
    cell_rows[[s]] <- res$cells; heat_rows[[s]] <- res$heatmap
    if (!is.null(res$fit)) fits[[s]] <- res$fit
  }
  if (!length(cell_rows)) return(status_table("NOT_EVALUABLE", "no sample passed the CNA gate"))
  cells <- do.call(rbind, cell_rows); heat <- if (length(heat_rows)) do.call(rbind, heat_rows) else data.frame()
  write_csv_safe(cells, file.path(ctx$dirs$advanced, "cnv_cell_scores.csv"))
  if (nrow(heat)) write_source_table(heat, file.path(ctx$dirs$advanced, "cnv_heatmap_long.csv.gz"))
  saveRDS(fits, file.path(ctx$dirs$advanced, "cnv_backend_objects.rds"))
  write_yaml(list(method = method, reference_labels = refs, candidate_labels = candidates,
                  interpretation = "expression-inferred aneuploidy evidence; not a definitive malignancy call"),
             file.path(ctx$dirs$advanced, "cnv_parameters.yml"))
  status_table("COMPLETED", paste(method, "expression-inferred CNA completed; orthogonal validation required"), nrow(cells))
}
