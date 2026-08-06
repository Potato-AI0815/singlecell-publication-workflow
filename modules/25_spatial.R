# Spatial transcriptomics workflow: Seurat preprocessing, spatial-variable features, reference transfer,
# optional RCTD deconvolution, and coordinate-neighborhood analysis. Every slice is retained as a distinct sample.

load_spatial_object <- function(cfg) {
  if (!is.null(cfg$object_path) && file.exists(cfg$object_path)) return(readRDS(cfg$object_path))
  if (!is.null(cfg$data_dir) && dir.exists(cfg$data_dir)) {
    return(Seurat::Load10X_Spatial(data.dir = cfg$data_dir, filename = cfg$matrix_filename %||% "filtered_feature_bc_matrix.h5",
                                   assay = cfg$assay %||% "Spatial", slice = cfg$slice_name %||% "slice1", filter.matrix = TRUE))
  }
  stop("spatial.object_path or spatial.data_dir is required")
}

spatial_coordinates <- function(obj, image_name = NULL) {
  imgs <- Seurat::Images(obj)
  if (!length(imgs)) stop("Seurat spatial object has no image/FOV")
  selected <- if (!is.null(image_name)) intersect(image_name, imgs) else imgs
  if (!length(selected)) stop("Requested spatial image/FOV was not found")
  rows <- lapply(selected, function(im) {
    co <- as.data.frame(Seurat::GetTissueCoordinates(obj, image = im), stringsAsFactors = FALSE)
    if (!"cell" %in% names(co)) co$cell <- rownames(co)
    xcol <- names(co)[tolower(names(co)) %in% c("x", "imagecol", "col", "pxl_col_in_fullres")][1]
    ycol <- names(co)[tolower(names(co)) %in% c("y", "imagerow", "row", "pxl_row_in_fullres")][1]
    if (is.na(xcol) || is.na(ycol)) {
      nums <- names(co)[vapply(co, is.numeric, logical(1))]
      if (length(nums) < 2L) stop("could not identify spatial coordinate columns for image ", im)
      xcol <- nums[1]; ycol <- nums[2]
    }
    data.frame(cell_id=as.character(co$cell), x=as.numeric(co[[xcol]]), y=as.numeric(co[[ycol]]), image=im, stringsAsFactors=FALSE)
  })
  out <- do.call(rbind, rows)
  out[!duplicated(out$cell_id), , drop=FALSE]
}

run_rctd_contract <- function(spatial, reference, coords, label_col, spatial_assay, cfg) {
  if (!as_flag(cfg$run_rctd, TRUE)) return(list(status="NOT_RUN"))
  if (!requireNamespace("spacexr", quietly=TRUE)) return(list(status="NOT_RUN", reason="spacexr unavailable"))
  reference <- join_layers_safe(reference,"RNA"); spatial <- join_layers_safe(spatial,spatial_assay)
  rc <- get_layer_safe(reference,"RNA","counts")
  labs <- factor(as.character(reference[[label_col]][,1])); names(labs) <- colnames(reference)
  good <- !is.na(labs) & !(labs %in% unlist(cfg$exclude_reference_labels %||% c("Unresolved","Ambiguous","Doublet-like")))
  rc <- rc[,good,drop=FALSE]; labs <- droplevels(labs[good])
  min_ref <- cfg$rctd_minimum_reference_cells %||% 20L
  supported <- names(table(labs))[table(labs)>=min_ref]
  keep_ref <- labs %in% supported; rc <- rc[,keep_ref,drop=FALSE]; labs <- droplevels(labs[keep_ref])
  if(length(levels(labs))<2L) return(list(status="NOT_EVALUABLE",reason="fewer than two supported reference cell types"))
  ref_obj <- spacexr::Reference(as.matrix(rc),labs,Matrix::colSums(rc))
  sc <- get_layer_safe(spatial,spatial_assay,"counts")
  common <- intersect(colnames(sc),coords$cell_id); coords <- coords[match(common,coords$cell_id),,drop=FALSE]
  sample_ids <- if("spatial_sample_id" %in% names(coords)) as.character(coords$spatial_sample_id) else rep("slice1",length(common))
  objects <- list(); weights_list <- list(); statuses <- list()
  for(sn in unique(sample_ids)) {
    cells <- common[sample_ids==sn]
    if(length(cells)<(cfg$minimum_locations %||% 100L)) {statuses[[sn]]<-data.frame(spatial_sample_id=sn,status="NOT_EVALUABLE",reason="too few locations",stringsAsFactors=FALSE);next}
    cc <- coords[match(cells,coords$cell_id),c("x","y"),drop=FALSE]; rownames(cc)<-cells
    spatial_obj <- spacexr::SpatialRNA(cc,as.matrix(sc[,cells,drop=FALSE]),Matrix::colSums(sc[,cells,drop=FALSE]))
    fit <- tryCatch({r<-spacexr::create.RCTD(spatial_obj,ref_obj,max_cores=cfg$rctd_cores %||% 1L,CELL_MIN_INSTANCE=min_ref); spacexr::run.RCTD(r,doublet_mode=cfg$rctd_doublet_mode %||% "full")},error=function(e)e)
    if(inherits(fit,"error")){statuses[[sn]]<-data.frame(spatial_sample_id=sn,status="FAILED",reason=conditionMessage(fit),stringsAsFactors=FALSE);next}
    w <- fit@results$weights %||% fit@results$weights_doublet %||% NULL
    if(is.null(w)){statuses[[sn]]<-data.frame(spatial_sample_id=sn,status="FAILED",reason="no weights",stringsAsFactors=FALSE);next}
    w<-as.matrix(w); objects[[sn]]<-fit; weights_list[[sn]]<-w
    statuses[[sn]]<-data.frame(spatial_sample_id=sn,status="COMPLETED",reason=paste(nrow(w),"locations"),stringsAsFactors=FALSE)
  }
  status_df<-do.call(rbind,statuses)
  if(!length(weights_list)) return(list(status="NOT_EVALUABLE",reason="no spatial sample passed RCTD",sample_status=status_df))
  all_types<-sort(unique(unlist(lapply(weights_list,colnames))))
  mats<-lapply(weights_list,function(w){out<-matrix(0,nrow(w),length(all_types),dimnames=list(rownames(w),all_types));out[,colnames(w)]<-w;out})
  list(status="COMPLETED",object=objects,weights=do.call(rbind,mats),sample_status=status_df)
}

spatial_knn_enrichment <- function(coords, labels, k = 6L, permutations = 200L, seed = 20260804) {
  ok <- complete.cases(coords) & !is.na(labels); coords <- as.matrix(coords[ok, , drop = FALSE]); labels <- as.character(labels[ok])
  if (nrow(coords) < 20L || length(unique(labels)) < 2L) return(data.frame())
  if (requireNamespace("FNN", quietly = TRUE)) nn <- FNN::get.knn(coords, k = min(k, nrow(coords) - 1L))$nn.index
  else {
    d <- as.matrix(stats::dist(coords)); diag(d) <- Inf; nn <- t(apply(d, 1, order)[seq_len(min(k, nrow(coords) - 1L)), , drop = FALSE])
  }
  lev <- sort(unique(labels)); count_edges <- function(lab) {
    tab <- matrix(0, length(lev), length(lev), dimnames = list(lev, lev))
    for (i in seq_len(nrow(nn))) for (j in nn[i, ]) tab[lab[i], lab[j]] <- tab[lab[i], lab[j]] + 1
    tab
  }
  obs <- count_edges(labels); set.seed(seed)
  null_sum <- null_sq <- matrix(0, nrow(obs), ncol(obs), dimnames = dimnames(obs))
  for (b in seq_len(permutations)) {x <- count_edges(sample(labels)); null_sum <- null_sum + x; null_sq <- null_sq + x^2}
  mu <- null_sum / permutations; sdv <- sqrt(pmax(null_sq / permutations - mu^2, 0)); z <- (obs - mu) / ifelse(sdv > 0, sdv, NA_real_)
  data.frame(source = rep(rownames(obs), times = ncol(obs)), target = rep(colnames(obs), each = nrow(obs)), observed = as.vector(obs), expected = as.vector(mu), z_score = as.vector(z), stringsAsFactors = FALSE)
}

run_spatial_module <- function(ctx, object_path = file.path(ctx$dirs$objects, "06_annotated.rds")) {
  cfg <- ctx$cfg$advanced_modules$spatial
  if (!as_flag(cfg$enabled, FALSE)) return(status_table("NOT_RUN", "disabled"))
  spatial <- tryCatch(load_spatial_object(cfg), error = function(e) e)
  if (inherits(spatial, "error")) return(status_table("NOT_EVALUABLE", conditionMessage(spatial)))
  assay <- cfg$assay %||% if ("Spatial" %in% names(spatial@assays)) "Spatial" else Seurat::DefaultAssay(spatial)
  Seurat::DefaultAssay(spatial) <- assay
  if (ncol(spatial) < (cfg$minimum_locations %||% 100L)) return(status_table("NOT_EVALUABLE", "too few spatial locations"))
  coords <- spatial_coordinates(spatial, cfg$image_name %||% NULL)
  md <- spatial[[]]
  spatial$spatial_sample_id <- md[[cfg$sample_column %||% "orig.ident"]] %||% (cfg$slice_name %||% "slice1")
  coords$spatial_sample_id <- as.character(spatial$spatial_sample_id[match(coords$cell_id, colnames(spatial))])
  coords$spatial_sample_id[is.na(coords$spatial_sample_id)] <- coords$image[is.na(coords$spatial_sample_id)]

  norm <- tolower(cfg$normalization %||% "sct")
  if (norm == "sct") spatial <- Seurat::SCTransform(spatial, assay = assay, new.assay.name = "SCT", verbose = FALSE, ncells = min(cfg$sct_model_cells %||% 3000L, ncol(spatial)))
  else spatial <- Seurat::NormalizeData(spatial, assay = assay, verbose = FALSE) |> Seurat::FindVariableFeatures(assay = assay, nfeatures = cfg$variable_features %||% 3000L, verbose = FALSE)
  active <- if (norm == "sct" && "SCT" %in% names(spatial@assays)) "SCT" else assay
  Seurat::DefaultAssay(spatial) <- active
  if (!"pca" %in% names(spatial@reductions)) spatial <- Seurat::RunPCA(spatial, assay = active, npcs = cfg$dimensions %||% 30L, verbose = FALSE)
  spatial <- Seurat::FindNeighbors(spatial, reduction = "pca", dims = seq_len(cfg$dimensions %||% 30L), verbose = FALSE)
  spatial <- Seurat::FindClusters(spatial, resolution = cfg$clustering_resolution %||% 0.4, verbose = FALSE)
  if (!"umap" %in% names(spatial@reductions)) spatial <- Seurat::RunUMAP(spatial, reduction = "pca", dims = seq_len(cfg$dimensions %||% 30L), seed.use = ctx$cfg$runtime$seed %||% 20260804, verbose = FALSE)
  spatial$spatial_cluster <- as.character(Seurat::Idents(spatial))

  spatial_features <- character()
  if (as_flag(cfg$find_spatially_variable_features, TRUE)) {
    feats <- head(Seurat::VariableFeatures(spatial), cfg$maximum_spatial_test_features %||% 2000L)
    spatial <- tryCatch(Seurat::FindSpatiallyVariableFeatures(spatial, assay = active, features = feats,
                                                               selection.method = cfg$spatial_variable_method %||% "moransi", verbose = FALSE),
                        error = function(e) {append_issue(ctx, "SPATIAL", "WARNING", conditionMessage(e)); spatial})
    spatial_features <- tryCatch(Seurat::SpatiallyVariableFeatures(spatial, method = cfg$spatial_variable_method %||% "moransi"), error = function(e) character())
    if (length(spatial_features)) write_csv_safe(data.frame(rank = seq_along(spatial_features), gene = spatial_features), file.path(ctx$dirs$advanced, "spatial_variable_features.csv"))
  }

  reference <- if (file.exists(cfg$reference_object_path %||% "")) readRDS(cfg$reference_object_path) else if (file.exists(object_path)) readRDS(object_path) else NULL
  label_col <- cfg$reference_label_column %||% "cell_type_l1"
  transfer_status <- "NOT_RUN"
  if (!is.null(reference) && label_col %in% names(reference[[]]) && as_flag(cfg$run_label_transfer, TRUE)) {
    reference <- join_layers_safe(reference, "RNA")
    ref_assay <- if ("SCT" %in% names(reference@assays)) "SCT" else "RNA"
    if (!"pca" %in% names(reference@reductions)) reference <- Seurat::RunPCA(reference, assay = ref_assay, npcs = cfg$dimensions %||% 30L, verbose = FALSE)
    anchors <- tryCatch(Seurat::FindTransferAnchors(reference = reference, query = spatial,
                                                     normalization.method = if (ref_assay == "SCT" && active == "SCT") "SCT" else "LogNormalize",
                                                     reference.assay = ref_assay, query.assay = active, dims = seq_len(cfg$dimensions %||% 30L)), error = function(e) e)
    if (!inherits(anchors, "error")) {
      pred <- Seurat::TransferData(anchorset = anchors, refdata = reference[[label_col]][, 1], prediction.assay = TRUE,
                                  weight.reduction = spatial[["pca"]], dims = seq_len(cfg$dimensions %||% 30L))
      spatial[["predictions"]] <- pred
      score_mat <- SeuratObject::LayerData(spatial, assay = "predictions", layer = "data")
      max_i <- max.col(t(as.matrix(score_mat)), ties.method = "first")
      spatial$predicted_cell_type <- rownames(score_mat)[max_i]
      spatial$prediction_score_max <- apply(as.matrix(score_mat), 2, max, na.rm = TRUE)
      transfer_status <- "COMPLETED"
      ps <- data.frame(cell_id = colnames(score_mat), predicted_cell_type = spatial$predicted_cell_type, prediction_score_max = spatial$prediction_score_max, stringsAsFactors = FALSE)
      write_csv_safe(ps, file.path(ctx$dirs$advanced, "spatial_label_transfer_predictions.csv"))
      score_long <- data.frame(cell_type = rep(rownames(score_mat), times = ncol(score_mat)), cell_id = rep(colnames(score_mat), each = nrow(score_mat)), score = as.vector(score_mat), stringsAsFactors = FALSE)
      write_source_table(score_long, file.path(ctx$dirs$advanced, "spatial_label_transfer_scores.csv.gz"))
    } else append_issue(ctx, "SPATIAL", "WARNING", paste("label transfer failed:", conditionMessage(anchors)))
  }

  rctd_status <- "NOT_RUN"
  if (!is.null(reference) && label_col %in% names(reference[[]])) {
    rr <- tryCatch(run_rctd_contract(spatial, reference, coords, label_col, assay, cfg), error = function(e) e)
    if (!inherits(rr, "error")) {
      rctd_status <- rr$status
      if (identical(rr$status, "COMPLETED")) {
        saveRDS(rr$object, file.path(ctx$dirs$advanced, "spatial_RCTD_object.rds"), compress = FALSE)
        if (!is.null(rr$sample_status)) write_csv_safe(rr$sample_status, file.path(ctx$dirs$advanced, "spatial_RCTD_sample_status.csv"))
        w <- rr$weights
        wt <- data.frame(cell_id = rep(rownames(w), times = ncol(w)), cell_type = rep(colnames(w), each = nrow(w)), weight = as.vector(w), stringsAsFactors = FALSE)
        write_source_table(wt, file.path(ctx$dirs$advanced, "spatial_RCTD_weights.csv.gz"))
      }
    } else append_issue(ctx, "SPATIAL", "WARNING", paste("RCTD failed:", conditionMessage(rr)))
  }

  md_export <- data.frame(cell_id=rownames(spatial[[]]),spatial[[]],check.names=FALSE)
  md_export <- md_export[, c("cell_id", setdiff(names(md_export), names(coords))), drop=FALSE]
  coords <- merge(coords, md_export, by="cell_id", all.x=TRUE, sort=FALSE)
  write_source_table(coords, file.path(ctx$dirs$advanced, "spatial_coordinates_metadata.csv.gz"))
  label_for_neighbors <- if ("predicted_cell_type" %in% names(coords)) coords$predicted_cell_type else coords$spatial_cluster
  enrich_rows <- lapply(split(seq_len(nrow(coords)), coords$spatial_sample_id), function(ii) {
    z <- spatial_knn_enrichment(coords[ii, c("x","y"), drop=FALSE], label_for_neighbors[ii], cfg$neighborhood_k %||% 6L,
                                cfg$neighborhood_permutations %||% 200L, ctx$cfg$runtime$seed %||% 20260804)
    if(nrow(z)) z$spatial_sample_id <- coords$spatial_sample_id[ii[1]]
    z
  })
  enrich_rows <- Filter(function(x) nrow(x)>0L,enrich_rows)
  enrich <- if(length(enrich_rows)) do.call(rbind,enrich_rows) else data.frame()
  if (nrow(enrich)) write_csv_safe(enrich, file.path(ctx$dirs$advanced, "spatial_neighborhood_enrichment.csv"))
  saveRDS(spatial, file.path(ctx$dirs$advanced, "spatial_processed_object.rds"), compress = FALSE)
  write_yaml(list(assay = assay, active_assay = active, normalization = norm, label_transfer = transfer_status, RCTD = rctd_status,
                  coordinate_image = unique(coords$image), spatial_variable_features = length(spatial_features),
                  interpretation = "spatial patterns and deconvolution are platform/model dependent; preserve slice-level replication"),
             file.path(ctx$dirs$advanced, "spatial_parameters.yml"))
  status_table("COMPLETED", paste("spatial workflow completed; transfer:", transfer_status, "RCTD:", rctd_status), ncol(spatial))
}
