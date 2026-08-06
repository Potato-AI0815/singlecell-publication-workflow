options(stringsAsFactors = FALSE, width = 180)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x

skill_script_path <- function() {
  x <- grep("^--file=", commandArgs(), value = TRUE)
  if (!length(x)) return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  normalizePath(sub("^--file=", "", x[[1]]), winslash = "/", mustWork = TRUE)
}

now_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
now_iso <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

safe_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

read_yaml_strict <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) stop("Package 'yaml' is required.")
  if (!file.exists(path)) stop("Config not found: ", path)
  yaml::read_yaml(path)
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(default)
  isTRUE(x) || identical(tolower(as.character(x)), "true")
}

write_yaml <- function(x, path) yaml::write_yaml(x, path)
write_csv_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fwrite(as.data.frame(x), path, na = "NA") else utils::write.csv(x, path, row.names = FALSE, na = "NA")
}
read_csv_safe <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) data.table::fread(path, data.table = FALSE) else utils::read.csv(path, check.names = FALSE)
}

append_issue <- function(ctx, stage, severity, message, evidence = NA_character_) {
  row <- data.frame(time = now_iso(), stage = stage, severity = severity, message = as.character(message), evidence = as.character(evidence), stringsAsFactors = FALSE)
  path <- file.path(ctx$dirs$audit, "issue_ledger.csv")
  old <- if (file.exists(path)) read_csv_safe(path) else row[0,]
  write_csv_safe(rbind(old, row), path)
  invisible(row)
}

log_msg <- function(...) cat(now_iso(), "|", paste(..., collapse = " "), "\n")

required_packages <- function() c("yaml","digest","data.table","Matrix","ggplot2","Seurat","SeuratObject","scales")
optional_packages <- function() c(
  "SingleCellExperiment","scDblFinder","SingleR","celldex","edgeR","limma","fgsea","msigdbr",
  "zellkonverter","sctransform","glmGamPoi","svglite","ragg","ggrepel","ggrastr","ggforce",
  "concaveman","dbscan","png","clustree","CellChat","slingshot","tradeSeq","BiocParallel",
  "GeneNMF","RcppML","NMF","hdWGCNA","WGCNA","UCell","copykat","infercna","spacexr","FNN",
  "oncoPredict","ComplexHeatmap","circlize","scTenifoldKnk","scTenifoldNet","igraph"
)

package_table <- function() {
  pkgs <- unique(c(required_packages(), optional_packages()))
  data.frame(package = pkgs, available = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE), version = vapply(pkgs, function(p) if (requireNamespace(p, quietly=TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1)), stringsAsFactors = FALSE)
}

sha256_file <- function(path, max_gb = Inf) {
  size <- file.info(path)$size
  if (is.na(size) || size > max_gb * 1024^3) return(NA_character_)
  if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

inventory_path <- function(path, hash = TRUE, max_gb = 20) {
  paths <- if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE) else path
  paths <- paths[file.exists(paths) & !dir.exists(paths)]
  if (!length(paths)) return(data.frame())
  info <- file.info(paths)
  data.frame(path = normalizePath(paths, winslash = "/", mustWork = TRUE), bytes = info$size, modified = format(info$mtime, "%Y-%m-%d %H:%M:%S %Z"), sha256 = if (hash) vapply(paths, sha256_file, character(1), max_gb=max_gb) else NA_character_, stringsAsFactors = FALSE)
}

init_context <- function(config_path, skill_root, mode_override = NULL) {
  cfg <- read_yaml_strict(config_path)
  if (!is.null(mode_override) && nzchar(mode_override)) cfg$analysis$mode <- mode_override
  project_root <- normalizePath(cfg$project_root, winslash = "/", mustWork = FALSE)
  dir.create(project_root, recursive = TRUE, showWarnings = FALSE)
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  mode <- cfg$analysis$mode %||% "full"
  output_base <- cfg$output$root %||% "auto"
  if (identical(output_base, "auto")) output_base <- file.path(project_root, "single_cell_results")
  output_base <- safe_dir(output_base)
  latest_pointer <- file.path(output_base, "LATEST_RUN.txt")
  existing_root <- cfg$output$existing_result_root %||% NULL
  if (identical(mode, "resume") && file.exists(latest_pointer)) {
    result_root <- trimws(readLines(latest_pointer, warn = FALSE)[1])
    if (!dir.exists(result_root)) stop("LATEST_RUN points to missing directory: ", result_root)
    result_root <- normalizePath(result_root, winslash = "/", mustWork = TRUE)
  } else {
    run_id <- cfg$run_id %||% "auto"
    if (identical(run_id, "auto")) run_id <- now_stamp()
    result_root <- safe_dir(file.path(output_base, paste0(cfg$project_name %||% "scRNA", "_", run_id)))
    writeLines(result_root, latest_pointer)
  }
  keys <- c("input","objects","qc","cluster","annotation","markers","composition","de","pathways","main","extended","source","text","logs","qa","audit","manifests","plot_rds","advanced")
  # The issue ledger is an audit artifact stored with the run logs; retain the
  # documented numbered output contract while exposing a semantic audit alias.
  vals <- c("00_input_audit","01_objects","02_qc","03_reduction_clustering","04_annotation","05_markers","06_composition","07_differential_expression","08_pathways","09_main_figures","10_extended_figures","11_source_data","12_methods_legends","13_logs","14_qa","13_logs","15_manifests","16_plot_rds","17_advanced_modules")
  dirs <- setNames(lapply(file.path(result_root, vals), safe_dir), keys)
  target_root <- if(!is.null(existing_root) && dir.exists(existing_root)) normalizePath(existing_root,winslash="/",mustWork=TRUE) else NULL
  if(mode%in%c("audit","figures_only") && is.null(target_root)) stop("output.existing_result_root is required for mode ",mode)
  ctx <- list(cfg=cfg, config_path=normalizePath(config_path,winslash="/",mustWork=TRUE), skill_root=normalizePath(skill_root,winslash="/",mustWork=TRUE), project_root=project_root, result_root=result_root, target_result_root=target_root, dirs=dirs)
  write_yaml(cfg, file.path(dirs$input, "config_resolved.yml"))
  manifest <- list(project_name=cfg$project_name, mode=mode, start_time=now_iso(), config_path=ctx$config_path, skill_root=ctx$skill_root, result_root=result_root)
  write_yaml(manifest, file.path(dirs$input, "run_manifest.yml"))
  if (!file.exists(file.path(dirs$audit,"issue_ledger.csv"))) write_csv_safe(data.frame(time=character(),stage=character(),severity=character(),message=character(),evidence=character()), file.path(dirs$audit,"issue_ledger.csv"))
  accepted <- cfg$input$accepted_object_path %||% NULL
  if (!is.null(accepted) && file.exists(accepted) && mode %in% c("reannotate","figures_only")) {
    target <- file.path(dirs$objects, if (mode=="reannotate") "05_clustered.rds" else "06_annotated.rds")
    if (!file.exists(target)) file.copy(accepted,target,overwrite=FALSE)
  }
  if(mode=="figures_only" && !is.null(target_root)) {
    copies <- list(
      c(file.path(target_root,"01_objects","06_annotated.rds"),file.path(dirs$objects,"06_annotated.rds")),
      c(file.path(target_root,"04_annotation","annotation_evidence.csv"),file.path(dirs$annotation,"annotation_evidence.csv")),
      c(file.path(target_root,"11_source_data","sample_celltype_composition.csv"),file.path(dirs$source,"sample_celltype_composition.csv")),
      c(file.path(target_root,"07_differential_expression","pseudobulk_edgeR_all_results.csv"),file.path(dirs$de,"pseudobulk_edgeR_all_results.csv")),
      c(file.path(target_root,"08_pathways","fgsea_all_results.csv"),file.path(dirs$pathways,"fgsea_all_results.csv"))
    )
    for(z in copies) if(file.exists(z[1])&&!file.exists(z[2])) file.copy(z[1],z[2],overwrite=FALSE)
    src_adv <- file.path(target_root,"17_advanced_modules")
    if (dir.exists(src_adv)) {
      files <- list.files(src_adv, recursive=TRUE, full.names=TRUE, all.files=TRUE, no..=TRUE)
      files <- files[file.exists(files) & !dir.exists(files)]
      for (f in files) {
        rel <- substring(normalizePath(f,winslash="/",mustWork=TRUE), nchar(normalizePath(src_adv,winslash="/",mustWork=TRUE))+2L)
        dest <- file.path(dirs$advanced, rel); dir.create(dirname(dest), recursive=TRUE, showWarnings=FALSE)
        if (!file.exists(dest)) file.copy(f,dest,overwrite=FALSE)
      }
    }
  }
  ctx
}

stage_done <- function(ctx, id) file.exists(file.path(ctx$dirs$logs, paste0(id, ".done")))

run_stage <- function(ctx, id, script) {
  force <- as_flag(ctx$cfg$analysis$force, FALSE)
  if (stage_done(ctx,id) && !force) { log_msg("SKIP", id, "already complete"); return(invisible(TRUE)) }
  path <- file.path(ctx$skill_root, "R", script)
  if (!file.exists(path)) stop("Missing stage script: ", path)
  log_path <- file.path(ctx$dirs$logs, paste0(id, "_", sub("\\.R$","",script), ".log"))
  con <- file(log_path, open="wt")
  sink(con, split=TRUE); sink(con, type="message")
  ok <- FALSE
  on.exit({sink(type="message"); sink(); close(con)}, add=TRUE)
  log_msg("START", id, script)
  tryCatch({
    env <- new.env(parent=globalenv())
    sys.source(path, envir=env)
    if (!exists("stage_main", envir=env, inherits=FALSE)) stop("stage_main(ctx) missing in ", script)
    env$stage_main(ctx)
    ok <- TRUE
  }, error=function(e){
    append_issue(ctx,id,"FAIL",conditionMessage(e),script)
    cat("ERROR:",conditionMessage(e),"\n")
  })
  if (ok) {
    writeLines(paste("completed",now_iso()), file.path(ctx$dirs$logs,paste0(id,".done")))
    log_msg("PASS",id,script)
  } else stop("Stage failed: ",id," ",script,". See ",log_path)
  invisible(ok)
}

resolve_alias <- function(nms, configured, aliases) {
  if (!is.null(configured) && !identical(configured,"auto")) return(if (configured %in% nms) configured else NA_character_)
  low <- tolower(nms)
  for (a in aliases) {
    hit <- nms[low==tolower(a)]
    if (length(hit)==1L) return(hit)
    if (length(hit)>1L) stop("Ambiguous duplicated metadata name for alias: ",a)
  }
  NA_character_
}

standardize_metadata <- function(obj, cfg) {
  md <- obj[[]]
  aliases <- list(sample=c("sample_id","sample","orig.ident","library","library_id","dataset","specimen"), patient=c("patient_id","patient","donor","donor_id","subject","case_id","individual"), condition=c("condition","group","status","treatment","response","disease"), batch=c("batch","run","lane","chemistry","library_batch","seq_batch"), tissue=c("tissue","site","organ","compartment"))
  mapping <- c(sample_id=resolve_alias(names(md),cfg$metadata$sample_column,aliases$sample), patient_id=resolve_alias(names(md),cfg$metadata$patient_column,aliases$patient), condition=resolve_alias(names(md),cfg$metadata$condition_column,aliases$condition), batch=resolve_alias(names(md),cfg$metadata$batch_column,aliases$batch), tissue=resolve_alias(names(md),cfg$metadata$tissue_column,aliases$tissue))
  for (target in names(mapping)) {
    src <- mapping[[target]]
    obj[[target]] <- if (!is.na(src)) as.character(md[[src]]) else if (target=="sample_id") "sample1" else NA_character_
  }
  obj$cell_id <- colnames(obj)
  obj$species <- cfg$input$species %||% NA_character_
  obj$input_source <- cfg$input$path
  for (nm in c("qc_pass","doublet_class","cluster_raw","cell_type_l1","cell_type_l2","annotation_confidence","annotation_reason","analysis_inclusion")) if (!nm %in% colnames(obj[[]])) obj[[nm]] <- NA
  attr(obj,"metadata_mapping") <- mapping
  obj
}

get_layer_safe <- function(obj, assay="RNA", layer="counts") {
  if (!assay %in% names(obj@assays)) stop("Assay missing: ",assay)
  layers <- SeuratObject::Layers(obj[[assay]])
  exact <- layers[layers==layer]
  if (length(exact)==1L) return(SeuratObject::LayerData(obj,assay=assay,layer=exact))
  hits <- layers[startsWith(layers,paste0(layer,"."))]
  if (length(hits)>1L) {
    tmp <- SeuratObject::JoinLayers(obj,assay=assay)
    return(SeuratObject::LayerData(tmp,assay=assay,layer=layer))
  }
  if (length(hits)==1L) return(SeuratObject::LayerData(obj,assay=assay,layer=hits))
  stop("Layer missing: ",assay,"/",layer)
}

join_layers_safe <- function(obj, assay="RNA") {
  if (!assay %in% names(obj@assays)) return(obj)
  ly <- SeuratObject::Layers(obj[[assay]])
  if (sum(startsWith(ly,"counts."))>1L || sum(startsWith(ly,"data."))>1L) obj <- SeuratObject::JoinLayers(obj,assay=assay)
  obj
}

choose_font <- function(requested="Arial") {
  if (requireNamespace("systemfonts",quietly=TRUE)) {
    fam <- unique(systemfonts::system_fonts()$family)
    if (requested %in% fam) return(requested)
    for (x in c("Arial","Liberation Sans","DejaVu Sans","sans")) if (x %in% fam || x=="sans") return(x)
  }
  "sans"
}

pub_theme <- function(ctx, base_size=NULL) {
  cfg <- ctx$cfg$figures
  fam <- choose_font(cfg$font_family %||% "Arial")
  bs <- base_size %||% cfg$base_size_pt %||% 10
  ggplot2::theme_classic(base_size=bs, base_family=fam) + ggplot2::theme(panel.grid=ggplot2::element_blank(), plot.background=ggplot2::element_rect(fill="white",colour=NA), panel.background=ggplot2::element_rect(fill="white",colour=NA), axis.text=ggplot2::element_text(size=max(8,bs-2)), strip.background=ggplot2::element_rect(fill="white",colour="#222222",linewidth=.3), strip.text=ggplot2::element_text(face="bold",size=max(8,bs-1)), legend.title=ggplot2::element_text(face="bold",size=max(8,bs-1)), legend.text=ggplot2::element_text(size=max(8,bs-2)), plot.title=ggplot2::element_text(face="bold",size=bs+2,hjust=0), plot.tag=ggplot2::element_text(face="bold",size=cfg$panel_tag_pt %||% 16), plot.margin=ggplot2::margin(4,5,4,5))
}

categorical_palette <- function(n) {
  base <- c("#4C78A8","#F58518","#54A24B","#E45756","#72B7B2","#B279A2","#FF9DA6","#9D755D","#BAB0AC","#2F4B7C","#A05195","#D45087","#F95D6A","#FF7C43","#FFA600","#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E")
  if (n<=length(base)) base[seq_len(n)] else grDevices::hcl.colors(n,"Dark 3")
}

# export_plot_set is defined only in R/lib/figure_engine.R so every figure follows the V4 sidecar and QA contract.

save_session <- function(ctx) {
  capture.output(sessionInfo(), file=file.path(ctx$dirs$manifests,"sessionInfo.txt"))
  write_csv_safe(package_table(),file.path(ctx$dirs$manifests,"package_versions.csv"))
}

score_program <- function(mat, positive, negative=character()) {
  pos <- intersect(positive,rownames(mat)); neg <- intersect(negative,rownames(mat))
  if (!length(pos)) return(rep(NA_real_,ncol(mat)))
  ps <- Matrix::colMeans(mat[pos,,drop=FALSE])
  ns <- if(length(neg)) Matrix::colMeans(mat[neg,,drop=FALSE]) else 0
  as.numeric(ps - 0.25*ns)
}

split_genes <- function(x) {
  if (is.na(x) || !nzchar(x)) character() else unique(trimws(unlist(strsplit(x,"[;,|]"))))
}

status_table <- function(status, reason, n=NA_integer_) data.frame(status=status,reason=reason,n=n,stringsAsFactors=FALSE)
