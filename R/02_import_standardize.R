read_table_counts <- function(path) {
  x <- data.table::fread(path,data.table=FALSE,check.names=FALSE)
  if(ncol(x)<3 || nrow(x)<100) stop("Count table is too small or malformed")
  genes <- as.character(x[[1]]); x[[1]] <- NULL
  m <- as.matrix(x); storage.mode(m) <- "numeric"; rownames(m) <- make.unique(genes)
  if(any(m<0,na.rm=TRUE)) stop("Negative values detected")
  if(nrow(m)<100 || ncol(m)<20) stop("Matrix requires >=100 genes and >=20 cells")
  Matrix::Matrix(m,sparse=TRUE)
}

read_metadata_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  else data.table::fread(path, data.table = FALSE, check.names = FALSE)
}

is_complete_10x_dir <- function(path) {
  if (!dir.exists(path)) return(FALSE)
  f <- tolower(basename(list.files(path, full.names = FALSE)))
  any(grepl("^matrix\\.mtx(\\.gz)?$", f)) &&
    any(grepl("^barcodes\\.tsv(\\.gz)?$", f)) &&
    any(grepl("^(features|genes)\\.tsv(\\.gz)?$", f))
}

list_complete_10x_dirs <- function(path) {
  if (!dir.exists(path)) return(character())
  subs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
  root_norm <- normalizePath(path, winslash = "/", mustWork = TRUE)
  subs <- subs[normalizePath(subs, winslash = "/", mustWork = TRUE) != root_norm]
  sort(subs[vapply(subs, is_complete_10x_dir, logical(1))])
}

detect_type <- function(path) {
  if(dir.exists(path)) {
    if (is_complete_10x_dir(path)) return("10x_mtx")
    if (length(list_complete_10x_dirs(path))) return("multi_10x")
  }
  ext <- tolower(tools::file_ext(path))
  if(ext=="h5") return("10x_h5")
  if(ext=="h5ad") return("h5ad")
  if(ext=="rds") return("rds")
  if(ext %in% c("csv","tsv","txt")) return("table")
  stop("Unable to detect input type")
}

stage_main <- function(ctx) {
  cfg <- ctx$cfg; path <- cfg$input$path; typ <- cfg$input$type %||% "auto"; if(typ=="auto") typ <- detect_type(path)
  obj <- switch(typ,
    "10x_mtx"={m<-Seurat::Read10X(path); if(is.list(m))m<-m[[1]]; Seurat::CreateSeuratObject(m,project=cfg$project_name,min.cells=0,min.features=0)},
    "10x_h5"={m<-Seurat::Read10X_h5(path); if(is.list(m))m<-m[[1]]; Seurat::CreateSeuratObject(m,project=cfg$project_name,min.cells=0,min.features=0)},
    "multi_10x"={
      subs <- list_complete_10x_dirs(path)
      if (!length(subs)) stop("No complete immediate 10x sample directories were found under: ", path)
      ids <- basename(subs)
      if (anyDuplicated(ids)) stop("Duplicated sample directory names in multi_10x input")
      objs <- lapply(seq_along(subs),function(i){
        m<-Seurat::Read10X(subs[[i]]);if(is.list(m))m<-m[[1]]
        z<-Seurat::CreateSeuratObject(m,project=ids[[i]],min.cells=0,min.features=0)
        z$sample_id<-ids[[i]]
        z
      })
      if(length(objs)==1)objs[[1]] else merge(objs[[1]],y=objs[-1],add.cell.ids=ids)
    },
    "table"={m<-read_table_counts(path);Seurat::CreateSeuratObject(m,project=cfg$project_name,min.cells=0,min.features=0)},
    "h5ad"={if(!requireNamespace("zellkonverter",quietly=TRUE))stop("zellkonverter required for h5ad");sce<-zellkonverter::readH5AD(path);Seurat::as.Seurat(sce,counts=if("counts"%in%SummarizedExperiment::assayNames(sce))"counts" else SummarizedExperiment::assayNames(sce)[1])},
    "seurat_rds"=readRDS(path),
    "sce_rds"={sce<-readRDS(path);Seurat::as.Seurat(sce,counts=if("counts"%in%SummarizedExperiment::assayNames(sce))"counts" else SummarizedExperiment::assayNames(sce)[1])},
    "rds"={x<-readRDS(path);if(inherits(x,"Seurat"))x else if(inherits(x,"SingleCellExperiment"))Seurat::as.Seurat(x,counts=if("counts"%in%SummarizedExperiment::assayNames(x))"counts" else SummarizedExperiment::assayNames(x)[1]) else stop("Unsupported RDS class")},
    stop("Unsupported input type: ",typ)
  )
  if(!inherits(obj,"Seurat")) stop("Import did not create a Seurat object")
  if(nrow(obj)<100 || ncol(obj)<20) stop("Object requires >=100 genes and >=20 cells")
  obj <- standardize_metadata(obj,cfg)

  sample_meta_path <- cfg$input$sample_metadata_path %||% NULL
  if(!is.null(sample_meta_path) && file.exists(sample_meta_path)) {
    sm <- read_metadata_table(sample_meta_path)
    key <- cfg$input$sample_metadata_key %||% "sample_id"
    if(!key %in% names(sm)) stop("Sample metadata key missing: ",key)
    if(anyDuplicated(as.character(sm[[key]]))) stop("Sample metadata key is duplicated: ",key)
    idx <- match(as.character(obj$sample_id),as.character(sm[[key]]))
    if(anyNA(idx)) stop("Sample metadata is missing sample_id values: ",paste(unique(obj$sample_id[is.na(idx)]),collapse=","))
    add <- sm[idx,setdiff(names(sm),key),drop=FALSE]
    rownames(add) <- colnames(obj)
    obj <- Seurat::AddMetaData(obj,add)
    obj <- standardize_metadata(obj,cfg)
    write_csv_safe(sm,file.path(ctx$dirs$input,"sample_metadata_input.csv"))
    write_csv_safe(data.frame(sample_id=unique(obj$sample_id),matched=unique(obj$sample_id)%in%as.character(sm[[key]])),file.path(ctx$dirs$input,"sample_metadata_match_audit.csv"))
  } else if(!is.null(sample_meta_path)) stop("Configured sample_metadata_path does not exist: ",sample_meta_path)

  if(!is.null(cfg$input$metadata_path) && file.exists(cfg$input$metadata_path)) {
    md <- read_metadata_table(cfg$input$metadata_path)
    key <- cfg$input$metadata_barcode_column %||% "cell_id"; if(!key %in% names(md)) stop("Metadata barcode key missing: ",key)
    idx <- match(colnames(obj),as.character(md[[key]])); match_rate<-mean(!is.na(idx));if(match_rate==0)stop("No cell barcodes matched the external metadata");if(match_rate<0.9)append_issue(ctx,"02","WARNING","Fewer than 90% of cells matched external metadata",round(match_rate,3));add <- md[idx,setdiff(names(md),key),drop=FALSE]; rownames(add)<-colnames(obj); obj<-Seurat::AddMetaData(obj,add)
    obj <- standardize_metadata(obj,cfg)
  }
  counts <- get_layer_safe(obj,cfg$input$assay %||% "RNA","counts")
  if(any(counts<0)) stop("Negative count values detected")
  vals<-if(inherits(counts,"sparseMatrix"))counts@x else as.numeric(counts);if(length(vals)>1000000){set.seed(cfg$runtime$seed %||% 20260804);vals<-sample(vals,1000000)};integer_like<-all(abs(vals-round(vals))<1e-6);write_csv_safe(data.frame(nonnegative=TRUE,integer_like=integer_like,n_genes=nrow(counts),n_cells=ncol(counts)),file.path(ctx$dirs$input,"count_layer_audit.csv"));if(!integer_like)append_issue(ctx,"02","WARNING","The count layer is not integer-like; count-based DE will be disabled")
  mapping <- attr(obj,"metadata_mapping")
  write_csv_safe(data.frame(role=names(mapping),source_column=unname(mapping)),file.path(ctx$dirs$input,"metadata_mapping.csv"))
  fmeta <- data.frame(feature_id=rownames(obj),feature_symbol=rownames(obj),stringsAsFactors=FALSE)
  write_csv_safe(fmeta,file.path(ctx$dirs$input,"feature_metadata.csv"))
  roles<-c("patient_id","condition","batch","tissue");for(r in roles){nu<-tapply(obj[[r]][,1],obj$sample_id,function(x)length(unique(x[!is.na(x)])));if(any(nu>1))stop("Metadata role ",r," has multiple values within sample(s): ",paste(names(nu)[nu>1],collapse=","))};sm <- unique(obj[[]][,c("sample_id","patient_id","condition","batch","tissue"),drop=FALSE]); sm$n_cells <- as.integer(table(obj$sample_id)[sm$sample_id]); write_csv_safe(sm,file.path(ctx$dirs$input,"sample_design_table.csv"))
  write_csv_safe(obj[[]],file.path(ctx$dirs$source,"all_cells_import_metadata.csv"))
  saveRDS(obj,file.path(ctx$dirs$objects,"02_raw_standardized_seurat.rds"),compress=FALSE)
  writeLines(typ,file.path(ctx$dirs$input,"detected_input_type.txt"))
}
