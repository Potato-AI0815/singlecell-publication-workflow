resolve_annotation_resource <- function(ctx, configured, default = NULL, required = FALSE) {
  x <- configured %||% "auto"
  if (identical(x, "auto")) x <- default
  if (is.null(x) || identical(x, "none")) return(NULL)
  candidates <- unique(c(as.character(x), file.path(ctx$skill_root, as.character(x)), file.path(dirname(ctx$config_path), as.character(x))))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    if (required) stop("Annotation resource does not exist: ", x)
    return(NULL)
  }
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

map_reference_l1 <- function(x, species = "human") {
  z <- tolower(as.character(x))
  out <- rep(NA_character_, length(z))
  common <- list(
    `T cell` = "t cell|t_cells|cd4|cd8|regulatory t|memory t|naive t",
    NK = "natural killer|\bnk\b|nk cell",
    `B cell` = "b cell|b_cells|memory b|naive b",
    `Plasma cell` = "plasma|plasmablast",
    `Dendritic cell` = "dendritic|\bdc\b|conventional dc",
    Neutrophil = "neutrophil",
    `Mast cell` = "mast",
    Erythroid = "eryth|red blood",
    Endothelial = "endothelial"
  )
  human <- list(
    Monocyte = "monocyte",
    Macrophage = "macrophage|kupffer",
    `pDC-like` = "plasmacytoid|\bpdc\b",
    `Megakaryocyte/platelet` = "megakaryocyte|platelet",
    Fibroblast = "fibroblast|stromal",
    `Pericyte/smooth muscle` = "pericyte|smooth muscle|vascular muscle",
    Epithelial = "epithelial|enterocyte|colonocyte|goblet|tuft|paneth"
  )
  mouse_cns <- list(
    Microglia = "microglia",
    `Peripheral myeloid-like` = "monocyte|macrophage|myeloid",
    Astrocyte = "astrocyte",
    Oligodendrocyte = "oligodendrocyte",
    OPC = "oligodendrocyte precursor|\bopc\b",
    Neuron = "neuron|neuronal",
    Pericyte = "pericyte",
    `Vascular smooth muscle` = "smooth muscle|vascular muscle",
    `Fibroblast/VLMC` = "fibroblast|vlmc|leptomeningeal",
    Ependymal = "ependymal",
    `Choroid plexus epithelial` = "choroid|epithelial"
  )
  rules <- c(common, if (tolower(species) == "human") human else mouse_cns)
  for (n in names(rules)) out[grepl(rules[[n]], z)] <- n
  out
}


run_singleR_clusters <- function(obj,cfg) {
  if(!as_flag(cfg$annotation$run_singleR,TRUE) || !all(vapply(c("SingleR","SingleCellExperiment","celldex"),requireNamespace,logical(1),quietly=TRUE))) return(NULL)
  tryCatch({
    obj2<-join_layers_safe(obj,"RNA");sce<-Seurat::as.SingleCellExperiment(obj2,assay="RNA")
    sp<-tolower(cfg$input$species %||% "human")
    ref<-if(sp=="mouse")celldex::MouseRNAseqData() else celldex::HumanPrimaryCellAtlasData()
    pr<-SingleR::SingleR(test=sce,ref=ref,labels=ref$label.main,clusters=obj2$cluster_raw)
    data.frame(cluster=rownames(pr),reference_label=pr$labels,reference_pruned=pr$pruned.labels,stringsAsFactors=FALSE)
  },error=function(e)NULL)
}

read_marker_dictionary <- function(path, level = "l1") {
  d <- read_csv_safe(path)
  need <- if (level == "l1") c("label","positive","negative") else c("parent","label","positive","negative")
  miss <- setdiff(need,names(d)); if(length(miss)) stop("Marker dictionary missing columns: ",paste(miss,collapse=","))
  d
}

cluster_score_summary <- function(md, clusters, score_names, labels) {
  do.call(rbind,lapply(clusters,function(k){
    idx<-md$cluster_raw==k
    s<-vapply(score_names,function(nm){v<-md[idx,nm,drop=TRUE];x<-mean(v,na.rm=TRUE);if(is.finite(x))x else -Inf},numeric(1))
    if(all(!is.finite(s) | s == -Inf)) return(data.frame(cluster=k,marker_label=NA_character_,marker_score=NA_real_,runner_up=NA_character_,runner_score=NA_real_,marker_margin=NA_real_,cluster_size=sum(idx),samples_present=length(unique(md$sample_id[idx])),stringsAsFactors=FALSE))
    ord<-order(s,decreasing=TRUE); top<-ord[[1]]; second<-if(length(ord)>1)ord[[2]] else ord[[1]]
    data.frame(cluster=k,marker_label=labels[[top]],marker_score=s[[top]],runner_up=labels[[second]],runner_score=s[[second]],marker_margin=if(length(ord)>1)s[[top]]-s[[second]] else Inf,cluster_size=sum(idx),samples_present=length(unique(md$sample_id[idx])),stringsAsFactors=FALSE)
  }))
}

stage_main <- function(ctx) {
  obj<-readRDS(file.path(ctx$dirs$objects,"05_clustered.rds"));cfg<-ctx$cfg;obj<-join_layers_safe(obj,"RNA");Seurat::DefaultAssay(obj)<-"RNA"
  if(!"data"%in%SeuratObject::Layers(obj[["RNA"]]))obj<-Seurat::NormalizeData(obj,verbose=FALSE)
  mat<-get_layer_safe(obj,"RNA","data");sp<-tolower(cfg$input$species %||% "human")
  default_l1<-if(sp=="mouse")file.path(ctx$skill_root,"resources","marker_dictionary_mouse_l1.csv") else file.path(ctx$skill_root,"resources","marker_dictionary_human_l1.csv")
  d1path<-resolve_annotation_resource(ctx,cfg$annotation$dictionary_l1,default_l1,required=TRUE)
  d1<-read_marker_dictionary(d1path,"l1")
  score_names<-make.unique(paste0("score_l1_",make.names(d1$label)))
  for(i in seq_len(nrow(d1)))obj[[score_names[[i]]]]<-score_program(mat,split_genes(d1$positive[i]),split_genes(d1$negative[i]))
  md<-obj[[]];cl<-sort(unique(as.character(md$cluster_raw)));agg<-cluster_score_summary(md,cl,score_names,as.character(d1$label))
  sr<-run_singleR_clusters(obj,cfg);if(is.null(sr)){sr<-data.frame(cluster=cl,reference_label=NA_character_,reference_pruned=NA_character_);append_issue(ctx,"06","WARNING","SingleR unavailable, disabled or reference retrieval failed; marker-program annotation used")}
  ev<-merge(agg,sr,by="cluster",all.x=TRUE);ev$reference_l1<-map_reference_l1(ifelse(is.na(ev$reference_pruned),ev$reference_label,ev$reference_pruned), species = sp);ev$agreement<-ifelse(is.na(ev$reference_l1),"reference_unavailable",ifelse(!is.na(ev$marker_label)&ev$reference_l1==ev$marker_label,"agree","conflict"))
  minscore<-cfg$annotation$marker_score_min %||% .05;minmargin<-cfg$annotation$marker_margin_min %||% .02
  unresolved <- cfg$annotation$unresolved_label %||% "Unresolved"
  ev$cell_type_l1<-ifelse(!is.finite(ev$marker_score) | !is.finite(ev$marker_margin) | ev$marker_score<minscore | ev$marker_margin<minmargin,unresolved,ev$marker_label)
  ev$annotation_confidence<-ifelse(ev$cell_type_l1==unresolved,"Unresolved",ifelse(ev$agreement=="agree","High",ifelse(ev$agreement=="conflict","Ambiguous","Medium")))
  ev$annotation_reason<-paste0("marker=",ifelse(is.na(ev$marker_label),"NA",ev$marker_label),"; score=",round(ev$marker_score,3),"; margin=",round(ev$marker_margin,3),"; reference=",ifelse(is.na(ev$reference_l1),"NA",ev$reference_l1),"; ",ev$agreement)
  map<-setNames(ev$cell_type_l1,ev$cluster);conf<-setNames(ev$annotation_confidence,ev$cluster);reason<-setNames(ev$annotation_reason,ev$cluster);obj$cell_type_l1<-unname(map[as.character(obj$cluster_raw)]);obj$annotation_confidence<-unname(conf[as.character(obj$cluster_raw)]);obj$annotation_reason<-unname(reason[as.character(obj$cluster_raw)])
  obj$cell_type_l2<-NA_character_
  default_l2<-if(sp=="human")file.path(ctx$skill_root,"resources","marker_dictionary_human_l2.csv") else NULL
  d2path<-resolve_annotation_resource(ctx,cfg$annotation$dictionary_l2,default_l2,required=FALSE)
  if(!is.null(d2path)){
    d2<-read_marker_dictionary(d2path,"l2");sn<-make.unique(paste0("score_l2_",make.names(d2$label)))
    for(i in seq_len(nrow(d2)))obj[[sn[[i]]]]<-score_program(mat,split_genes(d2$positive[i]),split_genes(d2$negative[i]))
    md2<-obj[[]]
    for(k in cl){
      idx<-which(md2$cluster_raw==k);parent<-unique(md2$cell_type_l1[idx]);cand<-which(d2$parent%in%parent)
      if(length(cand) && length(idx)>=(cfg$annotation$minimum_cluster_cells_fine %||% 50)){
        s<-vapply(sn[cand],function(nm){x<-mean(md2[idx,nm,drop=TRUE],na.rm=TRUE);if(is.finite(x))x else -Inf},numeric(1));ord<-order(s,decreasing=TRUE)
        if(length(ord)){margin<-s[ord[1]]-if(length(ord)>1)s[ord[2]] else Inf;if(is.finite(s[ord[1]])&&s[ord[1]]>=minscore&&margin>=minmargin)obj$cell_type_l2[idx]<-d2$label[cand[ord[1]]]}
      }
    }
  }
  obj$cell_type_l2[is.na(obj$cell_type_l2)]<-obj$cell_type_l1[is.na(obj$cell_type_l2)]
  ev$dictionary_l1<-d1path;ev$dictionary_l2<-d2path %||% NA_character_
  write_csv_safe(ev,file.path(ctx$dirs$annotation,"annotation_evidence.csv"));write_csv_safe(obj[[]][,c("cell_id","sample_id","cluster_raw","cell_type_l1","cell_type_l2","annotation_confidence","annotation_reason")],file.path(ctx$dirs$source,"cell_annotation.csv"))
  saveRDS(obj,file.path(ctx$dirs$objects,"06_annotated.rds"),compress=FALSE)
}
