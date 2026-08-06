stage_main <- function(ctx) {
  obj<-readRDS(file.path(ctx$dirs$objects,"04_normalized_pca.rds"));cfg<-ctx$cfg;set.seed(cfg$runtime$seed %||% 20260804)
  npcs<-ncol(Seurat::Embeddings(obj,"pca"));dims<-seq_len(min(30,npcs))
  obj<-Seurat::RunUMAP(obj,reduction="pca",dims=dims,reduction.name="umap.unintegrated",seed.use=cfg$runtime$seed %||% 20260804,verbose=FALSE)
  sm<-unique(obj[[]][,c("sample_id","condition","batch"),drop=FALSE]);sm<-sm[!duplicated(sm$sample_id),];if(all(!is.na(sm$condition))&&all(!is.na(sm$batch))){tb<-table(sm$batch,sm$condition);if(all(rowSums(tb>0)==1) && all(colSums(tb>0)==1))append_issue(ctx,"05","WARNING","Batch is perfectly confounded with condition; integration cannot identify a unique batch correction")}
  use_red<-"pca";integration_status<-"not_requested"
  multi<-length(unique(obj$sample_id))>1L; method<-cfg$analysis$integration %||% "auto"
  if(multi && method!="none") {
    integration_status<-"attempted_rpca"
    assay<-Seurat::DefaultAssay(obj); normmethod<-if(assay=="SCT")"SCT" else "LogNormalize"
    res<-try(Seurat::IntegrateLayers(object=obj,method=Seurat::RPCAIntegration,orig.reduction="pca",new.reduction="integrated.rpca",assay=assay,normalization.method=normmethod,verbose=FALSE),silent=TRUE)
    if(!inherits(res,"try-error")){obj<-res;use_red<-"integrated.rpca";integration_status<-"rpca_success"} else {integration_status<-"rpca_failed_unintegrated_used";append_issue(ctx,"05","WARNING","RPCA integration failed; unintegrated PCA used",as.character(res))}
  }
  obj<-Seurat::FindNeighbors(obj,reduction=use_red,dims=dims,graph.name=c("analysis_nn","analysis_snn"),verbose=FALSE)
  grid<-unlist(cfg$analysis$resolution_grid %||% c(.2,.4,.6,.8,1,1.2));if((cfg$analysis$mode %||% "full")=="fast")grid<-intersect(grid,c(.4,.6,.8))
  metrics<-list()
  for(r in grid){cn<-paste0("cluster_res_",gsub("\\.","_",r));obj<-Seurat::FindClusters(obj,graph.name="analysis_snn",resolution=r,cluster.name=cn,random.seed=cfg$runtime$seed %||% 20260804,verbose=FALSE);tab<-table(obj[[cn]][,1]);dom<-tapply(obj$sample_id,obj[[cn]][,1],function(x)max(prop.table(table(x))));metrics[[as.character(r)]]<-data.frame(resolution=r,n_clusters=length(tab),fraction_small=mean(tab<(cfg$analysis$minimum_cluster_cells %||% 20)),median_dominant_sample=median(dom,na.rm=TRUE))}
  met<-do.call(rbind,metrics);write_csv_safe(met,file.path(ctx$dirs$cluster,"resolution_sweep_metrics.csv"))
  sel<-cfg$analysis$clustering_resolution %||% "auto";if(identical(sel,"auto")){valid<-met$resolution[met$n_clusters>=2 & met$n_clusters<=50 & met$fraction_small<=.2 & met$median_dominant_sample<.95];sel<-if(length(valid))min(valid) else 0.6;if(!length(valid))append_issue(ctx,"05","WARNING","No resolution passed automatic structure checks; 0.6 selected provisionally")};sel<-as.numeric(sel)
  chosen<-paste0("cluster_res_",gsub("\\.","_",sel));if(!chosen%in%colnames(obj[[]])){obj<-Seurat::FindClusters(obj,graph.name="analysis_snn",resolution=sel,cluster.name=chosen,random.seed=cfg$runtime$seed %||% 20260804,verbose=FALSE)}
  obj$cluster_raw<-as.character(obj[[chosen]][,1]);obj$seurat_clusters<-obj$cluster_raw;Seurat::Idents(obj)<-"seurat_clusters"
  obj<-Seurat::RunUMAP(obj,reduction=use_red,dims=dims,reduction.name="umap",seed.use=cfg$runtime$seed %||% 20260804,verbose=FALSE)
  write_csv_safe(data.frame(integration_status=integration_status,reduction=use_red,selected_resolution=sel,pcs=paste(dims,collapse=",")),file.path(ctx$dirs$cluster,"clustering_parameters.csv"))
  saveRDS(obj,file.path(ctx$dirs$objects,"05_clustered.rds"),compress=FALSE)
  # The authoritative integrated cluster figure is generated once in stage 11.
  # Retain one independent pre-integration diagnostic using the same independent-figure engine.
  pre_fig<-plot_embedding_discrete(ctx,obj,"sample_id","umap.unintegrated",
    "Unintegrated embedding by sample",labels=FALSE,hulls=FALSE)
  export_publication_figure(ctx,pre_fig,"Embedding_preintegration_by_sample","extended",183,155)
}
