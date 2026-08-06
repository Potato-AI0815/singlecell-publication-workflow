stage_main <- function(ctx) {
  obj<-readRDS(file.path(ctx$dirs$objects,"03_primary_singlets.rds"));cfg<-ctx$cfg;set.seed(cfg$runtime$seed %||% 20260804)
  if(length(unique(obj$sample_id))>1L && length(SeuratObject::Layers(obj[["RNA"]]))<=2L) obj[["RNA"]]<-split(obj[["RNA"]],f=obj$sample_id)
  obj<-Seurat::NormalizeData(obj,assay="RNA",normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE)
  obj<-Seurat::FindVariableFeatures(obj,assay="RNA",selection.method="vst",nfeatures=3000,verbose=FALSE)
  norm<-cfg$analysis$normalization %||% "sct_v2"
  if(norm=="sct_v2" && requireNamespace("sctransform",quietly=TRUE)) {
    regress<-unlist(cfg$analysis$regress_covariates %||% list()); regress<-intersect(regress,colnames(obj[[]]))
    sct_try<-try(Seurat::SCTransform(obj,assay="RNA",new.assay.name="SCT",vst.flavor="v2",vars.to.regress=if(length(regress))regress else NULL,verbose=FALSE),silent=TRUE)
    if(!inherits(sct_try,"try-error")){obj<-sct_try;Seurat::DefaultAssay(obj)<-"SCT"} else {append_issue(ctx,"04","WARNING","SCTransform failed; log normalization used",as.character(sct_try));Seurat::DefaultAssay(obj)<-"RNA";obj<-Seurat::ScaleData(obj,features=Seurat::VariableFeatures(obj),verbose=FALSE)}
  } else {
    if(norm=="sct_v2") append_issue(ctx,"04","WARNING","sctransform unavailable; fell back to log normalization")
    Seurat::DefaultAssay(obj)<-"RNA";obj<-Seurat::ScaleData(obj,features=Seurat::VariableFeatures(obj),verbose=FALSE)
  }
  ncell<-ncol(obj);npcs<-cfg$analysis$pcs %||% "auto";if(identical(npcs,"auto"))npcs<-if(ncell<2000)20 else if(ncell<=50000)30 else 50;npcs<-max(15,min(50,as.integer(npcs)))
  obj<-Seurat::RunPCA(obj,npcs=npcs,verbose=FALSE,seed.use=cfg$runtime$seed %||% 20260804)
  stdev<-Seurat::Stdev(obj,"pca");pc<-data.frame(PC=seq_along(stdev),stdev=stdev,variance=stdev^2/sum(stdev^2),cumulative=cumsum(stdev^2/sum(stdev^2)))
  write_csv_safe(pc,file.path(ctx$dirs$cluster,"pca_variance.csv"));write_csv_safe(data.frame(feature=Seurat::VariableFeatures(obj)),file.path(ctx$dirs$cluster,"variable_features.csv"))
  p<-ggplot2::ggplot(pc,ggplot2::aes(PC,variance))+ggplot2::geom_line()+ggplot2::geom_point()+pub_theme(ctx)+ggplot2::labs(y="Variance fraction",title="PCA variance profile")
  export_plot_set(ctx,p,"PCA_variance_profile","extended",140,100,
    source_data=pc,parameters=list(type="pca_variance_profile",n_pcs=npcs))
  saveRDS(obj,file.path(ctx$dirs$objects,"04_normalized_pca.rds"),compress=FALSE)
}
