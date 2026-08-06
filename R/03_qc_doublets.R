robust_bounds <- function(x,k=3) {
  z <- log10(pmax(x,1)); med<-median(z,na.rm=TRUE); madv<-mad(z,constant=1,na.rm=TRUE)
  if(!is.finite(madv) || madv==0)return(c(lower=0,upper=Inf))
  c(lower=10^(med-k*madv),upper=10^(med+k*madv))
}

stage_main <- function(ctx) {
  obj <- readRDS(file.path(ctx$dirs$objects,"02_raw_standardized_seurat.rds")); cfg<-ctx$cfg
  species<-tolower(cfg$input$species %||% "human"); modality<-cfg$input$modality %||% "scRNA"
  mtpat<-if(species=="mouse")"^mt-" else "^MT-"; rbpat<-if(species=="mouse")"^Rp[sl]" else "^RP[SL]"; hbpat<-if(species=="mouse")"^Hb[ab]" else "^HB[AB]"
  obj[["percent.mt"]]<-Seurat::PercentageFeatureSet(obj,pattern=mtpat);obj[["percent.ribo"]]<-Seurat::PercentageFeatureSet(obj,pattern=rbpat);obj[["percent.hb"]]<-Seurat::PercentageFeatureSet(obj,pattern=hbpat)
  obj$complexity<-log10(pmax(obj$nFeature_RNA,1))/log10(pmax(obj$nCount_RNA,2))
  samples<-unique(as.character(obj$sample_id)); decisions<-list(); thresholds<-list(); k<-cfg$qc$mad_multiplier %||% 3
  obj$qc_pass_basic<-FALSE
  for(s in samples){idx<-which(obj$sample_id==s);bf<-robust_bounds(obj$nFeature_RNA[idx],k);bc<-robust_bounds(obj$nCount_RNA[idx],k);maxmt<-if(modality=="snRNA")cfg$qc$max_mito_snRNA %||% 10 else cfg$qc$max_mito_scRNA %||% 20
    pass<-obj$nFeature_RNA[idx]>=max(cfg$qc$min_features %||% 200,bf["lower"]) & obj$nFeature_RNA[idx]<=min(cfg$qc$max_features_guardrail %||% 8000,bf["upper"]) & obj$nCount_RNA[idx]>=max(cfg$qc$min_counts %||% 500,bc["lower"]) & obj$nCount_RNA[idx]<=bc["upper"] & obj$percent.mt[idx]<=maxmt
    obj$qc_pass_basic[idx]<-pass
    thresholds[[s]]<-data.frame(sample_id=s,min_features=max(cfg$qc$min_features %||% 200,bf["lower"]),max_features=min(cfg$qc$max_features_guardrail %||% 8000,bf["upper"]),min_counts=max(cfg$qc$min_counts %||% 500,bc["lower"]),max_counts=bc["upper"],max_percent_mt=maxmt)
  }
  obj$doublet_class<-"not_run";obj$doublet_score<-NA_real_
  if(as_flag(cfg$qc$run_doublets,TRUE) && requireNamespace("scDblFinder",quietly=TRUE) && requireNamespace("SingleCellExperiment",quietly=TRUE)) {
    dbl <- try({basic<-subset(obj,cells=colnames(obj)[obj$qc_pass_basic]); sce<-Seurat::as.SingleCellExperiment(basic,assay="RNA"); scDblFinder::scDblFinder(sce,samples=as.character(SummarizedExperiment::colData(sce)$sample_id))},silent=TRUE)
    if(!inherits(dbl,"try-error")){sce<-dbl;m<-match(colnames(obj),colnames(sce)); obj$doublet_class[!is.na(m)]<-as.character(sce$scDblFinder.class[m[!is.na(m)]]); obj$doublet_score[!is.na(m)]<-as.numeric(sce$scDblFinder.score[m[!is.na(m)]])} else append_issue(ctx,"03","WARNING","scDblFinder failed; QC-passing cells retained with doublet status not_run",as.character(dbl))
  } else append_issue(ctx,"03","WARNING","scDblFinder not run; primary analysis retains QC-passing cells as unclassified singlets")
  obj$qc_pass<-obj$qc_pass_basic & !obj$doublet_class %in% c("doublet")
  obj$analysis_inclusion<-ifelse(obj$qc_pass,"primary","excluded")
  qctab<-obj[[]][,c("cell_id","sample_id","nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.hb","complexity","qc_pass_basic","doublet_score","doublet_class","qc_pass","analysis_inclusion")]
  write_csv_safe(qctab,file.path(ctx$dirs$source,"cell_qc_decisions.csv"));write_csv_safe(do.call(rbind,thresholds),file.path(ctx$dirs$qc,"sample_qc_thresholds.csv"))
  sm<-do.call(rbind,lapply(split(obj$qc_pass,obj$sample_id),function(x)data.frame(retained=sum(x),total=length(x),fraction=mean(x))));sm$sample_id<-rownames(sm);rownames(sm)<-NULL;sm<-sm[,c("sample_id","retained","total","fraction")];sm$analysis_class<-ifelse(sm$retained>=(cfg$qc$minimum_cells_standard %||% 200),"standard",ifelse(sm$retained>=(cfg$qc$minimum_cells_descriptive %||% 50),"descriptive","audit_only"));write_csv_safe(sm,file.path(ctx$dirs$qc,"sample_qc_summary.csv"))
  if(any(sm$fraction<0.4)) append_issue(ctx,"03","WARNING","One or more samples lost >60% cells after QC",paste(sm$sample_id[sm$fraction<0.4],collapse=","));if(any(sm$analysis_class=="audit_only"))append_issue(ctx,"03","WARNING","Samples with fewer than the descriptive minimum were excluded from clustering",paste(sm$sample_id[sm$analysis_class=="audit_only"],collapse=","))
  eligible<-sm$sample_id[sm$analysis_class!="audit_only"];obj$analysis_inclusion<-ifelse(obj$qc_pass & obj$sample_id%in%eligible,"primary",ifelse(obj$qc_pass,"audit_only","excluded"));qctab$analysis_inclusion<-obj$analysis_inclusion;write_csv_safe(qctab,file.path(ctx$dirs$source,"cell_qc_decisions.csv"))
  saveRDS(obj,file.path(ctx$dirs$objects,"03_full_qc_annotated.rds"),compress=FALSE)
  primary<-subset(obj,cells=colnames(obj)[obj$analysis_inclusion=="primary"]); if(ncol(primary)<50)stop("Fewer than 50 eligible cells remain after QC")
  saveRDS(primary,file.path(ctx$dirs$objects,"03_primary_singlets.rds"),compress=FALSE)
  df<-obj[[]];df$status<-factor(ifelse(df$qc_pass,"Retained","Excluded"),levels=c("Retained","Excluded"))
  qc_pal<-c(Retained="#6FA3D2",Excluded="#D0D0D0")
  qc_plot<-function(metric,ylabel,title,log_scale=FALSE){
    p<-ggplot2::ggplot(df,ggplot2::aes(status,.data[[metric]],fill=status))+
      ggplot2::geom_violin(scale="width",trim=TRUE,linewidth=.3,alpha=.72)+
      ggplot2::geom_boxplot(width=.16,outlier.shape=NA,linewidth=.35,fill="white",alpha=.72)+
      ggplot2::scale_fill_manual(values=qc_pal)+pub_theme(ctx)+ggplot2::guides(fill="none")+
      ggplot2::labs(x=NULL,y=ylabel,title=title)
    if(log_scale)p<-p+ggplot2::scale_y_log10()
    p
  }
  p1<-qc_plot("nFeature_RNA","Detected genes","QC: detected genes",TRUE)
  p2<-qc_plot("nCount_RNA","UMI counts","QC: library size",TRUE)
  p3<-qc_plot("percent.mt","Mitochondrial %","QC: mitochondrial fraction",FALSE)
  export_plot_set(ctx,p1,"QC_01_detected_genes","extended",150,120,
    source_data=data.frame(cell_id=colnames(obj),sample_id=obj$sample_id,status=as.character(df$status),nFeature_RNA=obj$nFeature_RNA,qc_pass=obj$qc_pass),
    parameters=list(type="qc_violin",metric="nFeature_RNA"))
  export_plot_set(ctx,p2,"QC_02_library_size","extended",150,120,
    source_data=data.frame(cell_id=colnames(obj),sample_id=obj$sample_id,status=as.character(df$status),nCount_RNA=obj$nCount_RNA,qc_pass=obj$qc_pass),
    parameters=list(type="qc_violin",metric="nCount_RNA"))
  export_plot_set(ctx,p3,"QC_03_mitochondrial_fraction","extended",150,120,
    source_data=data.frame(cell_id=colnames(obj),sample_id=obj$sample_id,status=as.character(df$status),percent_mt=obj$percent.mt,qc_pass=obj$qc_pass),
    parameters=list(type="qc_violin",metric="percent.mt"))
}
