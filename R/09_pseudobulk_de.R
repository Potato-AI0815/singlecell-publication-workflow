stage_main <- function(ctx) {
  cfg<-ctx$cfg;if(!as_flag(cfg$analysis$run_pseudobulk_de,TRUE)){write_csv_safe(status_table("NOT_RUN","disabled"),file.path(ctx$dirs$de,"de_status.csv"));return(invisible())}
  if(!requireNamespace("edgeR",quietly=TRUE)){write_csv_safe(status_table("NOT_RUN","edgeR unavailable"),file.path(ctx$dirs$de,"de_status.csv"));append_issue(ctx,"09","WARNING","edgeR unavailable");return(invisible())}
  audit_path<-file.path(ctx$dirs$input,"count_layer_audit.csv");if(file.exists(audit_path)&&!isTRUE(read_csv_safe(audit_path)$integer_like[1])){write_csv_safe(status_table("NOT_EVALUABLE","count layer is not integer-like"),file.path(ctx$dirs$de,"de_status.csv"));return(invisible())}
  obj<-readRDS(file.path(ctx$dirs$objects,"06_annotated.rds"));obj<-join_layers_safe(obj,"RNA");md<-obj[[]];allowed<-unlist(cfg$comparisons$allowed_annotation_confidence %||% c("High","Medium"));excluded<-unlist(cfg$comparisons$exclude_celltypes %||% "Unresolved");keep_cells<-md$annotation_confidence%in%allowed & !md$cell_type_l1%in%excluded;if(sum(keep_cells)<50){write_csv_safe(status_table("NOT_EVALUABLE","fewer than 50 confidently annotated cells"),file.path(ctx$dirs$de,"de_status.csv"));return(invisible())};obj<-subset(obj,cells=colnames(obj)[keep_cells]);md<-obj[[]];counts<-get_layer_safe(obj,"RNA","counts")
  key<-paste(md$sample_id,md$cell_type_l1,sep="|||");lev<-unique(key);mm<-Matrix::sparse.model.matrix(~0+factor(key,levels=lev));colnames(mm)<-lev;pb<-counts%*%mm
  meta<-do.call(rbind,lapply(lev,function(k){idx<-which(key==k);data.frame(pb_id=k,sample_id=md$sample_id[idx[1]],cell_type=md$cell_type_l1[idx[1]],patient_id=md$patient_id[idx[1]],condition=md$condition[idx[1]],batch=md$batch[idx[1]],tissue=md$tissue[idx[1]],cell_count=length(idx),stringsAsFactors=FALSE)}));rownames(meta)<-meta$pb_id
  colnames(pb)<-lev;write_csv_safe(meta,file.path(ctx$dirs$de,"pseudobulk_sample_metadata.csv"));saveRDS(pb,file.path(ctx$dirs$de,"pseudobulk_raw_counts.rds"),compress=FALSE)
  ref<-cfg$comparisons$reference;tar<-cfg$comparisons$target;mincells<-cfg$comparisons$minimum_cells_per_pseudobulk %||% 20;minrep<-cfg$comparisons$minimum_replicates_per_condition %||% 3;outs<-list();stat<-list()
  for(tp in unique(meta$cell_type)){
    m<-meta[meta$cell_type==tp & meta$cell_count>=mincells & meta$condition%in%c(ref,tar),,drop=FALSE];nr<-sum(m$condition==ref);nt<-sum(m$condition==tar)
    if(nr<minrep || nt<minrep){stat[[tp]]<-data.frame(cell_type=tp,status="NOT_EVALUABLE",reason=paste0("replicates ",ref,"=",nr,", ",tar,"=",nt));next}
    paired<-as_flag(cfg$metadata$paired,FALSE) && all(!is.na(m$patient_id))
    if(paired){pt<-table(m$patient_id,m$condition);complete_pairs<-sum(rowSums(pt[,intersect(c(ref,tar),colnames(pt)),drop=FALSE]>0)==2);if(complete_pairs<minrep){stat[[tp]]<-data.frame(cell_type=tp,status="NOT_EVALUABLE",reason=paste("complete pairs",complete_pairs,"<",minrep));next}}
    ans<-try({
      y<-edgeR::DGEList(counts=pb[,m$pb_id,drop=FALSE]);cond<-relevel(factor(m$condition),ref=ref);df<-data.frame(condition=cond,patient_id=factor(m$patient_id),batch=factor(m$batch),tissue=factor(m$tissue),row.names=m$pb_id)
      covs<-intersect(unlist(cfg$comparisons$covariates %||% list()),names(df));terms<-c(if(paired)"patient_id", "condition", covs);form<-as.formula(paste("~",paste(unique(terms),collapse=" + ")));design<-model.matrix(form,df)
      if(qr(design)$rank<ncol(design))stop("design matrix is rank deficient; condition may be confounded")
      keep<-edgeR::filterByExpr(y,design=design);if(sum(keep)<50)stop("fewer than 50 expressed genes after filtering")
      y<-y[keep,,keep.lib.sizes=FALSE];y<-edgeR::calcNormFactors(y);y<-edgeR::estimateDisp(y,design);fit<-edgeR::glmQLFit(y,design,robust=TRUE)
      coefn<-grep("^condition",colnames(design),value=TRUE);if(length(coefn)!=1)stop("condition coefficient ambiguous")
      qlf<-edgeR::glmQLFTest(fit,coef=match(coefn,colnames(design)));tt<-edgeR::topTags(qlf,n=Inf,sort.by="none")$table;list(table=tt,formula=form)
    },silent=TRUE)
    if(inherits(ans,"try-error")){stat[[tp]]<-data.frame(cell_type=tp,status="NOT_EVALUABLE",reason=as.character(ans));next}
    tt<-ans$table;tt$gene<-rownames(tt);tt$cell_type<-tp;tt$contrast<-paste(tar,"vs",ref);tt$n_reference<-nr;tt$n_target<-nt;outs[[tp]]<-tt;stat[[tp]]<-data.frame(cell_type=tp,status="EVALUABLE",reason=paste("design",deparse(ans$formula)))
  }
  write_csv_safe(do.call(rbind,stat),file.path(ctx$dirs$de,"de_status.csv"));if(!length(outs)){append_issue(ctx,"09","WARNING","No cell type had sufficient pseudobulk replication");return(invisible())}
  all<-do.call(rbind,outs);write_csv_safe(all,file.path(ctx$dirs$de,"pseudobulk_edgeR_all_results.csv"));write_csv_safe(all,file.path(ctx$dirs$source,"pseudobulk_DE_source_data.csv"))
  all$neglog10FDR <- -log10(pmax(all$FDR, 1e-300))
  all$display <- all$FDR < (cfg$comparisons$fdr_threshold %||% .05) & abs(all$logFC) >= (cfg$comparisons$logfc_display_threshold %||% .5)
  write_csv_safe(all, file.path(ctx$dirs$source, "pseudobulk_DE_plot_source_data.csv"))
  # Publication volcano figures are generated centrally in stage 11, one per cell type.
}
