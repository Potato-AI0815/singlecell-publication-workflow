stage_main <- function(ctx) {
  obj<-readRDS(file.path(ctx$dirs$objects,"06_annotated.rds"));cfg<-ctx$cfg
  if(!as_flag(cfg$analysis$run_composition,TRUE)){write_csv_safe(status_table("NOT_RUN","disabled"),file.path(ctx$dirs$composition,"composition_status.csv"));return(invisible())}
  md<-obj[[]];obs<-as.data.frame(table(sample_id=md$sample_id,cell_type=md$cell_type_l1),stringsAsFactors=FALSE);names(obs)[3]<-"cell_count";grid<-expand.grid(sample_id=unique(md$sample_id),cell_type=unique(md$cell_type_l1),stringsAsFactors=FALSE);ct<-merge(grid,obs,by=c("sample_id","cell_type"),all.x=TRUE);ct$cell_count[is.na(ct$cell_count)]<-0
  totals<-aggregate(cell_count~sample_id,ct,sum);names(totals)[2]<-"sample_total";ct<-merge(ct,totals,by="sample_id");ct$percentage<-100*ct$cell_count/ct$sample_total
  cols<-c("sample_id","patient_id","condition","batch","tissue");sm0<-unique(md[,cols,drop=FALSE]);if(any(duplicated(sm0$sample_id)))stop("Contradictory sample-level metadata: multiple rows per sample_id");sm<-sm0;ct<-merge(ct,sm,by="sample_id",all.x=TRUE);write_csv_safe(ct,file.path(ctx$dirs$source,"sample_celltype_composition.csv"))
  ref<-cfg$comparisons$reference;tar<-cfg$comparisons$target;stats<-list()
  if(!all(c(ref,tar)%in%unique(ct$condition))){write_csv_safe(status_table("NOT_EVALUABLE","configured conditions absent"),file.path(ctx$dirs$composition,"composition_status.csv"))} else {
    for(tp in unique(ct$cell_type)){d<-ct[ct$cell_type==tp & ct$condition%in%c(ref,tar),];nref<-length(unique(d$sample_id[d$condition==ref]));ntar<-length(unique(d$sample_id[d$condition==tar]));p<-NA_real_;effect<-median(d$percentage[d$condition==tar])-median(d$percentage[d$condition==ref]);method<-"descriptive"
      if(as_flag(cfg$metadata$paired,FALSE) && all(!is.na(d$patient_id))){w<-reshape(d[,c("patient_id","condition","percentage")],idvar="patient_id",timevar="condition",direction="wide");a<-w[[paste0("percentage.",tar)]];b<-w[[paste0("percentage.",ref)]];ok<-complete.cases(a,b);if(sum(ok)>=3){p<-wilcox.test(a[ok],b[ok],paired=TRUE,exact=FALSE)$p.value;method<-"paired Wilcoxon";effect<-median(a[ok]-b[ok])}}
      else if(nref>=3 && ntar>=3){p<-wilcox.test(percentage~condition,data=d,exact=FALSE)$p.value;method<-"sample-level Wilcoxon"}
      stats[[tp]]<-data.frame(cell_type=tp,reference=ref,target=tar,n_reference=nref,n_target=ntar,effect_percentage_points=effect,p_value=p,method=method)
    }
    st<-do.call(rbind,stats);st$fdr<-p.adjust(st$p_value,"BH");write_csv_safe(st,file.path(ctx$dirs$composition,"composition_statistics.csv"))
  }
  # Publication composition figures are generated centrally in stage 11.

}
