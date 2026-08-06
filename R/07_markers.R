stage_main <- function(ctx) {
  obj<-readRDS(file.path(ctx$dirs$objects,"06_annotated.rds"));obj<-join_layers_safe(obj,"RNA");Seurat::DefaultAssay(obj)<-"RNA";Seurat::Idents(obj)<-"cluster_raw"
  if(!as_flag(ctx$cfg$analysis$run_markers,TRUE)){write_csv_safe(status_table("NOT_RUN","disabled"),file.path(ctx$dirs$markers,"marker_status.csv"));return(invisible())}
  if(length(unique(obj$cluster_raw))<2){write_csv_safe(status_table("NOT_EVALUABLE","fewer than two clusters"),file.path(ctx$dirs$markers,"marker_status.csv"));return(invisible())}
  markers<-Seurat::FindAllMarkers(obj,only.pos=TRUE,min.pct=.1,logfc.threshold=.25,return.thresh=.05,verbose=FALSE)
  if(is.null(markers)||!nrow(markers)){write_csv_safe(status_table("NOT_EVALUABLE","no markers passed filters"),file.path(ctx$dirs$markers,"marker_status.csv"));return(invisible())}
  write_csv_safe(markers,file.path(ctx$dirs$markers,"all_cluster_markers.csv"))
  fcname<-intersect(c("avg_log2FC","avg_logFC"),names(markers))[1];if(is.na(fcname))fcname<-names(markers)[grep("log.*FC",names(markers),ignore.case=TRUE)[1]]
  ord<-order(markers$cluster,-markers[[fcname]]);mk<-markers[ord,];top<-do.call(rbind,lapply(split(mk,mk$cluster),head,10));write_csv_safe(top,file.path(ctx$dirs$source,"top_cluster_markers.csv"))
  genes<-unique(top$gene);av<-Seurat::AverageExpression(obj,assays="RNA",features=genes,group.by="cluster_raw",layer="data",verbose=FALSE)$RNA
  z<-t(scale(t(as.matrix(av))));z[!is.finite(z)]<-0;df<-as.data.frame(as.table(z));names(df)<-c("gene","cluster","zscore");write_csv_safe(df,file.path(ctx$dirs$source,"marker_heatmap_source.csv"))
  # Publication marker figures are generated centrally in stage 11.
}
