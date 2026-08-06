load_pathways <- function(ctx) {
  gmt<-ctx$cfg$pathways$gmt_path
  if(!is.null(gmt) && file.exists(gmt)){if(!requireNamespace("fgsea",quietly=TRUE))return(NULL);return(fgsea::gmtPathways(gmt))}
  if(!requireNamespace("msigdbr",quietly=TRUE))return(NULL)
  sp<-if(tolower(ctx$cfg$input$species %||% "human")=="mouse")"Mus musculus" else "Homo sapiens";fm<-names(formals(msigdbr::msigdbr));args<-list(species=sp);if("collection"%in%fm)args$collection<-ctx$cfg$pathways$collection %||% "H" else args$category<-ctx$cfg$pathways$collection %||% "H";ms<-do.call(msigdbr::msigdbr,args);split(ms$gene_symbol,ms$gs_name)
}

stage_main <- function(ctx) {
  if(!as_flag(ctx$cfg$analysis$run_pathway_analysis,TRUE)){write_csv_safe(status_table("NOT_RUN","disabled"),file.path(ctx$dirs$pathways,"pathway_status.csv"));return(invisible())}
  dep<-file.path(ctx$dirs$de,"pseudobulk_edgeR_all_results.csv");if(!file.exists(dep) || !requireNamespace("fgsea",quietly=TRUE)){write_csv_safe(status_table("NOT_RUN","DE results or fgsea unavailable"),file.path(ctx$dirs$pathways,"pathway_status.csv"));return(invisible())}
  pathways<-load_pathways(ctx);if(is.null(pathways)){write_csv_safe(status_table("NOT_RUN","gene sets unavailable"),file.path(ctx$dirs$pathways,"pathway_status.csv"));return(invisible())}
  de<-read_csv_safe(dep);outs<-list();for(tp in unique(de$cell_type)){d<-de[de$cell_type==tp,];rank<-sign(d$logFC)*sqrt(pmax(d$F,0));names(rank)<-d$gene;rank<-sort(rank[is.finite(rank)],decreasing=TRUE);rank<-rank[!duplicated(names(rank))];if(length(rank)<100)next;fg<-fgsea::fgseaMultilevel(pathways=pathways,stats=rank,minSize=ctx$cfg$pathways$minimum_size %||% 10,maxSize=ctx$cfg$pathways$maximum_size %||% 500);fg$cell_type<-tp;outs[[tp]]<-as.data.frame(fg)}
  if(!length(outs)){write_csv_safe(status_table("NOT_EVALUABLE","no ranked cell type"),file.path(ctx$dirs$pathways,"pathway_status.csv"));return(invisible())}
  all<-do.call(rbind,outs);all$leadingEdge<-vapply(all$leadingEdge,paste,character(1),collapse=";");write_csv_safe(all,file.path(ctx$dirs$pathways,"fgsea_all_results.csv"));write_csv_safe(all,file.path(ctx$dirs$source,"pathway_enrichment_source_data.csv"))
  top <- do.call(rbind, lapply(split(all, all$cell_type), function(x) head(x[order(x$padj, -abs(x$NES)), ], 20)))
  write_csv_safe(top, file.path(ctx$dirs$source, "pathway_enrichment_top_source_data.csv"))
  # Publication enrichment figures are generated centrally in stage 11, one per cell type.
}
