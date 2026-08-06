args <- commandArgs(trailingOnly = TRUE)
project_root <- if(length(args)>=1) args[[1]] else "C:/Users/YHN/Desktop/Qoder/GSE160763_analysis"
raw_tar <- if(length(args)>=2) args[[2]] else "C:/Users/YHN/Desktop/Qoder/GSE160763_RAW.tar"
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
profile_root <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash="/", mustWork=TRUE))
if(!file.exists(raw_tar)) stop("Raw TAR not found: ", raw_tar)
dir.create(project_root, recursive=TRUE, showWarnings=FALSE)
input_root <- file.path(project_root,"input")
raw_flat <- file.path(input_root,"raw_flat")
by_sample <- file.path(input_root,"10x_by_sample")
dir.create(raw_flat,recursive=TRUE,showWarnings=FALSE);dir.create(by_sample,recursive=TRUE,showWarnings=FALSE)
if(!length(list.files(raw_flat,recursive=TRUE,all.files=FALSE))){message("Extracting TAR...");utils::untar(raw_tar,exdir=raw_flat)} else message("raw_flat is not empty; extraction skipped")
sm <- read.csv(file.path(profile_root,"sample_metadata.csv"),check.names=FALSE,stringsAsFactors=FALSE)
find_one <- function(gsm,kind){
  pat <- switch(kind,matrix=paste0("^",gsm,".*_matrix\\.mtx\\.gz$"),barcodes=paste0("^",gsm,".*_barcodes\\.tsv\\.gz$"),features=paste0("^",gsm,".*_(features|genes)\\.tsv\\.gz$"))
  x<-list.files(raw_flat,pattern=pat,recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
  if(length(x)!=1)stop("Expected exactly one ",kind," file for ",gsm,"; found ",length(x),": ",paste(x,collapse=","))
  x
}
count_gz_lines <- function(path){con<-gzfile(path,"rt");on.exit(close(con));n<-0L;repeat{z<-readLines(con,n=100000L,warn=FALSE);if(!length(z))break;n<-n+length(z)};n}
audit<-list()
for(i in seq_len(nrow(sm))){
  gsm<-sm$geo_accession[i];sid<-sm$sample_id[i];m<-find_one(gsm,"matrix");b<-find_one(gsm,"barcodes");g<-find_one(gsm,"features")
  dest<-file.path(by_sample,sid);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  ftype<-if(grepl("_genes\\.tsv\\.gz$",basename(g),ignore.case=TRUE))"genes.tsv.gz" else "features.tsv.gz"
  outm<-file.path(dest,"matrix.mtx.gz");outb<-file.path(dest,"barcodes.tsv.gz");outg<-file.path(dest,ftype)
  if(!file.exists(outm))file.copy(m,outm,overwrite=FALSE);if(!file.exists(outb))file.copy(b,outb,overwrite=FALSE);if(!file.exists(outg))file.copy(g,outg,overwrite=FALSE)
  nc<-count_gz_lines(outb)
  audit[[i]]<-data.frame(sample_id=sid,geo_accession=gsm,condition=sm$condition[i],matrix_source=normalizePath(m,winslash="/"),barcode_source=normalizePath(b,winslash="/"),feature_source=normalizePath(g,winslash="/"),feature_filename=ftype,raw_barcodes=nc,stringsAsFactors=FALSE)
}
aud<-do.call(rbind,audit);write.csv(aud,file.path(input_root,"preparation_audit.csv"),row.names=FALSE,quote=TRUE)
file.copy(file.path(profile_root,"sample_metadata.csv"),file.path(input_root,"sample_metadata.csv"),overwrite=TRUE)
obs<-aggregate(raw_barcodes~condition,aud,sum);exp<-read.csv(file.path(profile_root,"expected_group_counts.csv"),stringsAsFactors=FALSE);chk<-merge(exp,obs,by="condition",all=TRUE);chk$matches_expected<-chk$expected_raw_barcodes==chk$raw_barcodes
write.csv(chk,file.path(input_root,"group_count_audit.csv"),row.names=FALSE)
if(sum(aud$raw_barcodes)!=40666L)stop("Raw barcode total mismatch: observed ",sum(aud$raw_barcodes),", expected 40666")
if(any(!chk$matches_expected))stop("One or more group barcode totals differ from the supplied expected totals")
writeLines(c(paste0("prepared_at=",format(Sys.time(),tz="UTC")),paste0("raw_cells=",sum(aud$raw_barcodes)),paste0("samples=",nrow(aud))),file.path(input_root,"PREPARATION_COMPLETE.txt"))
message("Preparation PASS: 8 samples and 40,666 raw barcodes")
