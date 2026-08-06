args <- commandArgs(trailingOnly=TRUE); do_install <- "--install" %in% args
cran <- c("yaml","digest","data.table","Matrix","ggplot2","Seurat","SeuratObject","scales","svglite","ragg","ggrepel","ggrastr","ggforce","concaveman","dbscan","png","sctransform")
bioc <- c("SingleCellExperiment","scDblFinder")
missing <- function(x)x[!vapply(x,requireNamespace,logical(1),quietly=TRUE)]
cat("Missing CRAN:",paste(missing(cran),collapse=", "),"\n");cat("Missing Bioconductor:",paste(missing(bioc),collapse=", "),"\n")
if(do_install){cm<-missing(cran);if(length(cm))install.packages(cm,repos="https://cloud.r-project.org");bm<-missing(bioc);if(length(bm)){if(!requireNamespace("BiocManager",quietly=TRUE))install.packages("BiocManager",repos="https://cloud.r-project.org");BiocManager::install(bm,ask=FALSE,update=FALSE)}}
req<-c("yaml","digest","data.table","Matrix","ggplot2","Seurat","SeuratObject","scales");quit(status=if(length(missing(req)))2 else 0)
