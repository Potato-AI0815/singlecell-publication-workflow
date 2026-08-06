#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
install_core <- "--install" %in% args || "--install-all" %in% args
install_external <- "--install-github" %in% args || "--install-all" %in% args
repair_copykat <- "--repair-copykat" %in% args

cran <- c(
  "yaml","digest","data.table","Matrix","ggplot2","Seurat","SeuratObject","svglite","ragg","scales",
  "ggrepel","ggrastr","ggforce","concaveman","dbscan","png","sctransform","msigdbr","NMF","WGCNA",
  "FNN","igraph","scTenifoldNet","remotes","jsonlite"
)
bioc <- c(
  "SingleCellExperiment","scDblFinder","SingleR","celldex","edgeR","limma","fgsea","zellkonverter",
  "glmGamPoi","slingshot","tradeSeq","BiocParallel","ComplexHeatmap","UCell","bluster","impute",
  "preprocessCore","GeneOverlap"
)
external <- c(
  CellChat = "jinworks/CellChat",
  GeneNMF = "carmonalab/GeneNMF",
  hdWGCNA = "smorabit/hdWGCNA",
  spacexr = "dmcable/spacexr"
)
manual_optional <- c(
  copykat = "Restricted optional backend. It is not auto-installed. Review upstream academic/non-profit and commercial-use terms before using --repair-copykat.",
  infercna = "Install a compatible infercna implementation manually.",
  oncoPredict = "Install oncoPredict from its official distribution after checking its license and training-resource terms.",
  scTenifoldKnk = "Optional external backend. Install from the upstream repository only after reviewing its non-standard licensing notice; it is not installed or redistributed by this project."
)

missing <- function(x) x[!vapply(x, requireNamespace, logical(1), quietly = TRUE)]
load_status <- function(x) {
  data.frame(package = x, loadable = vapply(x, requireNamespace, logical(1), quietly = TRUE), stringsAsFactors = FALSE)
}

cat("Dependency mode:", if (install_core || install_external || repair_copykat) "install/repair" else "check", "\n")
cat("Missing CRAN packages:", paste(missing(cran), collapse = ", "), "\n")
cat("Missing Bioconductor packages:", paste(missing(bioc), collapse = ", "), "\n")
cat("Missing external packages:", paste(names(external)[!vapply(names(external), requireNamespace, logical(1), quietly = TRUE)], collapse = ", "), "\n")

if (install_core) {
  cm <- missing(cran)
  if (length(cm)) install.packages(cm, repos = "https://cloud.r-project.org")
  bm <- missing(bioc)
  if (length(bm)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(bm, ask = FALSE, update = FALSE)
  }
}

if (install_external) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
  for (pkg in names(external)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing ", pkg, " from ", external[[pkg]])
      tryCatch(remotes::install_github(external[[pkg]], upgrade = "never", dependencies = TRUE),
               error = function(e) warning(pkg, " installation failed: ", conditionMessage(e)))
    }
  }
}

repair_copykat_source <- function() {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
  tmp <- tempfile("copykat_repair_")
  dir.create(tmp, recursive = TRUE)
  archive <- file.path(tmp, "copykat.tar.gz")
  url <- "https://github.com/navinlabcode/copykat/archive/refs/heads/master.tar.gz"
  message("Downloading CopyKAT source for explicit namespace repair...")
  utils::download.file(url, archive, mode = "wb", quiet = FALSE)
  utils::untar(archive, exdir = tmp)
  src <- list.dirs(tmp, recursive = FALSE, full.names = TRUE)
  src <- src[grepl("copykat", basename(src), ignore.case = TRUE)][1]
  if (is.na(src) || !dir.exists(src)) stop("Could not locate extracted CopyKAT source directory.")
  from <- file.path(src, "data", "sysdata.rda")
  to_dir <- file.path(src, "R")
  to <- file.path(to_dir, "sysdata.rda")
  if (!file.exists(from) && !file.exists(to)) stop("CopyKAT source did not contain sysdata.rda in data/ or R/.")
  if (file.exists(from) && !file.exists(to)) {
    dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(from, to)) stop("Failed to move sysdata.rda into CopyKAT R/ directory.")
  }
  remotes::install_local(src, upgrade = "never", dependencies = TRUE, force = TRUE)
  if (!requireNamespace("copykat", quietly = TRUE) ||
      !exists("full.anno", envir = asNamespace("copykat"), inherits = FALSE)) {
    stop("CopyKAT repair completed installation but full.anno is still absent from the namespace.")
  }
  message("CopyKAT namespace repair succeeded.")
}

if (repair_copykat) repair_copykat_source()

# Runtime-load audit catches packages that are installed but unusable because dependencies are absent.
audit <- unique(c("Seurat","SeuratObject","WGCNA","scDblFinder","SingleR","edgeR","fgsea", names(external), "copykat"))
status <- load_status(audit)
cat("\nRuntime load audit:\n")
print(status, row.names = FALSE)
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  message("Hint: WGCNA commonly requires Bioconductor packages impute and preprocessCore.")
}
if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  message("Hint: install/update scDblFinder and bluster with BiocManager.")
}
if (requireNamespace("copykat", quietly = TRUE) &&
    !exists("full.anno", envir = asNamespace("copykat"), inherits = FALSE)) {
  message("CopyKAT is installed but full.anno is absent from its namespace. Re-run with --repair-copykat.")
}
cat("\nManual optional backends (not auto-installed):\n")
for (nm in names(manual_optional)) cat("- ", nm, ": ", manual_optional[[nm]], "\n", sep = "")

all_required <- c("yaml","digest","data.table","Matrix","ggplot2","Seurat","SeuratObject","scales")
quit(status = if (length(missing(all_required))) 2 else 0)
