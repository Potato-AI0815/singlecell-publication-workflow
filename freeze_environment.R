#!/usr/bin/env Rscript
if(!requireNamespace("renv",quietly=TRUE))install.packages("renv",repos="https://cloud.r-project.org")
renv::init(bare=TRUE)
renv::snapshot(prompt=FALSE)
cat("renv.lock created for the tested target environment.
")
