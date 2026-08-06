#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else "installed_dependency_licenses.csv"
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "scripts", "audit_installed_licenses.R")
skill_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
source(file.path(skill_root, "R", "00_utils.R"))
pkgs <- unique(c(required_packages(), optional_packages()))
rows <- lapply(pkgs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    return(data.frame(package = pkg, installed = FALSE, version = NA_character_, license = NA_character_, url = NA_character_, stringsAsFactors = FALSE))
  }
  d <- utils::packageDescription(pkg)
  data.frame(package = pkg, installed = TRUE,
             version = as.character(d$Version %||% NA_character_),
             license = as.character(d$License %||% NA_character_),
             url = as.character(d$URL %||% NA_character_), stringsAsFactors = FALSE)
})
write_csv_safe(do.call(rbind, rows), out)
cat("WROTE=", normalizePath(out, winslash = "/", mustWork = FALSE), "\n", sep = "")
