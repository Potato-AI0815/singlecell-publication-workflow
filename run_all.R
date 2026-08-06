#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript run_all.R /absolute/path/config.yml [--mode=full|fast|resume|reannotate|figures_only|audit|project_profile]")
config_args <- args[!grepl("^--", args)]
if (!length(config_args)) stop("A configuration YAML path is required.")
config_path <- config_args[[1]]
mode_args <- sub("^--mode=", "", grep("^--mode=", args, value = TRUE))
mode_override <- if (length(mode_args)) tail(mode_args, 1L) else NULL
valid_modes <- c("full","fast","resume","reannotate","figures_only","audit","project_profile")
if (!is.null(mode_override) && !mode_override %in% valid_modes) stop("Unsupported --mode value: ", mode_override)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
skill_root <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE)) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(skill_root, "R", "00_utils.R"))
source(file.path(skill_root, "R", "lib", "figure_engine.R"))
source(file.path(skill_root, "R", "lib", "advanced_figure_engine.R"))
ctx <- init_context(config_path, skill_root, mode_override = mode_override)
options(future.globals.maxSize = (ctx$cfg$runtime$future_max_gb %||% 16) * 1024^3)
set.seed(ctx$cfg$runtime$seed %||% 20260804)

mode <- ctx$cfg$analysis$mode %||% "full"
stages <- list(
  c("01", "01_preflight.R"), c("02", "02_import_standardize.R"), c("03", "03_qc_doublets.R"),
  c("04", "04_normalize_reduce.R"), c("05", "05_integrate_cluster.R"), c("06", "06_annotate.R"),
  c("07", "07_markers.R"), c("08", "08_composition.R"), c("09", "09_pseudobulk_de.R"),
  c("10", "10_enrichment.R"), c("11", "11_figures.R"), c("12", "12_methods_legends.R"), c("13", "13_final_qa.R")
)
if (mode == "fast") stages <- stages[!vapply(stages, function(x) x[[1]] %in% c("09", "10"), logical(1))]
if (mode == "figures_only") stages <- stages[vapply(stages, function(x) x[[1]] %in% c("11", "12", "13"), logical(1))]
if (mode == "audit") stages <- stages[vapply(stages, function(x) x[[1]] %in% c("01", "13"), logical(1))]
if (mode == "reannotate") stages <- stages[vapply(stages, function(x) x[[1]] %in% c("06", "07", "08", "09", "10", "11", "12", "13"), logical(1))]

run_enabled_advanced <- function(ctx) {
  mods <- list(
    cnv = c("20_cnv.R", "run_cnv_module"),
    communication = c("21_communication_cellchat.R", "run_communication_module"),
    trajectory = c("22_trajectory_slingshot.R", "run_trajectory_module"),
    nmf = c("23_nmf_programs.R", "run_nmf_module"),
    hdWGCNA = c("24_hdwgcna.R", "run_hdwgcna_module"),
    spatial = c("25_spatial.R", "run_spatial_module"),
    drug_response = c("26_drug_response.R", "run_drug_response_module"),
    virtual_knockout = c("27_virtual_knockout.R", "run_virtual_knockout_module")
  )
  status_path <- file.path(ctx$dirs$advanced, "advanced_module_status.csv")
  old_status <- if (file.exists(status_path)) read_csv_safe(status_path) else data.frame()
  rows <- list()
  force <- as_flag(ctx$cfg$analysis$force, FALSE)
  resume_mode <- identical(ctx$cfg$analysis$mode %||% "full", "resume")
  for (nm in names(mods)) {
    cfgm <- ctx$cfg$advanced_modules[[nm]]
    done_path <- file.path(ctx$dirs$logs, paste0("ADV_", sanitize_stem(nm), ".done"))
    if (!is.null(cfgm) && as_flag(cfgm$enabled, FALSE)) {
      prior <- if (nrow(old_status) && "module" %in% names(old_status)) old_status[as.character(old_status$module) == nm, , drop = FALSE] else data.frame()
      if (resume_mode && !force && file.exists(done_path)) {
        if (nrow(prior)) {
          res <- prior[1, setdiff(names(prior), "module"), drop = FALSE]
          res$reason <- paste0(res$reason, " [resume: reused completed module status]")
        } else res <- status_table("NOT_RUN", "resume marker exists; prior module status row was unavailable")
      } else {
        env <- new.env(parent = globalenv())
        script <- file.path(ctx$skill_root, "modules", mods[[nm]][1])
        if (!file.exists(script)) {
          res <- status_table("FAILED", paste("missing module script", script))
        } else {
          sys.source(script, envir = env)
          res <- try(env[[mods[[nm]][2]]](ctx), silent = TRUE)
          if (inherits(res, "try-error")) {
            append_issue(ctx, "ADVANCED", "WARNING", paste(nm, "module failed"), as.character(res))
            res <- status_table("FAILED", as.character(res))
          }
        }
        if (nrow(res) && !identical(as.character(res$status[1]), "FAILED")) {
          writeLines(paste(as.character(res$status[1]), now_iso()), done_path)
        }
      }
      res$module <- nm
      rows[[nm]] <- res
    } else {
      rows[[nm]] <- transform(status_table("NOT_RUN", "disabled"), module = nm)
    }
  }
  tab <- do.call(rbind, rows)
  write_csv_safe(tab, status_path)
  invisible(tab)
}

advanced_ran <- FALSE
for (s in stages) {
  if (s[[1]] == "11" && !advanced_ran && !mode %in% c("figures_only", "audit")) {
    run_enabled_advanced(ctx)
    advanced_ran <- TRUE
  }
  run_stage(ctx, s[[1]], s[[2]])
}
cat("RUN_COMPLETE=", ctx$result_root, "\n", sep = "")
