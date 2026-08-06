stage_main <- function(ctx) {
  req <- required_packages()
  avail <- vapply(req, requireNamespace, logical(1), quietly = TRUE)
  write_csv_safe(package_table(), file.path(ctx$dirs$input, "package_availability.csv"))
  if (any(!avail)) stop("Missing required packages: ", paste(req[!avail], collapse = ", "))

  engine_files <- c(
    "run_all.R", "R/00_utils.R", "R/lib/figure_engine.R", "R/lib/advanced_figure_engine.R",
    "R/11_figures.R", "R/13_final_qa.R", "templates/config.full.yml", "templates/config.schema.json",
    "modules/20_cnv.R", "modules/21_communication_cellchat.R", "modules/22_trajectory_slingshot.R",
    "modules/23_nmf_programs.R", "modules/24_hdwgcna.R", "modules/25_spatial.R", "modules/26_drug_response.R",
    "modules/27_virtual_knockout.R"
  )
  missing_engine <- engine_files[!file.exists(file.path(ctx$skill_root, engine_files))]
  if (length(missing_engine)) stop("Skill engine is incomplete: ", paste(missing_engine, collapse = ", "))
  if (!as_flag(ctx$cfg$figures$export_individual_figures, TRUE)) stop("This release requires figures.export_individual_figures: true")
  if (!as_flag(ctx$cfg$figures$composite_figures_forbidden, TRUE)) stop("This release requires figures.composite_figures_forbidden: true")
  if (as_flag(ctx$cfg$figures$build_contact_sheet, FALSE)) stop("This release requires figures.build_contact_sheet: false")
  if ((ctx$cfg$figures$dpi %||% 600) < 300) append_issue(ctx, "01", "WARNING", "Configured figure DPI is below 300", ctx$cfg$figures$dpi)

  ipath <- if ((ctx$cfg$analysis$mode %||% "full") == "audit") ctx$target_result_root else ctx$cfg$input$path
  if (is.null(ipath) || !file.exists(ipath)) stop("Input path does not exist: ", ipath)
  inv <- inventory_path(ipath, hash = as_flag(ctx$cfg$runtime$hash_large_files, TRUE), max_gb = ctx$cfg$runtime$hash_max_gb %||% 20)
  write_csv_safe(inv, file.path(ctx$dirs$input, "input_file_manifest.csv"))

  requested_advanced <- vapply(ctx$cfg$advanced_modules, function(z) as_flag(z$enabled, FALSE), logical(1))
  cap <- data.frame(
    capability = c("descriptive", "doublets", "integration", "annotation", "composition", "pseudobulk_DE", "enrichment", "publication_figure_engine", "CNV", "communication", "trajectory", "NMF", "hdWGCNA", "spatial", "drug_response", "virtual_knockout"),
    requested = c(
      TRUE,
      as_flag(ctx$cfg$qc$run_doublets, TRUE),
      ctx$cfg$analysis$integration != "none",
      TRUE,
      as_flag(ctx$cfg$analysis$run_composition, TRUE),
      as_flag(ctx$cfg$analysis$run_pseudobulk_de, TRUE),
      as_flag(ctx$cfg$analysis$run_pathway_analysis, TRUE),
      TRUE,
      requested_advanced
    ),
    status = "PENDING",
    reason = "prerequisites evaluated by the corresponding stage",
    stringsAsFactors = FALSE
  )
  write_csv_safe(cap, file.path(ctx$dirs$input, "analysis_capability_table.csv"))
  write_csv_safe(data.frame(file = engine_files, exists = file.exists(file.path(ctx$skill_root, engine_files)), stringsAsFactors = FALSE), file.path(ctx$dirs$input, "engine_file_audit.csv"))

  sys <- data.frame(
    item = c("R.version", "platform", "OS", "working_directory"),
    value = c(R.version.string, R.version$platform, paste(Sys.info()[c("sysname", "release")], collapse = " "), getwd()),
    stringsAsFactors = FALSE
  )
  write_csv_safe(sys, file.path(ctx$dirs$input, "runtime_environment.csv"))
  save_session(ctx)
}
