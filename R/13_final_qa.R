add_check <- function(tab, id, requirement, status, evidence) {
  rbind(tab, data.frame(check_id=id, requirement=requirement, status=status,
                        evidence=as.character(evidence), stringsAsFactors=FALSE))
}

split_manifest_paths <- function(x) {
  y <- trimws(unlist(strsplit(as.character(x), ";", fixed=TRUE)))
  y[nzchar(y) & !is.na(y)]
}

module_expected_artifacts <- function(module, root) {
  rel <- switch(module,
    cnv=c("cnv_cell_scores.csv"),
    communication=c("cellchat_interactions_by_sample.csv","cellchat_by_sample.rds"),
    trajectory=c("trajectory_cell_data.csv","trajectory_curve_data.csv","slingshot_sce.rds"),
    nmf=c("nmf_meta_program_genes.csv","nmf_scored_object.rds"),
    hdWGCNA=c("hdWGCNA_modules.csv","hdWGCNA_object.rds"),
    spatial=c("spatial_coordinates_metadata.csv.gz","spatial_processed_object.rds"),
    drug_response=c("drug_response_parameters.yml"),
    virtual_knockout=c("virtual_knockout_parameters.yml","virtual_knockout_run_manifest.csv","virtual_knockout_consensus_genes.csv","virtual_knockout_target_expression_audit.csv","virtual_knockout_target_subset_evaluability.csv"), character())
  file.path(root, rel)
}

stage_main <- function(ctx) {
  checks <- data.frame(check_id=integer(), requirement=character(), status=character(), evidence=character(), stringsAsFactors=FALSE)
  id <- 0L
  add <- function(requirement,status,evidence) {id <<- id+1L; checks <<- add_check(checks,id,requirement,status,evidence)}
  mode <- ctx$cfg$analysis$mode %||% "full"

  add("Independent-figure export is enabled", if(as_flag(ctx$cfg$figures$export_individual_figures,FALSE)) "PASS" else "FAIL",
      ctx$cfg$figures$export_individual_figures %||% NA)
  add("Composite figures are forbidden", if(as_flag(ctx$cfg$figures$composite_figures_forbidden,FALSE)) "PASS" else "FAIL",
      ctx$cfg$figures$composite_figures_forbidden %||% NA)
  add("Contact-sheet generation is disabled", if(!as_flag(ctx$cfg$figures$build_contact_sheet,FALSE)) "PASS" else "FAIL",
      ctx$cfg$figures$build_contact_sheet %||% NA)

  if (mode=="audit") {
    tr <- ctx$target_result_root
    expected <- c("00_input_audit","01_objects","04_annotation","11_source_data","14_qa","15_manifests")
    add("Audit target exists",if(dir.exists(tr))"PASS" else "FAIL",tr)
    missing_dirs <- expected[!dir.exists(file.path(tr,expected))]
    add("Expected result directories exist",if(!length(missing_dirs))"PASS" else "FAIL",paste(missing_dirs,collapse=","))
    mf <- file.path(tr,"15_manifests","figure_export_manifest.csv")
    add("Figure export manifest is present",if(file.exists(mf))"PASS" else "WARNING",mf)
  } else {
    req <- c(file.path(ctx$dirs$objects,"06_annotated.rds"),file.path(ctx$dirs$text,"Methods.md"),
             file.path(ctx$dirs$manifests,"figure_export_manifest.csv"))
    if(!mode%in%c("figures_only","reannotate")) req <- c(file.path(ctx$dirs$input,"input_file_manifest.csv"),
      file.path(ctx$dirs$objects,"02_raw_standardized_seurat.rds"),file.path(ctx$dirs$source,"cell_qc_decisions.csv"),req)
    miss <- req[!file.exists(req)]
    add("Required mode-specific artifacts exist",if(!length(miss))"PASS" else "FAIL",paste(miss,collapse=";"))

    object_path <- file.path(ctx$dirs$objects,"06_annotated.rds")
    if(file.exists(object_path)) {
      obj <- readRDS(object_path)
      rawok <- tryCatch(nrow(get_layer_safe(join_layers_safe(obj,"RNA"),"RNA","counts"))>0,error=function(e)FALSE)
      add("Raw RNA counts remain available",if(rawok)"PASS" else "FAIL",paste(nrow(obj),"genes",ncol(obj),"cells"))
      needed <- c("sample_id","patient_id","condition","cluster_raw","cell_type_l1","annotation_confidence")
      missing_md <- setdiff(needed,colnames(obj[[]]))
      add("Standard metadata fields are present",if(!length(missing_md))"PASS" else "FAIL",paste(missing_md,collapse=","))
      valid_conf <- c("High","Medium","Low","Ambiguous","Doublet-like","Unresolved")
      bad <- setdiff(unique(as.character(obj$annotation_confidence)),valid_conf)
      add("Annotation uncertainty is retained",if(!length(bad))"PASS" else "FAIL",paste(bad,collapse=","))
    }

    de <- file.path(ctx$dirs$de,"pseudobulk_edgeR_all_results.csv")
    ds <- file.path(ctx$dirs$de,"de_status.csv")
    add("DE is pseudobulk or explicitly unevaluable",if(file.exists(de)||file.exists(ds)||mode=="fast")"PASS" else "WARNING",
        if(file.exists(de))"pseudobulk results present" else "status table only, unavailable, or fast mode")

    mpath <- file.path(ctx$dirs$manifests,"figure_export_manifest.csv")
    if(file.exists(mpath)) {
      m <- read_csv_safe(mpath)
      add("At least one independent figure was exported",if(nrow(m)>0L)"PASS" else "FAIL",nrow(m))
      required_cols <- c("figure_id","pdf","tiff","png","proof_png","plot_rds","parameters","source_data","visual_qa","visual_qa_status")
      mc <- setdiff(required_cols,names(m))
      add("Figure manifest follows the V5 contract",if(!length(mc))"PASS" else "FAIL",paste(mc,collapse=","))
      if(!length(mc)) {
        fcols <- c("pdf","tiff","png","proof_png","plot_rds","parameters","visual_qa")
        fok <- vapply(fcols,function(nm) all(file.exists(m[[nm]]) & file.info(m[[nm]])$size>100),logical(1))
        add("Required figure files and sidecars exist",if(all(fok))"PASS" else "FAIL",paste(names(fok),fok,sep="=",collapse=";"))
        forbidden_names <- grepl("composite|contact[_ -]?sheet|montage|overview",m$figure_id,ignore.case=TRUE)
        add("No aggregate-layout figure was exported",if(!any(forbidden_names))"PASS" else "FAIL",paste(m$figure_id[forbidden_names],collapse=","))
        unique_dirs <- dirname(normalizePath(m$pdf,winslash="/",mustWork=FALSE))
        add("Each figure has a unique directory",if(!anyDuplicated(unique_dirs))"PASS" else "FAIL",paste("duplicates",anyDuplicated(unique_dirs)))
        plot_classes <- vapply(m$plot_rds,function(x){p<-tryCatch(readRDS(x),error=function(e)NULL); if(is.null(p))"READ_ERROR" else paste(class(p),collapse="|")},character(1))
        bad_plot <- grepl("patchwork|ggarrange",plot_classes,ignore.case=TRUE) | plot_classes=="READ_ERROR" | !grepl("gg|ggplot",plot_classes,ignore.case=TRUE)
        add("Every plot RDS contains one non-composite ggplot",if(!any(bad_plot))"PASS" else "FAIL",paste(m$figure_id[bad_plot],plot_classes[bad_plot],collapse=";"))
        qa_bad <- unique(m$visual_qa_status[!m$visual_qa_status%in%c("PASS","PASS_WITH_WARNINGS")])
        add("No figure failed technical visual QA",if(!length(qa_bad))"PASS" else "FAIL",paste(qa_bad,collapse=","))
        add("Figure identifiers are unique",if(!anyDuplicated(m$figure_id))"PASS" else "FAIL",anyDuplicated(m$figure_id))
        per_figure_source_ok <- vapply(m$source_data,function(x){sp<-split_manifest_paths(x); length(sp)>0L && all(file.exists(sp) & file.info(sp)$size>0)},logical(1))
        add("Every figure has traceable figure-specific Source Data",if(all(per_figure_source_ok))"PASS" else "FAIL",paste(m$figure_id[!per_figure_source_ok],collapse=","))
      }
    } else add("Figure export manifest exists","FAIL",mpath)

    fstat <- file.path(ctx$dirs$qa,"figure_generation_status.csv")
    if(file.exists(fstat)) {
      fs <- read_csv_safe(fstat)
      add("Figure generation status table exists","PASS",paste(table(fs$status),collapse=";"))
      add("At least one figure was exported or safely reused",if(any(fs$status%in%c("EXPORTED","SKIPPED_EXISTING")))"PASS" else "FAIL",paste(table(fs$status),collapse=";"))
      failed_figures <- fs$figure_id[as.character(fs$status)=="FAILED"]
      add("No requested figure export failed",if(!length(failed_figures))"PASS" else "FAIL",paste(failed_figures,collapse=","))
    } else add("Figure generation status table exists","WARNING",fstat)

    astat <- file.path(ctx$dirs$advanced,"advanced_module_status.csv")
    if(file.exists(astat)) {
      a <- read_csv_safe(astat)
      enabled <- names(Filter(function(z) !is.null(z)&&as_flag(z$enabled,FALSE),ctx$cfg$advanced_modules))
      observed <- a[a$module%in%enabled,,drop=FALSE]
      add("Every enabled advanced module reported a status",if(all(enabled%in%a$module))"PASS" else "FAIL",paste(enabled,collapse=","))
      failed <- observed$module[observed$status=="FAILED"]
      add("No enabled advanced module crashed",if(!length(failed))"PASS" else "FAIL",paste(failed,collapse=","))
      for(mod in enabled) {
        z <- a[a$module==mod,,drop=FALSE]
        if(nrow(z)&&z$status[1]=="COMPLETED") {
          ex <- module_expected_artifacts(mod,ctx$dirs$advanced); missing_ex <- ex[!file.exists(ex)]
          add(paste("Completed module has core artifacts:",mod),if(!length(missing_ex))"PASS" else "FAIL",paste(missing_ex,collapse=";"))
        } else if(nrow(z)) add(paste("Enabled module evaluability:",mod),"WARNING",paste(z$status[1],z$reason[1]))
      }
    } else if(mode%in%c("full","resume","reannotate")) add("Advanced-module status table exists","WARNING",astat)
  }

  issues_path <- file.path(ctx$dirs$audit,"issue_ledger.csv")
  issues <- if(file.exists(issues_path)) read_csv_safe(issues_path) else data.frame(time=character(),severity=character(),stage=character(),message=character())
  run_manifest_path <- file.path(ctx$dirs$input, "run_manifest.yml")
  run_manifest <- if (file.exists(run_manifest_path)) yaml::read_yaml(run_manifest_path) else list(start_time = NA_character_)
  run_start <- suppressWarnings(as.POSIXct(run_manifest$start_time, tz = "UTC"))
  issue_time <- if ("time" %in% names(issues)) suppressWarnings(as.POSIXct(issues$time, tz = "UTC")) else rep(as.POSIXct(NA), nrow(issues))
  current_mask <- if (!is.na(run_start)) is.na(issue_time) | issue_time >= run_start else rep(TRUE, nrow(issues))
  current_fail <- issues$severity == "FAIL" & current_mask
  historical_fail <- issues$severity == "FAIL" & !current_mask
  add("No current-run blocking failure", if(any(current_fail)) "FAIL" else "PASS",
      paste0("current_FAIL=",sum(current_fail),"; historical_FAIL=",sum(historical_fail)))
  save_session(ctx)
  write_csv_safe(checks,file.path(ctx$dirs$qa,"final_QA_checks.csv"))
  overall <- if(any(checks$status=="FAIL"))"FAIL" else if(any(checks$status=="WARNING")||any(issues$severity=="WARNING"))"PASS_WITH_WARNINGS" else "PASS"
  writeLines(overall,file.path(ctx$dirs$qa,"FINAL_STATUS.txt"))
  report <- c("# Final QA report","",paste0("Overall status: **",overall,"**"),"","## Checks","",
              paste0("- ",checks$status," — ",checks$requirement," — ",checks$evidence),"","## Recorded issues","",
              if(nrow(issues))paste0("- ",issues$severity," [",issues$stage,"]: ",issues$message) else "- None")
  writeLines(report,file.path(ctx$dirs$qa,"final_QA_report.md"),useBytes=TRUE)
  if(overall=="FAIL") stop("Final QA failed. See final_QA_report.md")
}
