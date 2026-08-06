figure_type_description <- function(type) {
  map <- c(
    embedding_discrete = "Two-dimensional embedding coloured by the declared categorical metadata. Soft hulls, when present, are descriptive visual regions only.",
    embedding_continuous = "Two-dimensional embedding coloured by the declared continuous feature. Raw and display-capped values are retained in Source Data.",
    split_violin_box = "Split violin density with an overlaid compact boxplot. Cell-level distributions are descriptive; formal inference remains sample-level.",
    marker_dotplot = "Marker dot plot in which point size denotes the percentage of cells expressing the gene and colour denotes group-scaled average expression.",
    marker_heatmap = "Heatmap of group-averaged expression scaled within each gene for visualization.",
    composition_stacked = "Recovered-cell composition calculated separately for each biological sample.",
    composition_by_condition = "Sample-level recovered-cell proportions displayed by condition.",
    paired_composition = "Within-patient recovered-cell proportion changes connected across paired conditions.",
    volcano = "Raw-count sample-by-cell-type pseudobulk differential-expression result for one evaluable cell type.",
    enrichment_dotplot = "Gene-set enrichment derived from pseudobulk-ranked statistics for one evaluable cell type.",
    trajectory_embedding = "Exploratory Slingshot lineage geometry with curves and direction derived from the fitted object rather than hand drawing.",
    pseudotime_trend = "Exploratory feature dynamics along fitted pseudotime.",
    tradeseq_knot_diagnostics = "tradeSeq knot-number diagnostic based on median AIC across evaluated genes.",
    tradeseq_gene_ranking = "Genes ranked by adjusted significance from the declared tradeSeq trajectory test.",
    cnv_heatmap = "Expression-inferred copy-number landscape; this is not a substitute for DNA-level copy-number measurement.",
    cnv_scatter = "Cell-level diagnostic summary from expression-inferred CNA analysis.",
    communication_bubble = "Sample-aware CellChat ligand-receptor probabilities; model-derived communication scores are non-causal.",
    communication_network = "Sample-aware CellChat sender-receiver network, restricted to supported cell groups.",
    communication_contrast = "Condition contrast computed from sample-level CellChat summaries.",
    nmf_rank_diagnostics = "NMF rank diagnostic from the selected backend.",
    nmf_similarity_heatmap = "Similarity among recurrent NMF component programs used to derive meta-programs.",
    nmf_gene_heatmap = "Top genes defining recurrent NMF meta-programs.",
    nmf_program_dotplot = "NMF program activity summarized across the declared groups.",
    nmf_program_embedding = "Cell-level score for one recurrent NMF meta-program on the embedding.",
    hdwgcna_soft_power = "hdWGCNA scale-free topology diagnostic across soft-thresholding powers.",
    hdwgcna_connectivity = "Network connectivity diagnostic across hdWGCNA soft-thresholding powers.",
    hdwgcna_hub_genes = "Top intramodular hub genes ranked by module eigengene connectivity.",
    hdwgcna_trait_heatmap = "Sample-level correlations between module eigengenes and declared traits.",
    hdwgcna_dme = "Sample-level differential module eigengene summary.",
    hdwgcna_sample_heatmap = "Sample-by-module eigengene activity heatmap.",
    hdwgcna_dendrogram = "hdWGCNA gene dendrogram and module assignment from the fitted network.",
    hdwgcna_module_embedding = "Cell-level module eigengene score on the embedding.",
    spatial_discrete = "Spatial locations coloured by the declared domain or transferred label.",
    spatial_continuous = "Spatial locations coloured by expression, transfer score or deconvolution weight.",
    spatial_variable_gene_ranking = "Genes ranked by the configured spatial-variability method.",
    spatial_neighborhood_heatmap = "Permutation-based spatial neighborhood enrichment among declared labels.",
    drug_reversal_ranking = "Exploratory ranking based on reversal of a disease-associated expression signature by a drug perturbation signature.",
    drug_reversal_scatter = "Gene-wise disease and drug perturbation effects for one exploratory reversal comparison.",
    drug_score_heatmap = "Exploratory sample-by-cell-type drug-signature scores.",
    drug_score_distribution = "Exploratory distribution of a drug-signature score across declared cell groups.",
    drug_therapeutic_cluster_embedding = "Exploratory transcriptomic clusters derived from drug-signature score profiles.",
    oncopredict_heatmap = "Exploratory oncoPredict estimates applied to sample-by-cell-type pseudobulk expression; not a clinical recommendation.",
    virtual_knockout_ranking = "Genes ranked by replicate-aware, unsigned regulatory-state displacement after a predeclared virtual knockout.",
    virtual_knockout_sample_support = "Biological-sample support for recurrent virtual-knockout network displacement.",
    virtual_knockout_target_expression = "Pre-analysis target-expression and evaluability audit by biological sample.",
    virtual_knockout_manifold = "Representative aligned WT and virtual-KO network states; arrows indicate state displacement, not expression direction.",
    virtual_knockout_network = "Recurrent target-centred WT-network topology displayed for context; edges are model-derived and non-causal.",
    virtual_knockout_run_similarity = "Jaccard stability among repeated outer-subsampling virtual-knockout runs.",
    virtual_knockout_pathways = "Over-representation analysis among replicate-supported, unsigned virtual-knockout impact genes.",
    qc_violin = "Single quality-control metric displayed independently before downstream analysis.",
    legacy_plot_call = "Independent single figure exported through the V5 multi-format contract."
  )
  if (is.null(type) || !type %in% names(map)) "Publication figure generated from the executed analysis and linked Source Data." else unname(map[[type]])
}

stage_main <- function(ctx) {
  cfg <- ctx$cfg
  sm <- if (file.exists(file.path(ctx$dirs$qc, "sample_qc_summary.csv"))) read_csv_safe(file.path(ctx$dirs$qc, "sample_qc_summary.csv")) else NULL
  cp <- if (file.exists(file.path(ctx$dirs$cluster, "clustering_parameters.csv"))) read_csv_safe(file.path(ctx$dirs$cluster, "clustering_parameters.csv")) else NULL
  ev <- if (file.exists(file.path(ctx$dirs$annotation, "annotation_evidence.csv"))) read_csv_safe(file.path(ctx$dirs$annotation, "annotation_evidence.csv")) else NULL
  mf_path <- file.path(ctx$dirs$manifests, "figure_export_manifest.csv")
  mf <- if (file.exists(mf_path)) read_csv_safe(mf_path) else NULL

  methods <- c(
    "# Methods generated from the executed configuration", "",
    paste0("Input data were imported from `", cfg$input$path, "` and analysed with the bundled single-cell publication workflow V5."),
    paste0("Quality control was performed separately by sample using robust median/MAD bounds and configured guardrails; high-confidence doublets were ", if (as_flag(cfg$qc$run_doublets, TRUE)) "evaluated using scDblFinder when available" else "not requested", "."),
    paste0("Normalization used `", cfg$analysis$normalization, "`. Dimensionality reduction used PCA; multi-sample integration was requested as `", cfg$analysis$integration, "`, with integrated values restricted to alignment, clustering and visualization."),
    "Cell identities were provisional and required marker-program evidence, competing-lineage review and optional reference support. Ambiguous and unresolved clusters were retained rather than forced into named types.",
    "Between-condition differential expression, when evaluable, used raw-count sample-by-cell-type pseudobulk profiles and edgeR quasi-likelihood models. Cells were not treated as biological replicates.",
    paste0("Every V5 independent figure was exported independently as PDF, SVG when available, 600-dpi TIFF/PNG, a ", cfg$figures$proof_width_mm %||% 183, "-mm proof PNG, a thumbnail, plot RDS, parameters and panel-level Source Data."),
    "Dense embedding points could be rasterized while text, axes, labels and hulls remained vector. Soft embedding hulls were used only as descriptive visual regions.",
    "Continuous embedding displays retained raw values and recorded any quantile-based display capping. Split violins were generated only when exactly two split levels were observed.",
    "Trajectory curves and arrows were displayed only when exported from a fitted Slingshot model; no directional path was hand drawn.",
    "Enabled advanced modules were executed only after their explicit evaluability gates passed. NMF used sample-aware recurrent program discovery; hdWGCNA preserved sample identity during metacell construction; spatial analyses retained location and slice metadata; drug-response analyses were labelled exploratory.",
    "Virtual knockout, when enabled, used raw RNA counts from a predeclared coherent cell population and baseline state. scTenifoldKnk was run on repeated outer cell subsamples within each biological sample, followed by cross-sample consensus. Differential-regulation distances were interpreted as unsigned network-state displacement and not as predicted expression induction, repression or causal regulation.",
    "Only independent single figures were exported. Aggregate layouts, montages and contact sheets were not generated."
  )
  if (!is.null(sm)) methods <- c(methods, "", paste0("After QC, ", sum(sm$retained), " of ", sum(sm$total), " cells were retained across ", nrow(sm), " samples."))
  if (!is.null(cp)) methods <- c(methods, paste0("The selected clustering resolution was ", cp$selected_resolution[1], " using `", cp$reduction[1], "`."))
  if (!is.null(ev)) methods <- c(methods, paste0("Annotation produced ", sum(ev$annotation_confidence == "High"), " high-confidence, ", sum(ev$annotation_confidence == "Ambiguous"), " ambiguous and ", sum(ev$annotation_confidence == "Unresolved"), " unresolved clusters."))
  writeLines(methods, file.path(ctx$dirs$text, "Methods.md"), useBytes = TRUE)

  legends <- c("# Figure legends generated from the export manifest", "")
  if (!is.null(mf) && nrow(mf)) {
    for (i in seq_len(nrow(mf))) {
      params <- tryCatch(read_yaml_strict(mf$parameters[i]), error = function(e) list())
      type <- params$type %||% NA_character_
      legends <- c(
        legends,
        paste0("## ", mf$figure_id[i]),
        figure_type_description(type),
        paste0("Source Data and plotting parameters are stored in the figure directory. Technical visual QA: **", mf$visual_qa_status[i], "**."),
        ""
      )
    }
  } else legends <- c(legends, "No figure manifest was available.")
  writeLines(legends, file.path(ctx$dirs$text, "Figure_legends.md"), useBytes = TRUE)

  notes <- c(
    "# Restrained result notes", "",
    "Automatic annotation is provisional and requires domain review.",
    "UMAP/tSNE geometry, soft hulls and visual proximity are descriptive and are not interpreted as quantitative lineage evidence.",
    "Analyses without adequate sample-level replication are reported as NOT_EVALUABLE rather than assigned a cell-level formal P value.",
    "Split violin distributions and continuous embeddings are descriptive; formal between-condition inference remains sample-level.",
    "Fitted trajectory and pseudotime outputs are exploratory and depend on subset, root and topology assumptions.",
    "Virtual-knockout outputs are exploratory network hypotheses. They do not replace real CRISPR/siRNA/Perturb-seq experiments and do not establish expression direction, direct binding, viability effects or therapeutic benefit.",
    "Technical visual QA does not replace direct human inspection of each independent proof PNG or biological validation."
  )
  writeLines(notes, file.path(ctx$dirs$text, "Result_notes.md"), useBytes = TRUE)
}
