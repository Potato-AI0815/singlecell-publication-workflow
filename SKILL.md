---
name: single-cell-publication-workflow
version: 0.1.0-alpha
description: Run, audit, annotate, statistically analyse and export publication-ready independent figures for human or mouse scRNA-seq/snRNA-seq projects using Seurat v5. Supports sample-level metadata manifests, configurable tissue-specific annotation dictionaries, and gated CNV, CellChat, Slingshot/tradeSeq, recurrent NMF programs, hdWGCNA, spatial transcriptomics, exploratory drug-response analysis and sample-aware virtual knockout with scTenifoldKnk. Never treat cells as biological replicates, never force unsupported annotations or advanced analyses, and never generate composite figures.
---

# SingleCell Publication Workflow v0.1.0-alpha

## 1. Purpose

Use this skill to convert a supported single-cell dataset into a reproducible analysis directory containing:

- standardized Seurat v5 objects;
- sample-aware QC and doublet decisions;
- normalization, optional integration, clustering and provisional annotation;
- cluster markers, sample-level composition, pseudobulk differential expression and enrichment;
- gated advanced analyses;
- one publication-ready file set for every independent figure;
- figure-specific Source Data, parameters, logs, manifests and QA reports.

This skill is an execution contract. Do not improvise a new workflow when a bundled stage, module or recipe already covers the request.

## V5.1 validation hardening

V5.1 adds three generic runtime safeguards discovered while preparing the first public-data validation profile:

- `input.sample_metadata_path` joins one row per biological sample onto all cells without requiring a giant cell-level metadata table;
- multi-sample 10x import filters only immediate directories that contain a complete matrix/barcode/feature triplet;
- `annotation.dictionary_l1` and `annotation.dictionary_l2` are now honored for tissue-specific marker programs, including mouse CNS profiles.

The bundled `validation_profiles/GSE160763/` profile is for runtime validation, not a frozen biological result.

V5.1.1 adds a Windows one-click launcher, explicit R exit-code propagation, resumable `--mode=` overrides, mixed 10x feature-file compatibility, and QA handling for empty optional Source Data tables. V5.1.2 additionally stores a grob-backed, renderable ggplot proxy in each plot RDS sidecar when ggplot2 would otherwise retain a Seurat-containing environment, preventing multi-gigabyte plot files and final-QA memory spikes.


## Alpha release status

This repository is an experimental public alpha. The core workflow has completed public-data runtime validation on GSE160763 and GSE132465. Paired pseudobulk, expression-inferred CNA, sample-wise CellChat and recurrent NMF programs have completed targeted runtime validation on the six-pair GSE132465 subset. hdWGCNA is partially validated through metacell, soft-power and TOM construction; trajectory, spatial, drug-response and virtual-knockout modules remain experimental until their dedicated public-data benchmarks pass.

The software does not replace expert review, orthogonal validation, clinical interpretation or a qualified statistician.

## 2. Non-negotiable scientific rules

1. Never overwrite an input object, accepted result, prior output or validated archive.
2. Create a new timestamped result directory for every run.
3. Preserve raw RNA counts. Integrated values are for alignment, clustering and visualization, not formal differential expression.
4. Use biological sample or patient as the inferential replicate. Never use cells as independent biological replicates.
5. Use cell-level tests only for descriptive cluster-marker discovery.
6. Use raw-count `sample × cell type` pseudobulk profiles for formal between-condition differential expression.
7. When replication or required inputs are inadequate, return `NOT_EVALUABLE`; do not substitute a cell-level P value.
8. Keep missing combinations as `NA`, not zero.
9. Keep `Ambiguous`, `Doublet-like` and `Unresolved` labels. Never force every cluster into a precise subtype.
10. Separate formal, sensitivity and exploratory branches. Label all advanced model-based analyses as exploratory unless independently validated.
11. Freeze subsets, gene sets, roots, reference populations and model parameters before inspecting the desired outcome.
12. Never infer lineage direction from UMAP/tSNE geometry alone.
13. Never interpret expression-inferred CNA as DNA-confirmed copy-number alteration.
14. Never interpret CellChat probability as direct physical interaction or causality.
15. Never interpret transcriptome-derived drug scores as clinical treatment recommendations.
16. Never describe virtual-knockout distance as predicted expression up- or down-regulation, and never treat it as real perturbation or causal evidence.
17. Stop on contradictory metadata, duplicated cell identifiers, missing raw counts, failed required QA or output drift.

## 3. Figure policy: independent single figures only

The workflow must not create:

- multi-panel composite figures;
- overview figures;
- patchwork/cowplot/ggarrange layouts;
- montages;
- contact sheets;
- image collages.

Every scientific display is exported as an independent figure with a unique `figure_id` and directory. A plot object with class `patchwork` or `ggarrange` must be rejected.

Each exported figure directory must contain:

```text
<figure_id>.pdf
<figure_id>.svg                         # when svglite is available
<figure_id>_600dpi.tiff
<figure_id>_600dpi.png
<figure_id>_proof183mm.png
<figure_id>_thumbnail.png
<figure_id>_plot.rds
<figure_id>_parameters.yml
<figure_id>_source_<table>.csv[.gz]
visual_QA.csv
VISUAL_QA_STATUS.txt
```

When ggplot2/S7 captures the analysis environment, `<figure_id>_plot.rds` is a lightweight grob-backed ggplot proxy of the already exported figure. Full cell-level evidence remains in the figure-specific Source Data tables; PDF/SVG/TIFF/PNG remain the publication assets.

Dense point layers may be rasterized, but text, axes, labels, curves and hulls should remain vector whenever supported. Raw values and any display-capped values must both remain traceable.

## 4. Select exactly one execution mode

- `full`: complete core workflow, enabled advanced modules, figures and QA.
- `fast`: QC, normalization, clustering, annotation and core diagnostic figures; formal DE/enrichment may be skipped.
- `resume`: continue the latest result directory without silently repeating completed stages.
- `reannotate`: start from an accepted clustered object and rerun annotation/downstream stages.
- `figures_only`: read accepted objects and Source Data from `output.existing_result_root`; do not recompute scientific results.
- `audit`: inspect an existing result directory without altering it.
- `project_profile`: apply a named project profile such as `profiles/BHSC200121.md`.

State the selected mode in the run log and final report.

## 5. Natural-language routing

Map the request to a bundled recipe before running code:

| User intent | Recipe or module |
|---|---|
| One-click basic analysis and provisional annotation | `recipes/01_full_basic_analysis.md` |
| Rapid first-pass QC and cell types | `recipes/02_fast_first_pass.md` |
| Reannotate without changing clustering | `recipes/03_reannotation_only.md` |
| Compare a cell type between groups | `recipes/04_compare_celltype.md` |
| Redraw accepted results as independent figures | `recipes/05_figures_only.md` |
| Audit an existing result directory | `recipes/06_audit_existing.md` |
| Infer tumour-like CNA states | `recipes/20_cnv.md` |
| Cell-cell communication | `recipes/21_communication.md` |
| Slingshot/tradeSeq lineage analysis | `recipes/22_trajectory.md` |
| Recurrent NMF programs | `recipes/23_nmf_programs.md` |
| hdWGCNA co-expression modules | `recipes/24_hdwgcna.md` |
| Spatial transcriptomics | `recipes/25_spatial.md` |
| Exploratory drug-response hypotheses | `recipes/26_drug_response.md` |
| Virtual knockout or GRN perturbation | `recipes/27_virtual_knockout.md` |
| Generate or regenerate independent figures | `recipes/30_publication_figure_engine.md` |

Do not combine advanced modules merely because the user says “run everything.” Evaluate each module independently.

## 6. Input contract

Supported primary inputs:

- Cell Ranger/10x matrix directory;
- 10x HDF5 matrix;
- Seurat RDS;
- SingleCellExperiment RDS;
- AnnData/H5AD when `zellkonverter` is available;
- a declared multi-sample 10x structure;
- matrix plus metadata table.

Required metadata semantics:

- `cell_id`: globally unique;
- `sample_id`: biological specimen/library identifier;
- `patient_id`: patient/donor identifier when available;
- `condition`: comparison group when applicable;
- `batch`: technical batch when available;
- `tissue`: site/compartment when available.

Do not infer paired design unless patient identifiers support it. Do not claim FASTQ-level reproducibility when the input starts from a matrix or RDS.

## 7. Core execution order

Run only through `run_all.R` unless a project profile explicitly says otherwise:

```bash
Rscript run_all.R /absolute/path/config.yml [--mode=full|fast|resume|reannotate|figures_only|audit|project_profile]
```

On Windows, the bundled one-click entry performs the dependency check, runs the selected mode, propagates failures, and verifies `14_qa/FINAL_STATUS.txt`:

```text
run_one_click.cmd -Config C:\path\to\config.yml -Mode fast
```

Use `-RscriptPath` or `-RShortcutPath` when R is not on `PATH`. A successful run ends with `FINAL_STATUS=PASS` or `PASS_WITH_WARNINGS`; a failed stage or final QA exits non-zero.

Core stages:

1. `01_preflight.R` — validate input, config, software and hashes.
2. `02_import_standardize.R` — import without overwriting and standardize metadata.
3. `03_qc_doublets.R` — sample-aware robust QC and optional scDblFinder.
4. `04_normalize_reduce.R` — LogNormalize or SCTransform v2 and PCA.
5. `05_integrate_cluster.R` — optional Seurat v5 RPCA integration, resolution sweep and UMAP.
6. `06_annotate.R` — provisional annotation with positive, competing and optional reference evidence.
7. `07_markers.R` — descriptive cluster markers and annotation diagnostics.
8. `08_composition.R` — sample-level recovered-cell composition.
9. `09_pseudobulk_de.R` — raw-count pseudobulk edgeR analysis where evaluable.
10. `10_enrichment.R` — enrichment from pseudobulk-ranked statistics.
11. Enabled advanced modules in the fixed order below.
12. `11_figures.R` — independent figure generation only.
13. `12_methods_legends.R` — methods, legends and restrained result notes.
14. `13_final_qa.R` — scientific, file and figure QA.

Do not declare completion while final QA is `FAIL`.

## 8. Provisional annotation contract

Use three evidence classes:

1. positive lineage markers or multi-gene programs;
2. competing-lineage/negative evidence;
3. optional reference mapping such as SingleR.

For each cluster export:

- top and second candidate labels;
- supporting and competing genes;
- marker-program score and margin;
- optional reference label and concordance;
- cell count and sample coverage;
- final provisional label;
- confidence: `High`, `Medium`, `Low`, `Ambiguous`, `Doublet-like` or `Unresolved`;
- human-review flag and reason.

Fine annotation is not allowed when the cluster is too small or lacks discriminating markers.

## 9. Advanced module dispatcher

Advanced scripts are stored under `modules/` and run before the figure stage:

```text
20_cnv.R
21_communication_cellchat.R
22_trajectory_slingshot.R
23_nmf_programs.R
24_hdwgcna.R
25_spatial.R
26_drug_response.R
27_virtual_knockout.R
```

Each module returns one of:

- `COMPLETED`;
- `NOT_RUN`;
- `NOT_EVALUABLE`;
- `FAILED`.

Write all statuses to `17_advanced_modules/advanced_module_status.csv`. A disabled or scientifically unevaluable module is not a fabricated negative result.

## 10. Expression-inferred CNA

### Gate

Require:

- tumour-relevant tissue context;
- explicit normal/diploid `reference_labels`;
- adequate reference cells overall and by sample;
- raw RNA counts;
- sufficient cells per sample;
- human/mouse genome compatibility with the selected backend.

### Backends

- preferred configured backend: `copykat` or `infercna`;
- do not execute deprecated `inferCNV` in this workflow.

Run by sample. Do not pool all patients to create a single pseudo-replicate.

### Required outputs

- cell-level CNA class and scores;
- ordered long-form CNA heatmap data;
- backend objects and parameters;
- sample-specific CNA heatmap figures;
- CNA signal/correlation diagnostic figure.

### Interpretation

Use “expression-inferred CNA state” or “aneuploid-like profile.” Require DNA-based or orthogonal validation before claiming a genomic alteration.

## 11. CellChat communication

### Gate

Require:

- supported cell labels;
- at least two eligible groups per sample;
- minimum cells per `sample × group`;
- multiple biological samples for condition contrasts.

Run CellChat separately by sample. Aggregate interaction probabilities at sample level for between-condition summaries. Exclude unresolved/doublet-like populations by default.

### Required outputs

- per-sample CellChat objects and status;
- ligand-receptor interactions by sample;
- sample-level interaction summaries;
- sender-receiver network edges;
- pathway-level edges;
- condition contrasts when replicated;
- independent bubble, network and contrast figures.

### Interpretation

Report model-derived communication potential, not direct interaction, signaling flux or causality.

## 12. Slingshot and tradeSeq

### Gate

Require:

- explicit lineage subset;
- biologically justified `root_cluster`;
- adequate cells and cluster coverage;
- a declared reduction;
- no outcome-guided root or branch selection.

Slingshot supplies lineage curves, nodes, pseudotime and weights. Direction arrows must come from the fitted object.

When enabled, tradeSeq must:

- use raw counts;
- evaluate knot number or use a predeclared value;
- fit negative-binomial GAMs;
- export association, start-versus-end, pattern and endpoint tests as applicable;
- adjust P values;
- preserve all tested genes.

### Independent figures

- fitted lineage embedding;
- one pseudotime trend per requested feature;
- knot diagnostic;
- one gene-ranking figure per supported tradeSeq test.

All trajectory results are exploratory and conditional on subset, root and topology assumptions.

## 13. Recurrent NMF programs

### Gate

Require:

- a biologically coherent subset;
- at least two biological samples, preferably more;
- adequate cells per sample;
- sufficient variable genes;
- a predeclared rank grid and random seed.

### Execution

1. Split the discovery object by `sample_id`.
2. Cap cells per sample by stratified sampling.
3. Prefer GeneNMF `multiNMF` plus `getMetaPrograms`.
4. Use the bundled `NMF` fallback only when GeneNMF is unavailable or fails.
5. Identify recurrent meta-programs across sample-specific components.
6. Export top genes, recurrence metrics and program similarity.
7. Score every retained cell with the frozen meta-program gene sets.
8. Do not name programs until marker/pathway review.

### Required outputs

- backend objects and standardized contract;
- rank diagnostics where available;
- meta-program gene table and weights;
- component-program similarity matrix;
- cell-level program scores;
- sample/group summaries;
- scored Seurat object;
- rank, similarity, top-gene, activity and per-program embedding figures.

Do not treat a program discovered in one sample as recurrent.

## 14. hdWGCNA

### Gate

Require:

- an explicit coherent cell population;
- at least three biological samples;
- adequate cells per sample;
- retained sample identity during metacell construction;
- a declared gene-selection rule and soft-power grid.

### Execution

1. `SetupForWGCNA` on the prespecified subset.
2. Construct metacells grouped by biological sample and cell population.
3. Test soft powers and freeze the selected power.
4. Construct a signed/declared network.
5. Calculate module eigengenes and connectivity.
6. Identify non-grey modules and hub genes.
7. Summarize eigengenes per biological sample.
8. Perform sample-level differential eigengene analysis only when replicated.
9. Correlate module eigengenes with predeclared traits at sample level.

### Required outputs

- power table and selected power;
- network object and TOM location;
- module assignments;
- hub genes and kME;
- cell-level eigengenes;
- sample-level module summary;
- differential eigengenes and trait correlations when evaluable;
- independent soft-power, connectivity, hub-gene, trait, sample-activity, dendrogram and module-embedding figures.

Modules require preservation or external validation before being described as robust biological networks.

## 15. Spatial transcriptomics

### Gate

Require:

- a supported Seurat spatial object or Space Ranger output;
- valid spatial coordinates and image/slice identifiers;
- sufficient locations per slice;
- a compatible annotated single-cell reference for mapping/deconvolution;
- preservation of slice-level replication.

### Execution

1. Load and normalize each spatial object without discarding coordinates.
2. Perform PCA, clustering and spatial-domain annotation.
3. Identify spatially variable genes with the configured method.
4. Run Seurat label transfer when a compatible reference is available.
5. Run RCTD only when reference counts and labels pass the gate.
6. Export prediction scores and deconvolution weights, not just hard labels.
7. Calculate neighborhood enrichment by permutation within appropriate spatial units.
8. Never merge slices before preserving `spatial_sample_id`.

### Required independent figures

- cluster map per slice;
- transferred-label map per slice;
- label-transfer confidence map per slice;
- cell-type transfer-score maps per slice;
- RCTD weight maps per slice;
- spatial-expression maps for selected genes per slice;
- spatially variable gene ranking;
- neighborhood-enrichment heatmap per slice.

Spatial mapping and deconvolution are model/platform dependent. Avoid interpreting spot-level labels as pure cell identities on multicellular platforms.

## 16. Exploratory drug-response analysis

### Allowed backends

1. `signature_reversal` — correlate pseudobulk disease effects with drug perturbation effects.
2. `signature_score` — score drug-response signatures across cells, then aggregate by sample and cell type.
3. `oncopredict` — apply a declared training expression/response resource to sample-by-cell-type pseudobulk expression.

### Gate

Require:

- a documented signature or training resource;
- sufficient overlapping genes;
- sample-aware disease effects or pseudobulk expression;
- explicit scale/direction semantics;
- no claim of patient-level clinical efficacy.

### Required outputs

- resource and overlap audit;
- drug ranking and gene contributions for signature reversal;
- cell and sample-level scores for signature scoring;
- optional transcriptomic therapeutic clusters;
- pseudobulk predictions for oncoPredict;
- independent ranking, contribution, heatmap, distribution and cluster-embedding figures;
- parameters containing a mandatory exploratory-use disclaimer.

Never recommend a drug, dose or clinical action from these outputs. Require external datasets, pharmacology and clinical context for validation.

## 17. Exploratory virtual knockout

### Primary backend and scope

Use `scTenifoldKnk` as the fixed R backend. It constructs a WT single-cell gene-regulatory network, simulates knockout by removing the target gene's outgoing edges, aligns WT and KO network states, and reports differential-regulation distance.

The module accepts WT/control observational scRNA-seq as input. It is a hypothesis-generation analysis, not a replacement for CRISPR, siRNA, Perturb-seq or functional validation.

Do not silently substitute:

- CellOracle, because it requires a separate Python workflow, TF/prior-GRN construction and trajectory-aware simulation contract;
- GEARS, because it requires perturbational training data and is not a generic WT-only backend.

### Gate

Require all of the following:

- one or more predeclared target genes;
- an explicit biologically coherent cell subset;
- explicit WT/control `baseline_conditions` when multiple conditions are present;
- retained raw RNA counts;
- biological `sample_id` metadata;
- sufficient cells and target expression in each evaluable sample;
- a frozen common gene universe and frozen network parameters;
- a predeclared outer `subsampling_fraction` in `(0,1]` for stability repeats;
- High/Medium annotation confidence by default.

If the target is absent or inadequately expressed, return `NOT_EVALUABLE`. Do not lower thresholds after viewing results.

### Sample-aware execution

Default to `sample_stratified_consensus`:

1. run each biological sample independently;
2. before each backend call, draw an independent seeded outer cell subset according to `subsampling_fraction`;
3. estimate repeat-level stability without relying only on backend-internal random behaviour;
4. aggregate to sample-level support;
5. derive cross-sample recurrence and consensus impact;
6. label results as replicate-consistent only when the configured minimum number of biological samples is met.

`balanced_pooled` is allowed only as a descriptive sensitivity analysis. It does not create biological replication.

### Required outputs

- target-symbol, expression and target/subset evaluability audits;
- gene-universe audit;
- run manifest with `COMPLETED`, `FAILED` and `NOT_EVALUABLE` rows;
- all differential-regulation results;
- sample-level repeat consensus;
- cross-sample consensus genes and evidence class;
- WT/KO manifold coordinates;
- recurrent WT-network edges used only as topology context;
- run-similarity table;
- optional pathway over-representation results;
- frozen parameters and explicit validation requirements.

### Required independent figures

For every evaluable target × cell subset, export separately:

- unsigned regulatory-impact ranking;
- sample-support heatmap;
- target-expression eligibility plot;
- representative WT–KO manifold displacement;
- target-centred recurrent WT-network context;
- run-similarity heatmap;
- pathway enrichment plot when evaluable.

The WT→KO arrows in manifold space show network-state displacement, not expression direction. Network-edge signs describe the fitted WT network and do not prove direct regulation.

### Interpretation

Allowed wording:

> Virtual knockout of TARGET was associated with recurrent network-state displacement involving GENE SET in CELL TYPE across N of M evaluable samples.

Forbidden wording:

> TARGET knockout upregulated/downregulated GENE.

Do not infer causality, cell viability, drug sensitivity or therapeutic benefit without real perturbation and orthogonal validation.

## 18. Independent figure engine

Use the fixed functions in `R/lib/figure_engine.R` and `R/lib/advanced_figure_engine.R`. Do not write ad hoc plotting code when a template exists.

Core templates include:

```text
plot_embedding_discrete
plot_embedding_continuous
plot_split_violin_box
plot_marker_dotplot
plot_marker_heatmap
plot_composition_stacked
plot_composition_by_condition
plot_paired_composition
plot_volcano_single
plot_enrichment_single
plot_trajectory_embedding
plot_pseudotime_trend
```

Advanced templates include:

```text
plot_advanced_heatmap
plot_rank_diagnostics
plot_program_prevalence_dotplot
plot_lollipop_ranking
plot_cnv_heatmap
plot_cnv_scatter
plot_communication_bubble
plot_network_edges
plot_soft_power_metric
plot_hub_genes
plot_spatial_discrete
plot_spatial_continuous
plot_drug_score_heatmap
plot_drug_reversal_scatter
plot_metric_curve
plot_gene_test_ranking
plot_group_score_distribution
plot_virtual_knockout_ranking
plot_virtual_knockout_sample_support
plot_virtual_knockout_manifold
plot_virtual_knockout_network
plot_virtual_knockout_pathways
plot_virtual_knockout_target_expression
```

For every figure:

- use a stable palette and readable minimum text size;
- preserve the same category colour across figures;
- use arrows only for fitted trajectory direction;
- draw soft embedding hulls only as descriptive regions;
- avoid overplotting through rasterized cell points or stratified display sampling;
- retain all cells in Source Data even if display sampling is used;
- record quantile clipping and scale limits;
- export `PASS`, `PASS_WITH_WARNINGS` or `FAIL` visual QA.

## 19. Output directory contract

```text
00_input_audit/
01_objects/
02_qc/
03_reduction_clustering/
04_annotation/
05_markers/
06_composition/
07_differential_expression/
08_pathways/
09_main_figures/          # independent primary figures
10_extended_figures/      # independent supporting figures
11_source_data/
12_methods_legends/
13_logs/
14_qa/
15_manifests/
16_plot_rds/
17_advanced_modules/
```

Do not place a composite or contact-sheet image anywhere in this structure.

## 20. Final QA requirements

Final QA must verify:

- raw counts remain available;
- metadata fields and cell IDs are valid;
- annotation uncertainty is retained;
- formal DE is pseudobulk or explicitly unevaluable;
- all enabled advanced modules reported a status;
- completed modules produced their core artifacts;
- no enabled module crashed;
- every figure has a unique directory;
- no figure ID contains `composite`, `overview`, `montage` or `contact_sheet`;
- every plot RDS is one ggplot and not patchwork/ggarrange;
- PDF, TIFF, PNG, proof PNG, parameter file and QA file exist;
- Source Data paths resolve;
- no figure failed technical visual QA;
- input objects and prior results remain unchanged.

Report `PASS`, `PASS_WITH_WARNINGS` and `FAIL` separately. Do not suppress warnings that affect interpretation.

## 21. Low-capability model operating algorithm

A lower-capability model must follow this exact sequence:

1. Read this file, the selected recipe and the configuration schema.
2. Inspect input filenames and metadata only; do not infer biology from names alone.
3. Copy `templates/config.full.yml` and fill only supported fields.
4. Keep every advanced module disabled until its required fields and scientific gate are satisfied.
5. Validate the YAML against `templates/config.schema.json`.
6. Run `install_dependencies.R` in check mode.
7. Run `run_all.R`; do not manually execute downstream scripts out of order.
8. Read `advanced_module_status.csv`, `figure_generation_status.csv` and `FINAL_STATUS.txt`.
9. When a module is `NOT_EVALUABLE`, report the missing prerequisite rather than editing thresholds to force output.
10. When a stage is `FAILED`, stop and inspect its log; do not continue with fabricated objects.
11. Review each independent proof PNG directly.
12. Report methods and conclusions only from generated Source Data and QA-approved outputs.

## 22. Local validation workflow

Before calling the package production-ready:

1. install the declared R dependencies;
2. freeze the environment with `freeze_environment.R` or `renv`;
3. run the core smoke dataset;
4. follow `tests/PUBLIC_DATA_VALIDATION_PLAN.md` and run one public dataset per enabled advanced module;
5. compare output schemas and selected numerical checkpoints with expected results;
6. inspect every independent proof PNG;
7. record package versions, runtime, warnings and failed branches;
8. create golden-output checks only after scientific review.

Static validation alone is not runtime or visual validation.

## 23. Invocation examples

```text
Use $single-cell-publication-workflow to run a full sample-aware Seurat v5 analysis with provisional annotation. Export independent figures only.
```

```text
Use $single-cell-publication-workflow to analyse recurrent NMF programs in the malignant-cell subset across samples. Do not run hdWGCNA unless its sample gate passes.
```

```text
Use $single-cell-publication-workflow to run hdWGCNA in the declared cell population, preserve sample identity in metacells, and export every module display as a separate figure.
```

```text
Use $single-cell-publication-workflow to map the annotated scRNA reference to the spatial dataset, run RCTD when evaluable, and export one map per slice and cell type.
```

```text
Use $single-cell-publication-workflow to generate exploratory drug-signature reversal rankings from pseudobulk DE. Do not provide treatment recommendations.
```


```text
Use $single-cell-publication-workflow to run scTenifoldKnk virtual knockout for the predeclared target in the specified WT cell population. Run samples independently with repeated subsampling, export independent figures only, and report network effects as unsigned exploratory perturbation.
```
