# Public-data runtime validation plan

## Purpose

Use public datasets to validate execution, compatibility, output contracts and visual quality. This plan does not validate biological truth in every tissue. Test one module at a time before enabling several modules together.

## Global acceptance gate

For every run:

1. freeze R/package versions with `renv`;
2. keep an untouched input copy and record SHA-256;
3. validate the completed YAML before running;
4. require `FINAL_STATUS.txt` to be `PASS` or an explicitly reviewed `PASS_WITH_WARNINGS`;
5. require no requested figure to have `FAILED` status;
6. confirm every exported figure has PDF, 600-dpi TIFF/PNG, proof PNG, plot RDS, parameter YAML and Source Data;
7. verify that each figure directory contains exactly one scientific plot and that no aggregate-layout figure exists;
8. inspect each proof PNG individually at final display size;
9. rerun with the same seed and compare result manifests, dimensions, row counts and selected numerical checkpoints;
10. record failures without changing thresholds after observing the desired biological direction.

## Phase 0 — core engine

Dataset characteristics: a small human PBMC multi-sample count matrix with at least two biological samples and known broad immune lineages.

Validate import, metadata mapping, per-sample QC, doublet handling, normalization, integration, clustering, provisional annotation, composition, pseudobulk evaluability and all basic independent figures. Formal differential testing requires adequate biological replication; otherwise the correct result is `NOT_EVALUABLE`.

## Phase 1 — NMF programs

Dataset characteristics: multi-patient tumour or immune data containing one coherent cell population with at least two well-populated samples.

Run the NMF module alone. Review rank diagnostics when available, component-program similarity, recurrent meta-program genes, sample coverage, group activity and per-program embeddings. Reject a meta-program that is driven by one sample, technical stress, mitochondrial/ribosomal content or an obviously mixed lineage unless it is explicitly retained as an artefact program.

## Phase 2 — hdWGCNA

Dataset characteristics: at least three biological samples with at least 100 cells per sample in one prespecified cell population.

Run hdWGCNA alone. Validate metacells preserve sample identity, the soft-power table is finite, network construction returns non-grey modules, module eigengenes map back to cells, hub genes are members of their modules, and sample-level eigengene summaries are present. Treat differential module activity and trait correlations as exploratory unless replicate and design requirements are met.

## Phase 3 — trajectory and tradeSeq

Dataset characteristics: a biologically supported differentiation system with a defensible root state and sufficient intermediate cells.

Declare subset, cluster column and root before running. Validate Slingshot curves and arrows arise from the fitted object, not manual geometry. Check pseudotime ordering against known markers. For tradeSeq, inspect knot diagnostics, fit status and gene-test tables. Branch-specific tests are only evaluable when more than one fitted lineage exists.

## Phase 4 — sample-wise CellChat

Dataset characteristics: replicated conditions with multiple supported cell types in each biological sample.

Run CellChat independently per sample. Confirm no pooled pseudo-replicate object is used for condition inference. Inspect per-sample interaction tables, condition summaries, network edges and contrast tables. A model-derived communication probability is descriptive and does not establish physical interaction or causality.

## Phase 5 — expression-inferred CNA

Dataset characteristics: tumour scRNA-seq with explicit normal reference populations in each evaluable sample and, ideally, orthogonal CNA information.

Test CopyKAT first; test infercna separately. Confirm each sample has sufficient normal reference cells, tumour candidates are not used as references, and every CNA heatmap is sample-specific. Compare inferred classes with known tumour/normal labels or orthogonal CNA data. Do not describe RNA-inferred CNA as DNA-confirmed.

## Phase 6 — spatial transcriptomics

Dataset characteristics: one or more Visium/Visium-HD slices plus a compatible scRNA-seq reference.

Validate each slice independently. Check coordinate orientation, image/FOV identity, spatial clustering, label-transfer confidence, RCTD weights, spatial-variable genes and neighborhood enrichment. Confirm maps from different slices are never overlaid in one coordinate system. RCTD and label-transfer results require reference compatibility and should not be interpreted as direct cell counts without platform-specific validation.

## Phase 7 — drug-response hypotheses

Run three backends separately:

- signature reversal using a transparent perturbation signature table;
- signed gene-set scores using up/down signatures;
- oncoPredict using explicitly supplied training expression and response matrices.

Confirm gene identifiers overlap, orientation of training matrices is correct, prediction outputs are sample × cell-type pseudobulk for oncoPredict, and all figure titles say exploratory. Never translate a score into a clinical treatment recommendation.


## Phase 8 — virtual knockout

First run a small synthetic smoke dataset to verify package/API compatibility and all output schemas. Then use a public dataset containing a real WT and knockout comparison, with the virtual-knockout model fitted only to WT/control cells. A recommended first benchmark is GEO `GSE167595`, which contains Ahr WT and KO intestinal stem-cell data and has been used as a scTenifoldKnk validation example.

Predeclare the target, coherent cell population and WT condition. Run each WT biological sample independently with repeated seeded subsampling, then build the cross-sample consensus. Validate:

- target-expression and sample-evaluability gates;
- stable gene universe and deterministic run manifest;
- repeat-level and sample-level recurrence;
- agreement of impacted genes/pathways with the real KO contrast at the level of overlap/rank enrichment, not signed expression prediction;
- independent single-figure output and figure-specific Source Data;
- explicit wording that manifold distances are unsigned and exploratory.

Do not use the real KO cells to tune target, subset, gene universe, thresholds or network parameters. Compare with the real KO only after the virtual analysis is frozen.

## Phase 9 — integrated stress test

After every module passes alone, enable compatible modules together on a moderate public dataset. Confirm advanced modules execute before figure export, module status is complete, every figure remains independent, file names are unique, memory use is acceptable, and rerunning does not overwrite the prior result directory.

## Visual review checklist

For every proof PNG: readable text at final size; no clipped labels; no hidden legend; stable colours across plots; finite colour limits; no overplotting that erases structure; no artificial trajectory arrows; spatial slices separated; heatmap rows/columns legible; and the plotted values match Source Data spot checks.
