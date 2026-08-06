# V4 upgrade report

## User requirement implemented

V4 removes all comprehensive or assembled figure outputs and treats every scientific graph as an independent publication artifact. It also completes the advanced-analysis layer requested for local public-data testing: CNA, sample-wise communication, Slingshot/tradeSeq, recurrent NMF programs, hdWGCNA, spatial transcriptomics and exploratory drug-response prediction.

## Advanced analysis modules

### Expression-inferred CNA

- CopyKAT and infercna backends;
- explicit normal reference labels;
- sample-specific execution and heatmaps;
- cell-level CNA class/score table;
- no DNA-confirmed language.

### CellChat

- one model per biological sample;
- minimum-cell gates by sample and cell group;
- interaction, pathway and network-edge tables;
- replicate-aware condition contrasts where evaluable;
- descriptive, non-causal interpretation.

### Slingshot and tradeSeq

- explicit subset and root cluster;
- fitted curve geometry, nodes, lineage weights and pseudotime;
- optional knot evaluation and negative-binomial GAM fitting;
- association, start-versus-end, pattern and differential-end tests;
- no hand-drawn arrows.

### Recurrent NMF programs

- sample splitting and cell caps;
- GeneNMF `multiNMF/getMetaPrograms` preferred;
- deterministic NMF fallback;
- recurrent meta-program genes, similarities, scores and sample/group summaries;
- per-program embeddings and independent diagnostic figures.

### hdWGCNA

- prespecified coherent population;
- sample-preserving metacells;
- soft-power testing and selected-power record;
- signed network, non-grey modules, eigengenes, kME and hub genes;
- sample-level module activity, differential eigengenes and trait correlations where evaluable;
- independent network-diagnostic, module and hub figures.

### Spatial transcriptomics

- Seurat spatial object or Space Ranger input;
- multi-image/FOV coordinate extraction;
- slice-aware normalization, clustering and spatial-variable genes;
- reference label transfer and confidence scores;
- optional RCTD run independently by spatial sample;
- neighborhood enrichment calculated within slice;
- every cluster, score, weight and gene map exported separately by slice.

### Exploratory drug response

- transparent disease-signature reversal;
- signed drug-gene signature scoring;
- sample × cell-type aggregation;
- exploratory therapeutic transcriptomic clusters;
- oncoPredict using explicitly supplied training matrices and pseudobulk test expression;
- mandatory non-clinical interpretation boundary.

## Figure-engine upgrades

- independent-only exporter rejects patchwork and ggarrange classes;
- component-aware embedding hulls;
- rasterized dense cell layers with vector annotations;
- split violin/boxplot, dotplot, heatmap, trajectory and pseudotime templates;
- advanced CNA, communication, NMF, hdWGCNA, spatial and drug-response templates;
- robust heatmap duplicate aggregation and clustering;
- finite-value checks and recorded display limits;
- one manifest row and one Source Data set per graph.

## Lower-capability-model safeguards

- fixed execution order through `run_all.R`;
- fixed recipes and module scripts;
- full and example configuration templates;
- per-module configuration fragments;
- explicit `NOT_RUN`, `NOT_EVALUABLE` and `FAILED` states;
- final QA blocks completion after an enabled module crash or missing completed-module artifacts;
- aggregate-layout figures are prohibited in config, code and QA.

## Remaining validation boundary

The build container lacks R. V4 has passed static/configuration checks but has not been executed against real or synthetic R data in this environment. Local public-data validation is therefore the next required gate, not optional polish. Follow `tests/PUBLIC_DATA_VALIDATION_PLAN.md` and freeze the successful local environment before declaring the package runtime validated.
