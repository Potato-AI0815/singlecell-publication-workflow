# Validation matrix

Version: `v0.1.0-alpha`

Status labels:

- `PUBLIC_DATA_VALIDATED`: completed on a public dataset with output and QA review;
- `RUNTIME_VALIDATED`: completed in the target Windows/R environment;
- `PARTIALLY_VALIDATED`: important internal stages completed, but the module did not reach its final output contract;
- `IMPLEMENTED_NOT_VALIDATED`: adapter exists but has not completed its dedicated benchmark;
- `NOT_APPLICABLE`: dataset did not satisfy the scientific gate.

| Capability | Dataset / environment | Status | Evidence and boundary |
|---|---|---|---|
| Multi-sample 10x import | GSE160763, mouse cortex | PUBLIC_DATA_VALIDATED | 8 samples, mixed `genes.tsv.gz`/`features.tsv.gz`, 40,666 raw cells |
| Sample metadata manifest | GSE160763 | PUBLIC_DATA_VALIDATED | 4 conditions, 2 samples per condition |
| QC | GSE160763 | PUBLIC_DATA_VALIDATED | 38,215/40,666 cells retained |
| SCTransform v2 + RPCA | GSE160763 | PUBLIC_DATA_VALIDATED | End-to-end core run completed |
| Mouse cortex broad annotation | GSE160763 | PUBLIC_DATA_VALIDATED | Major cortex lineages recovered; expert review still required |
| Independent figure engine | GSE160763 | PUBLIC_DATA_VALIDATED | 18 independent figure directories, no composite figure |
| Pseudobulk replication gate | GSE160763 | PUBLIC_DATA_VALIDATED | `n=2` groups correctly returned `NOT_EVALUABLE` |
| Windows launcher and R exit propagation | Windows, R 4.5.3 | RUNTIME_VALIDATED | `PASS_WITH_WARNINGS` accepted only after final QA |
| Dense raw UMI matrix + annotation-table import | GSE132465 6-pair subset | PUBLIC_DATA_VALIDATED | 12 samples, 20,074 input cells |
| Core tumour annotation | GSE132465 | PUBLIC_DATA_VALIDATED | 17,391 QC cells, 14 clusters, broad tumour/microenvironment lineages |
| Paired pseudobulk | GSE132465 | PUBLIC_DATA_VALIDATED | Design `~ patient_id + condition`; 5 cell types evaluable and 5 correctly gated |
| Expression-inferred CNA / CopyKAT | GSE132465 | PUBLIC_DATA_VALIDATED | 12 samples completed; Windows forced to one CopyKAT core |
| Sample-wise CellChat | GSE132465 | PUBLIC_DATA_VALIDATED | 12 sample objects, 6,247 interactions, paired n=6 contrasts |
| Recurrent GeneNMF programs | GSE132465 tumour epithelial subset | PUBLIC_DATA_VALIDATED | 6 meta-programs; low-cell samples correctly excluded |
| hdWGCNA | GSE132465 tumour epithelial subset | PARTIALLY_VALIDATED | Metacells, soft-power diagnostics and 433 MB TOM completed; run stopped during module detection |
| Drug signature reversal | GSE132465 | IMPLEMENTED_NOT_VALIDATED | Directional MSigDB-derived resource was prepared; module not executed |
| scDblFinder | GSE160763/GSE132465 | IMPLEMENTED_NOT_VALIDATED | Graceful fallback observed when unavailable; successful doublet run not yet benchmarked |
| SingleR reference annotation | GSE160763/GSE132465 | IMPLEMENTED_NOT_VALIDATED | Offline/reference-retrieval fallback validated; successful reference run not yet benchmarked |
| Slingshot/tradeSeq | — | IMPLEMENTED_NOT_VALIDATED | Dedicated differentiation dataset required |
| Spatial/RCTD | — | IMPLEMENTED_NOT_VALIDATED | Dedicated Visium benchmark required |
| Virtual knockout | — | IMPLEMENTED_NOT_VALIDATED | WT-only discovery plus held-out real-KO validation required |
| H5AD/zellkonverter import | — | IMPLEMENTED_NOT_VALIDATED | Adapter exists; no public-data runtime test yet |
| Linux runtime | — | IMPLEMENTED_NOT_VALIDATED | Static code is platform-aware; end-to-end Linux run not yet recorded |

## Release interpretation

This matrix supports the description **core beta / project alpha**. It does not support a claim that every advanced module or input type is production-ready.
