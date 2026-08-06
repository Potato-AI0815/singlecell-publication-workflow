# Changelog

## v0.1.0-alpha — 2026-08-06

First GitHub-oriented alpha release.

### Validated core

- multi-sample 10x import, QC, SCTransform v2, RPCA, clustering and provisional annotation;
- sample-level metadata manifests;
- independent publication figure export and Source Data;
- pseudobulk replication safeguards;
- Windows one-click execution and R exit-code propagation;
- public-data runtime validation on GSE160763 and the six-pair GSE132465 subset.

### Validated advanced paths

- paired pseudobulk on GSE132465;
- sample-wise CopyKAT CNA;
- sample-wise CellChat with paired condition contrast;
- recurrent GeneNMF meta-programs.

### Partially validated

- hdWGCNA through metacells, soft-power diagnostics and TOM construction.

### Runtime hardening

- mixed `genes.tsv.gz` and `features.tsv.gz` compatibility;
- Seurat v5 `JoinLayers()` compatibility;
- tissue-specific annotation dictionaries are honored;
- human SingleR broad-label mapping rules added;
- lightweight grob-backed plot RDS proxy prevents Seurat environments from producing multi-gigabyte sidecars;
- resume mode safely reuses complete existing figures as `SKIPPED_EXISTING`;
- historical failures remain in the issue ledger but no longer block a repaired resume run;
- CopyKAT is forced to one core on Windows;
- explicit CopyKAT namespace repair utility added;
- CellChat/self-loop network edges render without identical-endpoint failure;
- dependency load audit includes WGCNA/scDblFinder transitive hints.

### Not yet validated

- completed hdWGCNA module output;
- Slingshot/tradeSeq;
- spatial/RCTD;
- drug-response execution;
- virtual-knockout benchmark;
- successful SingleR and scDblFinder benchmark runs;
- Linux end-to-end runtime;
- H5AD import runtime.
