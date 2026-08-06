# V5 output contract

## Objects

- `02_raw_standardized_seurat.rds`: imported count-level object.
- `03_full_qc_annotated.rds`: all cells with QC/doublet decisions.
- `03_primary_singlets.rds`: primary QC-passing analysis cells.
- `04_normalized_pca.rds`: normalized object with PCA.
- `05_clustered.rds`: accepted provisional clustering.
- `06_annotated.rds`: provisional labels and evidence fields.

## Core tables

Input hashes, package versions, design audit, QC decisions, resolution metrics, annotation evidence, marker results, sample-level composition, pseudobulk counts/results, enrichment, issue ledger and QA.

## Advanced-module directory

`17_advanced_modules/` contains module status, objects, parameters and Source Data for CNV, CellChat, Slingshot/tradeSeq, NMF, hdWGCNA, spatial analysis and exploratory drug response. A completed module must produce its declared core artifacts.

## Independent figure directories

Every figure has one dedicated folder under `09_main_figures/` or `10_extended_figures/` containing PDF, optional SVG, 600-dpi TIFF/PNG, 183-mm proof PNG, thumbnail, plot RDS, parameters YAML, Source Data, `visual_QA.csv` and `VISUAL_QA_STATUS.txt`.

The workflow generates no aggregate-layout figure, montage or contact sheet. `15_manifests/figure_export_manifest.csv` links each independent figure to its sidecars. `14_qa/figure_generation_status.csv` records `EXPORTED`, `NOT_EVALUABLE` or `FAILED`.

## Completion states

- `PASS`: all required checks pass.
- `PASS_WITH_WARNINGS`: no blocking failure, but documented limitations remain.
- `FAIL`: a required scientific-integrity, artifact or technical figure check failed.

Image generation alone does not constitute completion.


## Virtual-knockout completed-module contract

A completed virtual-knockout module must include target and subset evaluability audits, a frozen gene universe, run manifest, sample and cross-sample consensus tables, parameters, and all available manifold/network/stability/pathway Source Data. Each target × subset visualization is exported as an independent figure directory; absence of a non-evaluable optional pathway plot is allowed only when recorded in figure-generation status.
