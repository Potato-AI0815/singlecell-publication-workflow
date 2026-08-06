# Advanced-module configuration snippets

These fragments are copied into `advanced_modules:` in a full project configuration. They are not standalone pipeline configurations. A low-capability model must copy one fragment without renaming keys, fill every `REQUIRED` field, validate the complete configuration, and enable only the module being tested.

Recommended local order: NMF → hdWGCNA → trajectory/tradeSeq → CellChat → CNA → spatial → drug response → virtual knockout. Each module must first be tested alone so that failures are attributable.

## Platform notes added in v0.1.0-alpha

- CopyKAT is forced to `n.cores = 1` on Windows because its `mclapply` path does not support multiple cores.
- WGCNA may be installed but not loadable until Bioconductor packages `impute` and `preprocessCore` are available.
- scDblFinder may require an updated Bioconductor installation including `bluster`.
- `scTenifoldKnk` is not auto-installed or redistributed because a standard open-source license was not confirmed in the upstream repository.
- Resume runs reuse complete existing figures as `SKIPPED_EXISTING`; incomplete figure directories remain errors.
